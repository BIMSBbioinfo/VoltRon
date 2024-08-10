####
# Spatial Interactive Plot (VoltRon) ####
####

####
## Background Shiny App ####
####

#' vrSpatialPlotInteractive
#'
#' @inheritParams shiny::runApp
#' @param plot_g the ggplot plot
#' @importFrom rstudioapi viewer
#'
#' @noRd
vrSpatialPlotInteractive <- function(host = getOption("shiny.host", "127.0.0.1"),
                                     port = getOption("shiny.port"), plot_g = NULL){
  shinyjs::useShinyjs()

  # UI
  ui <- mod_app_ui("app")

  # Server
  server <- function(input, output, session) {
    mod_app_server("app", plot_g = plot_g)
    session$onSessionEnded(function() {
      stopApp()
    })
  }

  # Start Shiny Application
  shiny::shinyApp(ui, server, options = list(host = host, port = port, launch.browser = rstudioapi::viewer),
                  onStart = function() {
                    cat("Doing application setup\n")
                    onStop(function() {
                      cat("Doing application cleanup\n")
                    })
                  })

}

#' App UI
#'
#' @param id id of the module
#'
#' @import shiny
#'
#' @noRd
mod_app_ui <- function(id) {
  ns <- NS(id)
  plotOutput(ns("image_plot"),
             height = "1000px",
             dblclick = ns("plot_dblclick"),
             brush = brushOpts(
               id = ns("plot_brush"), fill = "green",
               resetOnNew = TRUE
             ))
}

#' App Server
#'
#' @param id id of the module
#' @param plot_g the ggplot plot
#'
#' @import ggplot2
#'
#' @noRd
mod_app_server <- function(id, plot_g = NULL) {
  moduleServer(id, function(input, output, session) {

    ranges <- reactiveValues(x = plot_g$coordinates$limits$x, y = plot_g$coordinates$limits$y)
    observeEvent(input$plot_dblclick, {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- plot_g$coordinates$limits$x
        ranges$y <- plot_g$coordinates$limits$y
      }
    })

    # image output
    output$image_plot <- renderPlot({
      plot_g +
        ggplot2::coord_equal(xlim = ranges$x, ylim = ranges$y, ratio = 1)
    })
  })
}

####
# Spatial Interactive Plot (Vitessce) ####
####

#' vrSpatialPlotInteractive
#'
#' Interactive Plotting identification of spatially resolved cells, spots, and ROI on associated images from multiple assays in a VoltRon object.
#'
#' @param zarr.file The zarr file of a VoltRon object
#' @param group.by a grouping label for the spatial entities
#' @param reduction The name of the reduction to visualize an embedding alongside with the spatial plot.
#'
#' @noRd
vrSpatialPlotVitessce <- function(zarr.file, group.by = "Sample", reduction = NULL) {

  # check package
  if (!requireNamespace('vitessceR'))
    stop("Please install vitessceR package for using interactive visualization")

  # check file
  if(!dir.exists(zarr.file))
    stop(paste0(zarr.file, " is not found at the specified location!"))
  
  # get embedding
  if(is.null(reduction)){
    obs_embedding_paths <- c("obsm/spatial")
  } else {
    obs_embedding_paths <- c(paste0("obsm/", reduction), "obsm/spatial")
  }

  # initiate vitessceR
  vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.15", name = "MBrain")
  dataset <- vc$add_dataset("My dataset")
  
  # add ome tiff if exists
  ometiff.file <- gsub("zarr[/]?$", "ome.tiff", zarr.file)
  if(file.exists(ometiff.file)){
    w_img <- vitessceR::MultiImageWrapper$new(
      image_wrappers = list(
        vitessceR::OmeTiffWrapper$new(name="Test1", img_path=ometiff.file)
      )
    )
    dataset <- dataset$add_object(w_img)
  } 
  
  # add anndata
  w_data <- vitessceR::AnnDataWrapper$new(
    adata_path=zarr.file,
    obs_set_paths = c(paste0("obs/", group.by)),
    obs_set_names = c(group.by),
    obs_locations_path = "obsm/spatial",
    obs_segmentations_path = "obsm/segmentation",
    obs_embedding_paths = obs_embedding_paths
  )
  dataset <- dataset$add_object(w_data)

  # set up vitessce pane  
  spatial <- vc$add_view(dataset, vitessceR::Component$SPATIAL)
  cell_sets <- vc$add_view(dataset, vitessceR::Component$OBS_SETS)
  spatial_segmentation_layer_value <- list(opacity = 1, radius = 0, visible = TRUE, stroked = FALSE)
  spatial_layers <- vc$add_view(dataset, vitessceR::Component$LAYER_CONTROLLER)

  if(is.null(reduction)){
    vc$layout(
      vitessceR::hconcat(spatial, 
                         vitessceR::hconcat(cell_sets, spatial_layers))
    )
  } else {
    umap <- vc$add_view(dataset, vitessceR::Component$SCATTERPLOT, mapping = reduction)
    vc$layout(
      vitessceR::hconcat(spatial, 
                         vitessceR::vconcat(umap, 
                                            vitessceR::hconcat(cell_sets, spatial_layers)))
    )
  }

  vc$widget(theme = "light")
}