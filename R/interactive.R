####
# Spatial Interactive Plot (VoltRon) ####
####

####
## Background Shiny App ####
####

#' vrSpatialPlotInteractive
#'
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
#' @param plot.segments plot segments instead of points
#' @param group.ids a subset of categories defined with in the grouping label \code{group.by}
#' @param assay the assay name
#' @param reduction The name of the reduction to visualize an embedding alongside with the spatial plot.
#' @param background the background of the plot, either "image" for overlaying the image of the assay, or "black" or "white" background (suitable for IF based assays)
#' @param reg if TRUE, the registered coordinates will be used
#' @param crop whether to crop an image of a spot assay
#'
#' @noRd
vrSpatialPlotVitessce <- function(zarr.file, group.by = "Sample", plot.segments = FALSE, group.ids = NULL, assay = NULL, reduction = "umap",
                                  background = NULL, reg = FALSE, crop = FALSE) {

  # check package
  if (!requireNamespace('vitessceR'))
    stop("Please install vitessceR package for using interactive visualization")

  # get embedding
  if(is.null(reduction)){
    obs_embedding_paths <- c("obsm/spatial")
  } else {
    obs_embedding_paths <- c(paste0("obsm/", reduction), "obsm/spatial")
  }

  w <- vitessceR::AnnDataWrapper$new(
    adata_path=zarr.file,
    obs_set_paths = c(paste0("obs/", group.by)),
    obs_set_names = c(group.by),
    obs_locations_path = "obsm/spatial",
    obs_embedding_paths=obs_embedding_paths
  )
  vc <- vitessceR::VitessceConfig$new(schema_version = "1.0.15", name = "MBrain")
  dataset <- vc$add_dataset("My dataset")$add_object(w)
  spatial <- vc$add_view(dataset, vitessceR::Component$SCATTERPLOT, mapping = "spatial")
  cell_sets <- vc$add_view(dataset, vitessceR::Component$OBS_SETS)

  if(is.null(reduction)){
    vc$layout(
      vitessceR::hconcat(spatial, cell_sets)
    )
  } else {
    umap <- vc$add_view(dataset, vitessceR::Component$SCATTERPLOT, mapping = reduction)
    vc$layout(
      vitessceR::hconcat(spatial, vitessceR::vconcat(umap, cell_sets))
    )
  }

  vc$widget(theme = "light")
}

####
## Basilisk Environment ####
####

#' The Python Basilisk environment
#'
#' Defines a conda environment via Basilisk, which is used to convert R objects to Zarr stores.
#'
#' @importFrom basilisk BasiliskEnvironment
#'
#' @keywords internal
#'
#' @noRd
py_env <- basilisk::BasiliskEnvironment(
  envname="VoltRon_basilisk_env",
  pkgname="VoltRon",
  packages=c(
    "numpy==1.*",
    "pandas==1.*",
    "anndata==0.7.*",
    "h5py==3.*",
    "hdf5==1.*",
    "natsort==7.*",
    "packaging==20.*",
    "scipy==1.*",
    "sqlite==3.*",
    "zarr==2.*",
    "numcodecs==0.*"
  ),
  pip=c(
    "ome-zarr==0.2.1"
  )
)

####
## Conversion into Zarr for Vitessce ####
####

#' Title
#'
#' @param vrimage VoltRon image
#' @param out_path output path to ome.zarr
#' @param image_id image name
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @noRd
vrImage_to_zarr <- function (vrimage, out_path, image_id = "image_1")
{
  img_arr <- apply(as.matrix(magick::image_raster(vrimage, tidy = FALSE)), c(1, 2), col2rgb)
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(img_arr, image_id, out_path) {
    zarr <- reticulate::import("zarr")
    ome_zarr <- reticulate::import("ome_zarr")
    z_root <- zarr$open_group(out_path, mode = "w")
    obj_list <- function(...) {
      retval <- stats::setNames(list(), character(0))
      param_list <- list(...)
      for (key in names(param_list)) {
        retval[[key]] = param_list[[key]]
      }
      retval
    }
    default_window <- obj_list(start = 0, min = 0, max = 255, end = 255)
    ome_zarr$writer$write_image(image = img_arr,
                                group = z_root,
                                axes = "cyx",
                                omero = obj_list(name = image_id, version = "0.3",
                                                 rdefs = obj_list(),
                                                 channels = list(obj_list(label = "r", color = "FF0000", window = default_window),
                                                                 obj_list(label = "g", color = "00FF00", window = default_window),
                                                                 obj_list(label = "b", color = "0000FF", window = default_window))))
    return(TRUE)
  }, img_arr = img_arr, image_id = image_id, out_path = out_path)
  return(success)
}

vr_to_anndata_zarr <- function (object, out_path)
{
  # object data
  datax <- vrData(object, norm = FALSE)
  metadata <- Metadata(object)
  # feature.metadata <- vrFeatureData(object)

  # obsm
  obsm <- list()
  coords <- vrCoordinates(object)
  obsm[["spatial"]] <- t(as.matrix(coords))
  if (length(vrEmbeddingNames(object)) > 0) {
    for (embed_name in vrEmbeddingNames(object)) {
      embedding <- vrEmbeddings(object, type = embed_name)
      obsm[[embed_name]] <- t(as.matrix(embedding))
    }
  }

  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(datax, metadata, obsm, out_path) {
    anndata <- reticulate::import("anndata")
    zarr <- reticulate::import("zarr")
    make_numpy_friendly <- function(x, transpose = TRUE) {
      if (transpose) {
        x <- Matrix::t(x)
      }
      if (DelayedArray::is_sparse(x)) {
        methods::as(x, "dgCMatrix")
      }
      else {
        as.matrix(x)
      }
    }
    X <- make_numpy_friendly(datax)
    adata <- anndata$AnnData(X = X, obs = metadata)
    # adata <- anndata$AnnData(X = X, obs = metadata, var = feature.metadata)
    if (length(obsm) > 0) {
      obsm <- lapply(obsm, make_numpy_friendly)
      adata$obsm <- obsm
    }
    adata$write_zarr(out_path)
    return(TRUE)
  }, datax = datax, metadata = metadata, obsm = obsm, out_path = out_path)
  return(success)
}
