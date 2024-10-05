####
# Main Shiny App ####
####

#' annotateSpatialData
#'
#' A mini shiny app to for annotating spatial points
#'
#' @param object a VoltRon object
#' @param label the name of the new metadata column (default: annotation) annotating spatial points by selected polygons
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param annotation_assay name of the annotation assay ()
#' @param use.image.only if TRUE, use only the image
#' @param shiny.options a list of shiny options (launch.browser, host, port etc.) passed \code{options} arguement of \link{shinyApp}. For more information, see \link{runApp}
#' @param image_name the name/key of the image
#' @param channel the name of the main channel
#' @param ... additional parameters passed to \link{vrSpatialPlot}.
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs show hide
#' @importFrom stats median
#' @importFrom sp point.in.polygon
#' @import ggplot2
#'
#' @export
#' 
#' @examples
#' # Annotate based on images
#' visium_data <- annotateSpatialData(visium_data, use.image.only = TRUE)
#' 
#' # Annotate based on spatial plot
#' xenium_data <- annotateSpatialData(xenium_data, group.by = "clusters")
annotateSpatialData_new <- function(object, label = "annotation", assay = NULL, annotation_assay = "ROIAnnotation", use.image.only = FALSE, 
                                shiny.options = list(launch.browser = getOption("shiny.launch.browser", interactive())), 
                                image_name = NULL, channel = NULL, ...) {
  
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")
  
  ## Importing images ####
  
  # get assay names, and always get a single assay
  assay_names <- vrAssayNames(object, assay = assay)
  if(length(assay_names) > 0)
    assay <- assay_names[1]
  
  # get image name and channel
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object[[assay]])
  }
  
  # get image
  img <- vrImages(object[[assay]], name = image_name, channel = channel, as.raster = TRUE)
  if(!inherits(img, "Image_Array")){
    img <- magick::image_read(img)
  }
  if(!use.image.only){
    # get spatial plot
    g_spatial <- vrSpatialPlot(object, assay = assay, background = c(image_name, channel), scale.image = FALSE, ...)
    g_spatial <- g_spatial$layers[[2]]
  }
  
  # get image info
  max.pixel.size <- 1200
  imginfo <- getImageInfo(img)
  
  ## Get previous annotation ####
  
  # set label names
  sample_metadata <- SampleMetadata(object)
  metadata <- Metadata(object, assay = sample_metadata[assay, "Assay"])
  if(label %in% colnames(metadata)){
    unique_names <- make.unique(c(colnames(metadata)[grepl(paste0("^", label), colnames(metadata))], label))
    label <- unique_names[length(unique_names)]
  }
  
  # get segmentations (if exists) from the same layer
  if(!is.null(annotation_assay)){
    layer_metadata <- sample_metadata[sample_metadata$Layer == sample_metadata[assay, "Layer"] & sample_metadata$Sample == sample_metadata[assay, "Sample"],]
    if(annotation_assay %in% layer_metadata$Assay){
      
      # get segments
      segments <- vrSegments(object, assay = annotation_assay)
      segments <- lapply(segments, function(seg) seg[,colnames(seg)[colnames(seg) != "id"]])
      segment_names <- gsub("_Assay[0-9]+$", "", names(segments)) 
      
      # remove the latest annotation
      all_assay_names <- vrAssayNames(object, assay = "all")
      object <- subset(object, assays = all_assay_names[!all_assay_names %in% rownames(layer_metadata)[layer_metadata$Assay == annotation_assay]])
      
    } else {
      segments <- list()
      segment_names <- c()
    }
  }
  
  ## UI and Server ####
  
  # Define UI for the application
  ui <- fluidPage(
    sidebarLayout(position = "left",
                  
                  sidebarPanel(
                    
                    # margin settings
                    tags$style(make_css(list('.well', 'margin', '7%'))),
                    
                    # # specific settings for dealing with simultaneous click and brush events
                    # # https://jokergoo.github.io/2021/02/20/differentiate-brush-and-click-event-in-shiny/
                    tags$script(HTML("
            $('#plot').mousedown(function(e) {
                var parentOffset = $(this).offset();
                var relX = e.pageX - parentOffset.left;
                var relY = e.pageY - parentOffset.top;
                Shiny.setInputValue('x1', relX);
                Shiny.setInputValue('y1', relY);
            }).mouseup(function(e) {
                var parentOffset = $(this).offset();
                var relX = e.pageX - parentOffset.left;
                var relY = e.pageY - parentOffset.top;
                Shiny.setInputValue('x2', relX);
                Shiny.setInputValue('y2', relY);
                Shiny.setInputValue('action', Math.random());
            });
          ")),
                    
                    # Interface
                    fluidRow(
                      column(12,h4("Spatial Annotation")),
                      column(12,shiny::actionButton("reset_btn",     "Reset Points     ")),
                      column(12,shiny::actionButton("rmvlast_btn",   "Remove Last Point")),
                      column(12,shiny::actionButton("addregion_btn", "Add Region       ")),
                    ),
                    br(),
                    
                    fluidRow(
                      column(6,shiny::selectInput("region_type", label = "Region Type", choices = c("Polygon", "Circle"), selected = "Polygon")),
                      column(6,sliderInput("alpha", "Transparency", min = 0, max = 1, value = 0.2)),
                    ),
                    
                    # instructions
                    h4("How to use"),
                    p(style="font-size: 12px;", strong("Single-L-click"), " to select polygon or circle points"),
                    p(style="font-size: 12px;", strong("Add Region"), " to set points as a new region"),
                    p(style="font-size: 12px;", strong("Circles"), " require only 2 points"),
                    p(style="font-size: 12px;", strong("Polygons"), " require at least 3 points"),
                    br(),
                    
                    # Subsets
                    fluidRow(
                      column(12,h4("Selected Regions")),
                      br(),
                      uiOutput("textbox_ui"),
                      br()  
                    ),
                    
                    # Subsets
                    fluidRow(
                      column(12,shiny::actionButton("done", "Done"))
                    ),
                    
                    width = 4
                  ),
                  mainPanel(
                    shinyjs::useShinyjs(),
                    plotOutput("image_plot",
                               height = "1000px",
                               click = "plot_click",
                               dblclick = "plot_dblclick",
                               brush = brushOpts(
                                 id = "plot_brush", fill = "green",
                                 resetOnNew = TRUE
                               )),
                    width = 8
                  )
    )
  )
  
  # Define server logic required to create, add, and remove textboxes
  server <- function(input, output, session) {
    
    # Reactive values ####
    selected_corners_list <- reactiveVal(segments)
    selected_corners <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
    ranges <- reactiveValues(x = c(0, imginfo$width), y = c(0, imginfo$height))
    
    # Point click, double click event and zoom ####
    manageImageBrushOptions(img, ranges, max.pixel.size, input, output, session)
    
    # Point buttons ####
    # reset corners and remove last point 
    observeEvent(input$reset_btn, {
      selected_corners(data.frame(x = numeric(0), y = numeric(0)))
    })
    observeEvent(input$rmvlast_btn, {
      selected_corners(selected_corners()[-nrow(selected_corners()),])
    })
    
    # image output ####
    output$image_plot <- renderPlot({
      
      # get image and plot
      # img <- ImageArray::crop(img, ind = list(ranges$x,
      #                                         sort(imginfo$height - ranges$y, decreasing = FALSE)))
      zoom_info <- FromBoxToCrop(cbind(ranges$x, ranges$y), imageinfo = imginfo)
      img <- cropImage(img, zoom_info)
      g <- plotImage(img, max.pixel.size = max.pixel.size) + labs(title = "")
      if(!use.image.only){
        # g_spatial <- g_spatial + coord_equal(xlim = ranges$x, ylim = ranges$y, ratio = 1)
        g <- g + g_spatial 
      }
      
      # return graph
      g
    })
  }
  
  # Run App ####
  shiny.options <- configure_shiny_options(shiny.options)
  shiny::runApp(
    shiny::shinyApp(ui, server, options = list(host = shiny.options[["host"]], port = shiny.options[["port"]], launch.browser = shiny.options[["launch.browser"]]),
                    onStart = function() {
                      cat("Doing application setup\n")
                      onStop(function() {
                        cat("Doing application cleanup\n")
                      })
                    })
  )
}

####
# Annotation Utilities ####
####

manageImageBrushOptions <- function(image, ranges, max.pixel.size, input, output, session){
  imginfo <- getImageInfo(image)
  observeEvent(input$plot_dblclick, {
    brush <- isolate(input$plot_brush)
    if (!is.null(brush)) {
      
      # get brush
      brush_mat <- data.frame(x = c(brush$xmin, brush$xmax), 
                              y = c(brush$ymin, brush$ymax))

      # if width is large, then correct the brush event for the downsize effect
      limits <- data.frame(x = ranges$x, y = ranges$y)
      width <- limits[2,1]-limits[1,1]
      height <- limits[2,2]-limits[1,2]
      if(max(height,width) > max.pixel.size){
        if(inherits(img, "Image_Array")){
          n.series <- ImageArray::len(img)
          cur_width <- width
          cur_height <- height
          for(ii in 2:n.series){
            cur_width <- width/(2^(ii-1))
            cur_height <- height/(2^(ii-1))
            if(max(cur_height, cur_width) <= max.pixel.size){
              break
            }
          }
          brush_mat <- brush_mat*width/ceiling(cur_width)
        } else {
          brush_mat <- brush_mat*width/max.pixel.size
        }
      }
      
      # correct brush for the zoom effect
      brush_mat[,1] <- brush_mat[,1] + limits[1,1]
      brush_mat[,2] <- brush_mat[,2] + limits[1,2]
      # brush_mat[1,1] <- floor(brush_mat[1,1])
      # brush_mat[1,2] <- floor(brush_mat[1,2])
      # brush_mat[2,1] <- ceiling(brush_mat[2,1])
      # brush_mat[2,2] <- ceiling(brush_mat[2,2])

      # update ranges
      ranges$x <- brush_mat[,1]
      ranges$y <- brush_mat[,2]
      
    } else {
      ranges$x <- c(0, imginfo$width)
      ranges$y <- c(0, imginfo$height)
    }
  })
}