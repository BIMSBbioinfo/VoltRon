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
#' \dontrun{
#' # Annotate based on images
#' visium_data <- annotateSpatialData(visium_data, use.image.only = TRUE)
#' 
#' # Annotate based on spatial plot
#' xenium_data <- annotateSpatialData(xenium_data, group.by = "clusters")
#' }
annotateSpatialData <- function(object, label = "annotation", assay = NULL, annotation_assay = "ROIAnnotation", use.image.only = FALSE, 
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
    if(!requireNamespace("ImageArray")){
      message("Please install ImageArray package to speed up visualization")
      img <- magick::image_read(img)
    } else{
      img <- ImageArray::createImageArray(img)
    }
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
  coords <- vrCoordinates(object, assay = assay)
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
                    p(style="font-size: 12px;", strong("Polygons"), " require at least 4 points"),
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
    
    # Zoom Event ####
    manageImageBrush(img, ranges, max.pixel.size, input, output, session)
    
    # Corner Events ####
    observeEvent(input$reset_btn, {
      selected_corners(data.frame(x = numeric(0), y = numeric(0)))
    })
    observeEvent(input$rmvlast_btn, {
      selected_corners(selected_corners()[-nrow(selected_corners()),])
    })
    manageSelectedCorners(selected_corners, img, ranges, max.pixel.size, input, output, session)

    # Region Events ####
    n <- length(segments)
    textboxes <- reactiveVal(if (n > 0) 1:n else numeric(0))
    if (n > 0) {
      segment_names <- as.list(segment_names)
      names(segment_names) <- paste0("region", 1:n)
      textbox_values <- do.call("reactiveValues", segment_names)
    } else {
      textbox_values <- reactiveValues()
    }
    
    # Dynamically generate UI for textboxes and remove buttons
    output$textbox_ui <- renderUI({
      lapply(textboxes(), function(i) {
        column(12,
               textInputwithButton(textinputId = paste0("region", i), label = paste0("Region ", i),
                                   buttoninputId = paste0("remove", i), value = isolate(textbox_values[[paste0("region", i)]]), 
                                   onclick = sprintf('Shiny.setInputValue("remove", %d)', i))
               
        )
      })
    })
  
    ### update Region text ####
    
    observe({
      lapply(textboxes(), function(i) {
        observeEvent(input[[paste0("region", i)]], {
          textbox_values[[paste0("region", i)]] <- isolate(input[[paste0("region", i)]])
        }, ignoreNULL = FALSE)
      })
    })
    
    ### remove a Region ####
    
    # Observe event to remove textbox when the button is clicked
    observeEvent(input$remove, {
      
      # remove one point
      selected_corners_list(selected_corners_list()[!(textboxes() == as.numeric(isolate(input$remove)))])
      
      # Update the reactive value to remove the textbox
      textboxes(setdiff(textboxes(), as.numeric(isolate(input$remove))))
      
      # Remove the value from textbox_values
      textbox_values[[paste0("region", as.numeric(input$remove))]] <- NULL
      
    }, ignoreInit = TRUE)
    
    ### add a Region ####
    
    # Observe event to add a new textbox
    observeEvent(input$addregion_btn, {
      
      # Polygon selection
      if(isolate(input$region_type == "Polygon")){
        if(nrow(selected_corners()) > 3){
          
          # add to region list
          selected_corners_list(c(selected_corners_list(), list(selected_corners())))
          
          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))
          
          # add buttons
          new_id <- if (length(textboxes()) == 0) 1 else max(textboxes()) + 1
          textboxes(c(textboxes(), new_id))
          textbox_values[[paste0("region", new_id)]] <- ""
        } else {
          showNotification("You must selected at least 4 points for each polygon!")
        }
      } 
      
      # Circle selection
      if(isolate(input$region_type == "Circle")){
        if(nrow(selected_corners()) == 2){
          
          # add to region list
          circle <- makeCircleData(selected_corners())
          selected_corners_list(c(selected_corners_list(), list(circle)))
          
          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))
          
          # add buttons
          new_id <- if (length(textboxes()) == 0) 1 else max(textboxes()) + 1
          textboxes(c(textboxes(), new_id))
          textbox_values[[paste0("region", new_id)]] <- ""
        } else {
          showNotification("You must selected only 2 points for each circle!")
        }
      }
    })
    # Image output ####
    output$image_plot <- renderPlot({
      
      ## get image ####
      start <- proc.time()
      zoom_info <- FromBoxToCrop(cbind(ranges$x, ranges$y), imageinfo = imginfo)
      print(proc.time() - start)
      start <- proc.time()
      img <- cropImage(img, zoom_info)
      print(proc.time() - start)
      start <- proc.time()
      g <- plotImage(img, max.pixel.size = max.pixel.size) + labs(title = "")
      if(!use.image.only){
        g_spatial_clone <- cloneLayer(g_spatial)
        g <- g + transformSpatialLayer(g_spatial_clone, img, ranges, max.pixel.size)
      }
      print(proc.time() - start)
      
      # visualize currently selected corners ####
      transformed_corners <- transformSelectedCorners(selected_corners(), img, ranges, max.pixel.size)
      if(nrow(transformed_corners) > 0){
        if(isolate(input$region_type == "Polygon")){
          g <- g +
            ggplot2::geom_polygon(aes(x = x, y = y), data = transformed_corners, alpha = input$alpha, color = "red")
        } else {
          circle <- makeCircleData(transformed_corners)
          g <- g +
            ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0), 
                                  data = circle, alpha = input$alpha, color = "red", fill = "red")         
        }
        g <- g +
          ggplot2::geom_point(aes(x = x, y = y), data = transformed_corners, color = "red", shape = 16) 
      }
      
      # visualize already selected regions ####
      transformed_corners_list <- lapply(selected_corners_list(), transformSelectedCorners, img, ranges, max.pixel.size)
      if(length(transformed_corners_list) > 0){
        for (i in 1:length(transformed_corners_list)){
          cur_corners <- transformed_corners_list[[i]]
          
          # visualize regions
          if(ncol(cur_corners) == 2){
            g <- g +
              ggplot2::geom_polygon(aes(x = x, y = y, group = "region"), data = cur_corners, alpha = input$alpha, color = "red") 
          } else {
            g <- g +
              ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0), 
                                    data = cur_corners, alpha = input$alpha, color = "red", fill = "red")
          }
          
          # visualize labels
          if(ncol(cur_corners) == 2){
            cur_corners <- data.frame(x = mean(cur_corners[,1]), y = max(cur_corners[,2]), region = paste("Region ", isolate(textboxes()[i])))
          } else {
            cur_corners <- data.frame(x = cur_corners[,1], y = cur_corners[,2] + cur_corners[,3], region = paste("Region ", isolate(textboxes()[i])))
          }
          g <- g +
            ggrepel::geom_label_repel(mapping = aes(x = x, y = y, label = region), data = cur_corners,
                                      size = 5, direction = "y", nudge_y = 6, box.padding = 0, label.padding = 1, seed = 1, color = "red")
          
        }
      }
      
      # return graph
      g
    })
    
    # Return values for the shiny app ####
    observe({
      if(length(selected_corners_list()) > 0){
        shinyjs::show(id = "done")
      } else {
        shinyjs::hide(id = "done")
      }
    })
    observeEvent(input$done, {
      
      # selected list
      selected_polygon_list <- selected_corners_list()
      
      # collect labels
      selected_label_list <- sapply(1:length(selected_polygon_list), function(i) input[[paste0("region",i)]])
      
      if(length(selected_corners_list()) == 0){
        showNotification("You have not annotated the data yet!")
      } else if(any(selected_label_list == "")) {
        showNotification("Some regions have blank annotations (empty labels!)")
      } else {
        
        ### annotate spatial points ####
        if(inherits(metadata, "data.table")){
          spatialpoints <- metadata$id
        } else {
          spatialpoints <- rownames(metadata)
        }
        
        new_label <- rep("undefined", length(spatialpoints))
        names(new_label) <- spatialpoints
        result_list <- list()
        for(i in 1:length(selected_polygon_list)){
          cur_poly <- selected_polygon_list[[i]]
          if(ncol(cur_poly) > 2){
            in.list <- point.in.circle(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2], cur_poly[,3])
          } else {
            in.list <- sp::point.in.polygon(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2])
          }
          new_label[rownames(coords)[!!in.list]] <- selected_label_list[i]
        }
        
        # place annotation to metadata
        metadata[[label]] <- new_label
        Metadata(object, assay = sample_metadata[assay, "Assay"]) <- metadata
        
        # add polygons to a new assay ####
        segments <- list()
        for(i in 1:length(selected_label_list)){
          segments[[selected_label_list[i]]] <- data.frame(id = i, selected_polygon_list[[i]])
        }
        coords <- t(sapply(segments, function(seg){
          apply(seg[,c("x", "y")], 2, mean)
        }, simplify = TRUE))
        new_assay <- formAssay(coords = coords, segments = segments,
                               type = "ROI",
                               image = vrImages(object, assay = assay),
                               main_image = vrMainImage(object[[assay]]),
                               name = assay)
        object <- addAssay.VoltRon(object,
                                   assay = new_assay,
                                   metadata = data.frame(check.rows = FALSE, row.names = rownames(coords)),
                                   assay_name = annotation_assay,
                                   sample = sample_metadata[assay, "Sample"],
                                   layer = sample_metadata[assay, "Layer"])
        
        # stop app and return
        stopApp(object)
      }
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

manageImageBrush <- function(image, ranges, max.pixel.size, input, output, session){
  imginfo <- getImageInfo(image)
  observeEvent(input$plot_dblclick, {
    brush <- isolate(input$plot_brush)
    if (!is.null(brush)) {
      
      # get brush
      brush_mat <- data.frame(x = c(brush$xmin, brush$xmax), 
                              y = c(brush$ymin, brush$ymax))
      
      # if width is large, then correct the brush event for the downsize (scaling) effect
      limits <- data.frame(x = ranges$x, y = ranges$y)
      width <- limits[2,1]-limits[1,1]
      height <- limits[2,2]-limits[1,2]
      if(max(height,width) > max.pixel.size){
        if(inherits(image, "Image_Array")){
          n.series <- ImageArray::len(image)
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
      brush_mat[1,1] <- floor(brush_mat[1,1])
      brush_mat[1,2] <- floor(brush_mat[1,2])
      brush_mat[2,1] <- ceiling(brush_mat[2,1])
      brush_mat[2,2] <- ceiling(brush_mat[2,2])

      # update ranges
      ranges$x <- brush_mat[,1]
      ranges$y <- brush_mat[,2]
      
    } else {
      ranges$x <- c(0, imginfo$width)
      ranges$y <- c(0, imginfo$height)
    }
  })
}

manageSelectedCorners <- function(selected_corners, image, ranges, max.pixel.size, input, output, session){
  observeEvent(input$plot_click, {
    brush <- isolate(input$plot_brush)
    if(is.null(brush)) {
      
      # if a circle, dont allow more than two points
      if(isolate(input$region_type == "Circle")){
        if(nrow(selected_corners()) == 2){
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))
        }
      }
      
      # get points
      click <- input$plot_click
      pts <- c(click$x, click$y)

      # correct for zoom effect
      limits <- data.frame(x = ranges$x, y = ranges$y)
      width <- limits[2,1]-limits[1,1]
      height <- limits[2,2]-limits[1,2]
      if(max(height,width) > max.pixel.size){
        if(inherits(image, c("Image_Array"))){
          n.series <- ImageArray::len(image)
          cur_width <- width
          cur_height <- height
          for(ii in 2:n.series){
            cur_width <- width/(2^(ii-1))
            cur_height <- height/(2^(ii-1))
            if(max(cur_height, cur_width) <= max.pixel.size){
              break
            }
          }
          pts <- pts*width/ceiling(cur_width)
        } else {
          pts <- pts*width/max.pixel.size
        }
      }
    
      # correct for offset effect
      pts[1] <- pts[1] + limits[1,1]
      pts[2] <- pts[2] + limits[1,2]
      
      # Append new point to the data frame
      new_point <- data.frame(x = pts[1], y = pts[2])
      
      # add new point
      selected_corners(rbind(selected_corners(), new_point))
    }
  })
}

transformSelectedCorners <- function(selected_corners, image, ranges, max.pixel.size){
  
  # if empty, return
  if(nrow(selected_corners) == 0){
    return(selected_corners)
  }
  selected_corners <- as.matrix(selected_corners)
  
  # circle settings
  circle <- FALSE
  if(ncol(selected_corners) > 2)
    circle <- TRUE
  
  # get circle radius 
  if(circle){
    selected_radius <- selected_corners[,3:4, drop = FALSE]
    selected_corners <- selected_corners[,1:2, drop = FALSE]
  }
  
  # get image info
  imginfo <- getImageInfo(image)
  
  # correct for offset effect
  limits <- cbind(ranges$x, ranges$y)
  selected_corners <- selected_corners - matrix(rep(limits[1,], nrow(selected_corners)), nrow = nrow(selected_corners), byrow = T)

  # correct for zoom effect
  width <- limits[2,1]-limits[1,1]
  height <- limits[2,2]-limits[1,2]
  if(max(height,width) > max.pixel.size){
    if(inherits(image, "Image_Array")){
      n.series <- ImageArray::len(image)
      cur_width <- width
      cur_height <- height
      for(ii in 2:n.series){
        cur_width <- width/(2^(ii-1))
        cur_height <- height/(2^(ii-1))
        if(max(cur_height, cur_width) <= max.pixel.size){
          break
        }
      }
      selected_corners <- selected_corners*ceiling(cur_width)/width
      if(circle){
        selected_radius <- selected_radius*ceiling(cur_width)/width
      }
    } else {
      selected_corners <- selected_corners*max.pixel.size/width
      if(circle){
        selected_radius <- selected_radius*max.pixel.size/width
      }
    }
  }
  
  # make data frame
  if(circle){
    selected_corners <- as.data.frame(cbind(selected_corners,selected_radius))
    colnames(selected_corners) <- c("x", "y", "rx", "ry")
  } else {
    selected_corners <- as.data.frame(selected_corners)
    colnames(selected_corners) <- c("x", "y") 
  }
  
  # return
  selected_corners
}

#' @importFrom stats na.omit
transformSpatialLayer <- function(g_spatial, image, ranges, max.pixel.size){
  
  # correct for zoom
  ind <- (g_spatial$data$x > ranges$x[1] & g_spatial$data$x < ranges$x[2]) & (g_spatial$data$y > ranges$y[1] & g_spatial$data$y < ranges$y[2])
  g_spatial$data <- g_spatial$data[ind,]
  
  # is.na
  g_spatial$data <- stats::na.omit(g_spatial$data)
  
  # correct for zoom effect
  limits <- data.frame(x = ranges$x, y = ranges$y)
  width <- limits[2,1]-limits[1,1]
  height <- limits[2,2]-limits[1,2]
  if(max(height,width) > max.pixel.size){
    if(inherits(image, "Image_Array")){
      n.series <- ImageArray::len(image)
      cur_width <- width
      cur_height <- height
      for(ii in 2:n.series){
        cur_width <- width/(2^(ii-1))
        cur_height <- height/(2^(ii-1))
        if(max(cur_height, cur_width) <= max.pixel.size){
          break
        }
      }
      g_spatial$data$x <- g_spatial$data$x*ceiling(cur_width)/width
      g_spatial$data$y <- g_spatial$data$y*ceiling(cur_width)/width
    } else {
      g_spatial$data$x <- g_spatial$data$x*max.pixel.size/width
      g_spatial$data$y <- g_spatial$data$y*max.pixel.size/width
    }
  }
  
  # correct for offset effect
  g_spatial$data$x <- g_spatial$data$x - min(g_spatial$data$x)
  g_spatial$data$y <- g_spatial$data$y - min(g_spatial$data$y)
  
  # return
  g_spatial
}

####
# Shiny Utilities ####
####

#' internal Text input with button
#' 
#' Reproduced since it is not exported in the Shiny namespace.
#' 
#' @importFrom shiny tags validateCssUnit
#' 
#' @keywords internal
textInputwithButton <- function (textinputId, label, buttoninputId, value = "", width = NULL, placeholder = NULL, ...) 
{
  textvalue <- restoreInput(id = textinputId, default = value)
  buttonvalue <- restoreInput(id = buttoninputId, default = NULL)
  div(class = "form-group shiny-input-container", 
      style =  css(width = shiny::validateCssUnit(width), display = "inline-block"),
      shinyInputLabel(textinputId, label), 
      shiny::tags$input(id = textinputId, 
                        style = css(width = "80%", float = "left"),
                        type = "text", class = "shiny-input-text form-control", 
                        value = textvalue, placeholder = placeholder),
      shiny::tags$button(id = buttoninputId, 
                         style = css(width = "20%", float = "left"),
                         type = "button", class = "btn btn-default action-button", 
                         `data-val` = buttonvalue, disabled = NULL, list(shiny::icon("trash")), ...)
  )
}

#' Shiny's internal \code{shinyInputLabel} function
#' 
#' Reproduced since it is not exported in the Shiny namespace.
#' 
#' @importFrom shiny tags
#' 
#' @keywords internal
shinyInputLabel <- function(inputId, label=NULL) {
  shiny::tags$label(label,
                    class = "control-label",
                    class = if (is.null(label)) "shiny-label-null",
                    `for` = inputId
  )
}

####
# Auxiliary ####
####

#' makeCircleData
#' 
#' creates circle ROI from selected corners
#' 
#' @keywords internal
makeCircleData <- function(selected_corners) {
  radius <- sqrt((selected_corners[1,1] - selected_corners[2,1])^2 + (selected_corners[1,2] - selected_corners[2,2])^2)
  data.frame(x = selected_corners[1,1], y = selected_corners[1,2], rx = radius, ry = radius)
}

#' point.in.circle
#' 
#' checks if points are in the defined circle ROI with a (x,y) center and radius
#' 
#' @keywords internal
point.in.circle <- function(x,y,circle.x, circle.y, radius) {
  x <- x - circle.x
  y <- y - circle.y
  distall <- sqrt((x - circle.x)^2 + (y - circle.y)^2)
  return(distall < radius)
}

####
# Main Shiny App (old) ####
####

#' annotateSpatialData_old
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
#' \dontrun{
#' # Annotate based on images
#' visium_data <- annotateSpatialData_old(visium_data, use.image.only = TRUE)
#' 
#' # Annotate based on spatial plot
#' xenium_data <- annotateSpatialData_old(xenium_data, group.by = "clusters")
#' }
annotateSpatialData_old <- function(object, label = "annotation", assay = NULL, annotation_assay = "ROIAnnotation", use.image.only = FALSE, 
                                shiny.options = list(launch.browser = getOption("shiny.launch.browser", interactive())), image_name = NULL, channel = NULL, ...) {
  
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")
  
  ## Importing images ####
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # get assay names, and always get a single assay
  assay_names <- vrAssayNames(object, assay = assay)
  if(length(assay_names) > 0)
    assay <- assay_names[1]
  
  # metadata and coordinates
  metadata <- Metadata(object, assay = sample_metadata[assay, "Assay"])
  coords <- vrCoordinates(object, assay = assay)
  
  # set label names
  if(label %in% colnames(metadata)){
    unique_names <- make.unique(c(colnames(metadata)[grepl(paste0("^", label), colnames(metadata))], label))
    label <- unique_names[length(unique_names)]
  }
  
  # get image name and channel
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object[[assay]])
  }
  
  # get image
  max.pixel.size <- 1200
  img <- vrImages(object[[assay]], name = image_name, channel = channel, as.raster = TRUE)
  if(!inherits(img, "Image_Array")){
    img <- magick::image_read(img)
  }
  imginfo <- getImageInfo(img)

  # get spatial plot
  if(!use.image.only){
    g_spatial <- vrSpatialPlot(object, assay = assay, background = c(image_name, channel), scale.image = FALSE, ...)
    g_spatial <- g_spatial$layers[[2]]
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
    
    # initialize annotation segments 
    selected_corners_list <- reactiveVal(segments)
    selected_corners <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
    
    ## Point click, double click event and zoom ####
    # ranges <- reactiveValues(x = g$coordinates$limits$x, y = g$coordinates$limits$y)
    ranges <- reactiveValues(x = c(0, imginfo$width), y = c(0, imginfo$height))
    # ranges <- reactiveValues(x = c(0, max.pixel.size), 
    #                          y = c(0, ceiling(imginfo$height*(max.pixel.size/imginfo$width))))
    observeEvent(input$plot_dblclick, {
      brush <- isolate(input$plot_brush)
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        # ranges$x <- g$coordinates$limits$x
        # ranges$y <- g$coordinates$limits$y
        ranges$x <- c(0, imginfo$width)
        ranges$y <- c(0, imginfo$height)
        # ranges$x <- c(0, max.pixel.size)
        # ranges$y <- c(0, ceiling(imginfo$height*(max.pixel.size/imginfo$width)))
      }
    })
    observeEvent(input$plot_click, {
      brush <- isolate(input$plot_brush)
      if (is.null(brush)) {
        click <- input$plot_click
        x <- click$x
        y <- click$y
        
        # Append new point to the data frame
        new_point <- data.frame(x = x, y = y)
        
        # if a circle, dont allow more than two points
        if(isolate(input$region_type == "Circle")){
          if(nrow(selected_corners()) == 2){
            selected_corners(data.frame(x = numeric(0), y = numeric(0)))
          }
        }
        
        # add new point
        selected_corners(rbind(selected_corners(), new_point))
      }
    })
    
    ## Point buttons ####
    
    # reset corners
    observeEvent(input$reset_btn, {
      selected_corners(data.frame(x = numeric(0), y = numeric(0)))
    })
    
    # remove last point 
    observeEvent(input$rmvlast_btn, {
      selected_corners(selected_corners()[-nrow(selected_corners()),])
    })
    
    ## Region Annotators ####
    
    # length of segments
    n <- length(segments)
    
    ### Text Box Management ####
    
    # Reactive value to store the number of textboxes
    textboxes <- reactiveVal(if (n > 0) 1:n else numeric(0))
    
    # Initialize textbox values if n > 0, get already existing segments
    if (n > 0) {
      segment_names <- as.list(segment_names)
      names(segment_names) <- paste0("region", 1:n)
      textbox_values <- do.call("reactiveValues", segment_names)
    } else {
      textbox_values <- reactiveValues()
    }
    
    # Dynamically generate UI for textboxes and remove buttons
    output$textbox_ui <- renderUI({
      lapply(textboxes(), function(i) {
        column(12,
               textInputwithButton(textinputId = paste0("region", i), label = paste0("Region ", i),
                                   buttoninputId = paste0("remove", i), value = isolate(textbox_values[[paste0("region", i)]]), 
                                   onclick = sprintf('Shiny.setInputValue("remove", %d)', i))
               
        )
      })
    })
    
    # Observe changes in each textbox to update their values
    observe({
      lapply(textboxes(), function(i) {
        observeEvent(input[[paste0("region", i)]], {
          textbox_values[[paste0("region", i)]] <- isolate(input[[paste0("region", i)]])
        }, ignoreNULL = FALSE)
      })
    })
    
    ### Remove a Region ####
    
    # Observe event to remove textbox when the button is clicked
    observeEvent(input$remove, {
      
      # remove one point
      selected_corners_list(selected_corners_list()[!(textboxes() == as.numeric(isolate(input$remove)))])
      
      # Update the reactive value to remove the textbox
      textboxes(setdiff(textboxes(), as.numeric(isolate(input$remove))))
      
      # Remove the value from textbox_values
      textbox_values[[paste0("region", as.numeric(input$remove))]] <- NULL
      
    }, ignoreInit = TRUE)
    
    ### Add a Region ####
    
    # Observe event to add a new textbox
    observeEvent(input$addregion_btn, {
      
      # Polygon selection
      if(isolate(input$region_type == "Polygon")){
        if(nrow(selected_corners()) > 3){
          # add to region list
          selected_corners_list(c(selected_corners_list(), list(selected_corners())))
          
          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))
          
          # add buttons
          new_id <- if (length(textboxes()) == 0) 1 else max(textboxes()) + 1
          textboxes(c(textboxes(), new_id))
          textbox_values[[paste0("region", new_id)]] <- ""
        } else {
          showNotification("You must selected at least 4 points for each polygon!")
        }
      } 
      
      # Circle selection
      if(isolate(input$region_type == "Circle")){
        if(nrow(selected_corners()) == 2){
          
          # add to region list
          circle <- makeCircleData(selected_corners())
          selected_corners_list(c(selected_corners_list(), list(circle)))
          
          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))
          
          # add buttons
          new_id <- if (length(textboxes()) == 0) 1 else max(textboxes()) + 1
          textboxes(c(textboxes(), new_id))
          textbox_values[[paste0("region", new_id)]] <- ""
        } else {
          showNotification("You must selected only 2 points for each circle!")
        }
      }
    })
    
    ## image output ####
    
    output$image_plot <- renderPlot({
      
      # get image and plot
      g <- plotImage(img, max.pixel.size = NULL) + labs(title = "")
      if(!use.image.only){
        g <- g + g_spatial 
      }
      
      # visualize already selected polygons
      if(length(selected_corners_list()) > 0){
        for (i in 1:length(selected_corners_list())){
          cur_corners <- selected_corners_list()[[i]]
          if(ncol(cur_corners) == 2){
            g <- g +
              ggplot2::geom_polygon(aes(x = x, y = y, group = "region"), data = cur_corners, alpha = input$alpha, color = "red") 
          } else {
            g <- g +
              ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0), data = cur_corners, alpha = input$alpha, color = "red", fill = "red")
          }
        }
      }
      
      # add currently selected points
      g <- g +
        ggplot2::geom_point(aes(x = x, y = y), data = selected_corners(), color = "red", shape = 16) +
        coord_equal(xlim = ranges$x, ylim = ranges$y, ratio = 1) 
      
      # add label to currently selected points
      if(isolate(input$region_type == "Polygon")){
        g <- g +
          ggplot2::geom_polygon(aes(x = x, y = y), data = selected_corners(), alpha = input$alpha, color = "red")
      } else {
        if(nrow(selected_corners()) == 2){
          circle <- makeCircleData(selected_corners())
          g <- g +
            ggforce::geom_ellipse(aes(x0 = as.numeric(x), y0 = as.numeric(y), a = as.numeric(rx), b = as.numeric(ry), angle = 0), 
                                  data = circle, alpha = input$alpha, color = "red", fill = "red") +
            coord_equal(xlim = ranges$x, ylim = ranges$y, ratio = 1)
        }
      }
      
      # put labels of the already selected polygons
      if(length(selected_corners_list()) > 0){
        for (i in 1:length(selected_corners_list())){
          cur_corners <- selected_corners_list()[[i]]
          
          # create label data
          if(ncol(cur_corners) == 2){
            cur_corners <- data.frame(x = mean(cur_corners[,1]), y = max(cur_corners[,2]), region = paste("Region ", isolate(textboxes()[i])))
          } else {
            cur_corners <- data.frame(x = cur_corners[,1], y = cur_corners[,2] + cur_corners[,3], region = paste("Region ", isolate(textboxes()[i])))
          }
          
          # add label
          g <- g +
            ggrepel::geom_label_repel(mapping = aes(x = x, y = y, label = region), data = cur_corners,
                                      size = 5, direction = "y", nudge_y = 6, box.padding = 0, label.padding = 1, seed = 1, color = "red")
          
        }
      }
      
      # return graph
      g
    })
    
    ## Return values for the shiny app ####
    observe({
      if(length(selected_corners_list()) > 0){
        shinyjs::show(id = "done")
      } else {
        shinyjs::hide(id = "done")
      }
    })
    observeEvent(input$done, {
      
      # selected list
      selected_polygon_list <- selected_corners_list()
      
      # collect labels
      selected_label_list <- sapply(1:length(selected_polygon_list), function(i) input[[paste0("region",i)]])
      
      if(length(selected_corners_list()) == 0){
        showNotification("You have not annotated the data yet!")
      } else if(any(selected_label_list == "")) {
        showNotification("Some regions have blank annotations (empty labels!)")
      } else {
        
        ### annotate spatial points ####
        if(inherits(metadata, "data.table")){
          spatialpoints <- metadata$id
        } else {
          spatialpoints <- rownames(metadata)
        }
        
        new_label <- rep("undefined", length(spatialpoints))
        names(new_label) <- spatialpoints
        result_list <- list()
        for(i in 1:length(selected_polygon_list)){
          cur_poly <- selected_polygon_list[[i]]
          if(ncol(cur_poly) > 2){
            in.list <- point.in.circle(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2], cur_poly[,3])
          } else {
            in.list <- sp::point.in.polygon(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2])
          }
          new_label[rownames(coords)[!!in.list]] <- selected_label_list[i]
        }
        
        # place annotation to metadata
        metadata[[label]] <- new_label
        Metadata(object, assay = sample_metadata[assay, "Assay"]) <- metadata
        
        ## add polygons to a new assay ####
        segments <- list()
        for(i in 1:length(selected_label_list)){
          segments[[selected_label_list[i]]] <- data.frame(id = i, selected_polygon_list[[i]])
        }
        coords <- t(sapply(segments, function(seg){
          apply(seg[,c("x", "y")], 2, mean)
        }, simplify = TRUE))
        new_assay <- formAssay(coords = coords, segments = segments,
                               type = "ROI",
                               image = vrImages(object, assay = assay),
                               main_image = vrMainImage(object[[assay]]),
                               name = assay)
        object <- addAssay.VoltRon(object,
                                   assay = new_assay,
                                   metadata = data.frame(check.rows = FALSE, row.names = rownames(coords)),
                                   assay_name = annotation_assay,
                                   sample = sample_metadata[assay, "Sample"],
                                   layer = sample_metadata[assay, "Layer"])
        
        # stop app and return
        stopApp(object)
      }
    })
  }
  
  # configure options
  shiny.options <- configure_shiny_options(shiny.options)
  
  # run app
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