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
#' @param use.image if TRUE, use only the image
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
#' visium_data <- annotateSpatialData(visium_data, use.image = TRUE)
#' 
#' # Annotate based on spatial plot
#' xenium_data <- annotateSpatialData(xenium_data, group.by = "clusters")
annotateSpatialData <- function(object, label = "annotation", assay = NULL, annotation_assay = "ROIAnnotation", use.image = FALSE, image_name = NULL, channel = NULL, ...) {

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
    # image_name <- vrMainImage(object[[assay]])
    image_name <- vrMainSpatial(object[[assay]])
  }

  # get image
  if(use.image){
    g <- magick::image_ggplot(vrImages(object, assay = assay, name = image_name, channel = channel)) + labs(title = "")
  } else{
    g <- vrSpatialPlot(object, assay = assay, background = c(image_name, channel), scale.image = FALSE, ...) + labs(title = "")
  }
  
  # get segmentations (if exists)
  if(!is.null(annotation_assay)){
    layer_metadata <- sample_metadata[,sample_metadata$Layer == sample_metadata[assay, "Layer"] & sample_metadata$Sample == sample_metadata[assay, "Sample"]]
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

  # get the ui and server
  if (interactive()){
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
            column(12,h4("Annotation Interface")),
            column(12,shiny::actionButton("reset_btn",     "Reset Points     ")),
            column(12,shiny::actionButton("rmvlast_btn",   "Remove Last Point")),
            column(12,shiny::actionButton("addregion_btn", "Add Region       ")),
          ),
          br(),
          
          # Subsets
          fluidRow(
            column(12,h4("Selected Regions")),
            column(6,shiny::selectInput("region_type", label = "Region Type", choices = c("Polygon", "Circle"), selected = "Polygon")),
            column(6,sliderInput("alpha", "Transparency", min = 0, max = 1, value = 0.2)),
            br(),
            # getRegionAnnotators(segments, segment_names),
            uiOutput("textbox_ui"),
            br(),
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

    server <- function(input, output) {

      # initialize annotation segments 
      selected_corners_list <- reactiveVal(segments)
      counter <- reactiveValues(n = if(length(segments) > 0) seq_len(length(segments)) else c(), 
                                max = if(length(segments) > 0) length(segments) else NULL)
      
      # Initialize data frame to store points
      selected_corners <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
      ranges <- reactiveValues(x = g$coordinates$limits$x, y = g$coordinates$limits$y)

      ## point double click event and zoom ####
      observeEvent(input$plot_dblclick, {
        brush <- input$plot_brush
        if (!is.null(brush)) {
          ranges$x <- c(brush$xmin, brush$xmax)
          ranges$y <- c(brush$ymin, brush$ymax)
        } else {
          ranges$x <- g$coordinates$limits$x
          ranges$y <- g$coordinates$limits$y
        }
      })

      ## point click event ####
      observeEvent(input$plot_click, {
        brush <- input$plot_brush
        if (is.null(brush)) {
          click <- input$plot_click
          x <- click$x
          y <- click$y

          # Append new point to the data frame
          new_point <- data.frame(x = x, y = y)
          selected_corners(rbind(selected_corners(), new_point))
        }
      })

      ## point buttons ####
      
      # reset corners
      observeEvent(input$reset_btn, {
        selected_corners(data.frame(x = numeric(0), y = numeric(0)))
      })
      
      # remove last point 
      observeEvent(input$rmvlast_btn, {
        selected_corners(selected_corners()[-nrow(selected_corners()),])
      })

      # add region
      observeEvent(input$addregion_btn, {
        if(nrow(selected_corners()) > 3){

          # add to region list
          selected_corners_list(c(selected_corners_list(), list(selected_corners())))

          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))

          # Track input boxes to render
          if(is.null(counter$max)){
            counter$n <- 1
          } else {
            counter$n <- c(counter$n, counter$max + 1) 
          }
          counter$max <- max(counter$n)
        } else {
          showNotification("You must selected at least 4 points!")
        }
      })
    
      ## Region Annotators ####

      # counter for regions
      if(length(segment_names)){
        output$textbox_ui <- renderUI({
          lapply(counter$n, function(i) {
            column(12,
                   textInputwithButton(textinputId = paste0("region", i), label = paste0("Region ", i),
                                       buttoninputId = paste0("region_button", i), value = segment_names[i]))
          })
        })
        segment_names <- NULL
      } else {
        output$textbox_ui <- renderUI({textboxes()})
        textboxes <- reactive({
          if(length(counter$n) > 0){
            lapply(counter$n, function(i) {
              column(12,
                     textInputwithButton(textinputId = paste0("region", i), label = paste0("Region ", i), 
                                         buttoninputId = paste0("region_button", i), value = input[[paste0("region", i)]]))
            })
          }
        })
      }
      
      
      # remove region annotators when clicked on the button
      observe({
        if(length(counter$n) > 0){
          for (i in counter$n){
            observeEvent(input[[paste0("region_button", i)]], {

              # remove one point
              selected_corners_list(selected_corners_list()[-which(counter$n == i)])
              
              # remove one button
              counter$n <- setdiff(counter$n, i)
              
            })
          }
        }
      })

      ## image output ####
      output$image_plot <- renderPlot({

        # visualize already selected polygons
        if(length(selected_corners_list()) > 0){
          for (i in 1:length(selected_corners_list())){
            cur_corners <- selected_corners_list()[[i]]
            g <- g +
              ggplot2::geom_polygon(aes(x = x, y = y, group = "region"), data = cur_corners, alpha = input$alpha, color = "red")
          }
        }

        # add currently selected points
        g <- g +
          ggplot2::geom_point(aes(x = x, y = y), data = selected_corners(), color = "red", shape = 16) +
          coord_equal(xlim = ranges$x, ylim = ranges$y, ratio = 1)

        # add label to currently selected points
        datax_label_ind <- length(selected_corners_list()) + 1
        g <- g +
          ggplot2::geom_polygon(aes(x = x, y = y), data = selected_corners(), alpha = input$alpha, color = "red")

        # put labels of the already selected polygons
        if(length(selected_corners_list()) > 0){
          for (i in 1:length(selected_corners_list())){
            cur_corners <- selected_corners_list()[[i]]
            cur_corners <- data.frame(x = mean(cur_corners[,1]), y = max(cur_corners[,2]), region = paste("Region ", counter$n[i]))
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
            in.list <- sp::point.in.polygon(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2])
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

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

#' internal Text input with button
#' 
#' Reproduced since it is not exported in the Shiny namespace.
#' 
#' @importFrom htmltools css tags validateCssUnit
#' 
#' @keywords internal
textInputwithButton <- function (textinputId, label, buttoninputId, value = "", width = NULL, placeholder = NULL) 
{
  textvalue <- restoreInput(id = textinputId, default = value)
  buttonvalue <- restoreInput(id = buttoninputId, default = NULL)
  div(class = "form-group shiny-input-container", 
      style =  htmltools::css(width = htmltools::validateCssUnit(width), display = "inline-block"),
      shinyInputLabel(textinputId, label), 
      htmltools::tags$input(id = textinputId, 
                            style = htmltools::css(width = "80%", float = "left"),
                            type = "text", class = "shiny-input-text form-control", 
                            value = textvalue, placeholder = placeholder),
      htmltools::tags$button(id = buttoninputId, 
                             style = htmltools::css(width = "20%", float = "left"),
                             type = "button", class = "btn btn-default action-button", 
                             `data-val` = buttonvalue, disabled = NULL, list(shiny::icon("trash")))
      )
}

#' Shiny's internal \code{shinyInputLabel} function
#' 
#' Reproduced since it is not exported in the Shiny namespace.
#' 
#' @importFrom htmltools css tags validateCssUnit
#' 
#' @keywords internal
shinyInputLabel <- function(inputId, label=NULL) {
  htmltools::tags$label(label,
             class = "control-label",
             class = if (is.null(label)) "shiny-label-null",
             `for` = inputId
  )
}
