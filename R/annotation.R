####
# Main Shiny App ####
####

#' annotateSpatialData
#'
#' A mini shiny app to for annotating spatial points
#'
#' @param object a list of VoltRon (or Seurat) objects
#' @param assay a reference spatial data set, used only if \code{object_list} is \code{NULL}
#' @param ... additional parameters passed to \code{vrSpatialPlot}
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs show hide
#' @importFrom stats median
#' @importFrom sp point.in.polygon
#' @import ggplot2
#'
#' @export
#'
annotateSpatialData <- function(object, assay = NULL, ...) {

  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  ## Importing images ####

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  if(length(assay_names) > 0)
    assay <- assay_names[1]

  # get image
  g <- vrSpatialPlot(object, assay = assay, ...) + labs(title = "")

  ## UI and Server ####

  # get the ui and server
  if (interactive()){
    ui <- fluidPage(
      sidebarLayout(position = "left",

        sidebarPanel(
          tags$style(make_css(list('.well', 'margin', '7%'))),

          # Interface
          fluidRow(
            column(12,h4("Annotation Interface")),
            br()
          ),

          # points
          fluidRow(
            column(12,shiny::actionButton("reset_btn", "Reset Points")),
            column(12,shiny::actionButton("rmvlast_btn", "Remove Last Point")),
            column(12,shiny::actionButton("addregion_btn", "Add Region")),
            br()
          ),

          # Subsets
          fluidRow(
            column(12,h4("Selected Regions")),
            uiOutput("textbox_ui"),
            br()
          ),

          # Subsets
          fluidRow(
            column(12,h4("Finished Selecting?")),
            column(12,shiny::actionButton("done", "Yes!"))
          ),
          width = 4
        ),
        mainPanel(
          shinyjs::useShinyjs(),
          plotOutput("image_plot", click = "plot_click", height = "1000px"),
          width = 8
        )
      )
    )

    server <- function(input, output) {

      # Initialize data frame to store points
      selected_corners <- reactiveVal(data.frame(x = numeric(0), y = numeric(0)))
      selected_corners_list <- reactiveVal(list())
      selected_corners_list_label <- reactiveVal(list())

      # update summary
      output[["summary"]] <- renderUI({
        if(length(selected_corners_list_label()) > 0){
          htmltools::HTML(paste(unlist(selected_corners_list_label()), collapse = '<br/>'))
        }
      })

      # point click event
      observeEvent(input$plot_click, {
        click <- input$plot_click
        x <- click$x
        y <- click$y

        # Append new point to the data frame
        new_point <- data.frame(x = x, y = y)
        selected_corners(rbind(selected_corners(), new_point))
      })

      # reset and remove buttons
      observeEvent(input$reset_btn, {
        selected_corners(data.frame(x = numeric(0), y = numeric(0)))
      })
      observeEvent(input$rmvlast_btn, {
        selected_corners(selected_corners()[-nrow(selected_corners()),])
      })

      # add region
      observeEvent(input$addregion_btn, {
        if(nrow(selected_corners()) > 3){

          # add to region list
          selected_corners_list(c(selected_corners_list(), list(selected_corners())))
          print(selected_corners_list())

          # add to region label
          if(length(selected_corners_list_label) == 0){
            new_label <- "Region 1"
          } else {
            new_label <- paste0("Region ", length(selected_corners_list_label()) + 1)
          }
          selected_corners_list_label(c(selected_corners_list_label(), list(new_label)))
          print(selected_corners_list_label())

          # remove selected points
          selected_corners(data.frame(x = numeric(0), y = numeric(0)))

          # Track the number of input boxes to render
          counter$n <- counter$n + 1
        }
      })

      # counter for regions
      counter <- reactiveValues(n = 0)
      output$textbox_ui <- renderUI({ textboxes() })
      textboxes <- reactive({
        n <- counter$n
        if (n > 0) {
          lapply(seq_len(n), function(i) {
            fluidRow(
              column(12,textInput(inputId = paste0("region", i),
                                  label = paste0("Region ", i), value = paste0("Region ", i)))
            )
          })
        }
      })

      output$image_plot <- renderPlot({
        g <- g +
          ggplot2::geom_point(aes(x = x, y = y), data = selected_corners(), color = "red", shape = 16) +
          ggplot2::geom_polygon(aes(x = x, y = y, group = "region"), data = selected_corners(), alpha = 0.4, color = "red")

        datax <- selected_corners()
        datax_label_ind <- length(selected_corners_list()) + 1
        g <- g + ggrepel::geom_label_repel(mapping = aes(x = mean(datax[,1]), y = max(datax[,2]), label = paste("Region ", datax_label_ind)),
                                           size = 5, direction = "y", nudge_y = 6, box.padding = 0, label.padding = 1, seed = 1, color = "red")

        g
      })

      ## Return values for the shiny app ####
      observeEvent(input$done, {

      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}
