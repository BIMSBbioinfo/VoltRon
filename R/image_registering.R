####
# Main Shiny App ####
####

#' register_slides
#'
#' A mini shiny app to calculate spatial cell(barcode) projection matrix between a visium slide and a xenium slide
#'
#' @param reference_visium visium data
#' @param query_images a list of xenium data as query
#'
#' @return projection matrix
#'
#' @import magick
#'
#' @export
#'
register_slides <- function(reference_visium, query_image_paths, keypoints = NULL) {

  # shiny
  require(shiny)

  # get the reference image
  orig_image_ref <- image_read(reference_visium@images$slice1@image)

  # get the query images
  orig_image_query_list <- lapply(query_image_paths, function(img_path) {
    image_read(img_path)
  })

  # get the ui and server
  if (interactive()){
    # ui <- tagList(
    ui <- fluidPage(
      # use javascript extensions for Shiny
      useShinyjs(),

      sidebarLayout(position = "left",

                    # Side bar
                    sidebarPanel(

                      h4("Image Registration"),
                      fluidRow(
                        column(12,shiny::actionButton("manualregister", "Manuel Registration")),
                        br(),
                        column(12,shiny::actionButton("done", "Done")),
                      ),
                      br(),
                      fluidRow(
                        column(12,
                               htmlOutput("summary")
                               )
                      ),

                      # panel options
                      width = 3,
                    ),

                    mainPanel(

                      # Interface for the reference image
                      br(),
                      column(6,
                             # The tab set for different images
                             sortableTabsetPanel(id = 'image_tab_panel_ref',

                                                 # Image 1 Tab
                                                 tabPanel("Reference Image",
                                                          br(),
                                                          fluidRow(
                                                            column(6,
                                                                   selectInput("rotate_ref", "Rotate (ClockWise):",
                                                                               choices = c(0,90,180,270), selected = 0),
                                                            ),
                                                            column(6,
                                                                   selectInput("flipflop_ref", "Transform:",
                                                                               choices = c("None", "Flip", "Flop"), selected = "None"),
                                                            )
                                                          ),
                                                          fluidRow(
                                                            imageOutput("plot_ref", click = "click_plot_ref"),
                                                          ),
                                                          br(),
                                                          fluidRow(
                                                            column(6,
                                                                   tableOutput("xy_Table_ref"),
                                                                   shiny::actionButton("remove_ref", "Remove Point")
                                                            )
                                                          )
                                                 )
                             )
                      ),

                      # Interface for the query images
                      column(6,QueryTabPanels(length(orig_image_query_list))),

                      width = 9
                    )
      )

      # # get main fluid page
      # fluidPage(
      # ....
      # )
    )

    server <- function(input, output, session) {

      ## Set up reference keypoints and output ####
      if(is.null(keypoints)){
        xyTable_ref <- reactiveVal(tibble(KeyPoint = numeric(), x = numeric(), y = numeric()))
      } else {
        xyTable_ref <- reactiveVal(keypoints[[1]])
      }
      output$xy_Table_ref <- renderTable(xyTable_ref()) # output the keypoint table
      observeEvent(input$click_plot_ref, { # click events for updating the keypoints
        keypoint <- data.frame(x = input$click_plot_ref$x, y = input$click_plot_ref$y)
        keypoint <- transform_keypoints(trans_ref(), keypoint, "ref", input, session)
        xyTable_ref() %>%
          add_row(KeyPoint = nrow(xyTable_ref())+1, x = keypoint$x, y = keypoint$y) %>%
          xyTable_ref()
      })
      observeEvent(input$remove_ref, { # remove points by click
        if(nrow(xyTable_ref()) > 0)
          xyTable_ref() %>% filter(KeyPoint != nrow(xyTable_ref())) %>% xyTable_ref()
      })

      ## Set up query keypoints and output ####
      xyTable_query_list <- QueryKeypoints(length(orig_image_query_list), keypoints[2:length(keypoints)])
      QueryKeypointTable(xyTable_query_list, trans_query_list, input, output, session)

      ## Get transformed reference and query images ####
      trans_ref <- reactive(transform_magick_image(orig_image_ref, "ref", input, session))
      trans_query_list <- transform_magick_image_query_list(orig_image_query_list, input, session)

      ## Return Registered keypoints ####
      ManualRegisterMatrix_list <- QueryMatrices(length(orig_image_query_list))
      getManualRegisterMatrix(ManualRegisterMatrix_list, length(orig_image_query_list),
                              xyTable_query_list, xyTable_ref(), input, output, session)

      ## Main observable ####
      observe({

        # output the reference image
        output$plot_ref <- renderPlot({
          img <- transform_magick_image_keypoints(orig_image_ref, xyTable_ref(), "ref", input, session)
          img <- image_ggplot_keypoint(image_ggplot(img$image), img$keypoints)
          return(img)
        })

        # output the list of query images
        QueryImageOutput(orig_image_query_list, xyTable_query_list, input, output, session)
      })

      ## Return values for the shiny app ####
      observeEvent(input$done, {
        keypoints <- list(xyTable_ref())
        keypoints <- c(keypoints, reactiveValuesToList(xyTable_query_list))
        TransMatrix <- reactiveValuesToList(ManualRegisterMatrix_list)
        stopApp(list(Keypoints = keypoints, TransMatrix = TransMatrix))
      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

####
# Managing Images ####
####

#' QueryTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#'
#' @return tabsetpanel
#'
#' @export
#'
QueryTabPanels <- function(len_images){
  do.call(sortableTabsetPanel, c(id='image_tab_panel_query',lapply(1:len_images, function(i) {
    tabPanel(paste0("Query ",i),
             br(),
             fluidRow(
               column(6, selectInput(paste0("rotate_image",i), "Rotate (ClockWise):", choices = c(0,90,180,270), selected = 0)),
               column(6, selectInput(paste0("flipflop_image",i), "Transform:", choices = c("None", "Flip", "Flop"), selected = "None"))
             ),
             fluidRow(imageOutput(paste0("plot_query",i), click = paste0("click_plot",i))),
             fluidRow(
               tableOutput(paste0("xy_Table_query",i)),
               shiny::actionButton(paste0("remove_query",i), "Remove Point")
             )
    )
  })))
}

#' QueryImageOutput
#'
#' Shiny outputs for a set of magick images with keypoints
#'
#' @param image_list a list of magick images
#' @param keypoints_list a list of data frames, each having a set of keypoints
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @export
#'
QueryImageOutput <- function(image_list, keypoints_list, input, output, session){

  # get the length of images
  len_images <- length(image_list)

  # output query images
  lapply(1:len_images, function(i){
    output[[paste0("plot_query",i)]] <- renderPlot({
      keypoints <- keypoints_list[[paste0(i)]]
      img <- transform_magick_image_keypoints(image_list[[i]], keypoints, paste0("image",i), input, session)
      img <- image_ggplot_keypoint(image_ggplot(img$image), img$keypoints)
      return(img)
    })
  })

}

#' transform_magick_image
#'
#' Apply given transformations to a magick image
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#' @param session shiny session
#'
#' @return magick image
#'
#' @export
#'
transform_magick_image <- function(image, extension, input, session){

  # rotate image and keypoints
  input_rotate <- as.numeric(input[[paste0("rotate_", extension)]])
  image <- image_rotate(image, input_rotate)

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    image <- image_flip(image)
  } else if(input_flipflop == "Flop"){
    image <- image_flop(image)
  }

  image
}

#' transform_magick_image
#'
#' Apply given transformations to a magick image
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#' @param session shiny session
#'
#' @return magick image
#'
#' @export
#'
transform_magick_image_query_list <- function(image_list, input, session){

  trans_query_list <- lapply(1:length(image_list), function(i){
    reactive({
      transform_magick_image(image_list[[i]], paste0("image",i), input, session)
    })
  })
  return(trans_query_list)
}

####
# Managing Keypoints ####
####

#' QueryKeypoints
#'
#' Initiate shiny reactive values for keypoint dataframes
#'
#' @param len_images the number of query images
#' @param keypoints_list the keypoints list of query images
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @return a shiny reactive values object
#'
#' @export
#'
QueryKeypoints <- function(len_images, keypoints_list, input, output, session){

  # initiate keypoints
  if(is.null(keypoints_list)){
    keypoints_list <- lapply(1:len_images, function(i) {
      tibble(KeyPoint = numeric(), x = numeric(), y = numeric())
    })
  }

  # set names for keypoints
  names(keypoints_list) <- 1:len_images

  # return keypoints as reactive values
  do.call("reactiveValues", keypoints_list)
}

#' QueryKeypointTable
#'
#' A list of shiny observe events for keypoints tables, updating keypoints and auxiliry operations
#'
#' @param xyTable_list a list of keypoints x,y coordinates for each magick image
#' @param image_list a lost of magick image
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @export
#'
QueryKeypointTable <- function(xyTable_list, image_list, input, output, session){

  # get the length of tables
  len_tables <- length(xyTable_list)

  # output query keypoint tables
  lapply(1:len_tables, function(i){
    output[[paste0("xy_Table_query",i)]] <- renderTable({xyTable_list[[paste0(i)]]})
  })

  # set click operations for the set of keypoints
  lapply(1:len_tables, function(i){
    observeEvent(input[[paste0("click_plot",i)]], {
      keypoint <- data.frame(x = input[[paste0("click_plot",i)]]$x, y = input[[paste0("click_plot",i)]]$y)
      if(is.reactive(image_list[[i]])){
        image <- image_list[[i]]()
      } else {
        image <- image_list[[i]]
      }
      keypoint <- transform_keypoints(image, keypoint, paste0("image",i), input, session)
      temp <- xyTable_list[[paste0(i)]]
      temp <- temp %>%
        add_row(KeyPoint = nrow(temp)+1, x = keypoint$x, y = keypoint$y)
      xyTable_list[[paste0(i)]] <- temp
    })
  })

  # remove keypoints from query images
  lapply(1:len_tables, function(i){
    observeEvent(input[[paste0("remove_query",i)]], {
      temp <- xyTable_list[[paste0(i)]]
      if(nrow(temp) > 0){
        temp <- temp %>% filter(KeyPoint != nrow(temp))
        xyTable_list[[paste0(i)]] <- temp
      }
    })
  })
}

#' QueryMatrices
#'
#' Initiate shiny reactive values for registeration matrices
#'
#' @param len_images the number of query images
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @return a shiny reactive values object
#'
#' @export
#'
QueryMatrices <- function(len_images, input, output, session){

  # initiate matrices
  matrix_list <- lapply(1:len_images, function(i) return(NULL))
  # matrix_list <- list()
  # length(matrix_list) <- len_images
  names(matrix_list) <- 1:len_images

  # return matrices as reactive values
  do.call("reactiveValues", matrix_list)
}

#' ManualRegisterMatrix
#'
#' Manuel registeration of images using manually entered keypoints
#'
#' @param ManualRegisterMatrix_list a list of registration matrices of each query image
#' @param len_images length of query images
#' @param xyTable_list a list of keypoints x,y coordinates for query image
#' @param xyTable_ref the keypoints x,y coordinates of the reference image
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @export Morpho
#'
getManualRegisterMatrix <- function(ManualRegisterMatrix_list, len_images, xyTable_list, xyTable_ref, input, output, session){

  observeEvent(input$manualregister, {
    reference_landmark <- as.matrix(xyTable_ref[,c("x","y")])
    for(i in 1:len_images){
      target_landmark <- as.matrix(xyTable_list[[paste0(i)]][,c("x","y")])
      ManualRegisterMatrix_list[[paste0(i)]] <-
        computeTransform(reference_landmark, target_landmark, type = "tps")
    }
    output[["summary"]] <- renderUI({
      str1 <- paste0("Summary:")
      str2 <- paste0("# of Query Images: ", len_images)
      str3 <- paste0("# of Keypoints: ", nrow(reference_landmark))
      all_str <- c(str1, str2, str3)
      # for(i in 1:len_images){
      #   target_landmark <- as.matrix(xyTable_list[[paste0(i)]][,c("x","y")])
      #   all_str <- c(all_str,
      #                paste0("# of Query ", i, " Keypoints: ", nrow(target_landmark)))
      # }
      HTML(paste(all_str, collapse = '<br/>'))
    })
  })
}

#' transform_magick_image_keypoints
#'
#' Apply given transformations to a magick image and keypoints simultaneously for plotting
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param keypoints a set of keypoints
#' @param input shiny input
#' @param session shiny session
#'
#' @return a list of magick image and keypoints
#'
#' @export
#'
transform_magick_image_keypoints <- function(image, keypoints, extension, input, session){

  # get unrotated image info
  image_limits <- unlist(image_info(image)[1,c("width", "height")])
  image_origin <- image_limits/2

  # rotate image and keypoints
  input_rotate <- as.numeric(input[[paste0("rotate_", extension)]])
  image <- image_rotate(image, input_rotate)

  # get rotated image info
  rotated_image_limits <- unlist(image_info(image)[1,c("width", "height")])
  rotated_image_origin <- rotated_image_limits/2

  # rotate keypoints
  keypoints <- rotate_keypoint(keypoints, input_rotate, image_origin, image_limits, rotated_image_origin, rotated_image_limits)

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    image <- image_flip(image)
  } else if(input_flipflop == "Flop"){
    image <- image_flop(image)
  }

  # flipflop keypoints
  keypoints <- flipflop_keypoint(keypoints, rotated_image_limits, input_flipflop)

  # return both the image and the keypoints
  return(list(image = image, keypoints = keypoints))
}

#' transform_keypoints
#'
#' Apply transformations to keypoints given transformed images to find the keypoints locations in the original image
#'
#' @param image magick image
#' @param keypoints keypoints visualized on image
#' @param extension name extension for the shiny input parameter
#' @param input shiny input
#' @param session shiny session
#'
#' @return magick image
#'
#' @export
#'
transform_keypoints <- function(image, keypoints, extension, input, session){

  # get unrotated image info
  image_limits <- unlist(image_info(image)[1,c("width", "height")])
  image_origin <- image_limits/2

  # flip flop image and keypoints
  input_flipflop <- input[[paste0("flipflop_", extension)]]
  if(input_flipflop == "Flip"){
    image <- image_flip(image)
  } else if(input_flipflop == "Flop"){
    image <- image_flop(image)
  }
  keypoints <- flipflop_keypoint(keypoints, image_limits, input_flipflop)

  # rotate image (reverse) and keypoints
  input_rotate <- 360 - as.numeric(input[[paste0("rotate_", extension)]])
  image <- image_rotate(image, input_rotate)

  # get rotated image info
  rotated_image_limits <- unlist(image_info(image)[1,c("width", "height")])
  rotated_image_origin <- rotated_image_limits/2

  # rotate keypoints
  keypoints <- rotate_keypoint(keypoints, input_rotate, image_origin, image_limits, rotated_image_origin, rotated_image_limits)

  return(keypoints)
}

#' rotate_keypoint
#'
#' Find transformations of keypoints under clockwise rotations of the image
#'
#' @param keypoints dataset of keypoints
#' @param angle angle of rotation [0,360]
#' @param origin center of the image
#' @param limits limits of the image
#' @param rotated_origin center of the rotated image
#' @param rotated_limits limits of the rotated image
#'
#' @return keypoints
#'
#' @export
#'
rotate_keypoint <- function(keypoints, angle, origin, limits, rotated_origin, rotated_limits){

  # if there are no points, return
  if(nrow(keypoints) == 0)
    return(keypoints)

  # get coordinates from the keypoints dataset
  points <- keypoints[,c("x","y")]

  # set rotation matrix for angles
  radii <- ((360-angle)*pi/180)
  s = sin(radii);
  c = cos(radii);
  rotation_mat <- matrix(c(c, s, -s, c), nrow = 2, byrow = F)

  # rotate point
  points <- t(apply(points, 1, function(x) return(x - origin)))
  points <- t(apply(points, 1, function(x) return(x/limits)))
  rotated_points <- t(rotation_mat %*% t(points))
  rotated_points <- t(apply(rotated_points, 1, function(x) return(x*rotated_limits)))
  rotated_points <- t(apply(rotated_points, 1, function(x) return(x + rotated_origin)))

  # put rotated points back to keypoints
  keypoints[,c("x","y")] <- rotated_points

  return(keypoints)
}

#' flipflop_keypoint
#'
#' Find transformed keypoints on image given any flip or flop action by magick
#'
#' @param keypoints dataset of keypoints
#' @param image_limits limits of the images
#' @param flipflop a flip or flop action as string
#'
#' @param keypoints
#'
#' @export
#'
flipflop_keypoint <- function(keypoints, image_limits, flipflop){

  if(nrow(keypoints) == 0)
    return(keypoints)

  if(grepl("Flop", flipflop))
    keypoints$x = image_limits[1] - keypoints$x

  if(grepl("Flip", flipflop))
    keypoints$y = image_limits[2] - keypoints$y

  return(keypoints)
}

#' image_ggplot_keypoint
#'
#' add keypoints as points on ggplot object
#'
#' @param image magick image
#' @param keypoints keypoints to draw on image
#'
#' @return ggplot object
#'
#' @export
#'
image_ggplot_keypoint <- function(image, keypoints){

  # select keypoints and texts on image
  image <- image +
    geom_point(mapping = aes(x = x, y = y), keypoints, size = 10, shape = 21, fill = "white") +
    geom_text(mapping = aes(x = x, y = y, label = KeyPoint), keypoints, size = 7)
}
