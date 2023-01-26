#' register_slides
#'
#' A mini shiny app to calculate spatial cell(barcode) projection matrix between a visium slide and a xenium slide
#'
#' @param Visium_data visium data
#' @param Xenium_data_image_path xenium data
#'
#' @return projection matrix
#'
#' @export
#'
register_slides <- function(Visium_data, Xenium_data_image_path) {

  # shiny
  require(shiny)

  # get images from the slides
  file_image1 <- Visium_data@images$slice1@image
  file_image2 <- Xenium_data_image_path

  # get the ui and server
  if (interactive()){
    ui <- tagList(

      # use javascript extensions for Shiny
      useShinyjs(),

      # get main fluid page
      fluidPage(

        # Main Application in sidebar format
        br(),
        sidebarLayout(position = "left",

                      # Side bar
                      sidebarPanel(

                        # keypoints of image 1
                        h4("Image 1"),
                        fluidRow(
                          tableOutput("xy_Table1"),
                          shiny::actionButton("remove1", "Remove Point")
                        ),

                        # keypoints of image 2
                        h4("Image 2"),
                        fluidRow(
                          tableOutput("xy_Table2"),
                          shiny::actionButton("remove2", "Remove Point")
                        ),

                        # done button for finishing the shiny app and return values
                        br(),
                        fluidRow(
                          shiny::actionButton("done", "Done")
                        ),

                        # panel options
                        width = 3,
                      ),

                      # Main Panel
                      mainPanel(

                        column(6,
                               # The tab set for different images
                               tabsetPanel(id = 'image_tab_panel',

                                           # Image 1 Tab
                                           tabPanel("Image1",
                                                    br(),
                                                    fluidRow(
                                                      column(6,
                                                             selectInput("rotate_image1", "Rotate (ClockWise):",
                                                                         choices = c(0,90,180,270), selected = 0),
                                                      ),
                                                      column(6,
                                                             selectInput("flipflop_image1", "Transform:",
                                                                         choices = c("None", "Flip", "Flop"), selected = "None"),
                                                      )
                                                    ),
                                                    fluidRow(
                                                      imageOutput("plot1", click = "click_plot1"),
                                                    )
                                           )
                               )),
                        column(6,
                               # The tab set for different images
                               tabsetPanel(id = 'image_tab_panel2',

                                           # Image 2 Tab
                                           tabPanel("Image2",
                                                    br(),
                                                    fluidRow(
                                                      column(6,
                                                             selectInput("rotate_image2", "Rotate (ClockWise):",
                                                                         choices = c(0,90,180,270), selected = 0),
                                                      ),
                                                      column(6,
                                                             selectInput("flipflop_image2", "Transform:",
                                                                         choices = c("None", "Flip", "Flop"), selected = "None"),
                                                      )
                                                    ),
                                                    fluidRow(
                                                      imageOutput("plot2", click = "click_plot2"),
                                                    )
                                           ),
                               )),

                        # main Panel options
                        width = 9
                      )
        )
      )
    )

    server <- function(input, output, session) {

      # get a database for collecting chosen points
      xyTable1 <- reactiveVal(tibble(KeyPoint = numeric(), x = numeric(), y = numeric()))
      xyTable2 <- reactiveVal(tibble(KeyPoint = numeric(), x = numeric(), y = numeric()))

      # output the keypoint database
      output$xy_Table1 <- renderTable(xyTable1())
      output$xy_Table2 <- renderTable(xyTable2())

      # click events for updating the keypoints
      observeEvent(input$click_plot1, {
        keypoint <- data.frame(x = input$click_plot1$x, y = input$click_plot1$y)
        keypoint <- transform_keypoints(trans_image1(), keypoint, "image1", input, session)
        xyTable1() %>%
          add_row(KeyPoint = nrow(xyTable1())+1, x = keypoint$x, y = keypoint$y) %>%
          xyTable1()
      })
      observeEvent(input$click_plot2, {
        keypoint <- data.frame(x = input$click_plot2$x, y = input$click_plot2$y)
        keypoint <- transform_keypoints(trans_image2(), keypoint, "image2", input, session)
        xyTable2() %>%
          add_row(KeyPoint = nrow(xyTable2())+1, x = keypoint$x, y = keypoint$y) %>%
          xyTable2()
      })

      # remove points by click
      observeEvent(input$remove1, {
        if(nrow(xyTable1()) > 0)
          xyTable1() %>% filter(KeyPoint != nrow(xyTable1())) %>% xyTable1()
      })
      observeEvent(input$remove2, {
        if(nrow(xyTable2()) > 0)
          xyTable2() %>% filter(KeyPoint != nrow(xyTable2())) %>% xyTable2()
      })

      # get original images
      orig_image1 <- reactiveVal(image_read(file_image1))
      orig_image2 <- reactiveVal(image_read(file_image2))

      # get transformed images
      trans_image1 <- reactive(transform_magick_image(orig_image1(), "image1", input, session))
      trans_image2 <- reactive(transform_magick_image(orig_image2(), "image2", input, session))

      # main observable
      observe({

        # output images
        output$plot1 <- renderPlot({
          img <- transform_magick_image_keypoints(orig_image1(), xyTable1(), "image1", input, session)
          img <- image_ggplot_keypoint(image_ggplot(img$image), img$keypoints)
          return(img)
        })
        output$plot2 <- renderPlot({
          img <- transform_magick_image_keypoints(orig_image2(), xyTable2(), "image2", input, session)
          img <- image_ggplot_keypoint(image_ggplot(img$image), img$keypoints)
          return(img)
        })
      })

      # return values for the shiny app
      observeEvent(input$done, {
        returnValue <- xyTable1() # change this later
        stopApp(returnValue)
      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

#' transform_magick_image
#'
#' Apply given transformations to a magick image
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param input input
#' @param session session
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

#' transform_magick_image_keypoints
#'
#' Apply given transformations to a magick image and keypoints simultaneously for plotting
#'
#' @param image magick image
#' @param extension name extension for the shiny input parameter
#' @param keypoints a set of keypoints
#' @param input input
#' @param session session
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
#' @param input input
#' @param session session
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
