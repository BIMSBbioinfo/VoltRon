####
# Main Shiny App ####
####

#' register_slides
#'
#' A mini shiny app to calculate spatial cell(barcode) projection matrix between a visium slide and a xenium slide
#'
#' @param reference_seu Seurat object with reference spatial assay
#' @param query_list_seu a list of Seurat Object with query images
#'
#' @return keypoints
#'
#' @import magick
#'
#' @export
#'
register_slides <- function(reference_seu, query_list_seu, keypoints = NULL) {

  # shiny
  require(shiny)

  images <- getRegisteringImages(reference_seu, query_list_seu)
  orig_image_ref <- images$ref_image
  orig_image_query_list <- images$query_image_list

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
                                                            column(4,
                                                                   selectInput("rotate_ref", "Rotate (ClockWise):",
                                                                               choices = c(0,90,180,270), selected = 0),
                                                            ),
                                                            column(4,
                                                                   selectInput("flipflop_ref", "Transform:",
                                                                               choices = c("None", "Flip", "Flop"), selected = "None"),
                                                            ),
                                                            column(4,
                                                                   selectInput("negate_ref", "Negate Image:",
                                                                               choices = c("No", "Yes"), selected = "No"),
                                                            ),
                                                          ),
                                                          fluidRow(
                                                            imageOutput("plot_ref", click = "click_plot_ref"),
                                                          ),
                                                          br(),
                                                          fluidRow(
                                                            column(6,
                                                                   # tableOutput("xy_Table_ref"),
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
                              orig_image_query_list, orig_image_ref,
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

#' getRegisteringImages
#'
#' get the images from Spatial assays that will be used for registration
#'
#' @param reference Seurat object with reference spatial assay
#' @param query_list the list of Seurat object with query images for registration
#'
getRegisteringImages <- function(reference, query_list){

  # Seurat list of both reference and query list
  seu_list <- c(reference, query_list)

  # get images from Seurat list
  image_list <- lapply(seu_list, getSeuratImage)
  ref_image <- image_list[[1]]
  query_image_list <- image_list[-1]

  return(list(ref_image = ref_image, query_image_list = query_image_list))
}

#' getSeuratImage
#'
#' get the image from a Spatial assay
#'
#' @param seu Seurat object
#'
#' @import magick
#'
getSeuratImage <- function(seu){

  image_classes <- sapply(seu@images, class)

  if(any(grepl("FOV",image_classes))){
    image <- seu@images[[names(seu@images)[which(grepl("FOVImage", image_classes))]]]
    image <- image@image
  } else if(any(grepl("Visium",image_classes))) {
    image <- seu@images[[names(seu@images)[which(grepl("Visium", image_classes))]]]
    image <- magick::image_read(image@image)
  }

  return(image)
}


#' QueryTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#'
#' @return tabsetpanel
#'
QueryTabPanels <- function(len_images){
  do.call(sortableTabsetPanel, c(id='image_tab_panel_query',lapply(1:len_images, function(i) {
    tabPanel(paste0("Query ",i),
             br(),
             fluidRow(
               column(4, selectInput(paste0("rotate_image",i), "Rotate (ClockWise):", choices = c(0,90,180,270), selected = 0)),
               column(4, selectInput(paste0("flipflop_image",i), "Transform:", choices = c("None", "Flip", "Flop"), selected = "None")),
               column(4, selectInput(paste0("negate_image",i), "Negate Image:", choices = c("No", "Yes"), selected = "No"))
             ),
             fluidRow(imageOutput(paste0("plot_query",i), click = paste0("click_plot",i))),
             fluidRow(
               # tableOutput(paste0("xy_Table_query",i)),
               shiny::actionButton(paste0("remove_query",i), "Remove Point")
             ),
             br(),
             h5("Registered Image:"),
             fluidRow(imageOutput(paste0("plot_query_reg",i)))
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
transform_magick_image_keypoints <- function(image, keypoints, extension, input, session){

  # negate image
  input_negate <- input[[paste0("negate_", extension)]]
  if(input_negate == "Yes"){
    image <- image_negate(image)
  }

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
image_ggplot_keypoint <- function(image, keypoints){

  # select keypoints and texts on image
  image <- image +
    geom_point(mapping = aes(x = x, y = y), keypoints, size = 5, shape = 21, fill = "white") +
    geom_text(mapping = aes(x = x, y = y, label = KeyPoint), keypoints, size = 3)
}

####
# Image Registration ####
####

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
#' @param image_list the list of query images
#' @param xyTable_list a list of keypoints x,y coordinates for query image
#' @param xyTable_ref the keypoints x,y coordinates of the reference image
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @import Morpho
#'
getManualRegisterMatrix <- function(ManualRegisterMatrix_list, len_images, image_list, ref_image, xyTable_list, xyTable_ref, input, output, session){

  observeEvent(input$manualregister, {

    # Register keypoints
    reference_landmark <- as.matrix(xyTable_ref[,c("x","y")])
    for(i in 1:len_images){
      target_landmark <- as.matrix(xyTable_list[[paste0(i)]][,c("x","y")])
      ManualRegisterMatrix_list[[paste0(i)]] <-
        computeTransform(reference_landmark, target_landmark, type = "tps")
    }

    # Output summary
    output[["summary"]] <- renderUI({
      str1 <- paste0(" Registration Summary:")
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

    # Registered Images
    lapply(1:len_images, function(i){
      output[[paste0("plot_query_reg",i)]] <- renderPlot({
        getRegisteredImage(image_list[[i]], ref_image, ManualRegisterMatrix_list[[paste0(i)]])
      })
    })
  })
}

getRegisteredImage <- function(query_image, ref_image, transmatrix){

  # plot with raster
  ref_image_raster <- as.raster(ref_image) |> as.matrix() |> rast()
  query_image_raster <- as.raster(query_image) |> as.matrix() |> rast() |> stack()

  # apply transformation to image
  imageEx <- raster::extent(stack(ref_image_raster))
  imageRes <- raster::res(stack(ref_image_raster))
  query_image_raster_1 <- raster::as.data.frame(query_image_raster[[1]], xy = TRUE)
  query_image_raster_1_t <- Morpho::applyTransform(as.matrix(query_image_raster_1)[,1:2], transmatrix)
  r <- raster::raster(nrow = dim(query_image_raster)[1], ncol = dim(query_image_raster)[2], resolution = c(1,1))
  raster::extent(r) <- imageEx
  raster::res(r) <- imageRes
  query_image_raster_1_tr <- raster::rasterize(query_image_raster_1_t, field = query_image_raster_1[,3], r, fun = mean)
  query_image_raster_1_trf <- focal(query_image_raster_1_tr,
                                   w = matrix(1, nrow = 3, ncol = 3),
                                   fun = fill.na, pad = TRUE, na.rm = FALSE)
  query_image_raster_1_trf <- terra::rast(query_image_raster_1_trf, crs = "")

  # plot
  p <- recordPlot
  terra::plot(ref_image_raster)
  raster::plot(query_image_raster_1_trf, alpha = 0.2, add = TRUE, legend = FALSE)
  p
}

getRegisteredImage_old <- function(query_image, ref_image, transmatrix){

  # query image
  ref_raster <- as.raster(ref_image) |> as.matrix() |> rast() |> stack()
  ref_raster_1 <- raster::as.data.frame(ref_raster[[1]], xy = TRUE)
  # ref_raster_2 <- raster::as.data.frame(ref_raster[[2]], xy = TRUE)
  # ref_raster_3 <- raster::as.data.frame(ref_raster[[3]], xy = TRUE)
  imageEx <- raster::extent(ref_raster)
  r <- raster::raster(nrow = dim(ref_raster)[1], ncol = dim(ref_raster)[2])
  raster::extent(r) <- imageEx
  ref_raster_1r <- raster::rasterize(ref_raster_1[,1:2], field = ref_raster_1[,3], r, fun = mean)
  # ref_raster_2r <- raster::rasterize(ref_raster_2[,1:2], field = ref_raster_2[,3], r, fun = mean)
  # ref_raster_3r <- raster::rasterize(ref_raster_3[,1:2], field = ref_raster_3[,3], r, fun = mean)
  # ref_raster <- raster::stack(ref_raster_1r, ref_raster_2r, ref_raster_3r)
  ref_raster <- raster::stack(ref_raster_1r)
  ref_raster_data <- raster::as.data.frame(ref_raster, xy = TRUE)

  # make raster from magick
  img <- query_image
  query_raster <- as.raster(img) |> as.matrix() |> rast() |> stack()

  # transform
  query_raster_data <- raster::as.data.frame(query_raster[[1]], xy = TRUE)
  query_raster_data_reg <- Morpho::applyTransform(as.matrix(query_raster_data)[,1:2], transmatrix)
  r <- raster::raster(nrow = dim(query_raster)[1], ncol = dim(query_raster)[2])
  raster::extent(r) <- imageEx
  query_raster_data_regr <- raster::rasterize(query_raster_data_reg, field = query_raster_data[,3], r, fun = mean)
  # query_raster_data_regr <- raster::stack(query_raster_data_regr)

  # # visualize
  # query_raster_data_regr_data <- raster::as.data.frame(query_raster_data_regr, xy = TRUE)
  # ggplot() +
  #   geom_raster(aes(x=x,y=y, fill = layer, alpha = 0.5), data = ref_raster_data) +
  #   geom_raster(aes(x=x,y=y, fill = layer), data = query_raster_data_regr_data) +
  #   coord_fixed()
  #   # theme_bw() +
  #   # theme(axis.line=element_blank(),axis.text.x=element_blank(),
  #   #       axis.text.y=element_blank(),axis.ticks=element_blank(),
  #   #       axis.title.x=element_blank(),
  #   #       axis.title.y=element_blank(),
  #   #       legend.position="none",
  #   #       panel.background=element_blank(),
  #   #       panel.border=element_blank(),
  #   #       panel.grid.major=element_blank(),
  #   #       panel.grid.minor=element_blank(),
  #   #       plot.background=element_blank())
  ref_raster <- as.raster(ref_image) |> as.matrix() |> rast()
  query_raster_data_regr <- rast(query_raster_data_regr)
  ggplot() +
    geom_spatraster_rgb(data = ref_raster) +
    geom_spatraster(data = query_raster_data_regr, alpha = 0.5) +
    coord_fixed()
}
