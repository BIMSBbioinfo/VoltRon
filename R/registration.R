####
# Main Shiny App ####
####

#' register_slides
#'
#' A mini shiny app to calculate spatial cell(barcode) projection matrix between a visium slide and a xenium slide
#'
#' @param data_list a list of VoltRon objects
#' @param reference_spatdata a reference spatial data set, used only if \code{spatial_data_list} is \code{NULL}
#' @param query_spatdata a query spatial data set, used only if \code{spatial_data_list} is \code{NULL}
#'
#' @return keypoints
#'
#' @import magick
#'
#' @export
#'
SpatialRegistration <- function(data_list = NULL, reference_spatdata = NULL, query_spatdata = NULL, keypoints = NULL) {

  # shiny
  require(shiny)

  ## Importing images ####

  # check object classes
  # if the input is not a list, switch to reference vs query mode
  if(!is.null(data_list)){

    # check if all elements of the list are VoltRon
    if(!all(sapply(data_list, class)=="VoltRon")){
      stop("Please make sure that all objects in the list are of VoltRon class")
    } else {
      spatdata_list <- data_list
    }

    # reference and query indices
    centre <- floor(median(1:length(spatdata_list)))
    register_ind <- setdiff(1:length(spatdata_list), centre)

  # reference vs query mode
  } else {

    # check if the reference and query data sets are either Seurat or Giotto
    spat_classes <- c(class(reference_spatdata), class(query_spatdata))
    if(!all(spat_classes %in% c("Seurat", "Giotto"))) {
      stop("Please make sure that reference and query data sets are of either Seurat or Giotto class")
    } else {
      spatdata_list <- list(reference_spatdata, query_spatdata)
    }

    # reference and query indices
    centre <- 1
    register_ind <- 2
  }

  # get images from the list of objects
  orig_image_query_list <- unlist(lapply(spatdata_list, Image))


  ## UI and Server ####

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
                        # column(12,shiny::actionButton("manualregister", "Manual Registration", width = "80%")),
                        column(12,shiny::actionButton("manualregister", "Manual Registration")),
                        br(),
                        br(),
                        # column(12,shiny::actionButton("automaticregister", "Auto Registration", width = "80%")),
                        column(12,shiny::actionButton("automaticregister", "Auto Registration")),
                        br(),
                        br(),
                        column(12,textInput("GOOD_MATCH_PERCENT", "Match %", value = "0.20", width = "80%", placeholder = NULL)),
                        column(12,textInput("MAX_FEATURES", "# of Features", value = "1000", width = "80%", placeholder = NULL)),
                        br(),
                        column(12,shiny::actionButton("done", "Done")),
                      ),
                      br(),
                      fluidRow(
                        column(12,htmlOutput("summary"))
                      ),

                      # panel options
                      width = 3,
                    ),

                    mainPanel(

                      # Interface for the reference image
                      br(),
                      column(6,

                             # Reference Images
                             ImageTabPanels(length(orig_image_query_list), type = "ref"),

                             br(),

                             # Matching Alignment
                             AlignmentTabPanel(length(orig_image_query_list), centre, register_ind),
                      ),

                      # Interface for the query images
                      column(6,

                             # Query Images
                             ImageTabPanels(length(orig_image_query_list), type = "query"),

                             br(),

                             # Registered Query Images
                             RegisteredImageTabPanels(length(orig_image_query_list), centre, register_ind)
                      ),

                      # panel options
                      width = 9
                    )
      )
    )

    server <- function(input, output, session) {

      ### Manage interface ####
      UpdateSequentialTabPanels(input, output, session, centre, register_ind)

      ### Transform images ####
      trans_image_query_list <- transform_magick_image_query_list(orig_image_query_list, input, session)

      ### Manage reference and query keypoints ####
      xyTable_list <- initateKeypoints(length(orig_image_query_list), keypoints)
      manageKeypoints(centre, register_ind, xyTable_list, trans_image_query_list, input, output, session)

      ### Image registration ####
      registered_spatdata_list <- QueryMatrices(length(spatdata_list))
      getManualRegisteration(registered_spatdata_list, spatdata_list, orig_image_query_list, xyTable_list,
                             centre, register_ind, input, output, session)
      getAutomatedRegisteration(registered_spatdata_list, spatdata_list, orig_image_query_list,
                                centre, register_ind, input, output, session)

      ### Main observable ####
      observe({

        # output the list of query images
        ImageOutput(orig_image_query_list, xyTable_list, centre, input, output, session)

      })

      ## Return values for the shiny app ####
      observeEvent(input$done, {

        # keypoints
        keypoints <- reactiveValuesToList(xyTable_list)

        # Transformation (Mapping) matrix
        RegisteredSpatialData <- reactiveValuesToList(registered_spatdata_list)

        # stop app and return
        stopApp(list(keypoints = keypoints, registered_spat = RegisteredSpatialData))
      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

####
# User Interface ####
####

#' ImageTabPanels
#'
#' The UI for a set of reference/query spatial slides
#'
#' @param len_images the number of query images
#'
#' @return tabsetpanel
#'
ImageTabPanels <- function(len_images, type){

  # get panel label
  label <- ifelse(type == "ref", "Ref. ", "Query ")

  # call panels
  do.call(tabsetPanel, c(id=paste0('image_tab_panel_',type),lapply(1:len_images, function(i) {
    tabPanel(paste0(label,i),
             br(),
             fluidRow(
               column(4, selectInput(paste0("rotate_", type, "_image",i), "Rotate (ClockWise):", choices = c(0,90,180,270), selected = 0)),
               column(4, selectInput(paste0("flipflop_", type, "_image",i), "Transform:", choices = c("None", "Flip", "Flop"), selected = "None")),
               column(4, selectInput(paste0("negate_", type, "_image",i), "Negate Image:", choices = c("No", "Yes"), selected = "No"))
             ),
             fluidRow(
               shiny::actionButton(paste0("remove_", type, i), "Remove Point")
             ),
             fluidRow(imageOutput(paste0("plot_", type, i), click = paste0("click_plot_", type, i))),
    )
  })))
}

#' RegisteredImageTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#' @param centre center image index
#' @param register_ind query image indices
#'
#' @return tabsetpanel
#'
AlignmentTabPanel <- function(len_images, centre, register_ind){

  # tab panels
  do.call(tabsetPanel, c(id='image_tab_panel_alignment',lapply(register_ind, function(i) {
    tabPanel(paste0("Ali. ",i, "->", centre),
             br(),
             fluidRow(imageOutput(paste0("plot_alignment",i)))
    )
  })))
}

#' RegisteredImageTabPanels
#'
#' The UI for a set of query spatial slides
#'
#' @param len_images the number of query images
#' @param centre center image index
#' @param register_ind query image indices
#'
#' @return tabsetpanel
#'
RegisteredImageTabPanels <- function(len_images, centre, register_ind){

  # tab panels
  do.call(tabsetPanel, c(id='image_tab_panel_reg_query',lapply(register_ind, function(i) {
    tabPanel(paste0("Reg. ",i, "->", centre),
             br(),
             column(6, sliderInput(paste0("plot_query_reg_alpha",i), label = "Alpha Level", min = 0, max = 1, value = 0.2)),
             fluidRow(imageOutput(paste0("plot_query_reg",i)))
    )
  })))
}

#' UpdateSequentialTabPanels
#'
#' A function for automatized selection of reference/query images
#'
#' @param input input
#' @param output output
#' @param session session
#' @param centre center image index
#' @param register_ind query image indices
#'
UpdateSequentialTabPanels <- function(input, output, session, centre, register_ind){

  # number of panels
  npanels <- length(register_ind) + 1

  # observe changes in the reference tab panel
  observeEvent(input$image_tab_panel_ref,{
    selected_panel <- input$image_tab_panel_ref
    selected_panel_ind <- as.numeric(strsplit(selected_panel, split = " ")[[1]][2])

    query_panel_ind <- (selected_panel_ind + 1)
    if(query_panel_ind == 1) query_panel_ind <- npanels
    updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", query_panel_ind))
    updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",selected_panel_ind, "->", centre))

    if(selected_panel_ind == npanels)
      updateTabsetPanel(session, "image_tab_panel_ref", paste0("Ref. ", selected_panel_ind-1))
  })

  # observe changes in the query tab panel
  observeEvent(input$image_tab_panel_query,{
    selected_panel <- input$image_tab_panel_query
    selected_panel_ind <- as.numeric(strsplit(selected_panel, split = " ")[[1]][2])

    query_panel_ind <- (selected_panel_ind - 1)
    if(query_panel_ind == 0) query_panel_ind <- 1
    updateTabsetPanel(session, "image_tab_panel_ref", paste0("Ref. ", query_panel_ind))

    if(selected_panel_ind == 1){
      updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", selected_panel_ind+1))
      updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",selected_panel_ind+1, "->", centre))
    } else {
      query_panel_ind <- selected_panel_ind
      updateTabsetPanel(session, "image_tab_panel_reg_query", paste0("Reg. ",query_panel_ind, "->", centre))
    }
  })

  # observe changes in the registered query tab panel
  observeEvent(input$image_tab_panel_reg_query,{
    selected_panel <- input$image_tab_panel_reg_query
    selected_panel_ind <- strsplit(selected_panel, split = " ")[[1]][2]
    selected_panel_ind <- as.numeric(strsplit(selected_panel_ind, split = "->")[[1]][1])
    updateTabsetPanel(session, "image_tab_panel_query", paste0("Query ", selected_panel_ind))
  })
}

####
# Managing Cells/Barcodes ####
####

#' getRegisteredObject
#'
#' get a registered Spatial data object
#'
#' @param obj_list a list of spatial data object
#' @param mapping_list a list of transformation matrices
#' @param register_ind the indices of query images/spatialdatasets
#' @param centre the index of the central referance image/spatialdata
#'
getRegisteredObject <- function(obj_list, mapping_list, register_ind, centre, ...) {

  # check if the elements are VoltRon
  if(all(sapply(obj_list, class) == "VoltRon")){
    sr <- getRegisteredObject.VoltRon(obj_list, mapping_list, register_ind, centre, ...)
    registered_sr <- list()
    for(i in 1:length(sr)){
      registered_sr[[i]] <- sr[[i]]
    }
    return(registered_sr)

  # check if elements are Seurat
  } else if(all(sapply(obj_list, class) == "Seurat")) {
    return_list <- mapply(function(o,t) {
      if(is.null(t))
        return(NULL)
      getRegisteredObject.Seurat(o, t, ...)
    }, obj_list, mapping_list)
    return(return_list)
  }
}

#' getRegisteredObject.VoltRon
#'
#' get registered and merged VoltRon object composed of several Samples
#'
#' @param sr a list of VoltRon objects
#' @param mapping_list a list of transformation matrices
#' @param register_ind the indices of query images/spatialdatasets
#' @param centre the index of the central reference image/spatialdata
#' @param reg_mode the registration mode, either "auto" or "manual"
#'
#' @importFrom Morpho applyTransform
#'
getRegisteredObject.VoltRon <- function(sr, mapping_list, register_ind, centre, reg_mode = "manual"){

  # initiate registered VoltRon objects
  ref_ind <- centre
  registered_sr <- sr
  for(i in register_ind){
    mapping <- mapping_list[[i]]
    coords <- Coordinates(registered_sr[[i]])
    entities <- rownames(coords)
    for(kk in 1:length(mapping)){
      cur_mapping <- mapping[[kk]]
      if(reg_mode == "manual"){
        coords <- Morpho::applyTransform(coords, cur_mapping)
      } else {
        info <- image_info(Image(registered_sr[[i]])[[1]])
        coords[,2] <- info$height - coords[,2]
        coords <- perspectiveTransform(coords, cur_mapping)
        coords[,2] <- info$height - coords[,2]
      }
    }
    rownames(coords) <- entities
    Coordinates(registered_sr[[i]], reg = TRUE) <- coords
  }
  return(registered_sr)
}

#' getRegisteredObject.Seurat
#'
#' get the registered cell data/locations from a Seurat Object
#'
#' @param seu Seurat object
#' @param mapping a list of transformation mapping matrices for the registration
#'
getRegisteredObject.Seurat <- function(seu, mapping){

  image_classes <- sapply(seu@images, class)

  # Xenium (FOV)
  if(any(grepl("FOV",image_classes))){

    # get image and FOV data
    imagedata <- seu@images[[names(seu@images)[image_classes == "FOV"]]]
    image <- seu@images[[names(seu@images)[image_classes == "FOVImage"]]]

    # flip y coords and adjust coordinates to image scale
    cells <- imagedata$centroids@coords
    cells_box <- imagedata$segmentation@bbox
    cells[,2] <- max(cells[,2]) - cells[,2] + min(cells[,2])
    cells <- rescaleXeniumCells(cells, t(cells_box), image@image)

    # apply transformation to cells
    registered_cells <- as.matrix(cells)
    for(kk in 1:length(mapping)){
      cur_mapping <- mapping[[kk]]
      registered_cells <- applyTransform(registered_cells, cur_mapping)
    }
    registered_cells <- data.frame(x = registered_cells[,1], y = registered_cells[,2],
                                   cell = imagedata$centroids@cells)
    registered_segmentation_data <- list(centroids = CreateCentroids(registered_cells))
    coords <- CreateFOV(coords = registered_segmentation_data, type = "centroids", molecules = NULL, assay = "Spatial")
    seu[["registered_FOV"]] <- coords

  } else if(any(grepl("Visium",image_classes))) {

    # IMPLEMENT THIS LATER!!!!
    imagedata <- seu@images[[names(seu@images)[which(grepl("Visium", image_classes))]]]
    cells <- imagedata@coordinates[,c("imagerow", "imagecol")]
  }

  return(seu)
}

#' rescaleXeniumCells
#'
#' rescale Xenium cells coordinates for image registration
#'
#' @param cells coordinates of the cells from the Xenium assays
#' @param bbox the surrounding box of the Xenium cell coordinates
#' @param image reference image
#'
rescaleXeniumCells <- function(cells, bbox, image){

  # get image scales
  scales <- unlist(image_info(image)[c("width","height")])

  # rescale cell locations
  cells[,1] <- (cells[,1] - bbox[1,1])/(bbox[2,1] - bbox[1,1])
  cells[,1] <- cells[,1] * scales[1]
  cells[,2] <- (cells[,2] - bbox[1,2])/(bbox[2,2] - bbox[1,2])
  cells[,2] <- cells[,2] * scales[2]

  # return
  return(cells)
}

####
# Managing Keypoints ####
####

#' initateKeypoints
#'
#' Initiate shiny reactive values for keypoint dataframes for pairwise reference and query images
#'
#' @param len_images the length of images
#' @param keypoints_list the list of keypoint pairs
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @return a shiny reactive values object
#'
initateKeypoints <- function(len_images, keypoints_list, input, output, session){

  # initiate keypoints
  if(is.null(keypoints_list)){
    keypoints_list <- lapply(1:(len_images-1), function(i) {
      list(ref = tibble(KeyPoint = numeric(), x = numeric(), y = numeric()),
           query = tibble(KeyPoint = numeric(), x = numeric(), y = numeric()))

    })

    # set names for keypoints
    names(keypoints_list) <- paste0(1:(len_images-1),"-",2:len_images)
  }

  # return keypoints as reactive values
  do.call("reactiveValues", keypoints_list)
}

#' manageKeypoints
#'
#' A list of shiny observe events for tables and auxiliary operations for pairwise reference and query image
#'
#' @param centre center image index
#' @param register_ind query image indices
#' @param xyTable_list a list of keypoints x,y coordinates for each magick image
#' @param image_list a lost of magick image
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
manageKeypoints <- function(centre, register_ind, xyTable_list, image_list, input, output, session){

  # get image types
  image_types <- c("ref","query")

  # get the length of tables
  len_tables <- length(xyTable_list)

  # set click operations for reference and query points
  lapply(1:len_tables, function(i){
    lapply(image_types, function(type){

      # listen to click operations for reference/query plots
      observeEvent(input[[paste0("click_plot_", type ,i)]], {

        # get and transform keypoints
        keypoint <- data.frame(x = input[[paste0("click_plot_",type,i)]]$x,
                               y = input[[paste0("click_plot_",type,i)]]$y)
        if(is.reactive(image_list[[i]])){
          image <- image_list[[i]]()
        } else {
          image <- image_list[[i]]
        }
        image <- image[[type]]
        keypoint <- transform_keypoints(image, keypoint, paste0(type, "_image",i), input, session)

        # insert keypoint to associated table
        ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
        temp <- xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]
        temp <- temp %>%
          add_row(KeyPoint = nrow(temp)+1, x = keypoint$x, y = keypoint$y)
        xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]] <- temp
      })
    })
  })

  # remove keypoints from images
  lapply(1:len_tables, function(i){
    lapply(image_types, function(type){
      observeEvent(input[[paste0("remove_", type, i)]], {
        ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
        temp <- xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]
        if(nrow(temp) > 0){
          temp <- temp %>% filter(KeyPoint != nrow(temp))
          xyTable_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]] <- temp
        }
      })
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

  if(is.null(keypoints))
    return(list(image = image, keypoints = keypoints))

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

  if(is.null(keypoints))
    return(image)

  # select keypoints and texts on image
  image <- image +
    geom_point(mapping = aes(x = x, y = y), keypoints, size = 5, shape = 21, fill = "white") +
    geom_text(mapping = aes(x = x, y = y, label = KeyPoint), keypoints, size = 3)
}

####
# Managing Images ####
####

#' ImageOutput
#'
#' Shiny outputs for a set of magick images with keypoints
#'
#' @param image_list a list of magick images
#' @param keypoints_list a list of data frames, each having a set of keypoints
#' @param centre the center image index
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
ImageOutput <- function(image_list, keypoints_list = NULL, centre, input, output, session){

  # get image types
  image_types <- c("ref","query")

  # get the length of images
  len_images <- length(image_list)

  # output query images
  lapply(1:len_images, function(i){
    lapply(image_types, function(type){
      output[[paste0("plot_", type, i)]] <- renderPlot({

        # select keypoints
        ref_ind <- ifelse(type == "ref", i, i-1) # select reference image
        keypoints <- keypoints_list[[paste0(ref_ind, "-", ref_ind+1)]][[type]]

        # transform image and keypoints
        img <- transform_magick_image_keypoints(image_list[[i]], keypoints, paste0(type, "_image",i), input, session)
        img <- image_ggplot_keypoint(image_ggplot(img$image), img$keypoints)
        return(img)
      })
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

#' transform_magick_image_query_list
#'
#' Apply given transformations to a list of magick image and return shiny reactive
#'
#' @param image_list magick image
#' @param input shiny input
#' @param session shiny session
#'
#' @return magick image
#'
transform_magick_image_query_list <- function(image_list, input, session){

  # length of images
  len_register <- length(image_list) - 1

  trans_query_list <- lapply(1:len_register, function(i){
    reactive({
      list(ref = transform_magick_image(image_list[[i]], paste0("ref_image",i), input, session),
           query = transform_magick_image(image_list[[i+1]], paste0("query_image",i), input, session))
    })
  })
  return(trans_query_list)
}

####
# Manual Image Registration ####
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
  names(matrix_list) <- 1:len_images

  # return matrices as reactive values
  do.call("reactiveValues", matrix_list)
}

#' getManualRegisteration
#'
#' Manual registeration of images using manually entered keypoints
#'
#' @param registered_spatdata_list a list of registered Spatial data object of the query images, updated with image registration
#' @param spatdata_list a list of Spatial data object of the query images
#' @param image_list the list of query images
#' @param keypoints_list a list of keypoints x,y coordinates for query image
#' @param centre center image index
#' @param register_ind query image indices
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @import Morpho
#'
getManualRegisteration <- function(registered_spatdata_list, spatdata_list, image_list, keypoints_list,
                                   centre, register_ind, input, output, session){

  # the number of registrations
  len_register <- length(image_list) - 1

  # Registration events
  observeEvent(input$manualregister, {

    # Check keypoints
    keypoints_check_flag <- sapply(keypoints_list, function(key_list){
      nrow(key_list$ref) == nrow(key_list$query)
    })
    if(!all(unlist(keypoints_check_flag))){
      showNotification("The number of reference and query keypoints should be equal for all pairwise spatial datasets \n")
      return(NULL)
    }

    # Register keypoints
    mapping_list <- list()
    for(i in register_ind){

      # get a sequential mapping between a query and reference image
      mapping <- computeManualPairwiseTransform(keypoints_list, query_ind = i, ref_ind = centre)

      # save transformation matrix
      mapping_list[[i]] <- mapping
    }

    # get registered spatial datasets
    temp_reg_list <- getRegisteredObject(spatdata_list, mapping_list, register_ind, centre)
    for(i in 1:length(temp_reg_list))
      registered_spatdata_list[[paste0(i)]] <- temp_reg_list[[i]]

    # Plot registered images
    lapply(register_ind, function(i){
      cur_mapping <- mapping_list[[i]]
      images <- getManualRegisteredImage(image_list, cur_mapping, query_ind = i, ref_ind = centre, input)
      output[[paste0("plot_query_reg",i)]] <- renderPlot({
        info <- image_info(image_list[[centre]])
        r2 <- as.raster(images$query)
        r2 <- rasterGrob(apply(r2,2,scales::alpha, alpha = input[[paste0("plot_query_reg_alpha",i)]]))
        ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes_string("x","y")) +
          ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
          ggplot2::annotation_raster(image_list[[centre]], 0, info$width, info$height, 0, interpolate = FALSE) +
          ggplot2::annotation_custom(r2, 0, info$width, 0, info$height)
      })
    })

    # Output summary
    output[["summary"]] <- renderUI({
      str1 <- paste0(" Registration Summary:")
      str2 <- paste0("# of Images: ", length(image_list))
      str3 <- paste0("# of Registrations: ", len_register)
      all_str <- c(str1, str2, str3)
      HTML(paste(all_str, collapse = '<br/>'))
    })
  })
}

computeManualPairwiseTransform <- function(keypoints_list, query_ind, ref_ind){

  # determine the number of transformation to map from query to the reference
  indices <- query_ind:ref_ind
  mapping <- rep(indices,c(1,rep(2,length(indices)-2),1))
  mapping <- matrix(mapping,ncol=2,byrow=TRUE)

  # reference and target landmarks/keypoints
  mapping_list <- list()
  for(kk in 1:nrow(mapping)){
    cur_map <- mapping[kk,]
    if(which.min(cur_map) == 1){
      key_ind <- paste0(cur_map[1], "-", cur_map[2])
      keypoints <- keypoints_list[[key_ind]]
      target_landmark <- as.matrix(keypoints[["ref"]][,c("x","y")])
      reference_landmark <- as.matrix(keypoints[["query"]][,c("x","y")])
    } else {
      key_ind <- paste0(cur_map[2], "-", cur_map[1])
      keypoints <- keypoints_list[[key_ind]]
      reference_landmark <- as.matrix(keypoints[["ref"]][,c("x","y")])
      target_landmark <- as.matrix(keypoints[["query"]][,c("x","y")])
    }
    # compute and get transformation matrix
    mapping_list[[kk]] <- computeTransform(reference_landmark, target_landmark, type = "tps")
  }

  return(mapping_list)
}

getManualRegisteredImage <- function(images, transmatrix, query_ind, ref_ind, input){

  # plot with raster
  ref_image_raster <- as.raster(images[[ref_ind]]) |> as.matrix() |> rast()
  query_image_raster <- as.raster(images[[query_ind]]) |> as.matrix() |> rast() |> stack()

  # prepare image
  imageEx <- raster::extent(stack(ref_image_raster))
  imageRes <- raster::res(stack(ref_image_raster))
  query_image_raster_1 <- raster::as.data.frame(query_image_raster[[1]], xy = TRUE)

  # apply transformation as many as it is needed
  query_image_raster_1_t <- as.matrix(query_image_raster_1)[,1:2]
  for(trans in transmatrix){
    query_image_raster_1_t <- Morpho::applyTransform(query_image_raster_1_t, trans)
  }

  # finalize image
  r <- raster::raster(nrow = dim(query_image_raster)[1], ncol = dim(query_image_raster)[2], resolution = c(1,1))
  raster::extent(r) <- imageEx
  raster::res(r) <- imageRes
  query_image_raster_1_tr <- raster::rasterize(query_image_raster_1_t, field = query_image_raster_1[,3], r, fun = mean)
  query_image_raster_1_trf <- focal(query_image_raster_1_tr,
                                    w = matrix(1, nrow = 3, ncol = 3),
                                    fun = fill.na, pad = TRUE, na.rm = FALSE)
  query_image_raster_1_trf <- terra::rast(query_image_raster_1_trf, crs = "")

  return(list(ref = ref_image_raster, query = query_image_raster_1_trf))
}

####
# Manual Image Registration ####
####

#' getManualRegisteration
#'
#' Manual registeration of images using manually entered keypoints
#'
#' @param registered_spatdata_list a list of registered Spatial data object of the query images, updated with image registration
#' @param spatdata_list a list of Spatial data object of the query images
#' @param image_list the list of query images
#' @param centre center image index
#' @param register_ind query image indices
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#'
#' @import Morpho
#'
getAutomatedRegisteration <- function(registered_spatdata_list, spatdata_list, image_list, centre, register_ind,
                                      input, output, session){

  # the number of registrations
  len_register <- length(image_list) - 1

  # Registration events
  observeEvent(input$automaticregister, {

    # Register keypoints
    mapping_list <- list()
    aligned_image_list <- list()
    alignment_image_list <- list()
    for(i in register_ind){

      # get a sequential mapping between a query and reference image
      results <- computeAutomatedPairwiseTransform(image_list, query_ind = i, ref_ind = centre, input)

      # save transformation matrix
      mapping_list[[i]] <- results$mapping

      # save alignment
      aligned_image_list[[i]] <- results$aligned_image

      # save matches
      alignment_image_list[[i]] <- results$alignment_image
    }

    # get registered spatial datasets
    temp_reg_list <- getRegisteredObject(spatdata_list, mapping_list, register_ind, centre, reg_mode = "auto")
    for(i in 1:length(temp_reg_list))
      registered_spatdata_list[[paste0(i)]] <- temp_reg_list[[i]]

    # Plot registered images
    lapply(register_ind, function(i){
      cur_mapping <- mapping_list[[i]]
      cur_aligned_image <- aligned_image_list[[i]]
      output[[paste0("plot_query_reg",i)]] <- renderPlot({
        info <- image_info(image_list[[centre]])
        r2 <- as.raster(cur_aligned_image)
        r2 <- rasterGrob(apply(r2,2,scales::alpha, alpha = input[[paste0("plot_query_reg_alpha",i)]]))
        p <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes_string("x","y")) +
          ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
          ggplot2::annotation_raster(image_list[[centre]], 0, info$width, info$height, 0, interpolate = FALSE) +
          ggplot2::annotation_custom(r2, 0, info$width, 0, info$height)
        ggsave(filename = paste0("plot_query_reg_alpha",i, "pdf"), plot = p, device = "pdf", width = 7, height = 10, bg = "transparent")
        ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes_string("x","y")) +
          ggplot2::coord_fixed(expand = FALSE, xlim = c(0, info$width), ylim = c(0, info$height)) +
          ggplot2::annotation_raster(image_list[[centre]], 0, info$width, info$height, 0, interpolate = FALSE) +
          ggplot2::annotation_custom(r2, 0, info$width, 0, info$height)
      })
    })

    # Plot Alignment
    lapply(register_ind, function(i){
      cur_alignment_image <- alignment_image_list[[i]]
      output[[paste0("plot_alignment",i)]] <- renderPlot({
        image_ggplot(cur_alignment_image)
      })
    })

    # Output summary
    output[["summary"]] <- renderUI({
      str1 <- paste0(" Registration Summary:")
      str2 <- paste0("# of Images: ", length(image_list))
      str3 <- paste0("# of Registrations: ", len_register)
      all_str <- c(str1, str2, str3)
      HTML(paste(all_str, collapse = '<br/>'))
    })
  })
}

computeAutomatedPairwiseTransform <- function(image_list, query_ind, ref_ind, input){

  # determine the number of transformation to map from query to the reference
  indices <- query_ind:ref_ind
  mapping_mat <- rep(indices,c(1,rep(2,length(indices)-2),1))
  mapping_mat <- matrix(mapping_mat,ncol=2,byrow=TRUE)

  # reference and target landmarks/keypoints
  mapping <- list()
  aligned_image <- image_list[[query_ind]]
  for(kk in 1:nrow(mapping_mat)){
    cur_map <- mapping_mat[kk,]
    ref_image <- image_list[[cur_map[2]]]

    # compute and get transformation matrix
    reg <- automated_registration_rcpp(ref_image = ref_image, query_image = aligned_image,
                                       as.numeric(input$GOOD_MATCH_PERCENT), as.numeric(input$MAX_FEATURES))
    # reg <- automated_registration_rcpp(ref_image = ref_image, query_image = aligned_image,  0.20, 1000)
    mapping[[kk]] <- reg$transmat
    aligned_image <- reg$aligned_image
    alignment_image <- reg$alignment_image
  }

  return(list(mapping = mapping, aligned_image = aligned_image, alignment_image = alignment_image))
}

#' automated_registration
#'
#' automated registration with Rcpp/C++
#'
#' @param ref_image reference image
#' @param query_image query image
#' @param GOOD_MATCH_PERCENT the percentage of good matching keypoints
#' @param MAX_FEATURES maximum number of detected features, i.e. keypoints
#'
automated_registration_rcpp <- function(ref_image, query_image, GOOD_MATCH_PERCENT = 0.15, MAX_FEATURES = 500) {
  ref_image_rast <- image_data(ref_image)
  query_image_rast <- image_data(query_image)
  reg <- automated_registeration_rawvector(ref_image = ref_image_rast, query_image = query_image_rast,
                                           width1 = dim(ref_image_rast)[2], height1 = dim(ref_image_rast)[3],
                                           width2 = dim(query_image_rast)[2], height2 = dim(query_image_rast)[3],
                                           GOOD_MATCH_PERCENT, MAX_FEATURES)
  aligned_image <- image_read(reg[[2]])
  alignment_image <- image_read(reg[[3]])
  return(list(transmat = reg[[1]], aligned_image = aligned_image, alignment_image = alignment_image))
}

getAutomatedRegisteredImage <- function(images, transmatrix, query_ind, ref_ind, input){

  # plot with raster
  ref_image_raster <- as.raster(images[[ref_ind]]) |> as.matrix() |> rast()
  query_image_raster <- as.raster(images[[query_ind]]) |> as.matrix() |> rast() |> stack()

  # prepare image
  imageEx <- raster::extent(stack(ref_image_raster))
  imageRes <- raster::res(stack(ref_image_raster))
  query_image_raster_1 <- raster::as.data.frame(query_image_raster[[1]], xy = TRUE)

  # apply transformation as many as it is needed
  query_image_raster_1_t <- as.matrix(query_image_raster_1)[,1:2]
  for(trans in transmatrix){
    query_image_raster_1_t <- Morpho::applyTransform(query_image_raster_1_t, trans)
  }

  # finalize image
  r <- raster::raster(nrow = dim(query_image_raster)[1], ncol = dim(query_image_raster)[2], resolution = c(1,1))
  raster::extent(r) <- imageEx
  raster::res(r) <- imageRes
  query_image_raster_1_tr <- raster::rasterize(query_image_raster_1_t, field = query_image_raster_1[,3], r, fun = mean)
  query_image_raster_1_trf <- focal(query_image_raster_1_tr,
                                    w = matrix(1, nrow = 3, ncol = 3),
                                    fun = fill.na, pad = TRUE, na.rm = FALSE)
  query_image_raster_1_trf <- terra::rast(query_image_raster_1_trf, crs = "")

  return(list(ref = ref_image_raster, query = query_image_raster_1_trf))
}

