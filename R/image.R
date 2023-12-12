####
# Get Images ####
####

#' @param object A VoltRon object
#'
#' @rdname vrImages
#' @method vrImages VoltRon
#'
#' @export
#'
vrImages.VoltRon <- function(object, ...){
  images <- sapply(object@samples, function(samp) vrImages(samp, ...), USE.NAMES = TRUE)
  if(length(images) == 1){
    return(images[[1]])
  } else {
    return(images)
  }
}

#' @param object A vrSample object
#'
#' @rdname vrImages
#' @method vrImages vrSample
#'
#' @export
#'
vrImages.vrSample <- function(object, ...){
  sapply(object@layer, function(lay) vrImages(lay, ...), USE.NAMES = TRUE)
}

#' @param object A vrLayer object
#'
#' @rdname vrImages
#' @method vrImages vrLayer
#'
#' @export
#'
vrImages.vrLayer <- function(object, ...){
  sapply(object@assay, function(assy) vrImages(assy, ...), USE.NAMES = TRUE)
}

#' @param object A vrAssay object
#' @param main_image the name of the main image
#'
#' @rdname vrImages
#' @method vrImages vrAssay
#'
#' @importFrom magick image_read
#'
#' @export
#'
vrImages.vrAssay <- function(object, main_image = NULL, as.raster = FALSE){
  if(!is.null(object@image)){
    if(is.null(main_image)) {
      main_image <- object@main_image
    }
    if(main_image %in% vrImageNames(object)){
      if(as.raster){
        return(object@image[[main_image]])
      } else {
        return(magick::image_read(object@image[[main_image]]))
      }
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' @param object A vrAssay object
#' @param main_image the name of the main image
#' @param reg TRUE if registered coordinates are assigned
#' @param value new image
#'
#' @rdname vrImages
#' @method vrImages<- vrAssay
#'
#' @importFrom magick image_data
#'
#' @export
#'
"vrImages<-.vrAssay" <- function(object, main_image = NULL, reg = FALSE, ..., value) {
  if(is.null(main_image)) {
    main_image <- object@main_image
  }
  if(inherits(value, "bitmap")){
    object@image[[main_image]] <- value
  } else if(inherits(value, "magick-image")){
    object@image[[main_image]] <- magick::image_data(value)
  } else {
    stop("Please provide either a magick-image or bitmap class image object!")
  }
  return(object)
}

#' @param assay assay
#'
#' @rdname vrImageNames
#' @method vrImageNames VoltRon
#'
#' @export
#'
vrImageNames.VoltRon <- function(object, assay = NULL){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assay types
  image_names <- unique(unlist(lapply(assay_names, function(x) vrImageNames(object[[x]]))))

  return(image_names)
}

#' @rdname vrImageNames
#' @method vrImageNames vrAssay
#'
#' @export
#'
vrImageNames.vrAssay <- function(object){
  return(names(object@image))
}

#'
#' @rdname vrMainImage
#' @method vrMainImage vrAssay
#'
#' @export
#'
vrMainImage.vrAssay <- function(object){
  return(object@main_image)
}

#'
#' @rdname vrMainImage
#' @method vrMainImage<- vrAssay
#'
#' @export
#'
"vrMainImage<-.vrAssay" <- function(object, value){
  object@main_image <- value
  return(object)
}

####
# Managing Images ####
####

#' @rdname resizeImage
#' @method resizeImage VoltRon
#'
#' @export
#'
resizeImage.VoltRon <- function(object, ...){
  sample.metadata <- SampleMetadata(object)
  assay_names <- rownames(sample.metadata)
  for(assy in assay_names){
    object[[assy]] <- resizeImage(object[[assy]], ...)
  }
  return(object)
}

#' @param size the width of the resized image
#'
#' @rdname resizeImage
#' @method resizeImage vrAssay
#'
#' @importFrom magick image_info image_resize
#' @export
#'
resizeImage.vrAssay <- function(object, size){

  # check size
  if(!is.numeric(size))
    stop("width size should be numeric")
  if(!all.equal(size, as.integer(size)) & size > 0)
    stop("width size should be a positive integer")
  if(size < 100)
    stop("width size cannot be less than 100px")

  # resize coordinates
  sizefactor <- image_info(vrImages(object))$width
  vrCoordinates(object) <- (vrCoordinates(object)*size)/sizefactor

  # resize segments
  vrSegments(object) <- lapply(vrSegments(object), function(x) x*size/sizefactor)

  # resize images
  size <- paste0(size,"x")
  image_names <- vrImageNames(object)
  for(img in image_names){
    vrImages(object, main_image = img) <- image_resize(vrImages(object, main_image = img), geometry = size)
  }

  # return
  return(object)
}

#' @rdname modulateImage
#' @method modulateImage VoltRon
#'
#' @export
#'
modulateImage.VoltRon <- function(object, ...){
  sample.metadata <- SampleMetadata(object)
  assay_names <- rownames(sample.metadata)
  for(assy in assay_names){
    object[[assy]] <- modulateImage(object[[assy]], ...)
  }
  return(object)
}

#' @param brightness modulation of brightness as percentage of the current value (100 for no change)
#' @param saturation modulation of saturation as percentage of the current value (100 for no change)
#' @param hue modulation of hue is an absolute rotation of -180 degrees to +180 degrees from the current position corresponding to an argument range of 0 to 200 (100 for no change)
#'
#' @rdname modulateImage
#' @method modulateImage vrAssay
#'
#' @importFrom magick image_info image_resize
#' @export
#'
modulateImage.vrAssay <- function(object, brightness = 100, saturation = 100, hue = 100){

  # modulate image
  image_names <- vrImageNames(object)
  for(img in image_names){
    vrImages(object, main_image = img) <- magick::image_modulate(vrImages(object, main_image = img),
                                                                 brightness = brightness, saturation = saturation, hue = hue)
  }

  # return
  return(object)
}

####
# Image File Manipulation ####
####

#' generateXeniumImage
#'
#' Generate a low resolution DAPI image of the Xenium experiment
#'
#' @param dir.path Xenium output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param resolution_level the level of resolution within Xenium OME-TIFF image. Default: 7 (553x402)
#' @param overwrite_resolution if TRUE, the image "file.name" will be generated again although it exists at "dir.path"
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param file.name the name of the lowred morphology image. Default: morphology_lowres.tif
#' @param ... additional parameters passed to the \code{EBImage::writeImage} function
#'
#' @importFrom EBImage writeImage
#'
#' @details
#' The Xenium morphology_mip.ome.tif file that is found under the outs folder comes is an hyperstack of different resolutions of the DAPI image.
#' \code{generateXeniumImage} allows extracting only one of these layers by specifying the \code{resolution} parameter (Default: 7 for 553x402) among 1 to 8.
#' Lower incides of resolutions have higher higher resolutions, e.g. 1 for 35416x25778. Note that you may need to allocate larger memory of Java to import
#' higher resolution images.
#'
#' @export
#'
generateXeniumImage <- function(dir.path, increase.contrast = TRUE, resolution_level = 7, overwrite_resolution = FALSE, output.path = NULL, file.name = "morphology_lowres.tif", ...) {

  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/", file.name)
  output.file <- paste0(output.path, "/", file.name)

  # check if the file exists in either Xenium output folder, or the specified location
  if((file.exists(file.path) | file.exists(paste0(output.file))) & !overwrite_resolution){
    message(paste0(file.name, " already exists!"))
  } else {
    message("Loading morphology_mip.ome.tif \n")
    if (!requireNamespace('RBioFormats'))
      stop("Please install RBioFormats package to read the ome.tiff file!")
    morphology_image_lowres <- RBioFormats::read.image(paste0(dir.path, "/morphology_mip.ome.tif"), resolution = resolution_level)

    # pick a resolution level
    image_info <- morphology_image_lowres@metadata$coreMetadata
    message(paste0("Image Resolution (X:", image_info$sizeX, " Y:", image_info$sizeY, ") \n"))

    # increase contrast using EBImage
    if(increase.contrast) {
      message("Increasing Contrast \n")
      morphology_image_lowres <- (morphology_image_lowres/max(morphology_image_lowres))
    }

    # write to the same folder
    message("Writing Tiff File")
    if(is.null(output.path)){
      EBImage::writeImage(morphology_image_lowres, file = file.path, ...)
    } else {
      EBImage::writeImage(morphology_image_lowres, file = output.file, ...)
    }
  }
  invisible()
}

#' generateCosMxImage
#'
#' Generates a low resolution Morphology image of the CosMx experiment
#'
#' @param dir.path CosMx output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param ... additional parameters passed to the EBImage::writeImage function
#'
#' @importFrom magick image_read image_contrast
#' @importFrom EBImage writeImage
#' @importFrom reshape2 acast
#'
#' @export
#'
generateCosMxImage <- function(dir.path, increase.contrast = TRUE, output.path = NULL, ...) {

  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/CellComposite_lowres.tif")
  output.file <- paste0(output.path, "/CellComposite_lowres.tif")

  # check if the file exists in either Xenium output folder, or the specified location
  if(file.exists(file.path) | file.exists(paste0(output.file))){
    cat("CellComposite_lowres.tif already exists! \n")
    return(NULL)
  }

  # FOV positions of CosMx
  list_of_files <- list.files(dir.path)
  fov_positions_path <- paste0(dir.path, "/", list_of_files[grepl("fov_positions_file.csv$",list_of_files)][1])
  fov_positions <- read.csv(fov_positions_path)

  # manipulate fov positions matrix
  cat("Getting FOV Positions \n")
  relative_fov_positions <- fov_positions
  x_min <- min(relative_fov_positions$x_global_px)
  y_min <- min(relative_fov_positions$y_global_px)
  x_gap <- diff(unique(fov_positions$x_global_px))[1]
  y_gap <- diff(unique(fov_positions$y_global_px))[1]
  relative_fov_positions[,c("x_global_px","y_global_px")] <- t(apply(relative_fov_positions[,c("x_global_px","y_global_px")], 1, function(cell){
    c((cell[1]-x_min)/x_gap,(cell[2]-y_min)/y_gap)
  }))
  relative_fov_positions <- relative_fov_positions[order(relative_fov_positions$y_global_px, decreasing = TRUE),]

  # Combine Images of the FOV grid
  cat("Loading FOV tif files \n")
  image.dir.path <- paste0(dir.path,"/CellComposite/")
  morphology_image_data <- NULL
  for(i in relative_fov_positions$fov){
    image_path <- paste0(image.dir.path, "CellComposite_F", str_pad(as.character(i), 3, pad = 0), ".jpg")
    image_data <- magick::image_read(image_path) %>% magick::image_resize("x500") %>% magick::image_raster()
    if(is.null(morphology_image_data))
      dim_image <- apply(image_data[,1:2], 2, max)
    scale_dim <- relative_fov_positions[i,2:3]*dim_image
    image_data[,1:2] <- image_data[,1:2] +
      rep(1, nrow(image_data)) %o% as.matrix(scale_dim)[1,]
    morphology_image_data <- rbind(morphology_image_data, image_data)
  }
  morphology_image_data_array <- reshape2::acast(morphology_image_data, y ~ x)
  morphology_image <- magick::image_read(morphology_image_data_array) %>% magick::image_resize("x800")

  # pick a resolution level
  morphology_image_info <- image_info(morphology_image)
  cat(paste0("Image Resolution (X:", morphology_image_info$width, " Y:", morphology_image_info$height, ") \n"))

  # increase contrast using EBImage
  if(increase.contrast) {
    cat("Increasing Contrast \n")
    morphology_image <- magick::image_contrast(morphology_image, sharpen = 1)
  }

  # write to the same folder
  cat("Writing Tiff File \n")
  if(is.null(output.path)){
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = file.path, ...)
  } else {
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = output.file, ...)
  }

  return(NULL)
}

####
# Demultiplex Images ####
####

#' demuxVoltRon
#'
#' Subsetting/demultiplexing of the VoltRon Object using interactive shiny app
#'
#' @param object a VoltRon object
#' @param scale_width the initial width of the object image
#' @param use_points use spatial points instead of the reference image
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom magick image_scale image_info image_ggplot
#' @importFrom ggplot2 geom_rect
#' @importFrom htmltools HTML
#' @importFrom dplyr filter add_row tibble
#'
demuxVoltRon <- function(object, scale_width = 800, use_points = FALSE)
{
  # get images
  images <- vrImages(object)

  # check if there are only one assay in the object
  if(nrow(SampleMetadata(object)) > 1)
    stop("You can only subset a VoltRon assay with one image")

  # scale
  imageinfo <- magick::image_info(images)
  scale_factor <- imageinfo$width/scale_width
  scale_width_char <- paste0(scale_width, "x")
  images <- magick::image_scale(images, scale_width_char)

  # get the ui and server
  if (interactive()){
    ui <- fluidPage(
      # use javascript extensions for Shiny
      shinyjs::useShinyjs(),

      sidebarLayout(position = "left",

                    # Side bar
                    sidebarPanel(
                      tags$style(make_css(list('.well', 'margin', '7%'))),

                      # Interface
                      fluidRow(
                        column(12,h4("Subsetting Interface"))
                      ),

                      # Buttons
                      fluidRow(
                        column(12,shiny::actionButton("resetpoints", "Reset Points")),
                        br(),
                        column(12,shiny::actionButton("addbox", "Add Box")),
                        br()
                      ),

                      # Subsets
                      fluidRow(
                        column(12,h4("Selected Sections")),
                        column(12,htmlOutput("summary")),
                        br()
                      ),

                      # Subsets
                      fluidRow(
                        column(12,h4("Finished Selecting?")),
                        column(12,shiny::actionButton("done", "Yes!"))
                      ),
                      br(),
                      # panel options
                      width = 3,
                    ),

                    mainPanel(

                      # main image
                      br(),
                      br(),
                      fluidRow(
                        imageOutput("cropped_image", click = "choose_corner", height = "1000px")
                      ),

                      # panel options
                      width = 9
                    )
      )
    )

    server <- function(input, output, session) {

      # selected corners
      selected_corners <- reactiveVal(dplyr::tibble(x = numeric(), y = numeric()))

      # selected corner list
      selected_corners_list <- reactiveVal(dplyr::tibble(box = character()))

      # the image
      if(use_points){
        object_small <- resizeImage(object, size = scale_width)
        image_info_small <- magick::image_info(vrImages(object_small))
        coords <- as.data.frame(vrCoordinates(object_small, reg = FALSE))
        # coords[,2] <- max(coords[,2]) - coords[,2] + min(coords[,2])
        pl <- ggplot() + geom_point(aes_string(x = "x", y = "y"), coords, size = 1.5, color = "black") +
          theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                axis.line=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
                legend.margin = margin(0,0,0,0), plot.margin = unit( c(0,0,0,0),"in")) +
          coord_fixed()
      } else {
        pl <- magick::image_ggplot(images)
      }

      ### Main observable ####
      observe({

        # update summary
        output[["summary"]] <- renderUI({
          if(nrow(selected_corners_list()) > 0){
            corners <- selected_corners_list()$box
            print_selected <- paste0("Subset ", 1:length(corners), ": ", corners)
            htmltools::HTML(paste(print_selected, collapse = '<br/>'))
          }
        })

        # output image
        output[["cropped_image"]] <- renderPlot({
          corners <- apply(as.matrix(selected_corners()),2,as.numeric)
          if(nrow(selected_corners()) > 1){
            pl <- pl +
              ggplot2::geom_rect(aes(xmin = corners[1,1], xmax = corners[2,1], ymin = corners[1,2], ymax = corners[2,2]),
                                 fill = "green", alpha = 0.3, color = "black")
          }
          pl
        })
      })

      ## reset points ####
      observeEvent(input$resetpoints, {
        selected_corners() %>%
          dplyr::filter(FALSE) %>% selected_corners()
      })

      ## add box ####
      observeEvent(input$addbox, {
        if(nrow(selected_corners()) == 2){
          next_ind <- length(selected_corners_list()) + 1
          corners <- selected_corners()

          # adjust corners
          corners <- corners*scale_factor
          corners <- apply(corners,2,ceiling)

          # fix for limits
          corners[1,1] <- ifelse(corners[1,1] < 0, 0, corners[1,1])
          corners[1,1] <- ifelse(corners[1,1] > imageinfo$width, imageinfo$width, corners[1,1])
          corners[2,1] <- ifelse(corners[2,1] < 0, 0, corners[2,1])
          corners[2,1] <- ifelse(corners[2,1] > imageinfo$width, imageinfo$width, corners[2,1])
          corners[1,2] <- ifelse(corners[1,2] < 0, 0, corners[1,2])
          corners[1,2] <- ifelse(corners[1,2] > imageinfo$height, imageinfo$height, corners[1,2])
          corners[2,2] <- ifelse(corners[2,2] < 0, 0, corners[2,2])
          corners[2,2] <- ifelse(corners[2,2] > imageinfo$height, imageinfo$height, corners[2,2])

          # get crop info
          corners <- paste0(abs(corners[2,1]-corners[1,1]), "x",
                            abs(corners[2,2]-corners[1,2]), "+",
                            min(corners[,1]), "+", imageinfo$height - max(corners[,2]))

          # add to box list
          selected_corners_list() %>%
            dplyr::add_row(box = corners) %>%
            selected_corners_list()
          selected_corners() %>%
            dplyr::filter(FALSE) %>% selected_corners()
        }
      })

      ## select points on the image ####
      observeEvent(input$choose_corner, {
        if(nrow(selected_corners()) < 2) {
          keypoint <- data.frame(x = input[["choose_corner"]]$x,
                                 y = input[["choose_corner"]]$y)
          selected_corners() %>%
            dplyr::add_row(x = keypoint$x, y = keypoint$y) %>%
            selected_corners()
        } else {
          selected_corners() %>%
            dplyr::filter(FALSE) %>% selected_corners()
        }
      })

      ## select points on the image ####
      observeEvent(input$done, {
        if(nrow(selected_corners_list()) > 0){
          subsets <- list()
          box_list <- selected_corners_list()
          sample_names <- paste0("Sample", 1:length(box_list$box))
          print(sample_names)
          for(i in 1:length(box_list$box)){
            temp <- subset(object, image = box_list$box[i])
            temp$Sample <- sample_names[i]
            subsets[[sample_names[i]]] <- temp
          }
          stopApp(list(subsets = subsets, subset_info_list = box_list))
        } else {
          showNotification("You have not selected a subset yet!")
        }
      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

