####
# Objects and Classes ####
####

## FOVImage ####

#' The FOVImage class
#'
#' The FOVImage stores the morphology image from the Xenium assay
#'
#' @slot image A magick-image class object from magick package
#'
#' @name FOVImage-class
#' @rdname FOVImage-class
#' @exportClass FOVImage
#'
FOVImage <- setClass(
  Class = 'FOVImage',
  slots = list(
    'image' = 'magick-image'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'FOVImage',
  definition = function(object) {
    cat("FOV Morphology Image \n")
    cat(paste(format(image_info(object@image)), collapse = " \n "))
    return(invisible(x = NULL))
  }
)

####
# Get Images ####
####

#' @rdname vrImages
#' @method vrImages Seurat
#'
#' @importFrom magick image_read
#'
#' @export
#'
vrImages.Seurat <- function(seu, ...){

  # image class from Seurat
  image_classes <- sapply(seu@images, class)

  # get image given class
  if(any(grepl("FOV",image_classes))){
    image <- seu@images[[names(seu@images)[which(grepl("FOVImage", image_classes))]]]
    image <- image@image
  } else if(any(grepl("Visium",image_classes))) {
    image <- seu@images[[names(seu@images)[which(grepl("Visium", image_classes))]]]
    image <- magick::image_read(image@image)
  }

  # return image
  return(image)
}

#' @rdname vrImages
#' @method vrImages VoltRon
#'
#' @export
#'
vrImages.VoltRon <- function(object, ...){
  sapply(object@samples, function(samp) vrImages(samp), USE.NAMES = TRUE)
}

#' @rdname vrImages
#' @method vrImages vrSample
#'
#' @export
#'
vrImages.vrSample <- function(object, ...){
  sapply(object@layer, function(lay) vrImages(lay), USE.NAMES = TRUE)
}

#' @rdname vrImages
#' @method vrImages vrLayer
#'
#' @export
#'
vrImages.vrLayer <- function(object, ...){
  sapply(object@assay, function(assy) vrImages(assy), USE.NAMES = TRUE)
}

#' @rdname vrImages
#' @method vrImages vrAssay
#'
#' @export
#'
vrImages.vrAssay <- function(object){
  object@image
}

#' @param value new image
#'
#' @rdname vrImages
#' @method vrImages<- vrAssay
#'
#' @export
#'
"vrImages<-.vrAssay" <- function(object, reg = FALSE, ..., value) {
  object@image <- value
  return(object)
}

####
# Managing Images ####
####

#' addFOVImage
#'
#' Adding the Xenium image to the Seurat Object. Please run \code{generateXeniumImage} first.
#'
#' @param seu Seurat Object with Xenium Data
#' @param file the morphology image file created by \code{generateXeniumImage}.
#' @param fov FOV name, preferably the name used in the Seurat Object
#' @param overwrite Overwrite the existing FOV image
#'
#' @export
#'
addFOVImage <- function(seu, file, fov = "fov", overwrite = FALSE) {

  # fov image name
  fov_image <- paste0(fov, "_image")

  # check if the image file exists
  if(!file.exists(file))
    stop("FOV image was not generated. Please run generateXeniumImage() first!")

  # check if the image exists in the Seurat Object
  if(!is.null(seu@images[[fov_image]])){
    if(!overwrite)
      stop("FOV image is already provided")
  }

  # get image in FOVImage class
  image <- new(Class = "FOVImage", image = magick::image_read(file))

  # insert the image to the Seurat Object
  seu@images[[fov_image]] <- image

  # return Seurat Object
  return(seu)
}

#' @rdname resizeImage
#' @method resizeImage VoltRon
#'
#' @export
#'
resizeImage.VoltRon <- function(object, ...){
  sample.metadata <- SampleMetadata(object)
  assay_names <- rownames(sample.metadata)
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- resizeImage(cur_assay, ...)
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

  # resize coordinates
  # sizefactor <- image_info(object@image)$width
  sizefactor <- image_info(vrImages(object))$width
  object@coords <- (object@coords)*size/sizefactor

  # resize segments
  object@segments <- lapply(object@segments, function(x) x*size/sizefactor)

  # resize images
  size <- paste0(size,"x")
  object@image <- image_resize(object@image, geometry = size)

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
#' @param resolutionLevel the level of resolution within Xenium OME-TIFF image. Default: 7 (553x402)
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param ... additional parameters passed to the EBImage::writeImage function
#'
#' @importFrom RBioFormats read.image
#' @importFrom EBImage writeImage
#'
#' @export
#'
generateXeniumImage <- function(dir.path, increase.contrast = TRUE, resolutionLevel = 7, output.path = NULL, ...) {

  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/morphology_lowres.tif")
  output.file <- paste0(output.path, "/morphology_lowres.tif")

  # check if the file exists in either Xenium output folder, or the specified location
  if(file.exists(file.path) | file.exists(paste0(output.file))){
    cat("morphology_lowres.tif already exists! \n")
    return(NULL)
  }

  # read ome tiff file
  # options(java.parameters = "-Xmx4g") update 4g if more memory needed, see RBioFormats
  cat("Loading morphology_mip.ome.tif \n")
  morphology_image <- RBioFormats::read.image(paste0(dir.path, "/morphology_mip.ome.tif"))

  # pick a resolution level
  morphology_image_lowres <- morphology_image[[resolutionLevel]]
  image_info <- morphology_image_lowres@metadata$coreMetadata
  cat(paste0("Image Resolution (X:", image_info$sizeX, " Y:", image_info$sizeY, ") \n"))

  # increase contrast using EBImage
  if(increase.contrast) {
    cat("Increasing Contrast \n")
    morphology_image_lowres <- (morphology_image_lowres/max(morphology_image_lowres))
  }

  # write to the same folder
  cat("Writing Tiff File \n")
  if(is.null(output.path)){
    EBImage::writeImage(morphology_image_lowres, file = file.path, ...)
  } else {
    EBImage::writeImage(morphology_image_lowres, file = output.file, ...)
  }

  return(NULL)
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
#'
#' @import shiny
#' @import shinyjs
#'
demuxVoltRon <- function(object, scale_width = 800)
{
  # get images
  images <- vrImages(object)

  # check if there are only one assay in the object
  if(length(images) > 1)
    stop("You can only subset a VoltRon assay with one image")

  # scale
  imageinfo <- image_info(images[[1]])
  scale_factor <- imageinfo$width/800
  scale_width <- paste0(scale_width, "x")
  images <- image_scale(images[[1]], scale_width)

  # get the ui and server
  if (interactive()){
    # ui <- tagList(
    ui <- fluidPage(
      # use javascript extensions for Shiny
      useShinyjs(),

      sidebarLayout(position = "left",

                    # Side bar
                    sidebarPanel(

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
      selected_corners <- reactiveVal(tibble(x = numeric(), y = numeric()))

      # selected corner list
      selected_corners_list <- reactiveVal(tibble(box = character()))

      ### Main observable ####
      observe({

        # update summary
        output[["summary"]] <- renderUI({
          if(nrow(selected_corners_list()) > 0){
            corners <- selected_corners_list()$box
            print_selected <- paste0("Subset ", 1:length(corners), ": ", corners)
            HTML(paste(print_selected, collapse = '<br/>'))
          }
        })


        # output image
        output[["cropped_image"]] <- renderPlot({
          if(nrow(selected_corners()) < 2){
            corners <- apply(as.matrix(selected_corners()),2,as.numeric)
            image_ggplot(images)
          } else {
            corners <- apply(as.matrix(selected_corners()),2,as.numeric)
            image_ggplot(images) +
              geom_rect(aes(xmin = corners[1,1], xmax = corners[2,1], ymin = corners[1,2], ymax = corners[2,2]),
                        fill = "green", alpha = 0.3, color = "black")
          }
        })
      })

      ## reset points ####
      observeEvent(input$resetpoints, {
        selected_corners() %>%
          filter(FALSE) %>% selected_corners()
      })

      ## add box ####
      observeEvent(input$addbox, {
        if(nrow(selected_corners()) == 2){
          next_ind <- length(selected_corners_list()) + 1
          corners <- selected_corners()
          corners <- corners*scale_factor
          corners <- apply(corners,2,ceiling)
          corners <- paste0(abs(corners[2,1]-corners[1,1]), "x",
                            abs(corners[2,2]-corners[1,2]), "+",
                            min(corners[,1]), "+", imageinfo$height - max(corners[,2]))
          selected_corners_list() %>%
            add_row(box = corners) %>%
            selected_corners_list()
          selected_corners() %>%
            filter(FALSE) %>% selected_corners()
        }
      })

      ## select points on the image ####
      observeEvent(input$choose_corner, {
        if(nrow(selected_corners()) < 2) {
          keypoint <- data.frame(x = input[["choose_corner"]]$x,
                                 y = input[["choose_corner"]]$y)
          selected_corners() %>%
            add_row(x = keypoint$x, y = keypoint$y) %>%
            selected_corners()
        } else {
          selected_corners() %>%
            filter(FALSE) %>% selected_corners()
        }
      })

      ## select points on the image ####
      observeEvent(input$done, {
        if(nrow(selected_corners_list()) > 0){
          registration <- list()
          box_list <- selected_corners_list()
          sample_names <- paste0("Sample", 1:length(box_list$box))
          for(i in 1:length(box_list$box)){
            temp <- subset(object, image = box_list$box[i])
            temp$Sample <- sample_names[i]
            registration[[sample_names[i]]] <- temp
          }
          stopApp(list(registration = registration, subset_info_list = box_list))
        }
      })
    }

    shiny::runApp(shiny::shinyApp(ui, server))
  }
}

