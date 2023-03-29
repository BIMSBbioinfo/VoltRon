#' GenerateXeniumDAPIImage
#'
#' Generates a low resolution DAPI image from
#'
#' @param dir.path Xenium output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param resolutionLevel the level of resolution within Xenium OME-TIFF image. Default: 7 (553x402)
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param ... additional parameters passed to the EBImage::writeImage function
#'
#' @import RBioFormats
#' @import EBImage
#'
#' @export
#'
GenerateXeniumDAPIImage <- function(dir.path, increase.contrast = TRUE, resolutionLevel = 7, output.path = NULL, ...) {

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

#' GenerateCosMxImage
#'
#' Generates a low resolution Morphology image from CosMx
#'
#' @param dir.path CosMx output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param ... additional parameters passed to the EBImage::writeImage function
#'
#' @import RBioFormats
#' @import EBImage
#'
#' @export
#'
GenerateCosMxImage <- function(dir.path, increase.contrast = TRUE, output.path = NULL, ...) {

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
    print(i)
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
    morphology_image <- image_contrast(morphology_image, sharpen = 1)
  }

  # write to the same folder
  cat("Writing Tiff File \n")
  # morphology_image <- reshape2::acast(magick::image_raster(morphology_image), y ~ x)
  if(is.null(output.path)){
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = file.path, ...)
  } else {
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = output.file, ...)
  }

  return(NULL)
}

#' addFOVImage
#'
#' Adding the Xenium image to the Seurat Object. Please run \code{GenerateXeniumDAPIImage} first.
#'
#' @param seu Seurat Object with Xenium Data
#' @param file the morphology image file created by \code{GenerateXeniumDAPIImage}.
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
    stop("FOV image was not generated. Please run GenerateXeniumDAPIImage() first!")

  # check if the image exists in the Seurat Object
  if(!is.null(seu@images[[fov_image]])){
    if(!overwrite)
      stop("FOV image is already provided")
  }

  # # get image in FOVImage class
  # image <- new(Class = "FOVImage", image = magick::image_read(file))
  #
  # # insert the image to the Seurat Object
  # seu@images[[fov_image]] <- image

  # return Seurat Object
  return(seu)
}
