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


#' addFOVImage
#'
#' Adding the Xenium image to the Seurat Object. Please run \code{GenerateXeniumDAPIImage} first.
#'
#' @param seu Seurat Object with Xenium Data
#' @param file the morphology image file created by \code{GenerateXeniumDAPIImage}.
#' @param fov FOV name, preferably the name used in the Seurat Object
#' @param overwrite Overwrite the existing FOV image
#'
#' @import
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

  # get image in FOVImage class
  image <- new(Class = "FOVImage", image = magick::image_read(file))

  # insert the image to the Seurat Object
  seu@images[[fov_image]] <- image

  # return Seurat Object
  return(seu)
}
