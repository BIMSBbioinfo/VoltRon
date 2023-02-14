# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))

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

# set method for FOVImage class
setMethod(
  f = 'show',
  signature = 'FOVImage',
  definition = function(object) {
    cat("FOV Morphology Image \n")
    cat(paste(format(image_info(object@image)), collapse = " \n "))
    return(invisible(x = NULL))
  }
)
