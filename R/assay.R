####
# Objects and Classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))

## srAssay ####

#' The srAssay (SpaceRover Assay) Class
#'
#' @slot rawdata raw count table
#' @slot normdata normalized count table
#' @slot coord spatial coordinates of the assay
#' @slot image image of the spatial assay
#'
#' @name srAssay-class
#' @rdname srAssay-class
#' @exportClass srAssay
#'
srAssay <- setClass(
  Class = 'srAssay',
  slots = c(
    rawdata = 'matrix',
    normdata = 'matrix',
    coords = 'matrix',
    image = 'magick-image',
    type = "character"
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'srAssay',
  definition = function(object) {
    cat("srAssay (SpaceRover Assay) of", ncol(object@rawdata), "cells and", nrow(object@rawdata), "features. \n")
    return(invisible(x = NULL))
  }
)

####
# Methods ####
####

