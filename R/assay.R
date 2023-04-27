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
#' @slot type the type of the assay (cell, spot, ROI)
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
    coords_reg = 'matrix',
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

#' @rdname Coordinates
#' @method Coordinates srAssay
#'
#' @export
#'
Coordinates.srAssay <- function(object, reg = FALSE, ...) {
  if(reg){
    if(nrow(object@coords_reg) < 1) {
      return(object@coords)
    } else {
      return(object@coords_reg)
    }
  } else {
    return(object@coords)
  }
}

#' @rdname Coordinates
#' @method Coordinates<- srAssay
#'
#' @export
#'
"Coordinates<-.srAssay" <- function(object, reg = FALSE, ..., value) {

  # get coordinates
  coords <- Coordinates(object, ...)

  # stop if the rownames are not matching
  if(any(sapply(rownames(values),is.null)))
    stop("Provided coordinates data does not have cell/spot/ROI names")

  if(!all(rownames(values) %in% rownames(coords)))
    stop("Cant overwrite coordinates, non-existing cells/spots/ROIs!")

  # stop if the colnames are not matching
  if(!all(colnames(values) %in% colnames(coords)))
    stop("Cant overwrite coordinates, only x or y coordinates should be provided!")

  if(reg){
    slot(object = object, name = 'coords_reg') <- value
  } else{
    slot(object = object, name = 'coords') <- value
  }
  return(object)
}
