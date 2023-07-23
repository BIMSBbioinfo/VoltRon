####
# Objects and Classes ####
####

## vrSample ####

#' The vrSample (VoltRon Sample) Class
#'
#' @slot layer A list of layers (vrLayer)
#'
#' @name vrSample-class
#' @rdname vrSample-class
#' @exportClass vrSample
#'
vrSample <- setClass(
  Class = 'vrSample',
  slots = c(
    layer = 'list'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrSample',
  definition = function(object) {
    cat(class(x = object), "(VoltRon Sample) Object \n")
    layers <- names(unlist(object@layer))
    cat("Layer(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset ####

#' @importFrom methods slot
#'
setMethod(
  f = '[[',
  signature = 'vrSample',
  definition = function(x, i, j){

    # sample names
    layer_names <- names(methods::slot(x, "layer"))

    # check query sample name
    if(!i %in% layer_names){
      stop("There are no layers named ", i, " in this sample")
    }

    # return samples
    return(x@layer[[i]])
  }
)

## vrLayer ####

#' The vrLayer (VoltRon Layer) Class
#'
#' @slot assay A list of assays (vrAssay)
#'
#' @name vrLayer-class
#' @rdname vrLayer-class
#' @exportClass vrLayer
#'
vrLayer <- setClass(
  Class = 'vrLayer',
  slots = c(
    assay = 'list'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrLayer',
  definition = function(object) {
    cat(class(x = object), "(VoltRon Layer) Object \n")
    # cat("This object includes", length(object@assay), "assays \n")
    layers <- names(unlist(object@assay))
    cat("Assay(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset of assays ####

#' @importFrom methods slot
#'
setMethod(
  f = '[[',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ...){

    # if no assay were found, check sample names
    assay_names <- names(methods::slot(x, "assay"))

    # check query sample name
    if(!i %in% assay_names){
      stop("There are no assays named ", i, " in this object")
    } else {
      return(x@assay[[i]])
    }
  }
)

#' @importFrom methods slot
#'
setMethod(
  f = '[[<-',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ..., value){

    # if no assay were found, check sample names
    assay_names <- names(methods::slot(x, "assay"))

    # check query sample name
    if(!i %in% assay_names){
      stop("There are no assays named ", i, " in this object")
    }

    x@assay[[i]] <- value
    return(x)
  }
)

####
# Methods ####
####

### Merge vrSample object ####

#' Merging vrSample objects
#'
#' Given a vrSample object, and a list of vrSample objects, merge all.
#'
#' @param object a vrSample object
#' @param object_list a list of vrSample objects
#' @param samples the sample names
#'
#' @method merge vrSample
#'
merge.vrSample <- function(object, object_list, samples = NULL){

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)
  names(object_list) <- samples

  # set VoltRon class
  return(object_list)
}

### Subset vrSample object ####

#' Subsetting vrSample objects
#'
#' Given a vrSample object, subset the object given one of the attributes
#'
#' @param object a vrSample object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \code{magick::image_crop}
#'
#' @method subset vrSample
#'
#' @importFrom rlang enquo
#'
subset.vrSample <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  layers <- object@layer
  if(!is.null(assays)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, assays = assays)
    }, USE.NAMES = TRUE)
  } else if(!is.null(spatialpoints)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE)
  } else if(!is.null(image)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, image = image)
    }, USE.NAMES = TRUE)
  }

  # set VoltRon class
  return(object)
}

### Get Spatial Points of vrSample ####

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrSample
#'
#' @export
#'
vrSpatialPoints.vrSample <- function(object, ...) {

  layers <- object@layer
  spatialpoints <- lapply(layers, function(lay) {
    vrSpatialPoints(lay)
  })
  do.call(c, spatialpoints)
}

### Subset vrLayer objects ####

#' Subsetting vrLayer objects
#'
#' Given a vrLayer object, subset the object given one of the attributes
#'
#' @param object a vrLayer object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \code{magick::image_crop}
#'
#' @method subset vrLayer
#'
#' @importFrom rlang enquo
#'
subset.vrLayer <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){
    object@assay  <- object@assay[names(object@assay) %in% assays]
  } else if(!is.null(spatialpoints)){
    assay_set <- object@assay
    object@assay <- sapply(assay_set, function(assy) {
      subset.vrAssay(assy, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE)
  } else if(!is.null(image)){
    assay_set <- object@assay
    object@assay <- sapply(assay_set, function(assy) {
      subset.vrAssay(assy, image = image)
    }, USE.NAMES = TRUE)
  }

  # set VoltRon class
  return(object)
}

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrLayer
#'
#' @export
#'
vrSpatialPoints.vrLayer <- function(object, ...) {
  assays <- object@assay
  spatialpoints <- lapply(assays, function(assy) {
    vrSpatialPoints(assy)
  })
  do.call(c, spatialpoints)
}
