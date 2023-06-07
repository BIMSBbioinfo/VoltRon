####
# Objects and Classes ####
####

## vrSample ####

#' The vrSample (VoltRon Sample) Class
#'
#' @slot samples A list of layers for the this vrSample object
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
setMethod(
  f = '[[',
  signature = 'vrSample',
  definition = function(x, i, j){

    # sample names
    layer_names <- names(slot(x, "layer"))

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
#' @slot assay
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
setMethod(
  f = '[[',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ...){

    # if no assay were found, check sample names
    assay_names <- names(slot(x, "assay"))

    # check query sample name
    if(!i %in% assay_names){
      stop("There are no assays named ", i, " in this object")
    } else {
      return(x@assay[[i]])
    }
  }
)

setMethod(
  f = '[[<-',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ..., value){

    # if no assay were found, check sample names
    assay_names <- names(slot(x, "assay"))

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

#' @method merge vrSample
#'
#' @export
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

#' @method subset vrSample
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.vrSample <- function(object, subset, assays = NULL, entities = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  layers <- object@layer
  if(!is.null(assays)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, assays = assays)
    }, USE.NAMES = TRUE)
  } else if(!is.null(entities)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, entities = entities)
    }, USE.NAMES = TRUE)
  } else if(!is.null(image)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, image = image)
    }, USE.NAMES = TRUE)
  }

  # set VoltRon class
  return(object)
}

### Get entities vrSample ####

#' @rdname Entities
#' @method Entities vrSample
#'
#' @export
#'
Entities.vrSample <- function(object, ...) {

  layers <- object@layer
  entities <- lapply(layers, function(lay) {
    Entities(lay)
  })
  do.call(c, entities)
}

### Subset vrLayer objects ####

#' @method subset vrLayer
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.vrLayer <- function(object, subset, assays = NULL, entities = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){
    object@assay  <- object@assay[names(object@assay) %in% assays]
  } else if(!is.null(entities)){
    assay_set <- object@assay
    object@assay <- sapply(assay_set, function(assy) {
      subset.vrAssay(assy, entities = entities)
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

#' @rdname Entities
#' @method Entities vrLayer
#'
#' @export
#'
Entities.vrLayer <- function(object, ...) {
  assays <- object@assay
  entities <- lapply(assays, function(assy) {
    Entities(assy)
  })
  do.call(c, entities)
}
