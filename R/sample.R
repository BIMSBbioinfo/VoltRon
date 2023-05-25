####
# Objects and Classes ####
####

## srSample ####

#' The srSample (SpaceRover Sample) Class
#'
#' @slot samples A list of layers for the this srSample object
#'
#' @name srSample-class
#' @rdname srSample-class
#' @exportClass srSample
#'
srSample <- setClass(
  Class = 'srSample',
  slots = c(
    layer = 'list'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'srSample',
  definition = function(object) {
    cat(class(x = object), "(SpaceRover Sample) Object \n")
    layers <- names(unlist(object@layer))
    cat("Layer(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset ####
setMethod(
  f = '[[',
  signature = 'srSample',
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

## srLayer ####

#' The srLayer (SpaceRover Layer) Class
#'
#' @slot assay
#'
#' @name srLayer-class
#' @rdname srLayer-class
#' @exportClass srLayer
#'
srLayer <- setClass(
  Class = 'srLayer',
  slots = c(
    assay = 'list'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'srLayer',
  definition = function(object) {
    cat(class(x = object), "(SpaceRover Layer) Object \n")
    # cat("This object includes", length(object@assay), "assays \n")
    layers <- names(unlist(object@assay))
    cat("Assay(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset of assays ####
setMethod(
  f = '[[',
  signature = c('srLayer', "character"),
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
  signature = c('srLayer', "character"),
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

### Merge srSample object ####

#' @method merge srSample
#'
#' @export
#'
merge.srSample <- function(object, object_list, samples = NULL){

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)
  names(object_list) <- samples

  # set SpaceRover class
  return(object_list)
}

### Subset srSample object ####

#' @method subset srSample
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.srSample <- function(object, subset, assays = NULL, entities = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  layers <- object@layer
  if(!is.null(assays)){
    object@layer <- sapply(layers, function(lay) {
      subset.srLayer(lay, assays = assays)
    }, USE.NAMES = TRUE)
  } else if(!is.null(entities)){
    object@layer <- sapply(layers, function(lay) {
      subset.srLayer(lay, entities = entities)
    }, USE.NAMES = TRUE)
  } else if(!is.null(image)){
    object@layer <- sapply(layers, function(lay) {
      subset.srLayer(lay, image = image)
    }, USE.NAMES = TRUE)
  }

  # set SpaceRover class
  return(object)
}

### Get entities srSample ####

#' @rdname Entities
#' @method Entities srSample
#'
#' @export
#'
Entities.srSample <- function(object, ...) {

  layers <- object@layer
  entities <- lapply(layers, function(lay) {
    Entities(lay)
  })
  do.call(c, entities)
}

### Subset srLayer objects ####

#' @method subset srLayer
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.srLayer <- function(object, subset, assays = NULL, entities = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){
    object@assay  <- object@assay[names(object@assay) %in% assays]
  } else if(!is.null(entities)){
    assay_set <- object@assay
    object@assay <- sapply(assay_set, function(assy) {
      subset.srAssay(assy, entities = entities)
    }, USE.NAMES = TRUE)
  } else if(!is.null(image)){
    assay_set <- object@assay
    object@assay <- sapply(assay_set, function(assy) {
      subset.srAssay(assy, image = image)
    }, USE.NAMES = TRUE)
  }

  # set SpaceRover class
  return(object)
}

#' @rdname Entities
#' @method Entities srLayer
#'
#' @export
#'
Entities.srLayer <- function(object, ...) {
  assays <- object@assay
  entities <- lapply(assays, function(assy) {
    Entities(assy)
  })
  do.call(c, entities)
}
