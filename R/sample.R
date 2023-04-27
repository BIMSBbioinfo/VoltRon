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

### Subset ####

#' @method subset srSample
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.srSample <- function(object, subset, assays = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){
    layers <- object@layer
    lapply(layers, function(lay) {
      subset.srLayer(lay, assays = assays)
    })
  }

  # set SpaceRover class
  return(object)
}

### Subset srLayer objects ####

#' @method subset srLayer
#'
#' @importFrom rlang enquo
#' @import igraph
#'
#' @export
#'
subset.srLayer <- function(object, subset, assays = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){
    object@assay  <- object@assay[names(object@assay) %in% assays]
  }

  # set SpaceRover class
  return(object)
}
