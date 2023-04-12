####
# SpaceRover classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('igraph'))

## SpaceRover ####

#' The SpaceRover Class
#'
#' @slot samples A list of samples for the this project
#' @slot integrated.datasets A list of integrated data objects that indicate integrated spatial layers
#' @slot meta.data Contains meta-information about each sample
#' @slot project Name of the project
#'
#' @name SpaceRover-class
#' @rdname SpaceRover-class
#' @exportClass SpaceRover
#'
SpaceRover <- setClass(
  Class = 'SpaceRover',
  slots = c(
    samples = 'list',
    metadata = "srMetadata",
    sample.metadata = "data.frame",
    zstack = "igraph",
    main.assay = "character",
    project = 'character'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'SpaceRover',
  definition = function(object) {

    # print class
    cat(class(x = object), "Object \n")

    # print samples and layers
    for(samp in names(object@samples)){
      cat(samp, ": \n", sep = "")
      layers <- names(unlist(object@samples[[samp]]@layer))
      cat("  Layers:", paste(layers, collapse = " "), "\n")
    }

    # return invisible
    return(invisible(x = NULL))
  }
)

### subset of samples ####
setMethod(
  f = '[[',
  signature = c('SpaceRover', "character", "missing"),
  definition = function(x, i, j, ...){

    # sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){
      stop("There are no samples named ", i, " in this object")
    }

    return(x@samples[[i]])
  }
)

### subset of samples and layers ####
setMethod(
  f = '[[',
  signature = c('SpaceRover', "character", "character"),
  definition = function(x, i, j, ...){
    return(x[[i]]@layer[[j]])
  }
)

####
# Methods ####
####

### create SpaceRover object ####

#' CreateSpaceRover
#'
#' Create a SpaceRover object
#'
#' @param samples a list of samples list of class \code{srSamples}
#' @param metadata a metadata object of class \code{srMetadata}
#' @param sample.metadata a data frame of the sample metadata
#' @param zstack the zstack graph to determine the adjacency of spatial entities across layers
#' @param main.assay the name of the main assay of the object
#' @param project project name
#'
#' @export
#' @import igraph
#'
CreateSpaceRover <- function(samples, metadata = NULL, sample.metadata = NULL, zstack = NULL, main.assay = NULL, project = NULL){

  # set project name
  if(is.null(project))
    project <- "SpaceRover"

  # check for samples
  if(!is.list(samples))
    stop("Please introduce a list of samples")

  # set meta data
  if(is.null(metadata)){
    metadata <- setSRMetadata(cell = data.frame(), spot = data.frame(), ROI = data.frame())
  }

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setSRSampleMetadata(samples)
  }

  # set zgraph
  if(is.null(zstack)){
    spatial_entities <- Entities(metadata)
    zstack <- igraph::make_empty_graph(n = length(spatial_entities), directed = FALSE)
    igraph::V(zstack)$name <- spatial_entities
  }

  # set SpaceRover class
  new("SpaceRover", samples = samples, metadata = metadata, sample.metadata = sample.metadata, zstack = zstack, main.assay = main.assay, project = project)
}

### get main assay ####

#' @rdname MainAssay
#' @method MainAssay SpaceRover
#'
#' @export
#'
MainAssay.SpaceRover <- function(object, ...) {

  # get first sample and first layer with the main assay
  main.assay <- object@main.assay
  sample.info <- object@sample.metadata
  sample.info <- sample.info[sample.info$Assay == main.assay,, drop = FALSE]
  sample.info <- sample.info[1,] # GET FIRST FOR NOW!!!!

  # return assay
  return(object[[sample.info$Sample, sample.info$Layer]]@assay[[main.assay]])
}

### get entities ####

#' @rdname Entities
#' @method Entities SpaceRover
#'
#' @export
#'
Entities.SpaceRover <- function(object, ...) {
  return(Entities(object@metadata))
}
