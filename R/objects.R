# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))

####
# SpaceRover classes ####
####

## srMetadata ####

#' The srMetadata (SpaceRover Metadata) Class
#'
#' @slot samples A list of layers for the this srSample object
#'
#' @name srMetadata-class
#' @rdname srMetadata-class
#' @exportClass srMetadata
#'
srMetadata <- setClass(
  Class = 'srMetadata',
  slots = c(
    cell = 'data.frame',
    spot = 'data.frame',
    ROI = 'data.frame'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'srMetadata',
  definition = function(object) {
    cat("SpaceRover Metadata \n")
    cat("This object includes: \n")
    cat("  ", nrow(object@cell), "cells \n")
    cat("  ", nrow(object@spot), "spots \n")
    cat("  ", nrow(object@ROI), "ROIs \n")
    return(invisible(x = NULL))
  }
)

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

    # all_layers <- lapply(object@samples, function(x){
    #   return(unlist(x@layer))
    # })
    # all_layers <- unlist(all_layers)
    # cat("  ", length(object@samples), "sample(s) \n")
    # cat("  ", length(all_layers), "layer(s) \n")

    # return invisible
    return(invisible(x = NULL))
  }
)

### subset ####
setMethod(
  f = '[[',
  signature = 'SpaceRover',
  definition = function(x, i, j){

    # sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){
      stop("There are no samples named ", i, " in this object")
    }

    # return samples
    return(x@samples[[i]])
  }
)

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
    cat("This object includes", length(object@assay), "assays \n")
    return(invisible(x = NULL))
  }
)

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

## FOVImage ####

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

### show ####
setMethod(
  f = 'show',
  signature = 'FOVImage',
  definition = function(object) {
    cat("FOV Morphology Image \n")
    cat(paste(format(image_info(object@image)), collapse = " \n "))
    return(invisible(x = NULL))
  }
)

####
# SpaceRover Methods ####
####

CreateSpaceRover <- function(samples, metadata = NULL, sample.metadata = NULL, project = NULL){

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

  # set SpaceRover class
  new("SpaceRover", samples = samples, metadata = metadata, sample.metadata = sample.metadata, project = project)
}

setSRMetadata <- function(cell, spot, ROI){
  new("srMetadata", cell = cell, spot = spot, ROI = ROI)
}

setSRSampleMetadata <- function(samples){

  # imput missing sample names
  sample_name_ind <- sapply(names(samples), is.null)
  if(length(sample_name_ind) > 0){
    names_samples <- names(samples)
    if(any(sample_name_ind)){
      null_samples_ind <- which(sample_name_ind)
      names_samples[null_samples_ind] <- paste0("Sample", null_samples_ind)
    }
  } else {
    names_samples <- paste0("Sample", 1:length(samples))
  }

  # get sample metadata
  sample_list <- names(samples)
  sample.metadata <- NULL
  for(i in 1:length(sample_list)){
    layer_list <- samples[[sample_list[i]]]@layer
    layer_data <- NULL
    for(j in 1:length(layer_list)){
      assay_list <- layer_list[[j]]@assay
      layer_data <- rbind(layer_data, cbind(names(assay_list), names(layer_list)[j]))
    }
    sample.metadata <- rbind(sample.metadata, cbind(layer_data, sample_list[i]))
  }
  sample.metadata <- data.frame(sample.metadata, row.names = paste0("Assay", 1:nrow(sample.metadata)))
  colnames(sample.metadata) <- c("Assay", "Layer", "Sample")

  sample.metadata
}



