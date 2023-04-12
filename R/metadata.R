####
# Objects and Classes ####
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
    cat("SpaceRover Metadata Object \n")
    cat("This object includes: \n")
    if(nrow(object@cell) > 0)
      cat("  ", nrow(object@cell), "cells \n")
    if(nrow(object@spot) > 0)
      cat("  ", nrow(object@spot), "spots \n")
    if(nrow(object@ROI) > 0)
      cat("  ", nrow(object@ROI), "ROIs \n")
    return(invisible(x = NULL))
  }
)

####
# Methods ####
####

#' @rdname Entities
#' @method Entities srMetadata
#'
#' @export
#'
Entities.srMetadata <- function(object, ...) {

  # get the combination of cells, spots and ROIs
  points <- c(rownames(object@cell),
                rownames(object@spot),
                rownames(object@ROI))

  return(points)
}

####
# Functions ####
####

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
