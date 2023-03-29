####
# SpaceRover classes ####
####

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
    integrated.datasets = "list",
    meta.data = 'data.frame',
    project = 'character'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'SpaceRover',
  definition = function(object) {

    # print samples
    cat(class(x = object), "Class \n")
    cat("This object includes", length(object@samples), "sample(s). \n")

    # layers
    all_layers <- lapply(object@samples, function(x){
      return(unlist(x@layers))
    })
    all_layers <- unlist(all_layers)
    cat("This object includes", length(all_layers), "layer(s) in total. \n")

    # return invisible
    return(invisible(x = NULL))
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
    layers = 'list'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'srSample',
  definition = function(object) {
    cat(class(x = object), "(SpaceRover Sample) Class \n")
    cat("This object includes", length(object@layers), "layers. \n")
    return(invisible(x = NULL))
  }
)

####
# Auxiliary methods ####
####

makeSpaceRover <- function(samples, meta.data = NULL, project = NULL){

  # set project name
  if(is.null(project))
    project <- "SpaceRover"

  # check for samples
  if(!is.list(samples))
    stop("Please introduce a list of samples")

  # set meta data
  if(is.null(meta.data)){
    if(any(lapply(names(samples), is.null))){
      names_samples <- paste0("Sample", 1:length(samples))
    } else {
      names_samples <- names(samples)
    }
    meta.data <- data.frame(Sample = names_samples)
  }

  # set integrated datasets
  integrated.datasets <- list()

  # set SpaceRover class
  new("SpaceRover", samples = samples, integrated.datasets = integrated.datasets, meta.data = meta.data, project = project)
}





