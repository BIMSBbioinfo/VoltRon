#' @include generics.R
#'
NULL

####
# Nearest Neighbor graphs ####
####

#' @rdname findNeighbors
#' @concept clustering
#' @method findNeighbors SpaceRover
#'
#' @export
findNeighbors.SpaceRover <- function(object, assay = NULL, assay.type = NULL, ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    nngraph <- findNeighbors(cur_assay, ...)

  }
}

#' @rdname findNeighbors
#' @concept clustering
#' @method findNeighbors srAssay
#'
#' @export
findNeighbors.srAssay <- function(object, assay = NULL, assay.type = NULL, ...){

  # get data
  normdata <- object@normdata

  # find neighborhood
  FNN::get.knn(t(normdata))
}


