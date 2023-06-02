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
    nnedges <- findNeighbors(cur_assay, ...)
    object@zstack <- add_edges(object@zstack, edges = nnedges)
  }

  # return
  return(object)
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
  nnedges <- FNN::get.knn(t(normdata))
  nnedges <- nnedges$nn.index
  nnedges <- cbind(1:length(colnames(normdata)), nnedges)
  nnedges <- apply(nnedges, 1, function(x){
    do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
  })
  nnedges <- colnames(normdata)[nnedges]

  # return
  return(nnedges)
}

####
# Nearest Neighbor graphs ####
####

cluster <- function(object, resolution = 1, assay = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get assays
  object_subset <- subset(object, assays = assay_names)

  # graph
  object_graph <- Graph(object_subset)

  # clustering
  clusters <- cluster_leiden(object_graph, objective_function = "modularity", resolution_parameter = resolution)
  clusters <- clusters$membership

  # metadata
  metadata <- Metadata(object, type = AssayTypes(object))
  metadata_clusters <- rep(NA, nrow(metadata))
  metadata$clusters <- metadata_clusters
  entities <- Entities(object_subset)
  metadata[entities,]$clusters <- clusters
  Metadata(object, type = AssayTypes(object)) <- metadata

  # return
  return(object)
}


