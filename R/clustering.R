#' @include generics.R
#'
NULL

####
# Nearest Neighbor graphs ####
####

#' @rdname getNeighbors
#' @concept clustering
#' @method getNeighbors VoltRon
#'
#' @export
getNeighbors.VoltRon <- function(object, assay = NULL, data.type = "pca", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # list of plots
  gg <- list()

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    nnedges <- getNeighbors(cur_assay, data.type = data.type, ...)
    object@zstack <- add_edges(object@zstack, edges = nnedges)
  }

  # return
  return(object)
}

#' @rdname getNeighbors
#' @concept clustering
#' @method getNeighbors vrAssay
#'
#' @export
getNeighbors.vrAssay <- function(object, data.type = "pca", ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    nndata <- Data(object, norm = (data.type == "norm"))
    nndata <- t(nndata)
  } else {
    nndata <- Embeddings(object, type = data.type)
  }

  # find neighborhood
  nnedges <- FNN::get.knn(nndata)
  nnedges <- nnedges$nn.index
  nnedges <- cbind(1:length(Entities(object)), nnedges)
  nnedges <- apply(nnedges, 1, function(x){
    do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
  })
  nnedges <- Entities(object)[nnedges]

  # return
  return(nnedges)
}

####
# Nearest Neighbor graphs ####
####

getClusters <- function(object, resolution = 1, assay = NULL, label = "clusters"){

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
  metadata[[label]] <- metadata_clusters
  entities <- Entities(object_subset)
  metadata[entities,][[label]] <- clusters
  Metadata(object, type = AssayTypes(object)) <- metadata

  # return
  return(object)
}


