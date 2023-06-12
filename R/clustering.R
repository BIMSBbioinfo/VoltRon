#' @include generics.R
#'
NULL

####
# Nearest Neighbor graphs ####
####

#' @param assay assay
#' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#' @param ... additional parameters passed to \code{getNeighbors.vrAssay}
#'
#' @rdname getNeighbors
#' @method getNeighbors VoltRon
#'
#' @export
#'
getNeighbors.VoltRon <- function(object, assay = NULL, data.type = "pca", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    nnedges <- getNeighbors(cur_assay, data.type = data.type, ...)
    object@zstack <- add_edges(object@zstack, edges = nnedges)
  }

  # return
  return(object)
}

#' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#'
#' @rdname getNeighbors
#' @method getNeighbors vrAssay
#'
#' @export
#'
getNeighbors.vrAssay <- function(object, data.type = "pca", ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    nndata <- vrData(object, norm = (data.type == "norm"))
    nndata <- t(nndata)
  } else {
    nndata <- vrEmbeddings(object, type = data.type)
  }

  # find neighborhood
  nnedges <- FNN::get.knn(nndata)
  nnedges <- nnedges$nn.index
  nnedges <- cbind(1:length(vrSpatialPoints(object)), nnedges)
  nnedges <- apply(nnedges, 1, function(x){
    do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
  })
  nnedges <- vrSpatialPoints(object)[nnedges]

  # return
  return(nnedges)
}

####
# Nearest Neighbor graphs ####
####

#' getClusters
#'
#' Get clustering of the VoltRon object
#'
#' @param object A VoltRon object
#' @param resolution the resolution parameter for leiden clustering
#' @param assay assay
#' @param label the name for the newly created clustering column in the metadata
#'
#' @export
#'
getClusters <- function(object, resolution = 1, assay = NULL, label = "clusters"){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assays
  object_subset <- subset(object, assays = assay_names)

  # graph
  object_graph <- vrGraph(object_subset)

  # clustering
  clusters <- cluster_leiden(object_graph, objective_function = "modularity", resolution_parameter = resolution)
  clusters <- clusters$membership

  # metadata
  metadata <- Metadata(object, type = vrAssayTypes(object))
  metadata_clusters <- rep(NA, nrow(metadata))
  metadata[[label]] <- metadata_clusters
  entities <- vrSpatialPoints(object_subset)
  metadata[entities,][[label]] <- clusters
  Metadata(object, type = vrAssayTypes(object)) <- metadata

  # return
  return(object)
}


