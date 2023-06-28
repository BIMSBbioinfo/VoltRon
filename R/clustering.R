#' @include generics.R
#'
NULL

####
# Nearest Neighbor graphs ####
####

#' @param assay assay
#' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#' @param dims the set of dimensions of the embedding data
#' @param k number of neighbors for kNN
#' @param method the method used for graph construction, SNN or kNN
#' @param ... additional parameters passed to \code{FNN:get.knn}
#'
#' @rdname getNeighbors
#' @method getNeighbors VoltRon
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom FNN get.knn
#' @importFrom proxy dist
#'
getNeighbors.VoltRon <- function(object, assay = NULL, data.type = "pca", dims = 1:30, k = 10, method = "SNN", ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    nndata <- vrData(object, assay = assay, norm = (data.type == "norm"))
    nndata <- t(nndata)
  } else {
    nndata <- vrEmbeddings(object, assay = assay, type = data.type, dims = dims)
  }

  # find neighborhood
  nnedges <- FNN::get.knn(nndata, k = k + 1)
  if(method == "SNN"){
    nnedges <- apply(nnedges$nn.index, 1, function(x) {
      result <- rep(0, nrow(nndata))
      result[x] <- 1
      return(result)
    }, simplify = FALSE)
    nnedges <- do.call(rbind, nnedges)
    jaccard_matrix <- jaccard_similarity(nnedges)
    jaccard_matrix[jaccard_matrix < 0.5] <- 0
    nnedges <- apply(jaccard_matrix, 1, function(x) which(x > 0))
  } else {
    nnedges <- nnedges$nn.index
  }

  # add edges to the graph
  nnedges <- cbind(1:length(vrSpatialPoints(object)), nnedges)
  nnedges <- apply(nnedges, 1, function(x){
    do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
  })
  nnedges <- vrSpatialPoints(object)[nnedges]

  # make graph and add edges
  object@graph <- make_empty_graph(directed = FALSE) + vertices(V(object@graph)$name)
  object@graph <- add_edges(object@graph, edges = nnedges)

  # remove multiple edges
  object@graph <- simplify(object@graph, remove.multiple = TRUE, remove.loops = FALSE)

  # return
  return(object)
}

#' #' @param assay assay
#' #' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#' #' @param dims the set of dimensions of the embedding data
#' #' @param ... additional parameters passed to \code{getNeighbors.vrAssay} and \code{FNN:get.knn}
#' #'
#' #' @rdname getNeighbors_old
#' #' @method getNeighbors VoltRon
#' #'
#' #' @export
#' #'
#' getNeighbors.VoltRon <- function(object, assay = NULL, data.type = "pca", dims = 1:30, ...){
#'
#'   # sample metadata
#'   sample.metadata <- SampleMetadata(object)
#'
#'   # get assay names
#'   assay_names <- vrAssayNames(object, assay = assay)
#'
#'   # normalize assays
#'   for(assy in assay_names){
#'     cur_assay <- object[[assy]]
#'     nnedges <- getNeighbors(cur_assay, data.type = data.type, dims = dims, ...)
#'     object@graph <- add_edges(object@graph, edges = nnedges)
#'   }
#'
#'   # return
#'   return(object)
#' }

#' #' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#' #' @param dims the set of dimensions of the embedding data
#' #' @param ... additional parameters passed to \code{FNN::get.knn}
#' #'
#' #' @rdname getNeighbors
#' #' @method getNeighbors vrAssay
#' #'
#' #' @export
#' #'
#' getNeighbors.vrAssay <- function(object, data.type = "pca", dims = 1:30, ...){
#'
#'   # get data
#'   if(data.type %in% c("raw", "norm")){
#'     nndata <- vrData(object, norm = (data.type == "norm"))
#'     nndata <- t(nndata)
#'   } else {
#'     nndata <- vrEmbeddings(object, type = data.type, dims = dims)
#'   }
#'
#'   # find neighborhood
#'   nnedges <- FNN::get.knn(nndata, ...)
#'   nnedges <- nnedges$nn.index
#'   nnedges <- cbind(1:length(vrSpatialPoints(object)), nnedges)
#'   nnedges <- apply(nnedges, 1, function(x){
#'     do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
#'   })
#'   nnedges <- vrSpatialPoints(object)[nnedges]
#'
#'   # return
#'   return(nnedges)
#' }

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
  metadata <- Metadata(object)
  metadata_clusters <- rep(NA, nrow(metadata))
  metadata[[label]] <- metadata_clusters
  entities <- vrSpatialPoints(object_subset)
  metadata[entities,][[label]] <- clusters
  Metadata(object) <- metadata

  # return
  return(object)
}


