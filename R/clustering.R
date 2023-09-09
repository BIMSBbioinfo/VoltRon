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
#' @rdname getProfileNeighbors
#' @method getProfileNeighbors VoltRon
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom FNN get.knn
#'
#' @export
#'
getProfileNeighbors.VoltRon <- function(object, assay = NULL, data.type = "pca", dims = 1:30, k = 10, method = "kNN", ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    nndata <- vrData(object, assay = assay, norm = (data.type == "norm"))
    nndata <- t(nndata)
  } else {
    nndata <- vrEmbeddings(object, assay = assay, type = data.type, dims = dims)
  }

  # find neighborhood
  nnedges <- FNN::get.knn(nndata, k = k + 1)
  nnedges <-
    switch(method,
           SNN = {
             nnedges <- apply(nnedges$nn.index, 1, function(x) {
               result <- rep(0, nrow(nndata))
               result[x] <- 1
               return(result)
             }, simplify = FALSE)
             nnedges <- do.call(rbind, nnedges)
             jaccard_matrix <- jaccard_similarity(nnedges)
             jaccard_matrix[jaccard_matrix < 0.5] <- 0
             nnedges <- apply(jaccard_matrix, 1, function(x) which(x > 0))
             nnedges
           },
           kNN = {
             nnedges <- nnedges$nn.index
             nnedges
           })

  # add edges to the graph
  nnedges <- cbind(1:nrow(nndata), nnedges)
  nnedges <- apply(nnedges, 1, function(x){
    do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
  })
  nnedges <- rownames(nndata)[nnedges]

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(rownames(nndata))
  graph <- add_edges(graph, edges = nnedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, assay = assay, graph.type = "kNN") <- graph

  # return
  return(object)
}

####
# Clustering ####
####

#' getClusters
#'
#' Get clustering of the VoltRon object
#'
#' @param object A VoltRon object
#' @param resolution the resolution parameter for leiden clustering
#' @param assay assay
#' @param label the name for the newly created clustering column in the metadata
#' @param graph the graph type to be used
#'
#' @importFrom igraph cluster_leiden
#' @export
#'
getClusters <- function(object, resolution = 1, assay = NULL, label = "clusters", graph = "kNN"){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assays
  object_subset <- subset(object, assays = assay_names)

  # graph
  object_graph <- vrGraph(object_subset, assay = assay, graph.type = graph)

  # clustering
  clusters <- igraph::cluster_leiden(object_graph, objective_function = "modularity", resolution_parameter = resolution)
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


