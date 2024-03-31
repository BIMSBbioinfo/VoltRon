#' @include generics.R
#'
NULL

####
# Nearest Neighbor graphs ####
####

#' Get profile specific neighbors
#'
#' Get neighbors of spatial points
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param data.type the type of embedding used for neighborhood calculation, e.g. raw counts (raw), normalized counts (norm), PCA embeddings (pca), UMAP embeddings (umap) etc.
#' @param dims the set of dimensions of the embedding data
#' @param k number of neighbors for kNN
#' @param method the method used for graph construction, SNN or kNN
#' @param graph.key the name of the graph
#' @param ... additional parameters passed to \link{FNN::get.knn}
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices E<- E
#' @importFrom FNN get.knn
#'
#' @export
#'
getProfileNeighbors <- function(object, assay = NULL, data.type = "pca", dims = 1:30, k = 10, method = "kNN", graph.key = method, ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    nndata <- vrData(object, assay = assay, norm = (data.type == "norm"))
    nndata <- t(nndata)
  } else {
    embedding_names <- vrEmbeddingNames(object)
    if(data.type %in% vrEmbeddingNames(object)) {
      nndata <- vrEmbeddings(object, assay = assay, type = data.type, dims = dims)
    } else {
      stop("Please provide a data type from one of three choices: raw, norm and pca")
    }
  }

  # find profile neighbors
  nnedges <- FNN::get.knn(nndata, k = k + 1)
  nnedges <-
    switch(method,
           SNN = {
             g.out <- build_snn_number(nnedges$nn.index)
             nnedges <- g.out[[1]]
             weights <- g.out[[2]]
             weights <- weights/(2 * (k+2) - weights)
             nnedges
           },
           kNN = {
             nnedges <- nnedges$nn.index
             nnedges <- cbind(1:nrow(nndata), nnedges)
             nnedges <- apply(nnedges, 1, function(x){
               do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
             })
             nnedges
           })
  nnedges <- rownames(nndata)[nnedges]

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(rownames(nndata))
  graph <- add_edges(graph, edges = nnedges)
  if(!is.null(weights))
    igraph::E(graph)$weight <- weights
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, graph.type = graph.key) <- graph

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
#' @param object a VoltRon object
#' @param resolution the resolution parameter for leiden clustering
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param label the name for the newly created clustering column in the metadata
#' @param graph the graph type to be used
#' @param seed seed
#'
#' @importFrom igraph cluster_leiden
#' @export
#'
getClusters <- function(object, resolution = 1, assay = NULL, label = "clusters", graph = "kNN", seed = 1){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assays
  object_subset <- subset(object, assays = assay_names)

  # graph
  object_graph <- vrGraph(object_subset, assay = assay, graph.type = graph)

  # clustering
  set.seed(seed)
  clusters <- igraph::cluster_leiden(object_graph, objective_function = "modularity", resolution_parameter = resolution)
  clusters <- clusters$membership

  # metadata
  metadata <- Metadata(object)
  entities <- vrSpatialPoints(object_subset)
  if(inherits(metadata, "data.table")){
    metadata[[label]] <- as.numeric(NA)
    metadata[[label]][metadata$id %in% entities] <- clusters
  } else {
    metadata_clusters <- NA
    metadata[[label]] <- metadata_clusters
    metadata[entities,][[label]] <- clusters
  }
  Metadata(object) <- metadata

  # return
  return(object)
}


