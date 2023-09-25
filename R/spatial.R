#' @include generics.R
#'
NULL

####
# Spatial Neighbor graphs ####
####

#' @param assay assay
#' @param method the method spatial connectivity
#' @param ... additional parameters passed to \code{FNN:get.knn}
#'
#' @rdname getSpatialNeighbors
#' @method getSpatialNeighbors VoltRon
#'
#' @importFrom tripack tri.mesh neighbours
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom FNN get.knn
#'
#' @export
#'
getSpatialNeighbors.VoltRon <- function(object, assay = NULL, method = "delaunay", ...){

  # get coordinates
  coords <- vrCoordinates(object, assay = assay, reg = TRUE)

  # get spatial graph per assay
  assay_names <- vrAssayNames(object, assay = assay)
  spatialedges_list <- list()
  for(assy in assay_names){
    cur_coords <- coords[grepl(paste0(assy, "$"), rownames(coords)), ]
    spatialedges <-
      switch(method,
             delaunay = {
               tess <- tripack::tri.mesh(cur_coords[,1], cur_coords[,2])
               nnedges <- tripack::neighbours(tess)
               nnedges <- mapply(function(x,y) {
                 do.call(c,lapply(y, function(z) c(x,z)))
               }, 1:length(nnedges), nnedges)
               nnedges <- unlist(nnedges)
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             NN = {
               nnedges <- FNN::get.knn(cur_coords, k = 1)
               nnedges <- cbind(1:nrow(cur_coords), nnedges)
               nnedges <- apply(nnedges, 1, function(x){
                 do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
               })
               nnedges <- unlist(nnedges)
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             })
    spatialedges_list[[assy]] <- spatialedges
  }
  spatialedges <- unlist(spatialedges_list)

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(rownames(coords))
  graph <- add_edges(graph, edges = spatialedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, assay = assay, graph.type = "delaunay") <- graph

  # return
  return(object)
}

####
# Neighbor Enrichment test ####
####

#' vrNeighbourhoodEnrichment
#'
#' Conduct Neighborhood enrichment test for pairs of clusters for all assays
#'
#' @param object a VoltRon object
#' @param assay assay
#' @param group.by a column from metadata to seperate spatial points
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param num.sim the number of simulations
#' @param seed seed
#'
#' @export
#'
vrNeighbourhoodEnrichment <- function(object, assay = NULL, group.by = NULL, graph.type = "delaunay", num.sim = 1000, seed = 1){

  # set the seed
  set.seed(seed)

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # test for each assay
  neigh_results <- list()
  for(assy in assay_names){
    message("Testing Neighborhood Enrichment of '", group.by ,"' for '", assy, "'")
    object_subset <- subset(object, spatialpoints = vrSpatialPoints(object)[grepl(paste0(assy,"$"), vrSpatialPoints(object))])
    neigh_results[[assy]] <- vrNeighbourhoodEnrichmentSingle(object_subset, group.by = group.by, graph.type = graph.type,
                                                             num.sim = num.sim, seed = seed)
    neigh_results[[assy]] <- data.frame(neigh_results[[assy]], AssayID = assy, SampleMetadata(object_subset))
  }
  all_neigh_results <- do.call(rbind, neigh_results)

  # return
  return(all_neigh_results)
}

#' vrNeighbourhoodEnrichmentSingle
#'
#' Conduct Neighborhood enrichment test for pairs of clusters for a vrAssay
#'
#' @param object a VoltRon object
#' @param group.by a column from metadata to seperate spatial points
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param num.sim the number of simulations
#' @param seed seed
#'
#' @importFrom dplyr group_by bind_rows filter summarize mutate
#' @importFrom igraph neighborhood
#'
vrNeighbourhoodEnrichmentSingle <- function(object, group.by = NULL, graph.type = "delaunay", num.sim = 1000, seed = 1) {

  # set the seed
  set.seed(seed)

  # main object
  metadata <- Metadata(object)
  grp <- metadata[[group.by]]
  names(grp) <- rownames(metadata)

  # get graph and neighborhood
  graph <- vrGraph(object, graph.type = graph.type)
  neighbors_graph <- igraph::neighborhood(graph)
  neighbors_graph_data <- lapply(neighbors_graph, function(x) {
    cbind(x$name[1],x$name[-1])
  })
  neighbors_graph_data <- do.call(rbind, neighbors_graph_data)
  colnames(neighbors_graph_data) <- c("from", "to")

  # get simulations
  grp_sim <- sapply(1:1000, function(x) sample(grp))
  rownames(grp_sim) <- names(grp)

  # get adjacency for observed and simulated pairs
  neighbors_graph_data_list <- list(data.frame(neighbors_graph_data, from_value = grp[neighbors_graph_data[,1]], to_value = grp[neighbors_graph_data[,2]], type = "obs"))
  for(i in 2:(ncol(grp_sim)+1))
    neighbors_graph_data_list[[i]] <- data.frame(neighbors_graph_data, from_value = grp_sim[,i-1][neighbors_graph_data[,1]], to_value = grp_sim[,i-1][neighbors_graph_data[,2]], type = paste0("sim", i))
  neighbors_graph_data <- dplyr::bind_rows(neighbors_graph_data_list)

  # conduct randomized test
  `%notin%` <- Negate(`%in%`)
  neigh_results <- neighbors_graph_data %>%
    dplyr::group_by(from_value, to_value, type) %>%
    dplyr::summarize(mean_value = n()) %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::mutate(assoc_test = mean_value > ifelse("obs" %in% type, mean_value[type == "obs"], 0),
                  segreg_test = mean_value < ifelse("obs" %in% type, mean_value[type == "obs"], 0)) %>%
    dplyr::filter(type != "obs") %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::summarize(p_assoc = mean(assoc_test), p_segreg = mean(segreg_test)) %>%
    dplyr::mutate(p_assoc_adj = p.adjust(p_assoc, method = "fdr"),
                  p_segreg_adj = p.adjust(p_segreg, method = "fdr"))

  # return
  neigh_results
}
