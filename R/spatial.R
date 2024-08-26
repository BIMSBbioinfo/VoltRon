#' @include generics.R
#'
NULL

####
# Spatial Neighbor graphs ####
####

#' Get spatial neighbors
#'
#' get neighbors in an assay given spatial coordinates
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param method the method spatial connectivity: "delaunay", "spatialkNN", "radius"
#' @param k number of neighbors for kNN
#' @param radius When \code{method = "radius"} selected, determines the radius of a neighborhood ball around each spatial point
#' @param graph.key the name of the graph
#' @param ... additional parameters passed to \link{get.knn}
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom RCDT delaunay
#' @importFrom FNN get.knn
#' @importFrom RANN nn2
#' @importFrom reshape2 melt
#' @importFrom stats dist
#'
#' @export
getSpatialNeighbors <- function(object, assay = NULL, method = "delaunay", k = 10, radius = numeric(0), 
                                graph.key = method, ...){

  # get coordinates
  # coords <- vrCoordinates(object, assay = assay)
  spatialpoints <- vrSpatialPoints(object, assay = assay)

  # get spatial graph per assay
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get assay connectivity 
  assay_names <- getBlockConnectivity(object, assay = assay_names)
  
  # get spatial edges
  spatialedges_list <- list()
  for(assy in assay_names){
    # cur_coords <- coords[grepl(paste0(assy, "$"), rownames(coords)), ]
    cur_coords <- vrCoordinates(object, assay = assy)
    spatialedges <-
      switch(method,
             delaunay = {
               nnedges <- RCDT::delaunay(cur_coords)
               nnedges <- as.vector(t(nnedges$edges[,1:2]))
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             spatialkNN = {
               # nnedges <- FNN::get.knn(cur_coords, k = 1)
               nnedges <- RANN::nn2(cur_coords, k = k + 1)
               nnedges <- nnedges$nn.idx
               # nnedges <- cbind(1:nrow(cur_coords), nnedges)
               # nnedges <- apply(nnedges, 1, function(x){
               #   do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
               # })
               # nnedges <- unlist(nnedges)
               nnedges <- reshape2::melt(data.frame(nnedges), id.vars = "X1")
               nnedges <- subset(nnedges[,c("X1", "value")], value != 0 & X1 != 0)
               nnedges <- as.vector(t(as.matrix(nnedges)))
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             radius = {
               if(length(radius) == 0){
                 spot.radius <- vrAssayParams(object[[assy]], param = "nearestpost.distance")
                 radius <- ifelse(is.null(spot.radius), 1, spot.radius)
               }
               # distances <- as.matrix(stats::dist(cur_coords, method = "euclidean"))
               # nnedges <- apply(distances, 1, function(x){
               #   which(x < radius & x > .Machine$double.eps)
               # })
               nnedges <- RANN::nn2(cur_coords, searchtype = "radius", radius = radius, k = min(300, sqrt(nrow(cur_coords))/2))
               nnedges <- nnedges$nn.idx
               nnedges <- reshape2::melt(data.frame(nnedges), id.vars = "X1")
               nnedges <- subset(nnedges[,c("X1", "value")], value != 0 & X1 != 0)
               nnedges <- as.vector(t(as.matrix(nnedges)))
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             })
    spatialedges_list <- c(spatialedges_list, list(spatialedges))
  }
  spatialedges <- unlist(spatialedges_list)

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(spatialpoints)
  graph <- add_edges(graph, edges = spatialedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, graph.type = graph.key) <- graph

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
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
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
    # object_subset <- subset(object, spatialpoints = vrSpatialPoints(object)[grepl(paste0(assy,"$"), vrSpatialPoints(object))])
    object_subset <- subset(object, assays = assy)
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
#' @importFrom dplyr group_by bind_rows filter summarize mutate n
#' @importFrom igraph neighborhood
#'
#' @noRd
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

  neigh_results <- neighbors_graph_data %>%
    dplyr::group_by(from_value, to_value, type) %>%
    dplyr::summarize(mean_value = dplyr::n()) %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::mutate(assoc_test = mean_value > ifelse("obs" %in% type, mean_value[type == "obs"], 0),
                  segreg_test = mean_value < ifelse("obs" %in% type, mean_value[type == "obs"], 0)) %>%
    dplyr::mutate(majortype = ifelse(type == "obs", "obs", "sim")) %>% 
    dplyr::group_by(from_value, to_value) %>%
    dplyr::mutate(value = ifelse(sum(majortype == "obs") > 0, log(mean_value[majortype == "obs"]/mean(mean_value[majortype == "sim"])), 0)) %>% 
    dplyr::filter(type != "obs") %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::summarize(p_assoc = mean(assoc_test), p_segreg = mean(segreg_test), value = value[1]) %>%
    dplyr::mutate(p_assoc_adj = p.adjust(p_assoc, method = "fdr"),
                  p_segreg_adj = p.adjust(p_segreg, method = "fdr"))
  
  # number of samples
  grp_table <- table(grp)
  neigh_results$n_from <- grp_table[neigh_results$from_value]
  neigh_results$n_to <- grp_table[neigh_results$to_value]

  # return
  neigh_results
}
