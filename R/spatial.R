#' @include allgenerics.R
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
#' @param group.by a column of metadata from \link{Metadata} used as grouping label for the spatial entities.
#' @param group.ids a subset of categories defined in metadata column from \code{group.by}.
#' @param method the method spatial connectivity: "delaunay", "spatialkNN", "radius".
#' @param k number of neighbors for kNN.
#' @param radius When \code{method = "radius"} selected, determines the radius of a neighborhood ball around each spatial point.
#' @param graph.key the name of the graph.
#' @param verbose verbose
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom RCDT delaunay
#' @importFrom RANN nn2
#' @importFrom data.table data.table melt
#' @importFrom stats dist
#'
#' @export
getSpatialNeighbors <- function(object, 
                                assay = NULL, 
                                group.by = NULL,
                                group.ids = NULL,
                                method = "delaunay", 
                                k = 10, 
                                radius = numeric(0), 
                                graph.key = method, 
                                verbose = TRUE){

  # get coordinates
  spatialpoints <- vrSpatialPoints(object, assay = assay)

  # get spatial graph per assay
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get assay connectivity 
  assay_names_connected <- getBlockConnectivity(object, assay = assay_names)
  
  # get spatial edges
  spatialedges_list <- list()
  for(assy in assay_names_connected){
  
    coords <- NULL
    for(a in assy){
      
      # get coordinates and metadata
      cur_coords <- as.matrix(vrCoordinates(object, assay = a))
      metadata <- Metadata(object, assay = a)
      
      if(!is.null(group.by) && !is.null(group.ids)){
      
        # get groups
        if(is.list(group.by)){
          cur_group.by <- group.by[[a]]
          cur_group.ids <- group.ids[[a]]
        }
        
        if(verbose)
          message("Calculating Spatial Neighbors with group.by='", cur_group.by, "' and group.ids='", paste(cur_group.ids, collapse = ","), "'")
      
        if(!cur_group.by %in% colnames(metadata))
          stop("The column '", cur_group.by, "' was not found in the metadata!")
        if(inherits(metadata, "data.table")){
          column <- metadata[,get(names(metadata)[which(colnames(metadata) == cur_group.by)])]
          names(column) <- metadata$id
        } else {
          column <- metadata[,cur_group.by]
          if("id" %in% colnames(metadata)){
            names(column) <- as.vector(metadata$id)
          } else {
            names(column) <- rownames(metadata)
          }
        }
        if(a %in% names(group.ids)){
          len_set_diff <- length(setdiff(cur_group.ids,  column))
          if(len_set_diff > 0){
          } else if(len_set_diff == length(cur_group.ids)){ 
            stop("None of the groups defined in group.ids exist in cur_group.by!")
          } 
          column <- column[column %in% cur_group.ids] 
        }
        cur_coords <- cur_coords[names(column),]
        
      } else if(is.null(group.by) && is.null(group.ids)) {
        
      } else {
        stop("Either both 'group.by' and 'group.ids' should be specified or both should be null")
      }
      
      # merge coords
      coords <- rbind(coords, cur_coords)
    }
    
    # get edges
    spatialedges <-
      switch(method,
             delaunay = {
               nnedges <- RCDT::delaunay(coords)
               nnedges <- as.vector(t(nnedges$edges[,seq_len(2)]))
               nnedges <- rownames(coords)[nnedges]
               nnedges
             },
             spatialkNN = {
               # nnedges <- RANN::nn2(coords, k = k + 1)
               nnedges <- knn_annoy(coords, k = k + 1)
               names(nnedges) <- c("nn.index", "nn.dist")
               nnedges <- data.table::melt(data.table::data.table(nnedges$nn.index), id.vars = "V1")
               nnedges <- nnedges[,c("V1", "value")][V1 > 0 & value > 0]
               nnedges <- as.vector(t(as.matrix(nnedges)))
               nnedges <- rownames(coords)[nnedges]
               nnedges
             },
             radius = {
               if(length(radius) == 0){
                 spot.radius <- vrAssayParams(object[[assy]], param = "nearestpost.distance")
                 radius <- ifelse(is.null(spot.radius), 1, spot.radius)
               }
               nnedges <- suppressWarnings({RANN::nn2(coords, searchtype = "radius", radius = radius, k = min(300, sqrt(nrow(cur_coords))/2))})
               nnedges <- data.table::melt(data.table::data.table(nnedges$nn.idx), id.vars = "V1")
               nnedges <- nnedges[,c("V1", "value")][V1 > 0 & value > 0]
               nnedges <- as.vector(t(as.matrix(nnedges)))
               nnedges <- rownames(coords)[nnedges]
               nnedges
             })
    spatialedges_list <- c(spatialedges_list, list(spatialedges))
  }
  spatialedges <- unlist(spatialedges_list)

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + rownames(coords)
  graph <- add_edges(graph, edges = spatialedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, assay = assay_names, graph.type = graph.key) <- graph

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
#' @param group.by a column from metadata to separate spatial points
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param num.sim the number of simulations
#' @param seed seed
#' @param verbose verbose
#'
#' @export
vrNeighbourhoodEnrichment <- function(object, assay = NULL, group.by = NULL, graph.type = "delaunay", num.sim = 1000, seed = 1, verbose = TRUE){

  # set the seed
  set.seed(seed)

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # check groups
  if(is.list(group.by)){
    if(is.null(names(group.by)))
      stop("If 'group.by' is a list, then names should be Assay ids, e.g. Assay1, Assay4 etc.")
    if(!all(names(group.by) %in% assay_names))
      stop("Names of the 'group.by' list should of assays")
  }
    
  # get assay connectivity 
  assay_names_connected <- getBlockConnectivity(object, assay = assay_names)

  # test for each assay
  neigh_results <- list()
  for(assy in assay_names_connected){
    cur_group.by <- group.by
    if(is.list(group.by))
      cur_group.by <- group.by[assy]
    if(verbose)
      message("Testing Neighborhood Enrichment of '", 
              paste(cur_group.by, collapse = ", ") ,"' for '", 
              paste(assy, collapse = ", "), "'")
    object_subset <- subsetVoltRon(object, assays = assy)
    res <- vrNeighbourhoodEnrichmentSingle2(object_subset, 
                                           group.by = cur_group.by, 
                                           graph.type = graph.type,
                                           num.sim = num.sim, 
                                           seed = seed)
    res <- data.frame(res, AssayID = paste(assy, collapse = ", "), unique(object_subset$Sample))
    neigh_results <- c(neigh_results, list(res))
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
  
  # check groups
  if(!is.list(group.by)){
    group.by <- list(group.by)
    names(group.by) <- vrAssayNames(object)
  }

  # get groups
  grp <- NULL
  for(assy in names(group.by)){
    metadata <- Metadata(object, assay = assy)
    if(group.by[[assy]] %in% colnames(metadata)){
      cur_grp <- metadata[[group.by[[assy]]]]
      if("id" %in% colnames(metadata)){
        names(cur_grp) <- as.vector(metadata$id)
      } else {
        names(cur_grp) <- rownames(metadata) 
      }
    } else {
      stop("'", group.by[[assy]], "' is not available in metadata!")
    } 
    grp <- c(grp, cur_grp)
  }
  
  # get graph and neighborhood
  message("Getting Neighborhood graph ...")
  graphx <- vrGraph(object, graph.type = graph.type, assay = names(group.by))
  vertex_names <- igraph::V(graphx)$name
  neighbors_graph <- igraph::neighborhood(graphx, nodes = vertex_names)
  # graph <- vrGraph(object, graph.type = graph.type, assay = names(group.by))
  # neighbors_graph <- igraph::neighborhood(graph)
  # neighbors_graph_data <- lapply(neighbors_graph, function(x) {
  #   cbind(x$name[1],x$name)[-1,]
  # })
  # neighbors_graph_data <- do.call(rbind, neighbors_graph_data)
  # colnames(neighbors_graph_data) <- c("from", "to")
  
  # get neighborhood data
  node_names <- rep(vertex_names, rep(lapply(neighbors_graph, length)))
  target_names <- unlist(neighbors_graph)
  neighbors_graph_data <- cbind(node_names, names(target_names))
  colnames(neighbors_graph_data) <- c("from", "to")
  
  # get simulations
  message("Permute labels ...")
  grp <- grp[names(V(graphx))]
  # grp_sim <- vapply(seq_len(num.sim), function(x) sample(grp), grp)
  # rownames(grp_sim) <- names(grp)
  # 
  # # get adjacency for observed and simulated pairs
  # neighbors_graph_data_list <- list(data.frame(neighbors_graph_data, from_value = grp[neighbors_graph_data[,1]], to_value = grp[neighbors_graph_data[,2]], type = "obs"))
  # for(i in 2:(ncol(grp_sim)+1))
  #   neighbors_graph_data_list[[i]] <- data.frame(neighbors_graph_data, from_value = grp_sim[,i-1][neighbors_graph_data[,1]], to_value = grp_sim[,i-1][neighbors_graph_data[,2]], type = paste0("sim", i))
  # neighbors_graph_data <- dplyr::bind_rows(neighbors_graph_data_list)
  grp_sim <- replicate(num.sim, sample(grp), simplify = "matrix")
  rownames(grp_sim) <- names(grp)
  
  # get adjacency for observed and simulated pairs
  message("Simulate Edges ...")
  from_ids <- neighbors_graph_data[, 1]
  to_ids <- neighbors_graph_data[, 2]
  neighbors_graph_data_list <- vector("list", num.sim + 1)
  neighbors_graph_data_list[[1]] <- data.frame(
    neighbors_graph_data,
    from_value = grp[from_ids],
    to_value = grp[to_ids],
    type = "obs"
  )
  for(i in seq_len(num.sim)) {
    print(i)
    sim_grp <- grp_sim[, i]
    neighbors_graph_data_list[[i + 1]] <- data.frame(
      neighbors_graph_data,
      from_value = sim_grp[from_ids],
      to_value = sim_grp[to_ids],
    )
  }
  neighbors_graph_data <- dplyr::bind_rows(neighbors_graph_data_list)
  
  # get adjacency for observed and simulated pairs
  message("Calculate Statistics ...")
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
vrNeighbourhoodEnrichmentSingle2 <- function(object, group.by = NULL, graph.type = "delaunay", num.sim = 1000, seed = 1) {
  
  # set the seed
  set.seed(seed)
  
  # check groups
  if(!is.list(group.by)){
    group.by <- list(group.by)
    names(group.by) <- vrAssayNames(object)
  }
  
  # get groups
  grp <- NULL
  for(assy in names(group.by)){
    metadata <- Metadata(object, assay = assy)
    if(group.by[[assy]] %in% colnames(metadata)){
      cur_grp <- metadata[[group.by[[assy]]]]
      if("id" %in% colnames(metadata)){
        names(cur_grp) <- as.vector(metadata$id)
      } else {
        names(cur_grp) <- rownames(metadata) 
      }
    } else {
      stop("'", group.by[[assy]], "' is not available in metadata!")
    } 
    grp <- c(grp, cur_grp)
  }
  # spatialpoints <- names(grp)
  
  # get graph and neighborhood
  message("Getting neighborhood graph")
  graphx <- vrGraph(object, graph.type = graph.type, assay = names(group.by))
  vertex_names <- igraph::V(graphx)$name
  neighbors_graph <- igraph::neighborhood(graphx, nodes = vertex_names)
  spatialpoints <- vertex_names
  
  # get neighborhood data
  node_names <- rep(vertex_names, rep(lapply(neighbors_graph, length)))
  target_names <- unlist(neighbors_graph)
  neighbors_graph_data <- data.table::data.table(node_names, names(target_names))
  neighbors_graph_data <- setNames(neighbors_graph_data, c("from", "to"))

  # get simulations
  message("Permuting labels")
  grp <- grp[names(V(graphx))]
  grp_sim <- replicate(num.sim, sample(grp), simplify = "matrix")
  grp_sim <- cbind(grp, grp_sim)
  grp_sim <- data.table::as.data.table(grp_sim)

  # get adjacency for observed and simulated pairs
  message("Simulating edges")
  from_ids <- match(neighbors_graph_data$from, vertex_names)
  to_ids   <- match(neighbors_graph_data$to, vertex_names)
  # from_ids <- neighbors_graph_data$from
  # to_ids <- neighbors_graph_data$to
  message("art1")
  # idx_lookup <- setNames(seq_len(nrow(grp_sim)), rownames(grp_sim))
  # from_ids <- idx_lookup[neighbors_graph_data$from]
  # to_ids   <- idx_lookup[neighbors_graph_data$to]
  from_labels <- grp_sim[as.integer(from_ids), ]
  colnames(from_labels) <- c("obs", paste0("sim",1:(ncol(grp_sim)-1)))
  to_labels   <- grp_sim[as.integer(to_ids), ]
  colnames(to_labels) <- c("obs", paste0("sim",1:(ncol(grp_sim)-1)))
  message("art3")
  
  # statistics
  message("Calculating statistics")
  myfun <- function(x,y){
    print(head(x))
    print(head(y))
    as.data.frame(table(x, y))
  }
  # pair_counts <- data.table::copy(from_labels)
  tmp <- Map(myfun, from_labels, to_labels)
  # pair_counts[, names(from_labels) := Map(myfun, from_labels, to_labels)]
  
    
  # to_labels <- data.table::melt(data.table::data.table(to_labels),
  #                               measure.vars  = colnames(to_labels))
  # to_labels <- data.table::melt(to_labels,
  #                               measure.vars  = colnames(to_labels))
  # names(to_labels) <- c("type", "to")
  # from_labels <- data.table::melt(from_labels,
  #                               measure.vars  = colnames(from_labels))
  # names(from_labels) <- c("type", "from")
  # message("art1")
  # dt <- data.table(from_labels[,2], 
  #                  to_labels)
  # message("art2")
  # pair_counts <- dt[, .N, by = .(from, to, type)]
  # message("art3")
  # pair_counts[, `:=`(
  #   assoc_test = N > if ("obs" %in% type) N[type == "obs"] else 0,
  #   segreg_test = N < if ("obs" %in% type) N[type == "obs"] else 0,
  #   majortype = ifelse(type == "obs", "obs", "sim")
  # ), by = .(from, to)]
  # message("art4")
  # pair_counts[, value := {
  #   if (sum(majortype == "obs") > 0) {
  #     obs_val <- N[majortype == "obs"]
  #     sim_mean <- mean(N[majortype == "sim"])
  #     if (sim_mean > 0) log(obs_val / sim_mean) else 0
  #   } else {
  #     0
  #   }
  # }, by = .(from, to)]
  # pair_counts <- pair_counts[type != "obs"]
  # pair_counts <- pair_counts[, .(
  #   p_assoc = mean(assoc_test),
  #   p_segreg = mean(segreg_test),
  #   value = value[1]
  # ), by = .(from, to)]
  # pair_counts[, p_assoc_adj := p.adjust(p_assoc, method = "fdr")]
  # pair_counts[, p_segreg_adj := p.adjust(p_segreg, method = "fdr")]
  
  # return
  pair_counts
}

####
# Hot Spot Analysis ####
####

#' getHotSpotAnalysis
#'
#' Conduct hot spot detection
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param method the statistical method of conducting hot spot analysis. Default is "Getis-Ord"
#' @param features a set of features to be visualized, either from \link{vrFeatures} of raw or normalized data or columns of the \link{Metadata}.
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param alpha.value the alpha value for the hot spot analysis test. Default is 0.01
#' @param norm if TRUE, the normalized data is used
#' @param verbose verbose
#' 
#' @importFrom Matrix rowSums
#' @importFrom igraph as_adjacency_matrix
#' @importFrom stats pnorm
#'
#' @export
getHotSpotAnalysis <- function(object, assay = NULL, method = "Getis-Ord", features, graph.type = "delaunay", alpha.value = 0.01, norm = TRUE, verbose = TRUE){
  
  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")
  
  # metadata
  sample.metadata <- SampleMetadata(object)
  metadata <- Metadata(object, assay = assay)
  
  # initiate metadata columns
  for(feat in features){
    metadata[[paste0(feat,"_hotspot_flag")]] <- 
      metadata[[paste0(feat,"_hotspot_pvalue")]] <- 
        metadata[[paste0(feat,"_hotspot_stat")]] <- rep(NA, nrow(metadata))
  }
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # test for each assay
  neigh_results <- list()
  for(assy in assay_names){
    
    # verbose
    if(verbose)
      message("Running Hot Spot Analysis with '", method, "' for '", assy, "'")
    
    # get related data 
    graph <- vrGraph(object, assay = assy, graph.type = graph.type)
    adj_matrix <- igraph::as_adjacency_matrix(graph)
    cur_metadata <- subset_metadata(metadata, assays = assy)
    
    # get features
    data_features <- features[features %in% vrFeatures(object, assay = assy)]
    if(length(data_features) > 0){
      normdata <- vrData(object, assay = assy, features = data_features, norm = norm)
    }
    
    # for each feature
    for(feat in features){
      
      # Getis-Ord 
      if(method == "Getis-Ord"){

        # get feature
        if(feat %in% data_features){
          if(inherits(normdata, "IterableMatrix")){
            statistic <- as.matrix(normdata[feat,])[1,]
          } else {
            statistic <- normdata[feat,]
          }
        } else if(feat %in% colnames(cur_metadata)){
          if(inherits(cur_metadata, "data.table")){
            statistic <- cur_metadata[,get(names(cur_metadata)[which(colnames(cur_metadata) == feat)])]
          } else {
            statistic <- cur_metadata[,feat]
          }
        } else {
          stop("'", feat, "' is not found in either the data matrix or the metadata!")
        }
        
        # update statistics if not numeric
        if(!is.numeric(statistic)){
          statistic <- Matrix::rowSums(adj_matrix)
        }
        
        # initiate getis ord
        getisord <- list()
        length(getisord) <- 3
        names(getisord) <- c("hotspot_stat", "hotspot_pvalue", "hotspot_flag")
        
        # calculate getis ord
        n <- length(statistic)
        statistic <- statistic - (min(statistic)) ### correct for negative scores, correct for this later
        getisord_stat <- adj_matrix %*% statistic
        getisord_stat <- getisord_stat/(sum(statistic) - statistic)
        getisord[[1]] <- getisord_stat[,1]
          
        # calculate z score expectation and variance
        weight_sum <- Matrix::rowSums(adj_matrix)
        getisord_exp <- weight_sum/(n - 1)
        getisord_moment_1 <- (sum(statistic) - statistic)/(n - 1)
        getisord_moment_2 <- (sum(statistic^2) - statistic^2)/(n - 1) - getisord_moment_1^2
        getisord_var <- weight_sum*(n - 1 - weight_sum)*getisord_moment_2
        getisord_var <- getisord_var/((n-1)^2 *(n-2)*getisord_moment_1^2)
        
        # calculate z score 
        getisord_zscore <- (getisord[[1]] - getisord_exp)/sqrt(getisord_var)
        getisord_zscore[is.nan(getisord_zscore)] <- NA
        # getisord[[2]] <- 1-stats::pnorm(getisord_zscore)
        getisord[[2]] <- p.adjust(1-stats::pnorm(getisord_zscore), method = "bonferroni")
        getisord[[3]] <- ifelse(getisord[[2]] < alpha.value, "hot", "cold") 
        
        # get graph based hot spot filtering
        # at least some number of common neighbors should be hotspots
        
        # update metadata for features
        for(label in names(getisord))
          cur_metadata[[paste(feat, label, sep = "_")]] <-  getisord[[label]]
      }
    }
    
    # update metadata for assays
    if("id" %in% colnames(metadata)){
      ind <- match(as.vector(cur_metadata$id), as.vector(metadata$id))
      for(feat in features){
        for(label in names(getisord)){
          metadata_label <- paste(feat, label, sep = "_")
          metadata[[metadata_label]][ind] <- cur_metadata[[metadata_label]]
          object <- addMetadata(object, assay = assay, value = metadata[[metadata_label]], label = metadata_label)
        }
      }
    } else {
      for(feat in features){
        for(label in names(getisord)){
          metadata_label <- paste(feat, label, sep = "_")
          object <- addMetadata(object, assay = assay, value = metadata[[metadata_label]], label = metadata_label)
        }
      }
    }
  }
  
  # update metadata
  # Metadata(object, assay = assay) <- metadata
  
  # return
  return(object)
}

####
# Niche Analysis ####
####

#' getNicheAssay
#'
#' Create Niche Assays
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param label grouping label for Niche definition
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param new_feature_name the name of the new feature set created for the niche assay. Default: "Niche"
#'
#' @importFrom igraph V V<- neighborhood
#' @export
getNicheAssay <- function(object, assay = NULL, label = NULL, graph.type = "delaunay", new_feature_name = "Niche"){
  
  # get metadata
  sample.metadata <- SampleMetadata(object)
  metadata <- Metadata(object)
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get graph 
  graph <- vrGraph(object, assay = assay_names, graph.type = graph.type)
  
  # get label
  cur_metadata <- subset_metadata(metadata, assays = assay_names)
  if(label %in% colnames(cur_metadata)){
    label <- as.vector(cur_metadata[,label])
    if(!is.null(rownames(cur_metadata))){
      names(label) <- rownames(cur_metadata)
    } else {
      names(label) <- as.vector(cur_metadata$id)
    }
  } else {
    stop("'", label, "' is not found in the metadata!")
  }
  
  # get niche assay
  adj_matrix <- igraph::neighborhood(graph)
  unique_label <- unique(label)
  niche_counts <- vapply(adj_matrix, function(x){
    table(factor(label[x], levels = unique_label))
  }, numeric(length(na.omit(unique_label))))
  colnames(niche_counts) <- igraph::V(graph)$name
  
  # add cell type mixtures as new feature set
  for(assy in assay_names){
    cur_niche_counts <- niche_counts[,stringr::str_extract(colnames(niche_counts), "Assay[0-9]+") %in% assy]
    object <- addFeature(object, assay = assy, data = cur_niche_counts, feature_name = new_feature_name)
  }

  # return
  return(object)
}