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
    
    # get coordinates
    cur_coords <- as.matrix(vrCoordinates(object, assay = assy))
    
    # get groups
    if(!is.null(group.by) && !is.null(group.ids)){
      
      # metadata
      if(verbose)
        message("Calculating Spatial Neighbors with group.by='", group.by, "' and group.ids='", paste(group.ids, collapse = ","), "'")
      metadata = Metadata(object, assay = assy)
      if(!group.by %in% colnames(metadata))
        stop("The column '", group.by, "' was not found in the metadata!")
      if(inherits(metadata, "data.table")){
        cur_group.by <- metadata[,get(names(metadata)[which(colnames(metadata) == group.by)])]
        names(cur_group.by) <- metadata$id
      } else {
        cur_group.by <- metadata[,group.by]
        if("id" %in% colnames(metadata)){
          names(cur_group.by) <- as.vector(metadata$id)
        } else {
          names(cur_group.by) <- rownames(metadata)
        }
      }
      if(!is.null(group.ids)){
        len_set_diff <- length(setdiff(group.ids,  cur_group.by))
        if(len_set_diff > 0){
        } else if(len_set_diff == length(group.ids)){ 
          stop("None of the groups defined in group.ids exist in group.by!")
        } 
        cur_group.by <- cur_group.by[cur_group.by %in% group.ids]
        cur_coords <- cur_coords[names(cur_group.by),]
      }
      
    } else if(sum(is.null(group.by),is.null(group.ids)) == 2) {
      
    } else {
      stop("Either both 'group.by' and 'group.ids' should be specified or both should be null")
    }
    
    # get edges
    spatialedges <-
      switch(method,
             delaunay = {
               nnedges <- RCDT::delaunay(cur_coords)
               nnedges <- as.vector(t(nnedges$edges[,seq_len(2)]))
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             spatialkNN = {
               # nnedges <- RANN::nn2(cur_coords, k = k + 1)
               nnedges <- knn_annoy(cur_coords, k = k + 1)
               names(nnedges) <- c("nn.index", "nn.dist")
               nnedges <- data.table::melt(data.table::data.table(nnedges$nn.index), id.vars = "V1")
               nnedges <- nnedges[,c("V1", "value")][V1 > 0 & value > 0]
               nnedges <- as.vector(t(as.matrix(nnedges)))
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             radius = {
               if(length(radius) == 0){
                 spot.radius <- vrAssayParams(object[[assy]], param = "nearestpost.distance")
                 radius <- ifelse(is.null(spot.radius), 1, spot.radius)
               }
               nnedges <- suppressWarnings({RANN::nn2(cur_coords, searchtype = "radius", radius = radius, k = min(300, sqrt(nrow(cur_coords))/2))})
               nnedges <- data.table::melt(data.table::data.table(nnedges$nn.idx), id.vars = "V1")
               nnedges <- nnedges[,c("V1", "value")][V1 > 0 & value > 0]
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
#' @param group.by a column from metadata to seperate spatial points
#' @param graph.type the type of graph to determine spatial neighborhood
#' @param num.sim the number of simulations
#' @param seed seed
#' @param verbose verbose
#'
#' @export
#'
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

  # test for each assay
  neigh_results <- list()
  for(assy in assay_names){
    if(verbose)
      message("Testing Neighborhood Enrichment of '", group.by ,"' for '", assy, "'")
    object_subset <- subsetVoltRon(object, assays = assy)
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
  if(group.by %in% colnames(metadata)){
    grp <- metadata[[group.by]]
    names(grp) <- rownames(metadata) 
  } else {
    stop("'", group.by, "' is not available in metadata!")
  }

  # get graph and neighborhood
  graph <- vrGraph(object, graph.type = graph.type)
  neighbors_graph <- igraph::neighborhood(graph)
  neighbors_graph_data <- lapply(neighbors_graph, function(x) {
    # cbind(x$name[1],x$name[-1])
    cbind(x$name[1],x$name)[-1,]
    # if(dat)
  })
  neighbors_graph_data <- do.call(rbind, neighbors_graph_data)
  colnames(neighbors_graph_data) <- c("from", "to")
  
  # get simulations
  grp_sim <- vapply(seq_len(1000), function(x) sample(grp), grp)
  rownames(grp_sim) <- names(grp)
  
  # get adjacency for observed and simulated pairs
  neighbors_graph_data_list <- list(data.frame(neighbors_graph_data, from_value = grp[neighbors_graph_data[,1]], to_value = grp[neighbors_graph_data[,2]], type = "obs"))
  for(i in 2:(ncol(grp_sim)+1))
    neighbors_graph_data_list[[i]] <- data.frame(neighbors_graph_data, from_value = grp_sim[,i-1][neighbors_graph_data[,1]], to_value = grp_sim[,i-1][neighbors_graph_data[,2]], type = paste0("sim", i))
  neighbors_graph_data <- dplyr::bind_rows(neighbors_graph_data_list)

  # get adjacency for observed and simulated pairs
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
#' @param group.ids a subset of categories defined in metadata column from \code{features}, 
#' throws an error if the feature does not include 
#' @param alpha.value the alpha value for the hot spot analysis test. Default is 0.01
#' @param norm if TRUE, the normalized data is used
#' @param n.tile should points be rasterized along x-and y-axes, typically used for molecule assays to generate features for Getis-Ord statistics
#' @param verbose verbose
#' 
#' @importFrom igraph as_adjacency_matrix
#' @importFrom stats pnorm
#'
#' @export
getHotSpotAnalysis <- function(object, assay = NULL, method = "Getis-Ord", 
                               features, graph.type = NULL,
                               group.ids = NULL,
                               alpha.value = 0.01, norm = TRUE, 
                               n.tile = 0,
                               verbose = TRUE){
  
  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")
  
  # metadata
  sample.metadata <- SampleMetadata(object)
  metadata <- Metadata(object, assay = assay)
  
  # initiate metadata columns
  for(feat in features){
    metadata[[paste0(feat,"_hotspot_flag")]] <- 
      metadata[[paste0(feat,"_hotspot_padj")]] <- 
        metadata[[paste0(feat,"_hotspot_stat")]] <- rep(NA, nrow(metadata))
  }
  if(!is.null(group.ids)){
    if(length(features) > 1)
      stop("group.ids can only be specified for a single 'feautures' entry!")
  }
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  assay_types <- vrAssayTypes(object, assay = assay_names)
  if(length(unique(assay_types)) > 1)
    stop("Only one spatial entity type can be analyzed at time!")
  
  # test for each assay
  neigh_results <- list()
  for(assy in assay_names){
    
    # verbose
    if(verbose)
      message("Running Hot Spot Analysis with '", method, "' for '", assy, "'")
    
    # metadata
    cur_metadata <- subset_metadata(metadata, assays = assy)
    
    # get graph
    if(!is.null(graph.type)){
      graph <- vrGraph(object, assay = assy, graph.type = graph.type)
      adj_matrix <- igraph::as_adjacency_matrix(graph)
    } else {
      if(n.tile == 0)
        stop("If graph.type is not provided (NULL), then the data should be rasterized with argument n.tile > 0!")
    }

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
        
        # index
        if("id" %in% colnames(cur_metadata)){
          index <- cur_metadata$id
        } else {
          index <- rownames(index)
        }
        
        # update statistics if not numeric
        mapping <- NULL
        if(!is.numeric(statistic)){
          if(n.tile == 0){
            statistic <- rowSums(adj_matrix)
          } else {
            if(!is.null(graph.type))
              warning("n.tile > 0, hence 'graph.type'=", graph.type, " is ignored")
            if(!is.null(group.ids)){
              check <- group.ids %in% statistic
              if(any(!check))
                stop("some group.ids are not found in ", features)
              group.index <- index[statistic %in% group.ids]
            }
            coords <- vrCoordinates(object, assay = assy)
            rasterization <- .rasterize_points(coords, statistic, 
                                               group.ids = group.ids, n.tile)
            statistic <- rasterization$data$count
            adj_matrix <- rasterization$adj_matrix
            mapping <- rasterization$mapping
          }
        } else {
          if(!is.null(group.ids)){
            warning("group.ids is ignored since ", feat, " is not a categorical vector!")
            group.ids <- NULL
          }
        }
      
        # get getis ord stats
        getisord <- .calculate_getis_ord(statistic, adj_matrix, alpha.value)
        
        # if mapping is given for tiles, adjust for it 
        if(!is.null(mapping)){
           getisord[[1]] <- getisord[[1]][mapping]
           getisord[[2]] <- getisord[[2]][mapping]
           getisord[[3]] <- getisord[[3]][mapping]
           mapping <- NULL
           adj_matrix <- NULL
        }
        
        # if mapping is given for group.ids, adjust for it 
        if(!is.null(group.ids)){
          tmp <- setNames(rep(NA, length(index)), index)
          for(i in 1:3){
            cur_tmp <- tmp
            cur_tmp[group.index] <- getisord[[i]]
            getisord[[i]] <- cur_tmp 
          }
        } 
        
        # update metadata for features
        for(label in names(getisord))
          cur_metadata[[paste(feat, label, sep = "_")]] <-  getisord[[label]]
      }
    }
    
    # update metadata for assays
    if("id" %in% colnames(metadata)){
      ind <- match(as.vector(cur_metadata$id), as.vector(metadata$id))
    } else {
      ind <- match(as.vector(rownames(cur_metadata)), as.vector(rownames(metadata)))
    }
    
    for(feat in features){
      for(label in names(getisord)){
        metadata_label <- paste(feat, label, sep = "_")
        metadata[[metadata_label]][ind] <- cur_metadata[[metadata_label]]        
        object <- addMetadata(object, assay = assay, value = metadata[[metadata_label]], label = metadata_label)
      }
    }
  }
  
  # return
  return(object)
}

#' @importFrom Matrix rowSums
#' @noRd
.calculate_getis_ord <- function(statistic, adj_matrix, alpha.value){
  
  # initiate getis ord
  getisord <- list()
  length(getisord) <- 3
  names(getisord) <- c("hotspot_stat", "hotspot_padj", "hotspot_flag")
  
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
  getisord[[2]] <- p.adjust(1-stats::pnorm(getisord_zscore), method = "bonferroni")
  getisord[[3]] <- ifelse(getisord[[2]] < alpha.value, "hot", "cold") 
  
  # return
  getisord
}

#' @import ggplot2
#' @noRd
.rasterize_points <- function(coords, statistic, group.ids = NULL, n.tile = 0){
  
  # check tiles
  if(n.tile < 1){
    if((n.tile %% 1 != 0))
      stop("n.tile should be bigger than 1 and an integer!")  
  }
  
  # subset group.ids
  if(!is.null(group.ids)){
    ind <- statistic %in% group.ids
    coords <- coords[ind,]
  }
  
  # rasterize
  coords <- data.frame(coords[,c("x", "y")])
  ranges <- apply(coords, 2, range)
  ranges <- ranges[2,] - ranges[1,]
  binwidth <- min(ranges/n.tile)
  gplot <- ggplot() + 
    ggplot2::stat_bin_2d(mapping = aes(x = .data[["x"]], y = .data[["y"]]),
                         data = coords,
                         binwidth = c(binwidth, binwidth), 
                         drop = FALSE)
  hex_count_data <- suppressWarnings(ggplot_build(gplot)$data)[[1]]

  # adjacency matrix 
  adj_matrix <- .make_raster_adjacency_matrix(hex_count_data)
  
  # mapping between points and rasters
  coords_tile <- sweep(coords, 2, c(binwidth, binwidth), FUN = "/")
  coords_tile <- ceiling(coords_tile)
  width_ind <- max(coords_tile$x)
  coords_tile$id <- (coords_tile$y-1)*width_ind + coords_tile$x
  mapping <- setNames(coords_tile$id, rownames(coords))
  
  # return
  hex_count_data <- hex_count_data[,c("xbin", "ybin", "count")]
  return(list(data = hex_count_data, 
              adj_matrix = adj_matrix, 
              mapping = mapping))
}

#' @importFrom igraph add_edges simplify make_empty_graph get.adjacency
#' @importFrom RANN nn2
#' @importFrom data.table data.table melt
#' @noRd
.make_raster_adjacency_matrix <- function(hex_count_data){
  
  # find nnedges
  hex_count_data$id <- 1:nrow(hex_count_data)
  # distance <- (hex_count_data$x[2] - hex_count_data$x[1]) * 1.75
  distance <- 1.75
  nnedges <- suppressWarnings({RANN::nn2(hex_count_data[,c("xbin", "ybin")], 
                                         searchtype = "radius", 
                                         radius = distance, k = 9)})

  # make edges
  nnedges <- data.table::melt(data.table::data.table(nnedges$nn.idx), id.vars = "V1")
  nnedges <- nnedges[,c("V1", "value")][V1 > 0 & value > 0]
  nnedges <- as.vector(t(as.matrix(nnedges)))
  
  # make graph
  graph <- make_empty_graph(directed = FALSE) + vertices(hex_count_data$id)
  graph <- add_edges(graph, edges = nnedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  graph <- igraph::get.adjacency(graph)
    
  # return
  graph
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
  if(is.null(label)){
    vrdata <- vrData(object, assay = assay_names, norm = FALSE)
    
    # get niche assay
    adj_matrix <- igraph::get.adjacency(graph, type = "both")
    niche_counts <- vrdata %*% adj_matrix
    
    
  } else {
    cur_metadata <- subset_metadata(metadata, assays = assay_names)
    if(label %in% colnames(cur_metadata)){
      label <- as.vector(cur_metadata[,label])
      if("id" %in% colnames(cur_metadata)){
        names(label) <- as.vector(cur_metadata$id)
      } else {
        names(label) <- rownames(cur_metadata)
      }
    } else {
      stop("'", label, "' is not found in the metadata!")
    }
    
    # get niche assay
    unique_label <- unique(label)
    adj_matrix <- igraph::neighborhood(graph)
    niche_counts <- vapply(adj_matrix, function(x){
      table(factor(label[x], levels = unique_label))
    }, numeric(length(na.omit(unique_label))))
    colnames(niche_counts) <- igraph::V(graph)$name
  }
  niche_counts <- as.matrix(niche_counts)

  # add cell type mixtures as new feature set
  for(assy in assay_names){
    cur_niche_counts <- niche_counts[,stringr::str_extract(colnames(niche_counts), "Assay[0-9]+") %in% assy]
    object <- addFeature(object, assay = assy, data = cur_niche_counts, feature_name = new_feature_name)
  }

  # return
  return(object)
}