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
#' @param calculate.distances If TRUE, distances between neighbors are also stored in the graph.
#' @param verbose verbose
#'
#' @importFrom igraph add_edges simplify make_empty_graph vertices set_edge_attr
#' @importFrom RCDT delaunay
#' @importFrom RANN nn2
#' @importFrom data.table data.table melt
#' @importFrom stats dist
#'
#' @export
getSpatialNeighbors <- function(
  object,
  assay = NULL,
  group.by = NULL,
  group.ids = NULL,
  method = "delaunay",
  k = 10,
  radius = numeric(0),
  graph.key = method,
  calculate.distances = FALSE, 
  verbose = TRUE
) {
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
        } else {
          cur_group.by <- group.by
          cur_group.ids <- setNames(group.ids,a)
        }
        
        # find neighbors based on group.by and group.ids
        if(verbose)
          message("Finding Spatial Neighbors with group.by='", 
                  cur_group.by, 
                  if(!is.null(cur_group.ids)){
                    paste0("' and group.ids='", paste(cur_group.ids, 
                                                      collapse = ","))
                  } else {
                    NULL
                  }, "'")
        
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
               nnresults <- RCDT::delaunay(coords)
               nnedges <- data.table::data.table(nnresults$edges[,seq_len(2)])
               colnames(nnedges) <- c("V1", "value")
               # nnedges <- as.vector(t(nnresults$edges[,seq_len(2)]))
               # nnedges <- rownames(coords)[nnedges]
               # nnedges$V1 <- rownames(coords)[nnedges$V1]
               # nnedges$value <- rownames(coords)[nnedges$value]
               nnedges
             },
             spatialkNN = {
               # nnedges <- RANN::nn2(coords, k = k + 1)
               nnresults <- knn_annoy(coords, k = k + 1)
               names(nnresults) <- c("nn.index", "nn.dist")
               nnedges <- data.table::melt(data.table::data.table(nnresults$nn.index), id.vars = "V1")
               if(calculate.distances){
                 nndist <- data.table::data.table(nnresults$nn.dist)
                 nndist <- data.table::melt(nndist, id.vars = "V1")
                 nnedges <- data.table::data.table(nnedges, dist = nndist$value)
               }
               nnedges <- nnedges[V1 > 0 & value > 0, !"variable"]
               # nnedges <- as.vector(t(as.matrix(nnedges)))
               # nnedges$V1 <- rownames(coords)[nnedges$V1]
               # nnedges$value <- rownames(coords)[nnedges$value]
               nnedges
             },
             radius = {
               if(length(radius) == 0){
                 spot.radius <- vrAssayParams(object[[assy]], param = "nearestpost.distance")
                 radius <- ifelse(is.null(spot.radius), 1, spot.radius)
               }
               nnresults <- suppressWarnings({RANN::nn2(coords, searchtype = "radius", radius = radius, k = min(300, sqrt(nrow(cur_coords))/2))})
               nnedges <- data.table::melt(data.table::data.table(nnresults$nn.idx), id.vars = "V1")
               if(calculate.distances){
                 nndist <- data.table::data.table(nnresults$nn.dists)
                 nndist <- data.table::melt(nndist, id.vars = "V1")
                 nnedges <- data.table::data.table(nnedges, dist = nndist$value)
               }
               nnedges <- nnedges[V1 > 0 & value > 0, !"variable"]
               # nnedges <- as.vector(t(as.matrix(nnedges)))
               # nnedges$V1 <- rownames(coords)[nnedges$V1]
               # nnedges$value <- rownames(coords)[nnedges$value]
               nnedges
             })
    spatialedges$V1 <- rownames(coords)[spatialedges$V1]
    spatialedges$value <- rownames(coords)[spatialedges$value]
    spatialedges_list <- c(spatialedges_list, list(spatialedges))
  }
  # spatialedges <- unlist(spatialedges_list)
  spatialresults <- do.call(rbind, spatialedges_list)
  spatialedges <- as.vector(t(as.matrix(spatialresults[,c("V1", "value")])))
  spatialedges
  
  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(spatialpoints)
  graph <- add_edges(graph, edges = spatialedges)
  
  # add dist if available
  if(calculate.distances){
    spatialdist <- spatialresults$dist
    graph <- igraph::set_edge_attr(graph, name = "dist", value = spatialdist) 
  }
  
  # add to VoltRon
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE, 
                    edge.attr.comb = "first")
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
  if (!inherits(object, "VoltRon")) {
    stop("Please provide a VoltRon object!")
  }

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
    res <- vrNeighbourhoodEnrichmentSingle(
      object_subset,
      group.by = group.by,
      graph.type = graph.type,
      num.sim = num.sim,
      seed = seed,
      verbose = verbose
    )
    neigh_results <- c(neigh_results, 
                       list(
                         data.frame(
                           res,
                           AssayID = paste(assy, collapse = ","),
                           unique(object_subset$Sample),
                           row.names = NULL
                         )
                       )
    )
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
#' @param verbose verbose
#'
#' @importFrom dplyr group_by bind_rows filter summarize mutate n
#' @importFrom igraph neighborhood
#'
#' @noRd
vrNeighbourhoodEnrichmentSingle <- function(
    object,
    group.by = NULL,
    graph.type = "delaunay",
    num.sim = 1000,
    seed = 1, 
    verbose = FALSE
) {
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
  if(verbose)
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
  grp_sim <- replicate(num.sim, sample(grp), simplify = "matrix")
  rownames(grp_sim) <- names(grp)
  
  # optional
  baseDT  <- as.data.table(neighbors_graph_data)
  
  # get adjacency for observed and simulated pairs
  if(verbose)
    message("Simulate Edges ...")
  from_ids <- neighbors_graph_data[, 1]
  to_ids <- neighbors_graph_data[, 2]
  neighbors_graph_data_list <- vector("list", num.sim + 1)
  # neighbors_graph_data_list[[1]] <- data.frame(
  #   neighbors_graph_data,
  #   from_value = grp[from_ids],
  #   to_value = grp[to_ids],
  #   type = "obs"
  # ) %>% 
  #   dplyr::group_by(from_value, to_value) %>%
  #   dplyr::summarize(mean_value = dplyr::n(), type = type[1])
  # TODO: put this back, if above doesn't work
  for(i in seq_len(num.sim + 1)) {
    if(verbose)
      message("Iteration ", i, " ...")
    sim_grp <- if(i == 1) grp else grp_sim[, i-1]
    neighbors_graph_data_list[[i]] <- data.frame(
      neighbors_graph_data,
      from_value = sim_grp[from_ids],
      to_value = sim_grp[to_ids],
      type = if(i == 1) "obs" else paste0("sim", i)
    ) %>%
      dplyr::group_by(from_value, to_value) %>%
      dplyr::summarize(mean_value = dplyr::n(), type = type[1])
  }
  neighbors_graph_data <- dplyr::bind_rows(neighbors_graph_data_list)
  
  # get adjacency for observed and simulated pairs
  if(verbose)
    message("Calculate Statistics ...")
  neigh_results <- neighbors_graph_data %>%
    # dplyr::group_by(from_value, to_value, type) %>%
    # dplyr::summarize(mean_value = dplyr::n()) %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::mutate(
      assoc_test = mean_value >
        ifelse("obs" %in% type, mean_value[type == "obs"], 0),
      segreg_test = mean_value <
        ifelse("obs" %in% type, mean_value[type == "obs"], 0)
    ) %>%
    dplyr::mutate(majortype = ifelse(type == "obs", "obs", "sim")) %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::mutate(
      value = ifelse(
        sum(majortype == "obs") > 0,
        log(
          mean_value[majortype == "obs"] / mean(mean_value[majortype == "sim"])
        ),
        0
      )
    ) %>%
    dplyr::filter(type != "obs") %>%
    dplyr::group_by(from_value, to_value) %>%
    dplyr::summarize(
      p_assoc = mean(assoc_test),
      p_segreg = mean(segreg_test),
      value = value[1]
    ) %>%
    dplyr::mutate(
      p_assoc_adj = p.adjust(p_assoc, method = "fdr"),
      p_segreg_adj = p.adjust(p_segreg, method = "fdr")
    )
  
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
getHotSpotAnalysis <- function(
  object,
  assay = NULL,
  method = "Getis-Ord",
  features,
  graph.type = NULL,
  group.ids = NULL,
  alpha.value = 0.01,
  norm = TRUE,
  n.tile = 0,
  verbose = TRUE
) {
  # check object
  if (!inherits(object, "VoltRon")) {
    stop("Please provide a VoltRon object!")
  }

  # metadata
  sample.metadata <- SampleMetadata(object)
  metadata <- Metadata(object, assay = assay)

  # initiate metadata columns
  for (feat in features) {
    metadata[[paste0(feat, "_hotspot_flag")]] <-
      metadata[[paste0(feat, "_hotspot_padj")]] <-
        metadata[[paste0(feat, "_hotspot_stat")]] <- rep(NA, nrow(metadata))
  }
  if (!is.null(group.ids)) {
    if (length(features) > 1) {
      stop("group.ids can only be specified for a single 'feautures' entry!")
    }
  }

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  assay_types <- vrAssayTypes(object, assay = assay_names)
  if (length(unique(assay_types)) > 1) {
    stop("Only one spatial entity type can be analyzed at time!")
  }

  # test for each assay
  neigh_results <- list()
  for (assy in assay_names) {
    # verbose
    if (verbose) {
      message("Running Hot Spot Analysis with '", method, "' for '", assy, "'")
    }

    # metadata
    cur_metadata <- subset_metadata(metadata, assays = assy)

    # get graph
    if (!is.null(graph.type)) {
      graph <- vrGraph(object, assay = assy, graph.type = graph.type)
      adj_matrix <- igraph::as_adjacency_matrix(graph)
    } else {
      if (n.tile == 0) {
        stop(
          "If graph.type is not provided (NULL), then the data should be rasterized with argument n.tile > 0!"
        )
      }
    }

    # get features
    data_features <- features[features %in% vrFeatures(object, assay = assy)]
    if (length(data_features) > 0) {
      normdata <- vrData(
        object,
        assay = assy,
        features = data_features,
        norm = norm
      )
    }

    # for each feature
    for (feat in features) {
      # Getis-Ord
      if (method == "Getis-Ord") {
        # get feature
        if (feat %in% data_features) {
          if (inherits(normdata, "IterableMatrix")) {
            statistic <- as.matrix(normdata[feat, ])[1, ]
          } else {
            statistic <- normdata[feat, ]
          }
        } else if (feat %in% colnames(cur_metadata)) {
          if (inherits(cur_metadata, "data.table")) {
            statistic <- cur_metadata[, get(names(cur_metadata)[which(
              colnames(cur_metadata) == feat
            )])]
          } else {
            statistic <- cur_metadata[, feat]
          }
        } else {
          stop(
            "'",
            feat,
            "' is not found in either the data matrix or the metadata!"
          )
        }

        # index
        if ("id" %in% colnames(cur_metadata)) {
          index <- cur_metadata$id
        } else {
          index <- rownames(cur_metadata)
        }

        # update statistics if not numeric
        mapping <- NULL
        if (!is.numeric(statistic)) {
          if (n.tile == 0) {
            statistic <- rowSums(adj_matrix)
          } else {
            if (!is.null(graph.type)) {
              warning(
                "n.tile > 0, hence 'graph.type'=",
                graph.type,
                " is ignored"
              )
            }
            if (!is.null(group.ids)) {
              check <- group.ids %in% statistic
              if (any(!check)) {
                stop("some group.ids are not found in ", features)
              }
              group.index <- index[statistic %in% group.ids]
            }
            coords <- vrCoordinates(object, assay = assy)
            rasterization <- .rasterize_points(
              coords,
              statistic,
              group.ids = group.ids,
              n.tile
            )
            statistic <- rasterization$data$count
            adj_matrix <- rasterization$adj_matrix
            mapping <- rasterization$mapping
          }
        } else {
          if (!is.null(group.ids)) {
            warning(
              "group.ids is ignored since ",
              feat,
              " is not a categorical vector!"
            )
            group.ids <- NULL
          }
        }

        # get getis ord stats
        getisord <- .calculate_getis_ord(statistic, adj_matrix, alpha.value)

        # if mapping is given for tiles, adjust for it
        if (!is.null(mapping)) {
          getisord[[1]] <- getisord[[1]][mapping]
          getisord[[2]] <- getisord[[2]][mapping]
          getisord[[3]] <- getisord[[3]][mapping]
          mapping <- NULL
          adj_matrix <- NULL
        }

        # if mapping is given for group.ids, adjust for it
        if (!is.null(group.ids)) {
          tmp <- setNames(rep(NA, length(index)), index)
          for (i in 1:3) {
            cur_tmp <- tmp
            cur_tmp[group.index] <- getisord[[i]]
            getisord[[i]] <- cur_tmp
          }
        }

        # update metadata for features
        for (label in names(getisord)) {
          cur_metadata[[paste(feat, label, sep = "_")]] <- getisord[[label]]
        }
      }
    }

    # update metadata for assays
    if ("id" %in% colnames(metadata)) {
      ind <- match(as.vector(cur_metadata$id), as.vector(metadata$id))
    } else {
      ind <- match(
        as.vector(rownames(cur_metadata)),
        as.vector(rownames(metadata))
      )
    }

    for (feat in features) {
      for (label in names(getisord)) {
        metadata_label <- paste(feat, label, sep = "_")
        metadata[[metadata_label]][ind] <- cur_metadata[[metadata_label]]
        object <- addMetadata(
          object,
          assay = assay,
          value = metadata[[metadata_label]],
          label = metadata_label
        )
      }
    }
  }

  # return
  return(object)
}

#' @importFrom Matrix rowSums
#' @noRd
.calculate_getis_ord <- function(statistic, adj_matrix, alpha.value) {
  # initiate getis ord
  getisord <- list()
  length(getisord) <- 3
  names(getisord) <- c("hotspot_stat", "hotspot_padj", "hotspot_flag")

  # calculate getis ord
  n <- length(statistic)
  statistic <- statistic - (min(statistic)) ### correct for negative scores, correct for this later
  getisord_stat <- adj_matrix %*% statistic
  getisord_stat <- getisord_stat / (sum(statistic) - statistic)
  getisord[[1]] <- getisord_stat[, 1]

  # calculate z score expectation and variance
  weight_sum <- Matrix::rowSums(adj_matrix)
  getisord_exp <- weight_sum / (n - 1)
  getisord_moment_1 <- (sum(statistic) - statistic) / (n - 1)
  getisord_moment_2 <- (sum(statistic^2) - statistic^2) /
    (n - 1) -
    getisord_moment_1^2
  getisord_var <- weight_sum * (n - 1 - weight_sum) * getisord_moment_2
  getisord_var <- getisord_var / ((n - 1)^2 * (n - 2) * getisord_moment_1^2)

  # calculate z score
  getisord_zscore <- (getisord[[1]] - getisord_exp) / sqrt(getisord_var)
  getisord_zscore[is.nan(getisord_zscore)] <- NA
  getisord[[2]] <- p.adjust(
    1 - stats::pnorm(getisord_zscore),
    method = "bonferroni"
  )
  getisord[[3]] <- ifelse(getisord[[2]] < alpha.value, "hot", "cold")

  # return
  getisord
}

#' @import ggplot2
#' @noRd
.rasterize_points <- function(coords, statistic, group.ids = NULL, n.tile = 0) {
  # check tiles
  if (n.tile < 1) {
    if ((n.tile %% 1 != 0)) {
      stop("n.tile should be bigger than 1 and an integer!")
    }
  }

  # subset group.ids
  if (!is.null(group.ids)) {
    ind <- statistic %in% group.ids
    coords <- coords[ind, ]
  }

  # rasterize
  coords <- data.frame(coords[, c("x", "y")])
  ranges <- apply(coords, 2, range)
  ranges <- ranges[2, ] - ranges[1, ]
  binwidth <- min(ranges / n.tile)
  gplot <- ggplot() +
    ggplot2::stat_bin_2d(
      mapping = aes(x = .data[["x"]], y = .data[["y"]]),
      data = coords,
      binwidth = c(binwidth, binwidth),
      drop = FALSE
    )
  hex_count_data <- suppressWarnings(ggplot_build(gplot)$data)[[1]]

  # adjacency matrix
  adj_matrix <- .make_raster_adjacency_matrix(hex_count_data)

  # mapping between points and rasters
  coords_tile <- sweep(coords, 2, c(binwidth, binwidth), FUN = "/")
  coords_tile <- ceiling(coords_tile)
  width_ind <- max(coords_tile$x)
  coords_tile$id <- (coords_tile$y - 1) * width_ind + coords_tile$x
  mapping <- setNames(coords_tile$id, rownames(coords))

  # return
  hex_count_data <- hex_count_data[, c("xbin", "ybin", "count")]
  return(list(
    data = hex_count_data,
    adj_matrix = adj_matrix,
    mapping = mapping
  ))
}

#' @importFrom igraph add_edges simplify make_empty_graph as_adjacency_matrix
#' @importFrom RANN nn2
#' @importFrom data.table data.table melt
#' @noRd
.make_raster_adjacency_matrix <- function(hex_count_data) {
  # find nnedges
  hex_count_data$id <- 1:nrow(hex_count_data)
  # distance <- (hex_count_data$x[2] - hex_count_data$x[1]) * 1.75
  distance <- 1.75
  nnedges <- suppressWarnings({
    RANN::nn2(
      hex_count_data[, c("xbin", "ybin")],
      searchtype = "radius",
      radius = distance,
      k = 9
    )
  })

  # make edges
  nnedges <- data.table::melt(
    data.table::data.table(nnedges$nn.idx),
    id.vars = "V1"
  )
  nnedges <- nnedges[, c("V1", "value")][V1 > 0 & value > 0]
  nnedges <- as.vector(t(as.matrix(nnedges)))

  # make graph
  graph <- make_empty_graph(directed = FALSE) + vertices(hex_count_data$id)
  graph <- add_edges(graph, edges = nnedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  graph <- igraph::as_adjacency_matrix(graph)

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
getNicheAssay <- function(
  object,
  assay = NULL,
  label = NULL,
  graph.type = "delaunay",
  new_feature_name = "Niche"
) {
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get metadata
  sample.metadata <- SampleMetadata(object)
  metadata <- Metadata(object, assay = assay_names)

  # get graph
  graph <- vrGraph(object, assay = assay_names, graph.type = graph.type)

  # get label
  if (is.null(label)) {
    vrdata <- vrData(object, assay = assay_names, norm = FALSE)

    # get niche assay
    adj_matrix <- igraph::as_adjacency_matrix(graph, type = "both")
    adj_matrix <- adj_matrix[colnames(vrdata), colnames(vrdata)]
    niche_counts <- vrdata %*% adj_matrix
  } else {
    cur_metadata <- subset_metadata(metadata, assays = assay_names)
    if (label %in% colnames(cur_metadata)) {
      label <- as.vector(cur_metadata[, label])
      if ("id" %in% colnames(cur_metadata)) {
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
    niche_counts <- vapply(
      adj_matrix,
      function(x) {
        table(factor(label[x], levels = unique_label))
      },
      numeric(length(na.omit(unique_label)))
    )
    colnames(niche_counts) <- igraph::V(graph)$name
  }
  niche_counts <- as.matrix(niche_counts)

  # add cell type mixtures as new feature set
  for (assy in assay_names) {
    cur_niche_counts <- niche_counts[,
      stringr::str_extract(colnames(niche_counts), "Assay[0-9]+") %in% assy
    ]
    object <- addFeature(
      object,
      assay = assy,
      data = cur_niche_counts,
      feature_name = new_feature_name
    )
  }

  # return
  return(object)
}
