####
# Data Transfer ####
####

#' transferData
#'
#' transfer data across assays
#'
#' @param object a VoltRon object
#' @param from the name or class of assay whose data transfered to the second assay
#' @param to the name or class of target assay where data is transfered to
#' @param features the set of features from \link{vrFeatures} or metadata columns from \link{Metadata} that are transferred. 
#' Only one metadata feature can be transfered at a time.
#' @param expand if TRUE, metadata features will be transformed into
#' dummy features where each category in the feature will be a new feature. 
#' If FALSE, metadata features will not be transformed and transfered as
#' metadata columns, else the decision will be made automatically.
#' @param new_feature_name the name of the new feature set created 
#' from the source assay defined in \code{from} argument.
#' Only used when a new assay in created.
#'
#' @export
transferData <- function(object, 
                         from = NULL, 
                         to = NULL, 
                         features = NULL, 
                         expand = NULL,
                         new_feature_name = NULL){
  
  # assay list
  assaytypes <- c("ROI", "spot", "cell", "molecule", "tile")
  
  # get Assay IDs from Names and IDs
  from <- vrAssayNames(object, assay = from)
  to <- vrAssayNames(object, assay = to)
  
  # check assay names
  if(length(from) > 1 | length(to) > 1){
    stop("For now, label transfer can only be accomplished across two assays")
  }
  
  # check if assays are in the same block
  sample.metadata <- SampleMetadata(object)
  from_assayclass <- sample.metadata[from, "Assay"]
  to_assayclass <- sample.metadata[to, "Assay"]
  samples <- sample.metadata[c(from, to), "Sample"]
  
  if(length(unique(samples)) > 1)
    stop("Selected assays have to be within the same sample block!")
  
  # get assay types
  to_object_type <- vrAssayTypes(object[[to]])
  from_object_type <- vrAssayTypes(object[[from]])
  
  # get metadata transfer material and approach
  if(!is.null(expand)){
    expand_feature <- expand
  } else {
    from_object_ind <- which(assaytypes %in% from_object_type)
    to_object_ind <- which(assaytypes %in% to_object_type)
    expand_feature <- (from_object_ind > to_object_ind) & to_object_ind < 4
  }
  transfer_info <- .prepare_transfer(object, from, features, expand_feature)
  transfer_data <- transfer_info$data
  transfer_type <- transfer_info$type
  
  if(transfer_type == "metadata_feature"){
    return(transferLabels(object = object, 
                          from = from, 
                          to = to, 
                          transfer_data = transfer_data))
  } else {
    return(transferFeatureData(object = object, 
                               from = from, 
                               to = to, 
                               transfer_data = transfer_data, 
                               new_feature_name = new_feature_name))
  }
}

#' transferFeatureData
#'
#' transfer feature data across assays
#'
#' @param object a VoltRon object
#' @param from The ID of assay whose data transfer to the second assay
#' @param to The ID of the target assay where data is transfered to
#' @param features The set of data or metadata features that are transfered. Only one metadata feature can be transfered at a time.
#' @param new_feature_name the name of the new feature set created from the source assay defined in \code{from}.
#'
#' @noRd
transferFeatureData <- function(object, from = NULL, to = NULL, transfer_data = NULL, new_feature_name = NULL){
  
  # get assays and metadata
  from_object <- object[[from]]
  from_metadata <- Metadata(object, assay = from)
  to_object <- object[[to]]
  to_metadata <- Metadata(object, assay = to)
  
  # get assay types
  to_object_type <- vrAssayTypes(to_object)
  from_object_type <- vrAssayTypes(from_object)
  
  # get transfer data type
  if(to_object_type == "spot"){
    if(from_object_type == "cell"){
      new_assay <- getSpotsFromCells(from_object, from_metadata, 
                                     to_object, transfer_data = transfer_data)
    }
  } else if(to_object_type == "cell"){
    if(from_object_type == "tile"){
      stop("Tile to cell feature data transfer is currently not supported")
      # new_assay <- getCellsFromTiles(from_object, from_metadata, 
      #                                to_object, transfer_data = transfer_data)
    } else if(from_object_type == "spot"){
      new_assay <- getCellsFromSpots(from_object, from_metadata, 
                                     to_object, transfer_data = transfer_data)
    }
  } else if(to_object_type == "ROI"){
    if(from_object_type == "cell"){
      new_assay <- getROIsFromCells(from_object, from_metadata, 
                                    to_object, transfer_data = transfer_data)
    }
  }
  
  # add new feature set
  if(is.null(new_feature_name)){
    new_feature_name <- paste(vrMainFeatureType(object, assay = from)$Feature, 
                              "pseudo", sep = "_")
  }
  object <- addFeature(object, assay = to, 
                       data = new_assay, feature_name = new_feature_name)
  
  # return
  object
}

#' transferLabels
#'
#' transfer labels across assays
#'
#' @param object a VoltRon object
#' @param from The ID of assay whose data transfer to the second assay
#' @param to The ID of the target assay where data is transfered to
#' @param features The set of data or metadata features that are transferred. Only one metadata feature can be transferred at a time.
#'
#' @noRd
transferLabels <- function(object, from = NULL, to = NULL, transfer_data = NULL){
  
  # get assays and metadata
  from_object <- object[[from]]
  from_metadata <- Metadata(object, assay = from)
  to_object <- object[[to]]
  to_metadata <- Metadata(object, assay = to)
  
  # get assay types
  to_object_type <- vrAssayTypes(to_object)
  from_object_type <- vrAssayTypes(from_object)
  
  # get transfer data type from ROI to others
  if(from_object_type == "ROI" & to_object_type != "ROI"){
    
    # transfer labels
    transferedLabelsMetadata <- 
      transferLabelsFromROI(from_object, 
                            from_metadata, 
                            to_object, 
                            to_metadata, 
                            transfer_data = transfer_data)
    
    # update metadata
    object <- addMetadata(object, 
                          assay = to, 
                          value = transferedLabelsMetadata[[colnames(transfer_data)]], 
                          label = colnames(transfer_data))
  }
  
  # get transfer data type from tile to cell
  if(from_object_type == "tile" & to_object_type == "cell"){
    transferedLabelsMetadata <-
      transferLabelsFromTiles2Cells(from_object, 
                                    from_metadata, 
                                    to_object, 
                                    to_metadata, 
                                    transfer_data = transfer_data)
    # update metadata
    object <- addMetadata(object, 
                          assay = to, 
                          value = transferedLabelsMetadata[[colnames(transfer_data)]], 
                          label = colnames(transfer_data))
  }
  
  # return
  object
}


#' getSpotsFromCells
#'
#' Generate Psuedo counts per spots and insert to as a separate image assay to the Visium object
#'
#' @param from_object The vrAssay object whose data transfer to the second assay
#' @param from_metadata the metadata associated with \code{from_object}
#' @param to_object The ID of the target vrAssay object where data is transfered to
#' @param features the name of the metadata feature to transfer, if NULL, the rawdata will be transfered
#'
#' @importFrom dplyr %>% right_join
#' @importFrom stats aggregate
#' @importFrom magick image_data
#'
#' @noRd
#'
getSpotsFromCells <- function(from_object, from_metadata = NULL, to_object, transfer_data = NULL) {
  
  # get the spot radius of Visium spots
  Vis_spotradius <- vrAssayParams(to_object, param = "spot.radius")
  
  # get cell and spot coordinates
  message("Getting cell to spot distances ...")
  coords_spots <- vrCoordinates(to_object)
  coords_cells <- vrCoordinates(from_object)
  
  # get distances from cells to spots
  # alldist <- flexclust::dist2(coords_cells, coords_spots)
  # cell_to_spot <- FNN::get.knnx(coords_spots, coords_cells, k = 1)
  cell_to_spot <- knn_annoy(coords_spots, coords_cells, k = 1)
  names(cell_to_spot) <- c("nn.index", "nn.dist")
  cell_to_spot_nnid <- vrSpatialPoints(to_object)[cell_to_spot$nn.index[,1]]
  names(cell_to_spot_nnid) <- rownames(coords_cells)
  cell_to_spot_nndist <- cell_to_spot$nn.dist[,1]
  names(cell_to_spot_nndist) <- rownames(coords_cells)
  cell_to_spot_nnid <- cell_to_spot_nnid[cell_to_spot_nndist < Vis_spotradius]
  
  # find associated spot for each cell
  message("Find associated spots for each cell ...")
  cell_to_spot_id <- names(cell_to_spot_nnid)
  
  # get data
  # if(is.null(features)){
  #   raw_counts <- vrData(from_object, norm = FALSE)
  # } else {
  #   data_features <- features[features %in% vrFeatures(from_object)]
  #   metadata_features <- features[features %in% colnames(from_metadata)]
  #   if(length(data_features) > 0){
  #     if(length(metadata_features) > 0){
  #       stop("Data and metadata features cannot be transfered in the same time!")
  #     } else {
  #       raw_counts <- vrData(from_object, norm = FALSE)
  #       raw_counts <- raw_counts[features,]
  #       message("There are ", 
  #               length(setdiff(features, data_features)), 
  #               " unknown features!")
  #     }
  #   } else {
  #     if(length(metadata_features) > 1){
  #       stop("Only one metadata column can be transfered at a time")
  #     } else if(length(metadata_features) == 1) {
  #       raw_counts <- from_metadata[,metadata_features, drop = FALSE]
  #       rownames_raw_counts <- rownames(raw_counts)
  #       raw_counts <- dummy_cols(raw_counts, remove_first_dummy = FALSE)
  #       raw_counts <- raw_counts[,-1]
  #       raw_counts <- t(raw_counts)
  #       colnames(raw_counts) <- rownames_raw_counts
  #       rownames(raw_counts) <- gsub(paste0("^", metadata_features, "_"), "", rownames(raw_counts))
  #     } else {
  #       stop("Features cannot be found in data and metadata!")
  #     }
  #   }
  # }
  transfer_data <- transfer_data[,cell_to_spot_id, drop = FALSE]
  
  # pool cell counts to Spots
  message("Aggregating cell profiles into spots ...")
  aggregate_transfer_data <- stats::aggregate(t(as.matrix(transfer_data)), list(cell_to_spot_nnid), sum)
  aggregate_transfer_data <- data.frame(barcodes = vrSpatialPoints(to_object)) %>% dplyr::right_join(aggregate_transfer_data, by = c("barcodes" = "Group.1"))
  rownames(aggregate_transfer_data) <- aggregate_transfer_data$barcodes
  aggregate_transfer_data <- t(aggregate_transfer_data[,-1])
  aggregate_transfer_data[is.na(aggregate_transfer_data)] <- 0
  
  # return
  return(aggregate_transfer_data)
}

#' getSpotsFromCells
#'
#' Generate Psuedo counts per spots and insert to as a separate image assay to the Visium object
#'
#' @param from_object The vrAssay object whose data transfer to the second assay
#' @param from_metadata the metadata associated with \code{from_object}
#' @param to_object The ID of the target vrAssay object where data is transfered to
#' @param features the name of the metadata feature to transfer, if NULL, the rawdata will be transfered
#'
#' @importFrom dplyr %>% right_join
#' @importFrom stats aggregate
#' @importFrom magick image_data
#'
#' @noRd
#'
getCellsFromSpots <- function(from_object, from_metadata = NULL, to_object, transfer_data = NULL) {
  
  # get the spot radius of Visium spots
  radius <- vrAssayParams(from_object, param = "nearestpost.distance")/2
  
  # get cell and spot coordinates
  message("Getting spot to cell distances ...")
  coords_spots <- vrCoordinates(from_object)
  coords_cells <- vrCoordinates(to_object)
  
  # get distances from cells to spots
  spot_to_cell <- knn_annoy(coords_spots, coords_cells, k = 1)
  names(spot_to_cell) <- c("nn.index", "nn.dist")
  nnindex <- spot_to_cell$nn.index[,1]
  names_nnindex <- names(nnindex)
  nnindex <- vrSpatialPoints(from_object)[nnindex]
  names(nnindex) <- names_nnindex
  nndist <- spot_to_cell$nn.dist[,1]
  nnindex <- nnindex[nndist < radius]
  
  # find associated spot for each cell
  message("Find associated spot for each cell ...")
  
  # get data
  # if(is.null(features)){
  #   raw_counts <- vrData(from_object, norm = FALSE)
  # } else {
  #   data_features <- features[features %in% vrFeatures(from_object)]
  #   metadata_features <- features[features %in% colnames(from_metadata)]
  #   if(length(data_features) > 0){
  #     if(length(metadata_features) > 0){
  #       stop("Data and metadata features cannot be transfered in the same time!")
  #     } else {
  #       raw_counts <- vrData(from_object, norm = FALSE)
  #       raw_counts <- raw_counts[features,]
  #       message("There are ", length(setdiff(features, data_features)), " unknown features!")
  #     }
  #   } else {
  #     if(length(metadata_features) > 1){
  #       stop("Only one metadata column can be transfered at a time")
  #     } else if(length(metadata_features) == 1) {
  #       raw_counts <- from_metadata[,metadata_features, drop = FALSE]
  #       rownames_raw_counts <- rownames(raw_counts)
  #       raw_counts <- dummy_cols(raw_counts, remove_first_dummy = FALSE)
  #       raw_counts <- raw_counts[,-1]
  #       raw_counts <- t(raw_counts)
  #       colnames(raw_counts) <- rownames_raw_counts
  #       rownames(raw_counts) <- gsub(paste0("^", metadata_features, "_"), "", rownames(raw_counts))
  #     } else {
  #       stop("Features cannot be found in data and metadata!")
  #     }
  #   }
  # }
  transfer_data <- transfer_data[,nnindex, drop = FALSE]
  colnames(transfer_data) <- names(nnindex)
  
  # return
  return(transfer_data)
}

#' getROIsFromCells
#'
#' Generate Psuedo counts per ROIs and insert to as a separate image assay to the Visium object
#'
#' @param from_object The vrAssay object whose data transfer to the second assay
#' @param from_metadata the metadata associated with \code{from_object}
#' @param to_object The ID of the target vrAssay object where data is transfered to
#' @param features the name of the metadata feature to transfer, if NULL, the rawdata will be transfered
#'
#' @importFrom dplyr %>% right_join
#' @importFrom stats aggregate
#' @importFrom magick image_data
#'
#' @noRd
#'
getROIsFromCells <- function(from_object, from_metadata = NULL, to_object, transfer_data = NULL) {
  
  # get cell and ROIs coordinates
  message("Getting cell to ROI distances ...")
  segments_rois <- vrSegments(to_object)
  coords_cells <- vrCoordinates(from_object)
  
  # find associated spot for each cell
  message("Find associated ROIs for each cell ...")
  cell_to_roi_id <- NULL
  cell_to_roi_labelid <- NULL
  names_segments_rois <- names(segments_rois)
  for(i in seq_len(length(segments_rois))){
    cur_segt <- segments_rois[[i]]
    if(ncol(cur_segt) > 3){
      in.list <- point.in.circle(coords_cells[,1], coords_cells[,2], cur_segt$x, cur_segt$y, cur_segt$rx)
    } else {
      in.list <- sp::point.in.polygon(coords_cells[,1], coords_cells[,2], cur_segt$x, cur_segt$y)
    }
    in.list.cells <- rownames(coords_cells)[!!in.list]
    cell_to_roi_id <- c(cell_to_roi_id, in.list.cells)
    cell_to_roi_labelid <- c(cell_to_roi_labelid, rep(names_segments_rois[i], length(in.list.cells)))
  }
  
  # get data
  # if(is.null(features)){
  #   raw_counts <- vrData(from_object, norm = FALSE)
  # } else {
  #   data_features <- features[features %in% vrFeatures(from_object)]
  #   metadata_features <- features[features %in% colnames(from_metadata)]
  #   if(length(data_features) > 0){
  #     if(length(metadata_features) > 0){
  #       stop("Data and metadata features cannot be transfered in the same time!")
  #     } else {
  #       raw_counts <- vrData(from_object, norm = FALSE)
  #       raw_counts <- raw_counts[features,]
  #       message("There are ", length(setdiff(features, data_features)), " unknown features!")
  #     }
  #   } else {
  #     if(length(metadata_features) > 1){
  #       stop("Only one metadata column can be transfered at a time")
  #     } else if(length(metadata_features) == 1) {
  #       raw_counts <- from_metadata[,metadata_features, drop = FALSE]
  #       rownames_raw_counts <- rownames(raw_counts)
  #       raw_counts <- dummy_cols(raw_counts, remove_first_dummy = FALSE)
  #       raw_counts <- raw_counts[,-1]
  #       raw_counts <- t(raw_counts)
  #       colnames(raw_counts) <- rownames_raw_counts
  #       rownames(raw_counts) <- gsub(paste0("^", metadata_features, "_"), "", rownames(raw_counts))
  #     } else {
  #       stop("Features cannot be found in data and metadata!")
  #     }
  #   }
  # }
  transfer_data <- transfer_data[,cell_to_roi_id, drop = FALSE]
  
  # pool cell counts to Spots
  message("Aggregating cell profiles into spots ...")
  aggregate_transfer_data <- stats::aggregate(t(as.matrix(transfer_data)), list(cell_to_roi_labelid), sum)
  aggregate_transfer_data <- data.frame(barcodes = vrSpatialPoints(to_object)) %>% dplyr::right_join(aggregate_transfer_data, by = c("barcodes" = "Group.1"))
  rownames(aggregate_transfer_data) <- aggregate_transfer_data$barcodes
  aggregate_transfer_data <- t(aggregate_transfer_data[,-1])
  aggregate_transfer_data[is.na(aggregate_transfer_data)] <- 0
  
  # return
  return(aggregate_transfer_data)
}

# getCellsFromTiles <- function(from_object, from_metadata = NULL, to_object, features = NULL, k = 1) {
#   
#   # get cell and spot coordinates
#   message("Tile to Cell distances")
#   coords_cells <- vrCoordinates(to_object)
#   coords_tiles <- vrCoordinates(from_object)
#   
#   # get distances from cells to spots
#   tile_to_cell <- knn_annoy(coords_tiles, coords_cells, k = k)
#   names(tile_to_cell) <- c("nn.index", "nn.dist")
#   tile_to_cell_nnid <- data.frame(id = rownames(coords_cells), tile_to_cell$nn.index)
#   tile_to_cell_nnid <- data.table::melt(tile_to_cell_nnid, id.vars = "id")
#   tile_id <- vrSpatialPoints(from_object)[tile_to_cell_nnid$value]
#   tile_to_cell_nnid <- tile_to_cell_nnid$id
#   
#   # get data
#   raw_counts <- vrData(from_object, norm = FALSE)
#   raw_counts <- raw_counts[,tile_id]
#   
#   # pool cell counts to Spots
#   message("Aggregating tile profiles in cells")
#   aggregate_raw_counts <- stats::aggregate(t(as.matrix(raw_counts)), list(tile_to_cell_nnid), mean)
#   aggregate_raw_counts <- data.frame(barcodes = vrSpatialPoints(to_object)) %>% dplyr::right_join(aggregate_raw_counts, by = c("barcodes" = "Group.1"))
#   rownames(aggregate_raw_counts) <- aggregate_raw_counts$barcodes
#   aggregate_raw_counts <- t(aggregate_raw_counts[,-1])
#   aggregate_raw_counts[is.na(aggregate_raw_counts)] <- 0
#   
#   # return
#   return(aggregate_raw_counts)
# }

transferLabelsFromROI <- function(from_object, from_metadata = NULL, to_object, to_metadata = NULL, transfer_data = NULL) {
  
  # get ROI and other coordinates
  segments_roi <- vrSegments(from_object)
  coords <- vrCoordinates(to_object)
  spatialpoints <- rownames(coords)
  
  # # check if all features are in from_metadata
  # if(!all(features %in% colnames(from_metadata))){
  #   stop("Some features are not found in the ROI metadata!")
  # }
  
  # annotate points in the to object
  for(feat in colnames(transfer_data)){
    
    # get from metadata labels
    feat_labels <- transfer_data[,feat]
    
    # get to metadata
    new_label <- rep("undefined", length(spatialpoints))
    names(new_label) <- spatialpoints
    
    for(i in seq_len(length(segments_roi))){
      cur_poly <- segments_roi[[i]][,c("x","y")]
      in.list <- sp::point.in.polygon(coords[,1], coords[,2], cur_poly[,1], cur_poly[,2])
      new_label[rownames(coords)[!!in.list]] <- feat_labels[i]
    }
    
    to_metadata[[feat]] <- new_label
  }
  
  # return label
  return(to_metadata)
}

transferLabelsFromTiles2Cells <- function(from_object, from_metadata = NULL, 
                                          to_object, to_metadata = NULL, 
                                          transfer_data = NULL, k = 1){
  
  # get cell and spot coordinates
  message("Getting tile to cell distances ...")
  coords_cells <- vrCoordinates(to_object)
  if(!inherits(coords_cells, "IterableMatrix")){
    coords_cells <- as.data.frame(coords_cells)
  } 
  coords_cells <- as.matrix(coords_cells)
  coords_tiles <- vrCoordinates(from_object)
  if(!inherits(coords_tiles, "IterableMatrix")){
    coords_tiles <- as.data.frame(coords_tiles)
  } 
  coords_tiles <- as.matrix(coords_tiles)
  
  # spatial points of cells
  spatialpoints <- rownames(coords_cells)
  
  # get distances from cells to tile
  tile_to_cell <- knn_annoy(coords_tiles, coords_cells, k = k)
  names(tile_to_cell) <- c("nn.index", "nn.dist")
  tile_to_cell_nnid <- data.table::data.table(
    id = rownames(coords_cells), 
    tile_to_cell$nn.index)
  tile_to_cell_nnid <- data.table::melt(
    tile_to_cell_nnid, 
    id.vars = "id")
  
  # annotate points in the to object
  for(feat in colnames(transfer_data)){
    
    # get to metadata
    new_label <- rep("undefined", length(spatialpoints))
    names(new_label) <- spatialpoints
    
    # map labels
    feat_labels <- transfer_data[,feat]
    to_metadata[[feat]] <- feat_labels[tile_to_cell_nnid$value]
  }
  
  # return label
  return(to_metadata)
}

.prepare_transfer <- function(object, from, features, expand = FALSE){
  
  # get from data
  from_metadata <- Metadata(object, assay = from)
  from_object <- object[[from]]
  
  # check feature direction
  data_type <- "data_feature"
  if(is.null(features)){
    raw_counts <- vrData(from_object, norm = FALSE)
  } else {
    data_features <- features[features %in% vrFeatures(from_object)]
    metadata_features <- features[features %in% colnames(from_metadata)]
    if(length(data_features) > 0){
      if(length(metadata_features) > 0){
        stop("Data and metadata features cannot be transfered in the same time!")
      } else {
        raw_counts <- vrData(from_object, norm = FALSE)
        raw_counts <- raw_counts[features,]
        message("There are ", length(setdiff(features, data_features)), " unknown features!")
      }
    } else {
      if(length(metadata_features) > 1){
        stop("Only one metadata column can be transfered at a time")
      } else if(length(metadata_features) == 1) {
        if(inherits(from_metadata, "data.table")){
          raw_counts <- from_metadata[, get(names(from_metadata)[which(colnames(from_metadata) == metadata_features)])]
        } else{
          raw_counts <- from_metadata[,metadata_features, drop = FALSE]
        }
        raw_counts <- as.data.frame(raw_counts)
        colnames(raw_counts) <- metadata_features
        if(expand){
          rownames_raw_counts <- rownames(raw_counts)
          raw_counts <- dummy_cols(raw_counts, remove_first_dummy = FALSE)
          raw_counts <- raw_counts[,-1]
          raw_counts <- t(raw_counts)
          colnames(raw_counts) <- rownames_raw_counts
          rownames(raw_counts) <- gsub(paste0("^", metadata_features, "_"), "", rownames(raw_counts))
        } else {
          data_type <- "metadata_feature"
        }
      } else {
        stop("Features cannot be found in data and metadata!")
      }
    }
  }
  
  list(data = raw_counts, type = data_type)
}


####
# Embedding ####
####

getJointPCA <- function(object, assay.1 = NULL, assay.2 = NULL, dims = 30, seed = 1, ...){
  
  # get assay names
  assay <- assay.1
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if there are features of a VoltRon object, then get variable features too
  assay_features <- vrFeatures(object, assay = assay)
  if(length(assay_features) > 0) {
    features <- getVariableFeatures(object, assay = assay)
    vrMainAssay(object) <- object@sample.metadata[assay, "Assay"]
    object_subset <- subsetVoltRon(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)
    if(dims > length(features)){
      message("Requested more PC dimensions than existing features: dims = length(features) now!")
      dims <- length(features)
    }
  } else {
    object_subset <- object
  }
  
  # set seed
  set.seed(seed)
  
  # get PCA embedding
  normdata.1 <- vrData(object_subset, assay = assay, norm = TRUE)
  # scale.data <- apply(normdata.1, 1, scale)
  # pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  pr.data <- BiocSingular::runPCA(t(normdata.1),
                                  rank=dims,
                                  scale=TRUE,
                                  center=TRUE,
                                  BSPARAM=BiocSingular::FastAutoParam())
  pr.data.1 <- pr.data$x
  colnames(pr.data.1) <- paste0("PC", seq_len(dims))
  rownames(pr.data.1) <- colnames(normdata.1)
  normdata.1 <- normdata.1/(pr.data$sdev[1]^2)
  
  # get assay names
  assay <- assay.2
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if there are features of a VoltRon object, then get variable features too
  assay_features <- vrFeatures(object, assay = assay)
  if(length(assay_features) > 0) {
    features <- getVariableFeatures(object, assay = assay)
    vrMainAssay(object) <- object@sample.metadata[assay, "Assay"]
    object_subset <- subsetVoltRon(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)
    if(dims > length(features)){
      message("Requested more PC dimensions than existing features: dims = length(features) now!")
      dims <- length(features)
    }
  } else {
    object_subset <- object
  }
  
  # get PCA embedding
  normdata.2 <- vrData(object_subset, assay = assay, norm = TRUE)
  # scale.data <- apply(normdata.2, 1, scale)
  # pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  pr.data <- BiocSingular::runPCA(t(normdata.2),
                                  rank=dims,
                                  scale=TRUE,
                                  center=TRUE,
                                  BSPARAM=BiocSingular::FastAutoParam())
  pr.data.2 <- pr.data$x
  colnames(pr.data.2) <- paste0("PC", seq_len(dims))
  rownames(pr.data.2) <- colnames(normdata.2)
  normdata.2 <- normdata.2/(pr.data$sdev[1]^2)
  
  # get joint PCA
  normdata <- rbind(normdata.1, normdata.2)
  # scale.data <- apply(normdata, 1, scale)
  # pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  pr.data <- BiocSingular::runPCA(t(normdata),
                                  rank=dims,
                                  scale=TRUE,
                                  center=TRUE,
                                  BSPARAM=BiocSingular::FastAutoParam())
  prsdev <- pr.data$sdev
  pr.data <- pr.data$x
  colnames(pr.data) <- paste0("PC", seq_len(dims))
  rownames(pr.data) <- colnames(normdata)
  normdata <- normdata/(prsdev[1]^2)
  
  # set Embeddings
  vrEmbeddings(object, type = "pca_joint", assay = assay.1, ...) <- pr.data
  vrEmbeddings(object, type = "pca_joint", assay = assay.2, ...) <- pr.data
  
  # return
  return(object)
}