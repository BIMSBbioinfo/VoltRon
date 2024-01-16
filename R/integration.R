####
# Data Transfer ####
####

#' transferData
#'
#' transfer data across assays
#'
#' @param object a VoltRon object
#' @param from The ID of assay whose data transfer to the second assay
#' @param to The ID of the target assay where data is transfered to
#' @param features The set of data or metadata features that are transfered. Only one metadata feature can be transfered at a time.
#' @param new_assay_name the name of the new assay created from the source assay defined in \code{from}.
#'
#' @export
#'
transferData <- function(object, from = NULL, to = NULL, features = NULL, new_assay_name = NULL){

  # check if assays are in the same block
  samples <- SampleMetadata(object)[c(from, to), "Sample"]

  if(length(unique(samples)) > 1)
    stop("Selected assays have to be within the same sample block!")

  # get from assay
  from_object <- object[[from]]
  from_metadata <- Metadata(object, assay = from, type = vrAssayTypes(from_object))

  # get to assay
  to_object <- object[[to]]
  to_metadata <- Metadata(object, assay = to, type = vrAssayTypes(to_object))

  # get assay types
  to_object_type <- vrAssayTypes(to_object)
  from_object_type <- vrAssayTypes(from_object)

  # get transfer data type
  if(to_object_type == "spot"){
    if(from_object_type == "cell"){
      new_assay <- getSpotsFromCells(from_object, from_metadata, to_object, features = features)
    }
  } else if(to_object_type == "cell"){
    if(from_object_type == "tile"){
      new_assay <- getCellsFromTiles(from_object, from_metadata, to_object, features = features)
    }
  }

  # add new assay
  if(is.null(new_assay_name)){
    new_assay_name <- paste(SampleMetadata(object)[to, "Assay"], "pseudo", sep = "_")
  }
  cat("Adding cell type compositions as new assay:", new_assay_name, "\n")
  object <- addAssay(object,
                     assay = new_assay,
                     metadata = to_metadata[vrSpatialPoints(new_assay),],
                     assay_name = new_assay_name,
                     sample = SampleMetadata(object)[to, "Sample"],
                     layer = SampleMetadata(object)[to, "Layer"])

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
#' @importFrom FNN get.knnx
#' @importFrom magick image_data
#' @importFrom fastDummies dummy_cols
#'
#' @noRd
#'
getSpotsFromCells <- function(from_object, from_metadata = NULL, to_object, features = NULL) {

  # get the spot radius of Visium spots
  # Vis_spotradius <- to_object@params$spot.radius
  Vis_spotradius <- vrAssayParams(to_object, param = "spot.radius")

  # get cell and spot coordinates
  cat("Cell to Spot Distances \n")
  coords_spots <- vrCoordinates(to_object)
  coords_cells <- vrCoordinates(from_object, reg = TRUE)

  # get distances from cells to spots
  # alldist <- flexclust::dist2(coords_cells, coords_spots)
  cell_to_spot <- FNN::get.knnx(coords_spots, coords_cells, k = 1)
  cell_to_spot_nnid <- vrSpatialPoints(to_object)[cell_to_spot$nn.index[,1]]
  names(cell_to_spot_nnid) <- rownames(coords_cells)
  cell_to_spot_nndist <- cell_to_spot$nn.dist[,1]
  names(cell_to_spot_nndist) <- rownames(coords_cells)
  cell_to_spot_nnid <- cell_to_spot_nnid[cell_to_spot_nndist < Vis_spotradius]

  # find associated spot for each cell
  cat("Find associated spots for each cell \n")
  cell_to_spot_id <- names(cell_to_spot_nnid)

  # get data
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
        raw_counts <- from_metadata[,metadata_features, drop = FALSE]
        rownames_raw_counts <- rownames(raw_counts)
        raw_counts <- fastDummies::dummy_cols(raw_counts, remove_first_dummy = FALSE)
        raw_counts <- raw_counts[,-1]
        raw_counts <- t(raw_counts)
        colnames(raw_counts) <- rownames_raw_counts
        rownames(raw_counts) <- gsub(paste0("^", metadata_features, "_"), "", rownames(raw_counts))
      } else {
        stop("Features cannot be found in data and metadata!")
      }
    }
  }
  raw_counts <- raw_counts[,cell_to_spot_id]

  # pool cell counts to Spots
  cat("Aggregating cell profiles in spots \n")
  aggregate_raw_counts <- stats::aggregate(t(as.matrix(raw_counts)), list(cell_to_spot_nnid), sum)
  aggregate_raw_counts <- data.frame(barcodes = vrSpatialPoints(to_object)) %>% dplyr::right_join(aggregate_raw_counts, by = c("barcodes" = "Group.1"))
  rownames(aggregate_raw_counts) <- aggregate_raw_counts$barcodes
  aggregate_raw_counts <- t(aggregate_raw_counts[,-1])
  aggregate_raw_counts[is.na(aggregate_raw_counts)] <- 0

  # create new assay
  # images <- list()
  # for(img in vrImageNames(to_object)){
  #   images[[img]] <- magick::image_data(vrImages(to_object, name = img))
  # }
  new_assay <- formAssay(data = aggregate_raw_counts,
                         coords = vrCoordinates(to_object)[colnames(aggregate_raw_counts),],
                         image = vrImages(to_object),
                         type = vrAssayTypes(to_object),
                         main_image = to_object@main_image,
                         params = to_object@params)
  new_assay@image <- to_object@image
  new_assay <- subset(new_assay, spatialpoints = colnames(aggregate_raw_counts))

  # return
  return(new_assay)
}

getCellsFromTiles <- function(from_object, from_metadata = NULL, to_object, features = NULL, k = 5) {

  # get cell and spot coordinates
  cat("Tile to Cell Distances \n")
  coords_cells <- vrCoordinates(to_object)
  coords_tiles <- vrCoordinates(from_object, reg = TRUE)

  # get distances from cells to spots
  # tile_to_cell <- FNN::get.knnx(coords_cells, coords_tiles, k = 1)
  tile_to_cell <- FNN::get.knnx(coords_tiles, coords_cells, k = k)
  tile_to_cell_nnid <- data.frame(id = rownames(coords_cells), tile_to_cell$nn.index)
  tile_to_cell_nnid <- reshape2::melt(tile_to_cell_nnid, id.vars = "id")
  tile_id <- vrSpatialPoints(from_object)[tile_to_cell_nnid$value]
  tile_to_cell_nnid <- tile_to_cell_nnid$id

  # get data
  raw_counts <- vrData(from_object, norm = FALSE)
  raw_counts <- raw_counts[,tile_id]

  # pool cell counts to Spots
  cat("Aggregating tile profiles in cells \n")
  aggregate_raw_counts <- stats::aggregate(t(as.matrix(raw_counts)), list(tile_to_cell_nnid), mean)
  aggregate_raw_counts <- data.frame(barcodes = vrSpatialPoints(to_object)) %>% dplyr::right_join(aggregate_raw_counts, by = c("barcodes" = "Group.1"))
  rownames(aggregate_raw_counts) <- aggregate_raw_counts$barcodes
  aggregate_raw_counts <- t(aggregate_raw_counts[,-1])
  aggregate_raw_counts[is.na(aggregate_raw_counts)] <- 0

  # create new assay
  new_assay <- formAssay(data = aggregate_raw_counts,
                         coords = vrCoordinates(to_object)[colnames(aggregate_raw_counts),],
                         image = vrImages(to_object),
                         type = vrAssayTypes(to_object),
                         main_image = to_object@main_image,
                         params = to_object@params)
  new_assay@image <- to_object@image
  new_assay <- subset(new_assay, spatialpoints = colnames(aggregate_raw_counts))

  # return
  return(new_assay)
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
    object_subset <- subset(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)
    if(dims > length(features)){
      message("Requested more PC dimensions than existing features: dims = length(features) now!")
      dims <- length(features)
    }
  } else {
    object_subset <- object
  }

  # get PCA embedding
  set.seed(seed)
  normdata.1 <- vrData(object_subset, assay = assay, norm = TRUE)
  scale.data <- apply(normdata.1, 1, scale)
  pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  pr.data.1 <- pr.data$x
  colnames(pr.data.1) <- paste0("PC", 1:dims)
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
    object_subset <- subset(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)
    if(dims > length(features)){
      message("Requested more PC dimensions than existing features: dims = length(features) now!")
      dims <- length(features)
    }
  } else {
    object_subset <- object
  }

  # get PCA embedding
  set.seed(seed)
  normdata.2 <- vrData(object_subset, assay = assay, norm = TRUE)
  scale.data <- apply(normdata.2, 1, scale)
  pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  pr.data.2 <- pr.data$x
  colnames(pr.data.2) <- paste0("PC", 1:dims)
  rownames(pr.data.2) <- colnames(normdata.2)
  normdata.2 <- normdata.2/(pr.data$sdev[1]^2)

  # get joint PCA
  normdata <- rbind(normdata.1, normdata.2)
  scale.data <- apply(normdata, 1, scale)
  pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  prsdev <- pr.data$sdev
  pr.data <- pr.data$x
  colnames(pr.data) <- paste0("PC", 1:dims)
  rownames(pr.data) <- colnames(normdata)
  normdata <- normdata/(prsdev[1]^2)

  # set Embeddings
  vrEmbeddings(object, type = "pca_joint", assay = assay.1, ...) <- pr.data
  vrEmbeddings(object, type = "pca_joint", assay = assay.2, ...) <- pr.data

  # return
  return(object)
}
