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

  # get assays
  to_object <- object[[to]]
  from_object <- object[[from]]
  from_metadata <- Metadata(object, assay = from, type = from_object@type)
  from_metadata <- from_metadata[grepl(paste0(from, "$"), rownames(from_metadata)),]
  to_metadata <- Metadata(object, assay = to, type = to_object@type)
  to_metadata <- to_metadata[grepl(paste0(to, "$"), rownames(to_metadata)),]

  # get assay types
  to_object_type <- vrAssayTypes(to_object)
  from_object_type <- vrAssayTypes(from_object)

  # get transfer data type
  if(to_object_type == "spot"){
    if(from_object_type == "cell"){
      new_assay <- getSpotsFromCells(from_object, from_metadata, to_object, features = features)
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
#'
getSpotsFromCells <- function(from_object, from_metadata = NULL, to_object, features = NULL) {

  # get the spot radius of Visium spots
  Vis_spotradius <- to_object@params$spot.radius

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

  # Add as new assay
  images <- list()
  for(img in vrImageNames(to_object)){
    images[[img]] <- magick::image_data(vrImages(to_object, main_image = img))
  }
  new_assay <- new("vrAssay",
                   rawdata = aggregate_raw_counts, normdata = aggregate_raw_counts,
                   coords = vrCoordinates(to_object)[colnames(aggregate_raw_counts),],
                   image = images,
                   type =  vrAssayTypes(to_object),
                   main_image = to_object@main_image,
                   params = to_object@params)

  # return
  return(new_assay)
}
