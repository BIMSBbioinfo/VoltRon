

transferData <- function(object, from = NULL, to = NULL){

  # get assays
  to_object <- object[[to]]
  from_object <- object[[from]]

  # get assay types
  to_object_type <- vrAssayTypes(to_object)
  from_object_type <- vrAssayTypes(from_object)

  # get transfer data type
  if(to_object_type == "spot"){
    if(from_object_type == "cell"){
      new_assay <- getSpotsFromCells(from_object, to_object)
    }
  }

  # add new assay
  cat("Adding cell type compositions as new assay:", paste(SampleMetadata(object)[to, "Assay"], "pseudo", sep = "_"), "\n")
  object <- addAssay(object,
                     assay = new_assay,
                     assay_name = paste(SampleMetadata(object)[to, "Assay"], "pseudo", sep = "_"),
                     sample = SampleMetadata(object)[to, "Sample"],
                     layer = SampleMetadata(object)[to, "Layer"])
  object
}


#' getSpotsFromCells
#'
#' Generate Psuedo counts per spots and insert to as a separate image assay to the Visium object
#'
#' @param from_object The first vrAssay
#' @param to_object The second vrAssay
#'
#' @importFrom dplyr %>% right_join
#' @importFrom stats aggregate
#' @importFrom FNN get.knnx
#' @importFrom magick image_data
#'
#' @export
#'
getSpotsFromCells <- function(from_object, to_object) {

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

  # pool cell counts to Spots
  raw_counts <- vrData(from_object, norm = FALSE)
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
