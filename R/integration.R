#' getSpotsFromCells
#'
#' Generate Psuedo counts per spots and insert to as a separate image assay to the Visium object
#'
#' @param Visium Visium Seurat Object
#' @param Xenium Xenium Seurat Object
#'
#' @importFrom dplyr %>% left_join
#' @importFrom stats aggregate
#' @importFrom FNN get.knnx
#' @importFrom magick image_data
#'
#' @export
#'
getSpotsFromCells <- function(Visium, Xenium) {

  # get the spot radius of Visium spots
  Vis_spotradius <- Visium[["Assay1"]]@params$spot.radius

  # get cell and spot coordinates
  cat("Cell to Spot Distances \n")
  coords_spots <- vrCoordinates(Visium)
  coords_cells <- vrCoordinates(Xenium, reg = TRUE)

  # get distances from cells to spots
  # alldist <- flexclust::dist2(coords_cells, coords_spots)
  cell_to_spot <- FNN::get.knnx(coords_spots, coords_cells, k = 1)
  cell_to_spot_nnid <- vrSpatialPoints(Visium)[cell_to_spot$nn.index[,1]]
  names(cell_to_spot_nnid) <- rownames(coords_cells)
  cell_to_spot_nndist <- cell_to_spot$nn.dist[,1]
  names(cell_to_spot_nndist) <- rownames(coords_cells)
  cell_to_spot_nnid <- cell_to_spot_nnid[cell_to_spot_nndist < Vis_spotradius]

  # find associated spot for each cell
  cat("Find associated spots for each cell \n")
  cell_to_spot_id <- names(cell_to_spot_nnid)

  # pool cell counts to Spots
  raw_counts <- vrData(Xenium, norm = FALSE)
  raw_counts <- raw_counts[,cell_to_spot_id]

  # pool cell counts to Spots
  cat("Aggregating cell profiles in spots \n")
  aggregate_raw_counts <- stats::aggregate(t(as.matrix(raw_counts)), list(cell_to_spot_nnid), sum)
  aggregate_raw_counts <- data.frame(barcodes = vrSpatialPoints(Visium)) %>% dplyr::left_join(aggregate_raw_counts, by = c("barcodes" = "Group.1"))
  rownames(aggregate_raw_counts) <- aggregate_raw_counts$barcodes
  aggregate_raw_counts <- t(aggregate_raw_counts[,-1])
  aggregate_raw_counts[is.na(aggregate_raw_counts)] <- 0

  # Add as new assay
  cat("Adding cell type compositions as new assay:", paste(SampleMetadata(Visium)[["Assay"]], "pseudo", sep = "_"), "\n")
  images <- list()
  for(img in vrImageNames(Visium)){
    images[[img]] <- magick::image_data(vrImages(Visium, main_image = img))
  }
  assy <- new("vrAssay",
                 rawdata = aggregate_raw_counts, normdata = aggregate_raw_counts,
                 coords = vrCoordinates(Visium),
                 image = images,
                 type =  vrAssayTypes(Vis_seu_reg),
                 main_image = Vis_seu_reg[["Assay1"]]@main_image,
                 params = Vis_seu_reg[["Assay1"]]@params)
  Visium <- addAssay(Visium,
                     assay = assy,
                     assay_name = paste(SampleMetadata(Visium)[["Assay"]], "pseudo", sep = "_"),
                     sample = SampleMetadata(Visium)[["Sample"]],
                     layer = SampleMetadata(Visium)[["Layer"]])

  # return
  return(Visium)
}
