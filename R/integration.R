#' #' getSpotsFromCells
#' #'
#' #' Generate Psuedo counts per spots and insert to as a separate image assay to the Visium object
#' #'
#' #' @param Visium Visium Seurat Object
#' #' @param Xenium Xenium Seurat Object
#' #' @param scale_factors scale factors from Visium object
#' #' @param Vis_image_key the name of image with Visium spots
#' #' @param Xen_image_key the name of image with Xenium cells
#' #'
#' #' @importFrom dplyr %>% left_join
#' #' @importFrom stats aggregate
#' #'
#' #' @export
#' #'
#' getSpotsFromCells <- function(Visium, Xenium, scale_factors, Vis_image_key = "slice1", Xen_image_key = "registered_FOV") {
#'
#'   # get the spot radius of Visium spots
#'   Vis_spotradius <- scale_factors$tissue_lowres_scalef*scale_factors$spot_diameter_fullres
#'
#'   # get distances between cells to spots
#'   cat("Cell to Spot Distances \n")
#'   coords_spots <- GetTissueCoordinates(Visium, image = Vis_image_key)
#'   coords_spots[,1] = max(coords_spots[,1]) - coords_spots[,1] + min(coords_spots[,1])
#'   coords_spots <- coords_spots[,c(2,1)]
#'   coords_cells <- Xenium[[Xen_image_key]]$centroids@coords
#'   alldist <- flexclust::dist2(coords_cells, coords_spots)
#'
#'   # find associated spot for each cell
#'   cat("Find associated spots \n")
#'   spot_full <- apply(as.matrix(alldist), 1, function(x) {
#'     colnames(alldist)[order(x)[x[order(x)] < Vis_spotradius]][1]
#'   })
#'
#'   # pool cell counts to Spots
#'   cat("Cell counts to Spots \n")
#'   DataAssay <- DefaultAssay(Xenium)
#'   raw_counts <- GetAssayData(Xenium, slot = "counts", assay = DataAssay)
#'
#'   # aggregate mean
#'   raw_counts <- raw_counts[,!is.na(spot_full)]
#'   spot <- spot_full[!is.na(spot_full)]
#'   aggregate_raw_counts <- stats::aggregate(t(as.matrix(raw_counts)), list(spot), sum)
#'   aggregate_raw_counts <- data.frame(barcodes = colnames(Visium)) %>% left_join(aggregate_raw_counts, by = c("barcodes" = "Group.1"))
#'   rownames(aggregate_raw_counts) <- aggregate_raw_counts$barcodes
#'   aggregate_raw_counts <- t(aggregate_raw_counts[,-1])
#'
#'   # Pseudo Xenium
#'   cat("Create new assay in the Visium Seurat Object \n")
#'   Visium[["Pseudo_Xenium"]] <- CreateAssayObject(counts = aggregate_raw_counts)
#'
#'   # return
#'   return(list(Visium = Visium, spot = spot_full, coords_spots = coords_spots, coords_cells = coords_cells, alldist = alldist))
#' }
