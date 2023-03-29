#' ImportXenium
#'
#' ImportXenium
#'
#' @param dir.path path to Xenium
#' @param selected_assay selected assay
#'
#' @import hdf5r
#' @import magick
#'
#' @export
#'
ImportXenium <- function (dir.path, selected_assay = "Gene Expression", layer_name = NULL, sample_name = NULL, assay_name = "Xenium", assay_type = "cell", ...)
{
  # raw counts
  infile <- Read10X_h5(filename = paste0(dir.path, "/cell_feature_matrix.h5"))
  rawdata <- as.matrix(infile[[selected_assay]])
  colnames(rawdata) <- paste0("Cell",1:ncol(rawdata))

  # coordinates
  Xenium_coords <- read.csv(file = paste0(dir.path, "/cells.csv.gz"))
  coords <- as.matrix(Xenium_coords[,c("x_centroid", "y_centroid")])
  rownames(coords) <-  paste0("Cell",Xenium_coords$cell_id)

  # image
  morphology_image <- paste0(dir.path, "/morphology_lowres.tif")
  image <-  magick::image_read(morphology_image)

  # create srAssay
  Xenium_assay <- new("srAssay", rawdata = rawdata, normdata = rawdata, coords = coords, image = image, type = assay_type)
  listofAssays <- list(Xenium_assay)
  names(listofAssays) <- assay_name

  # create layers
  layer_name <- ifelse(is.null(layer_name), "slide1", layer_name)
  sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  listofLayers <- list(new("srLayer", assay = listofAssays))
  names(listofLayers) <- layer_name

  # create samples
  sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  listofSamples <- list(new("srSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # create SpaceRover
  CreateSpaceRover(listofSamples, ...)
}

#' ImportVisium
#'
#' @param dir.path path to Xenium
#'
#' @import hdf5r
#' @import magick
#'
#' @export
#'
ImportVisium <- function(dir.path, layer_name = NULL, sample_name = NULL, assay_name = "Visium", assay_type = "Spot", ...)
{
  # raw counts
  listoffiles <- list.files(dir.path)
  infile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  infile <- Read10X_h5(filename = paste0(dir.path, "/", infile))
  rawdata <- as.matrix(infile)

  # coordinates
  Visium_coords <- read.csv(file = paste0(dir.path, "/spatial/tissue_positions.csv"))
  print(head(Visium_coords))
  coords <- as.matrix(Visium_coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  rownames(coords) <-  Visium_coords$barcode

  # image
  image <- paste0(dir.path, "/spatial/tissue_lowres_image.png")
  image <-  magick::image_read(image)

  # create srAssay
  Visium_assay <- new("srAssay", rawdata = rawdata, normdata = rawdata, coords = coords, image = image, type = assay_type)
  listofAssays <- list(Visium_assay)
  names(listofAssays) <- assay_name

  # create layers
  layer_name <- ifelse(is.null(layer_name), "slide1", layer_name)
  sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  listofLayers <- list(new("srLayer", assay = listofAssays))
  names(listofLayers) <- layer_name

  # create samples
  sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  listofSamples <- list(new("srSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # create SpaceRover
  CreateSpaceRover(listofSamples, ...)
}
