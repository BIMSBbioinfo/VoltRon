####
# Xenium ####
####

#' ImportXenium
#'
#' Importing a Xenium data
#'
#' @param dir.path path to Xenium
#' @param selected_assay selected assay
#' @param assay_name the assay name of the SR object
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import magick
#'
#' @export
#'
ImportXenium <- function (dir.path, selected_assay = "Gene Expression", assay_name = "Xenium", ...)
{
  # raw counts
  datafile <- paste0(dir.path, "/cell_feature_matrix.h5")
  if(file.exists(datafile)){
    rawdata <- Read10X_h5(filename = datafile)
    rawdata <- as.matrix(rawdata[[selected_assay]])
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # image
  image_file <- paste0(dir.path, "/morphology_lowres.tif")
  if(file.exists(image_file)){
    image <-  magick::image_read(image_file)
  } else {
    stop("There are no spatial image files in the path")
  }

  # cell boundaries
  bound_file <- paste0(dir.path, "/cell_boundaries.csv.gz")
  if(file.exists(bound_file)){
    Xenium_boundaries <- read.csv(bound_file)
    Xenium_box <- apply(Xenium_boundaries[,-1], 2, range)
  } else {
    stop("There are no files named 'cell_boundaries.csv.gz' in the path")
  }

  # coordinates
  coord_file <- paste0(dir.path, "/cells.csv.gz")
  if(file.exists(coord_file)){
    Xenium_coords <- read.csv(file = coord_file)
    coords <- as.matrix(Xenium_coords[,c("x_centroid", "y_centroid")])
    colnames(coords) <- c("x","y")
    coords[,2] <- max(coords[,2]) - coords[,2] + min(coords[,2])
    coords <- rescaleXeniumCells(coords, Xenium_box, image)
  } else {
    stop("There are no files named 'cells.csv.gz' in the path")
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, assay.type = "cell", ...)
}

#' rescaleXeniumCells
#'
#' rescale Xenium cells coordinates for image registration
#'
#' @param cells coordinates of the cells from the Xenium assays
#' @param bbox the surrounding box of the Xenium cell coordinates
#' @param image reference image
#'
rescaleXeniumCells <- function(cells, bbox, image){

  # get image scales
  scales <- unlist(image_info(image)[c("width","height")])

  # rescale cell locations
  cells[,1] <- (cells[,1] - bbox[1,1])/(bbox[2,1] - bbox[1,1])
  cells[,1] <- cells[,1] * scales[1]
  cells[,2] <- (cells[,2] - bbox[1,2])/(bbox[2,2] - bbox[1,2])
  cells[,2] <- cells[,2] * scales[2]

  # return
  return(cells)
}

####
# Visium ####
####

#' ImportVisium
#'
#' @param dir.path path to Xenium
#' @param assay_name the assay name
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import hdf5r
#' @import magick
#' @import jsonlite
#'
#' @export
#'
ImportVisium <- function(dir.path, assay_name = "Visium", InTissue = TRUE, ...)
{
  # raw counts
  listoffiles <- list.files(dir.path)
  datafile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  datafile <- paste0(dir.path, "/", datafile)
  if(file.exists(datafile)){
    rawdata <- Read10X_h5(filename = datafile)
    rawdata <- as.matrix(rawdata)
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # image
  image_file <- paste0(dir.path, "/spatial/tissue_lowres_image.png")
  if(file.exists(image_file)){
    image <-  magick::image_read(image_file)
  } else {
    stop("There are no spatial image files in the path")
  }

  # coordinates
  coords_file <- paste0(dir.path, "/spatial/tissue_positions.csv")
  if(file.exists(coords_file)){
    coords <- read.csv(file = coords_file)
  } else {
    stop("There are no files named 'tissue_positions.csv' in the path")
  }
  coords$pxl_row_in_fullres <- max(coords$pxl_row_in_fullres) - coords$pxl_row_in_fullres + min(coords$pxl_row_in_fullres)
  if(InTissue){
    coords <- coords[coords$in_tissue==1,]
    rawdata <- rawdata[,colnames(rawdata) %in% coords$barcode]
  }
  spotID <- coords$barcode
  coords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  colnames(coords) <- c("x", "y")
  rownames(coords) <- spotID

  # scale coordinates
  scale_file <- paste0(dir.path, "/spatial/scalefactors_json.json")
  if(file.exists(scale_file)){
    scalefactors <- jsonlite::read_json(path = scale_file)
    scales <- scalefactors$tissue_lowres_scalef
    coords <- coords*scales
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, assay.type = "spot", ...)
}


#' ImportVisium_old
#'
#' @param dir.path path to Xenium
#' @param assay_name the assay name
#' @param layer_name the layer name
#' @param sample_name the sample name
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#'
#' @import hdf5r
#' @import magick
#' @import jsonlite
#'
#' @export
#'
ImportVisium_old <- function(dir.path, assay_name = "Visium", layer_name = NULL, sample_name = NULL, InTissue = TRUE, ...)
{
  # raw counts
  listoffiles <- list.files(dir.path)
  infile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  infile <- Read10X_h5(filename = paste0(dir.path, "/", infile))
  rawdata <- as.matrix(infile)

  # labels
  layer_name <- ifelse(is.null(layer_name), "slide1", layer_name)
  sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  colnames(rawdata) <- paste(colnames(rawdata),
                             paste(c(assay_name, layer_name, sample_name), collapse = "_"),
                             sep = "_")

  # coordinates
  coords <- read.csv(file = paste0(dir.path, "/spatial/tissue_positions.csv"))
  coords$pxl_row_in_fullres <- max(coords$pxl_row_in_fullres) - coords$pxl_row_in_fullres + min(coords$pxl_row_in_fullres)
  if(InTissue){
    coords <- coords[coords$in_tissue==1,]
    rawdata <- rawdata[rownames(rawdata) %in% coords$barcode,]
  }
  spotID <- coords$barcode
  coords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  rownames(coords) <- spotID

  # scale coordinates
  scalefactors <- jsonlite::read_json(path = paste0(dir.path, "/spatial/scalefactors_json.json"))
  scales <- scalefactors$tissue_lowres_scalef
  coords <- coords*scales

  # image
  image <- paste0(dir.path, "/spatial/tissue_lowres_image.png")
  image <-  magick::image_read(image)

  # create srAssay
  Visium_assay <- new("srAssay", rawdata = rawdata, normdata = rawdata, coords = coords, image = image, type = "spot")
  listofAssays <- list(Visium_assay)
  names(listofAssays) <- assay_name

  # create layers and samples
  listofLayers <- list(new("srLayer", assay = listofAssays))
  names(listofLayers) <- layer_name
  listofSamples <- list(new("srSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # create metadata
  metadata <- setSRMetadata(cell = data.frame(),
                            spot = data.frame(Count = colSums(rawdata), row.names = colnames(rawdata)),
                            ROI = data.frame())

  # create SpaceRover
  CreateSpaceRover(listofSamples, metadata = metadata, main.assay = "Visium", ...)
}
