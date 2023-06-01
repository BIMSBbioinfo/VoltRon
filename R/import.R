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
#' @importFrom magick image_read
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
    image <-  image_read(image_file)
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
    info <- image_info(image)
  } else {
    stop("There are no spatial image files in the path")
  }

  # coordinates
  coords_file <- list.files(paste0(dir.path, "/spatial/"), full.names = TRUE)
  coords_file <- coords_file[grepl("tissue_positions",coords_file)]
  if(length(coords_file) == 1){
    if(grepl("tissue_positions_list.csv", coords_file)) {
      coords <- read.csv(file = coords_file, header = FALSE)
      colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    } else {
      coords <- read.csv(file = coords_file, header = FALSE)
    }
  } else if(length(coords_file) > 1) {
    stop("There are more than 1 position files in the path")
  } else {
    stop("There are no files named 'tissue_positions.csv' in the path")
  }
  # coords$pxl_row_in_fullres <- max(coords$pxl_row_in_fullres) - coords$pxl_row_in_fullres + min(coords$pxl_row_in_fullres)
  if(InTissue){
    coords <- coords[coords$in_tissue==1,]
    rawdata <- rawdata[,colnames(rawdata) %in% coords$barcode]
  }
  coords <- coords[match(colnames(rawdata), coords$barcode),]
  spotID <- coords$barcode
  coords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")])
  colnames(coords) <- c("x", "y")
  rownames(coords) <- spotID

  # scale coordinates
  scale_file <- paste0(dir.path, "/spatial/scalefactors_json.json")
  if(file.exists(scale_file)){
    scalefactors <- jsonlite::read_json(path = scale_file)
    scales <- scalefactors$tissue_lowres_scalef
    # spot.radius is the half of the diameter, but we visualize by a factor of 1.5 larger
    params <- list(spot.radius = scalefactors$spot_diameter_fullres*scalefactors$tissue_lowres_scalef*1.5)
    coords <- coords*scales
    coords[,2] <- info$height - coords[,2]
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, main.assay = assay_name, params = params, assay.type = "spot", ...)
}

####
# GeoMx ####
####

#' ImportGeoMx
#'
#' @param dir.path path to GeoMx run
#' @param pkc_file path to the pkc file
#' @param summarySegment the annotation excel file
#' @param summarySegmentSheetName the sheet name of the excel file, \code{summarySegment}
#' @param assay_name the assay name, default: GeoMx
#' @param segment_polygons if TRUE, the ROI polygons are parsed from the OME.TIFF file
#' @param ome.tiff the OME.TIFF file of the GeoMx experiment if exists
#' @param ... additional parameters passed to \code{CreateSpaceRover}
#'
#' @import dplyr
#' @importFrom xlsx read.xlsx
#'
#' @export
#'
ImportGeoMx <- function(dir.path, pkc_file, summarySegment, summarySegmentSheetName, assay_name = "GeoMx", segment_polygons = FALSE, ome.tiff = NULL, ...)
{
  if (!requireNamespace('GeomxTools'))
    stop("Please install Seurat package for using Seurat objects")

  # Get pkc file
  pkcdata <- GeomxTools::readPKCFile(pkc_file)

  # Get dcc file
  dcc_files <- dir(dir.path, pattern = ".dcc$", full.names = TRUE)
  dcc_files <- dcc_files[!grepl("A01.dcc$", dcc_files)]
  dcc_filenames <- dir(dir.path, pattern = ".dcc$", full.names = FALSE)
  dcc_filenames <- dcc_filenames[!grepl("A01.dcc$", dcc_filenames)]
  dcc_filenames <- gsub(".dcc$", "", dcc_filenames)
  dcc_filenames <- gsub("-", "_", dcc_filenames)
  dccData <- sapply(dcc_files, GeomxTools::readDccFile, simplify = FALSE, USE.NAMES = FALSE)
  names(dccData) <- dcc_filenames

  # merge dcc files
  rawdata <- NULL
  for(i in 1:length(dccData)){
    cur_data <- dccData[[i]]$Code_Summary
    colnames(cur_data) <- c("RTS_ID", dcc_filenames[i])
    if(i == 1){
      rawdata <- cur_data
    } else {
      suppressMessages(rawdata <- rawdata %>% full_join(cur_data))
    }
  }
  rawdata[is.na(rawdata)] <- 0
  rownames(rawdata) <- rawdata$RTS_ID
  rawdata <- rawdata[,!colnames(rawdata) %in% "RTS_ID"]

  # get genes
  NegProbes <- pkcdata$RTS_ID[pkcdata$Target == "NegProbe-WTX"]
  rawdata <- rawdata[!rownames(rawdata) %in% NegProbes, ]
  rownames(rawdata) <- pkcdata$Target[match(rownames(rawdata), pkcdata$RTS_ID)]
  rawdata <- as.matrix(rawdata)

  # get segment summary
  segmentsummary <- read.xlsx(summarySegment, sheetName = summarySegmentSheetName)

  # get image
  image <- image_read(paste0(dir.path, "/morphology.tiff"))
  geomx_image_info <- image_info(image)

  # get coordinates
  coords <- segmentsummary[,c("X","Y")]
  colnames(coords) <- c("x", "y")
  rownames(coords) <- segmentsummary$ROI.name
  coords <- rescaleGeoMxPoints(coords, segmentsummary, geomx_image_info)

  # get ROI segments (polygons)
  segments <- list()
  if(segment_polygons){
    if(is.null(ome.tiff)){
      ome.tiff <- paste0(dir.path, "/geomx.ome.tiff")
    }
    segments <- ImportGeoMxSegments(ome.tiff, segmentsummary, geomx_image_info)
  }

  # create SpaceRover
  CreateSpaceRover(rawdata, metadata = NULL, image, coords, segments, main.assay = assay_name, assay.type = "ROI", ...)
}


#' ImportGeoMxSegments
#'
#' Import ROI polygons from the OME.TIFF file
#'
#' @param ome.tiff the OME.TIFF file of the GeoMx Experiment
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#'
#' @importFrom RBioFormats read.omexml
#' @importFrom XML xmlToList
#'
ImportGeoMxSegments <- function(ome.tiff, summary, imageinfo){

  # get the xml file
  # xmltemp <- xmlExtraction(ometiff = ome.tiff)
  omexml <- read.omexml(ome.tiff)
  omexml <- xmlToList(omexml, simplify = TRUE)

  # get ROIs
  ROIs <- omexml[which(names(omexml) == "ROI")]

  # Y-axis height
  sizeY <- as.numeric(omexml$Image$Pixels$.attrs['SizeY'])

  # get masks for each ROI
  mask_lists <- list()
  for(i in 1:length(ROIs)){
    cur_ROI <- ROIs[[i]]
    # if the shape is a polygon
    if("Polygon" %in% names(cur_ROI$Union)){
      coords <- strsplit(cur_ROI$Union$Polygon, split = "\\n")
      coords <- strsplit(coords$Points, split = " ")[[1]]
      coords <- sapply(coords, function(x) as.numeric(strsplit(x, split = ",")[[1]]), USE.NAMES = FALSE)
      coords <- as.data.frame(t(coords))
      colnames(coords) <- c("x", "y")
      coords <- rescaleGeoMxPoints(coords, summary, imageinfo)
      mask_lists[[cur_ROI$Union$Label[["Text"]]]] <- data.frame(coords)

    # if the shape is an ellipse
    } else if("Ellipse" %in% names(cur_ROI$Union)){
      coords <- as.numeric(cur_ROI$Union$Ellipse[c("X","Y", "RadiusX", "RadiusY")])
      coords <-  as.data.frame(matrix(coords, nrow = 1))
      colnames(coords) <- c("x", "y", "rx", "ry")
      coords[,c("x", "y")] <- rescaleGeoMxPoints(coords[,c("x", "y")], summary, imageinfo)
      coords$rx <- coords$rx * imageinfo$width/summary$Scan.Width[1]
      coords$ry <- coords$ry * imageinfo$height/summary$Scan.Height[1]
      mask_lists[[cur_ROI$Union$Label[["Text"]]]] <- coords
    }
  }

  mask_lists
}

#' rescaleGeoMxROIs
#'
#' rescale GeoMx point (center or polygon corners of ROI) coordinates for image registration
#'
#' @param pts coordinates of the cells from the Xenium assays
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#'
rescaleGeoMxPoints <- function(pts, summary, imageinfo){

  xRatio = imageinfo$width/summary$Scan.Width[1]
  yRatio = imageinfo$height/summary$Scan.Height[1]
  pts$x = (pts$x - summary$Scan.Offset.X[1]) * xRatio
  pts$y = (pts$y - summary$Scan.Offset.Y[1]) * yRatio
  pts$y = imageinfo$height - pts$y
  pts <- as.matrix(pts)

  # return
  return(pts)
}
