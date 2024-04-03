####
# 10X Genomics ####
####

####
## Xenium ####
####

#' importXenium
#'
#' Importing Xenium data
#'
#' @param dir.path path to Xenium output folder
#' @param selected_assay selected assay from Xenium
#' @param assay_name the assay name of the SR object
#' @param sample_name the name of the sample
#' @param use_image if TRUE, the DAPI image will be used.
#' @param morphology_image the name of the lowred morphology image. Default: morphology_lowres.tif
#' @param resolution_level the level of resolution within Xenium OME-TIFF image, see \link{generateXeniumImage}. Default: 7 (553x402)
#' @param overwrite_resolution if TRUE, the image "file.name" will be generated again although it exists at "dir.path"
#' @param image_name the image name of the Xenium assay, Default: main
#' @param channel_name the channel name of the image of the Xenium assay, Default: DAPI
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @importFrom ids random_id
#'
#' @export
#'
importXenium <- function (dir.path, selected_assay = "Gene Expression", assay_name = "Xenium", sample_name = NULL, use_image = TRUE, morphology_image = "morphology_lowres.tif", resolution_level = 7, overwrite_resolution = FALSE, image_name = "main", channel_name = "DAPI", import_molecules = FALSE, ...)
{
  # cell assay
  message("Creating cell level assay ...")

  # raw counts
  datafile <- paste0(dir.path, "/cell_feature_matrix.h5")
  if(file.exists(datafile)){
    rawdata <- import10Xh5(filename = datafile)
    if(any(names(rawdata) %in% selected_assay)){
      rawdata <- as.matrix(rawdata[[selected_assay]])
    } else {
      stop("There are no assays called ", selected_assay, " in the h5 file!")
    }
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # image
  if(use_image){
    suppressMessages(generateXeniumImage(dir.path, file.name = morphology_image, resolution_level = resolution_level, overwrite_resolution = overwrite_resolution))
    image_file <- paste0(dir.path, "/", morphology_image)
    if(file.exists(image_file)){
      image <-  image_read(image_file)
    } else {
      stop("There are no spatial image files in the path")
    }
    # scale the xenium image instructed by 10x Genomics help page
    scaleparam <- 0.2125*(2^(resolution_level-1))
  } else {
    image <- NULL
  }

  # coordinates
  coord_file <- paste0(dir.path, "/cells.csv.gz")
  if(file.exists(coord_file)){
    Xenium_coords <- utils::read.csv(file = coord_file)
    coords <- as.matrix(Xenium_coords[,c("x_centroid", "y_centroid")])
    colnames(coords) <- c("x","y")
    if(use_image) {
      coords <- coords/scaleparam
      imageinfo <- unlist(magick::image_info(image)[c("height")])
      range_coords <- c(0,imageinfo)
    } else {
      range_coords <- range(coords[,2])
    }
    coords[,2] <- range_coords[2] - coords[,2] + range_coords[1]
  } else {
    stop("There are no file named 'cells.csv.gz' in the path")
  }

  # segments
  segments_file <- paste0(dir.path, "/cell_boundaries.csv.gz")
  if(file.exists(segments_file)){
    segments <- as.data.frame(data.table::fread(segments_file))
    colnames(segments) <- c("cell_id", "x", "y")
    if(use_image){
      segments[,c("x","y")] <- segments[,c("x","y")]/scaleparam
    }
    segments[,"y"] <- range_coords[2] - segments[,"y"]  + range_coords[1]
    segments <- segments %>% dplyr::group_split(cell_id)
    segments <- as.list(segments)
    names(segments) <- rownames(coords)
  } else {
    stop("There are no file named 'cell_boundaries.csv.gz' in the path")
  }

  # create VoltRon object for cells
  cell_object <- formVoltRon(rawdata, metadata = NULL, image = image, coords, segments = segments, main.assay = assay_name, assay.type = "cell", image_name = image_name, main_channel = channel_name, sample_name = sample_name, ...)

  # molecule assay
  if(!import_molecules){
    return(cell_object)
  } else {
    message("Creating molecule level assay ...")
    # transcripts
    transcripts_file <- paste0(dir.path, "/transcripts.csv.gz")
    if(!file.exists(transcripts_file)){
      stop("There are no file named 'transcripts.csv.gz' in the path")
    } else {
      # get subcellur data components
      subcellular_data <- data.table::fread(transcripts_file)
      subcellular_data <- subcellular_data[,c("transcript_id", colnames(subcellular_data)[!colnames(subcellular_data) %in% "transcript_id"]), with = FALSE]
      colnames(subcellular_data)[colnames(subcellular_data)=="transcript_id"] <- "id"
      colnames(subcellular_data)[colnames(subcellular_data)=="feature_name"] <- "gene"
      subcellular_data <- subcellular_data[subcellular_data$qv >= 20, ]

      # coordinates
      coords <- as.matrix(subcellular_data[,c("x_location", "y_location")])
      colnames(coords) <- c("x", "y")
      if(use_image){
        coords <- coords/scaleparam
      }
      coords[,"y"] <- range_coords[2] - coords[,"y"]  + range_coords[1]

      # metadata
      mol_metadata <- subcellular_data[,colnames(subcellular_data)[!colnames(subcellular_data) %in% c("cell_id", "transcript_id", "x_location", "y_location")], with = FALSE]
      set.seed(nrow(mol_metadata$id))
      mol_metadata$postfix <- paste0("_", ids::random_id(bytes = 3, use_openssl = FALSE))
      mol_metadata$assay_id <- "Assay1"
      mol_metadata[, id:=do.call(paste0,.SD), .SDcols=c("id", "postfix")]

      # coord names
      rownames(coords) <- mol_metadata$id

      # create VoltRon assay for molecules
      mol_assay <- formAssay(coords = coords, image = image, type = "molecule", main_image = image_name, main_channel = channel_name)

      # merge assays in one section
      message("Merging assays ...")
      sample.metadata <- SampleMetadata(cell_object)
      object <- addAssay(cell_object,
                         assay = mol_assay,
                         metadata = mol_metadata,
                         assay_name = paste0(assay_name, "_mol"),
                         sample = sample.metadata["Assay1", "Sample"],
                         layer = sample.metadata["Assay1", "Layer"])

      # connectivity
      connectivity <- data.table::data.table(transcript_id = mol_metadata$id, cell_id = subcellular_data[["cell_id"]])
      if(is.numeric(connectivity$cell_id)){
        connectivity <- subset(connectivity, cell_id != -1)
        connectivity[["cell_id"]] <- vrSpatialPoints(cell_object)[connectivity[["cell_id"]]]
      } else {
        connectivity <- subset(connectivity, cell_id != "UNASSIGNED")
        connectivity$cell_assay_id <- "_Assay1"
        connectivity[, cell_id:=do.call(paste0,.SD), .SDcols=c("cell_id", "cell_assay_id")]
        connectivity$cell_assay_id <- NULL
      }

      # add connectivity of spatial points across assays
      object <- addConnectivity(object,
                                connectivity = connectivity,
                                sample = sample.metadata["Assay1", "Sample"],
                                layer = sample.metadata["Assay1", "Layer"])

      # return
      return(object)
    }
  }
}

#' generateXeniumImage
#'
#' Generate a low resolution DAPI image of the Xenium experiment
#'
#' @param dir.path Xenium output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param resolution_level the level of resolution within Xenium OME-TIFF image. Default: 7 (553x402)
#' @param overwrite_resolution if TRUE, the image "file.name" will be generated again although it exists at "dir.path"
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param file.name the name of the lowred morphology image. Default: morphology_lowres.tif
#' @param ... additional parameters passed to the \link{writeImage} function
#'
#' @importFrom EBImage writeImage
#'
#' @details
#' The Xenium morphology_mip.ome.tif file that is found under the outs folder comes is an hyperstack of different resolutions of the DAPI image.
#' \link{generateXeniumImage} allows extracting only one of these layers by specifying the \code{resolution} parameter (Default: 7 for 553x402) among 1 to 8.
#' Lower incides of resolutions have higher higher resolutions, e.g. 1 for 35416x25778. Note that you may need to allocate larger memory of Java to import
#' higher resolution images.
#'
#' @export
#'
generateXeniumImage <- function(dir.path, increase.contrast = TRUE, resolution_level = 7, overwrite_resolution = FALSE, output.path = NULL, file.name = "morphology_lowres.tif", ...) {

  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/", file.name)
  output.file <- paste0(output.path, "/", file.name)

  # check if the file exists in either Xenium output folder, or the specified location
  if((file.exists(file.path) | file.exists(paste0(output.file))) & !overwrite_resolution){
    message(paste0(file.name, " already exists!"))
  } else {
    message("Loading morphology_mip.ome.tif \n")
    if (!requireNamespace('RBioFormats'))
      stop("Please install RBioFormats package to read the ome.tiff file!")
    morphology_image_lowres <- RBioFormats::read.image(paste0(dir.path, "/morphology_mip.ome.tif"), resolution = resolution_level)

    # pick a resolution level
    image_info <- morphology_image_lowres@metadata$coreMetadata
    message(paste0("Image Resolution (X:", image_info$sizeX, " Y:", image_info$sizeY, ") \n"))

    # increase contrast using EBImage
    if(increase.contrast) {
      message("Increasing Contrast \n")
      morphology_image_lowres <- (morphology_image_lowres/max(morphology_image_lowres))
    }

    # write to the same folder
    message("Writing Tiff File")
    if(is.null(output.path)){
      EBImage::writeImage(morphology_image_lowres, file = file.path, ...)
    } else {
      EBImage::writeImage(morphology_image_lowres, file = output.file, ...)
    }
  }
  invisible()
}

####
## Visium ####
####

#' importVisium
#'
#' Importing Visium data
#'
#' @param dir.path path to Visium output folder
#' @param selected_assay selected assay from Visium
#' @param assay_name the assay name
#' @param sample_name the name of the sample
#' @param image_name the image name of the Visium assay, Default: main
#' @param channel_name the channel name of the image of the Visium assay, Default: H&E
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#' @param resolution_level the level of resolution of Visium image: "lowres" (default) or "hires"
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read
#' @importFrom rjson fromJSON
#' @importFrom utils read.csv
#'
#' @export
#'
importVisium <- function(dir.path, selected_assay = "Gene Expression", assay_name = "Visium", sample_name = NULL, image_name = "main", channel_name = "H&E", inTissue = TRUE, resolution_level = "lowres", ...)
{
  # raw counts
  listoffiles <- list.files(dir.path)
  datafile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  datafile <- paste0(dir.path, "/", datafile)
  if(file.exists(datafile)){
    rawdata <- import10Xh5(filename = datafile)
    if(any(names(rawdata) %in% selected_assay)){
      rawdata <- as.matrix(rawdata[[selected_assay]])
    } else {
      stop("There are no assays called ", selected_assay, " in the h5 file!")
    }
  } else {
    stop("There are no files named 'filtered_feature_bc_matrix.h5' in the path")
  }

  # resolution
  if(!resolution_level %in% c("lowres","hires"))
    stop("resolution_level should be either 'lowres' or 'hires'!")

  # image
  image_file <- paste0(dir.path, paste0("/spatial/tissue_", resolution_level, "_image.png"))
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
      coords <- utils::read.csv(file = coords_file, header = FALSE)
      colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    } else {
      coords <- utils::read.csv(file = coords_file, header = TRUE)
    }
  } else if(length(coords_file) > 1) {
    stop("There are more than 1 position files in the path")
  } else {
    stop("There are no files named 'tissue_positions.csv' in the path")
  }

  if(inTissue){
    coords <- coords[coords$in_tissue==1,]
    rawdata <- rawdata[,colnames(rawdata) %in% coords$barcode]
  }
  coords <- coords[match(colnames(rawdata), coords$barcode),]
  spotID <- coords$barcode
  coords <- as.matrix(coords[,c("pxl_col_in_fullres", "pxl_row_in_fullres")], )
  colnames(coords) <- c("x", "y")
  rownames(coords) <- spotID

  # scale coordinates
  scale_file <- paste0(dir.path, "/spatial/scalefactors_json.json")
  if(file.exists(scale_file)){
    scalefactors <- rjson::fromJSON(file = scale_file)
    scales <- scalefactors[[paste0("tissue_", resolution_level, "_scalef")]]
    # spot.radius is the half of the diameter, but we visualize by a factor of 1.5 larger
    # params <- list(spot.radius = scalefactors$spot_diameter_fullres*scalefactors$tissue_lowres_scalef*1.5)
    params <- list(
      nearestpost.distance = 200*scales, # distance to nearest spot
      spot.radius = scalefactors$spot_diameter_fullres*scales,
      vis.spot.radius = scalefactors$spot_diameter_fullres*scales*2)
    coords <- coords*scales
    coords[,2] <- info$height - coords[,2]
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }

  # create VoltRon
  formVoltRon(rawdata, metadata = NULL, image, coords, main.assay = assay_name, params = params, assay.type = "spot", image_name = image_name, main_channel = channel_name, sample_name = sample_name, ...)
}

####
## Auxiliary ####
####

#' import10Xh5
#'
#' import the sparse matrix from the H5 file
#'
#' @param filename the path tp h5 file
#'
#' @importFrom hdf5r H5File readDataSet
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
#'
import10Xh5 <- function(filename){

  # check file
  if(file.exists(filename)){
    input.file <- hdf5r::H5File$new(filename = filename, mode = "r")
    matrix.10X <- input.file[["matrix"]]
  } else {
    stop("There are no files named ", filename," in the path")
  }

  # get data, barcodes and feature
  features <- hdf5r::readDataSet(matrix.10X[["features"]][["name"]])
  feature_type <- hdf5r::readDataSet(matrix.10X[["features"]][["feature_type"]])
  cells <- hdf5r::readDataSet(matrix.10X[["barcodes"]])
  mat <- hdf5r::readDataSet(matrix.10X[["data"]])
  indices <- hdf5r::readDataSet(matrix.10X[["indices"]])
  indptr <- hdf5r::readDataSet(matrix.10X[["indptr"]])
  sparse.mat <- sparseMatrix(i = indices + 1, p = indptr,
                             x = as.numeric(mat), dims = c(length(features), length(cells)),
                             repr = "T")
  colnames(sparse.mat) <- cells
  rownames(sparse.mat) <- make.unique(features)

  # separate feature types
  matrix.10X <- list()
  for(feat in unique(feature_type)){
    cur_features <- features[feature_type %in% feat]
    cur_mat <- sparse.mat[features %in% cur_features,]
    matrix.10X[[feat]] <- cur_mat
  }

  return(matrix.10X)
}

####
# Nanostring ####
####

####
## GeoMx ####
####

#' importGeoMx
#'
#' Import GeoMx data
#'
#' @param dcc.path path to the folder where the dcc files are found
#' @param pkc.file path to the pkc file
#' @param summarySegment the metadata csv (sep = ";") or excel file, if the file is an excel file, \code{summarySegmentSheetName} should be provided as well.
#' @param summarySegmentSheetName the sheet name of the excel file, \code{summarySegment}
#' @param assay_name the assay name, default: GeoMx
#' @param image the reference morphology image of the GeoMx assay
#' @param segment_polygons if TRUE, the ROI polygons are parsed from the OME.TIFF file
#' @param ome.tiff the OME.TIFF file of the GeoMx experiment if exists
#' @param resolution_level the level of resolution within GeoMx OME-TIFF image, Default: 3
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom dplyr %>% full_join
#' @importFrom utils read.csv
#' @importFrom magick image_info image_read
#'
#' @export
#'
importGeoMx <- function(dcc.path, pkc.file, summarySegment, summarySegmentSheetName, assay_name = "GeoMx",
                        image = NULL, segment_polygons = FALSE, ome.tiff = NULL, resolution_level = 3, ...)
{
  # Get pkc file
  if(file.exists(pkc.file)){
    pkcdata <- readPKC(pkc.file)
  } else {
    stop("pkc file is not found!")
  }

  # Get dcc file
  if(file.exists(dcc.path)){
    dcc_files <- dir(dcc.path, pattern = ".dcc$", full.names = TRUE)
    if(length(dcc_files) == 0){
      stop("no dcc files are found under ", dcc.path)
    } else {
      dcc_files <- dcc_files[!grepl("A01.dcc$", dcc_files)]
      dcc_filenames <- dir(dcc.path, pattern = ".dcc$", full.names = FALSE)
      dcc_filenames <- dcc_filenames[!grepl("A01.dcc$", dcc_filenames)]
      dcc_filenames <- gsub(".dcc$", "", dcc_filenames)
      dcc_filenames <- gsub("-", "_", dcc_filenames)
      dccData <- sapply(dcc_files, readDCC, simplify = FALSE, USE.NAMES = FALSE)
      names(dccData) <- dcc_filenames
    }
  } else {
    stop("path to dcc files does not exist!")
  }

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

  # get negative probes and targets
  # NegProbes <- pkcdata$RTS_ID[pkcdata$Target == "NegProbe-WTX"]
  NegProbes <- pkcdata$RTS_ID[grepl("NegProbe", pkcdata$Target)]

  # negative probes
  rawdata_neg <- rawdata[rownames(rawdata) %in% NegProbes, ]
  rownames(rawdata_neg) <- paste0(rownames(rawdata_neg), "_",
                                  pkcdata$Target[match(rownames(rawdata_neg), pkcdata$RTS_ID)])
  rawdata_neg <- as.matrix(rawdata_neg)

  # other probes
  rawdata <- rawdata[!rownames(rawdata) %in% NegProbes, ]
  rownames(rawdata) <- pkcdata$Target[match(rownames(rawdata), pkcdata$RTS_ID)]
  rawdata <- as.matrix(rawdata)

  # get segment summary
  if(file.exists(summarySegment)){
    if(grepl(".xls$|.xlsx$", summarySegment)){
      if (!requireNamespace('xlsx'))
        stop("Please install xlsx package for using the read.xlsx function!")
      if(!is.null(summarySegmentSheetName)){
        segmentsummary <- xlsx::read.xlsx(summarySegment, sheetName = summarySegmentSheetName)
      } else {
        stop("Please provide 'summarySegmentSheetName' for the excel sheet name!")
      }
    } else if(grepl(".csv$", summarySegment)) {
      segmentsummary <- utils::read.csv(summarySegment, row.names = NULL, header = T, sep = ";")
    }
    rownames(segmentsummary) <- gsub(".dcc$", "", segmentsummary$Sample_ID)
    rownames(segmentsummary) <- gsub("-", "_", rownames(segmentsummary))
    if(all(dcc_filenames %in% rownames(segmentsummary))){
      segmentsummary <- segmentsummary[dcc_filenames, ]
    } else{
      stop("Some GeoMx dcc files are not represented in the segment summary file!")
    }
  } else {
    stop(summarySegment, " is not found!")
  }

  # get image
  if(!is.null(image)){
    if(file.exists(image)){
      image <- magick::image_read(image)
    } else {
      stop(image, " is not found!")
    }
  } else {
    stop("Please provide a image for the GeoMx data!")
  }
  geomx_image_info <- magick::image_info(image)

  # get coordinates
  coords <- segmentsummary[,c("X","Y")]
  colnames(coords) <- c("x", "y")
  rownames(coords) <- segmentsummary$ROI.name
  coords <- rescaleGeoMxPoints(coords, segmentsummary, geomx_image_info)

  # get ROI segments (polygons)
  segments <- list()
  if(!is.null(ome.tiff)){

    # get segments
    segments <- importGeoMxSegments(ome.tiff, segmentsummary, geomx_image_info)
    segments <- segments[rownames(coords)]

    # parse channels if ome.tiff is given
    channels <- importGeoMxChannels(ome.tiff, segmentsummary, geomx_image_info, resolution_level = resolution_level)
    image <- c(list(scanimage = image), channels)
  }

  # create VoltRon for non-negative probes
  object <- formVoltRon(rawdata, metadata = segmentsummary, image, coords, segments, main.assay = assay_name, assay.type = "ROI", ...)

  # add negative probe assay
  new_assay <- formAssay(data = rawdata_neg,
                         coords = coords, segments = segments,
                         image = image, type = "ROI")
  new_assay@image <- object[["Assay1"]]@image
  sample.metadata <- SampleMetadata(object)
  object <- addAssay(object,
                     assay = new_assay,
                     metadata = Metadata(object),
                     assay_name = paste(sample.metadata["Assay1", "Assay"], "NegProbe", sep = "_"),
                     sample = sample.metadata["Assay1", "Sample"],
                     layer = sample.metadata["Assay1", "Layer"])

  # add connectivity of spatial points across assays
  connectivity <- cbind(vrSpatialPoints(object, assay = assay_name),
                        vrSpatialPoints(object, assay = paste(assay_name, "NegProbe", sep = "_")))
  object <- addConnectivity(object,
                            connectivity = connectivity,
                            sample = sample.metadata["Assay1", "Sample"],
                            layer = sample.metadata["Assay1", "Layer"])

  # return
  return(object)
}

#' readPKC
#'
#' Read a NanoString Probe Kit Configuration (PKC) file. Adapted from \code{GeomxTools} package.
#'
#' @param file A character string containing the path to the PKC file.
#' @param default_pkc_vers Optional list of pkc file names to use as default if more than one pkc version of each module is provided.
#'
#' @importFrom utils capture.output
#' @importFrom rjson fromJSON
#' @importFrom S4Vectors metadata DataFrame
#'
#' @noRd
readPKC <- function (file, default_pkc_vers = NULL)
{
  pkc_json_list <- lapply(file, function(pkc_file) {
    rjson::fromJSON(file = pkc_file)
  })
  pkc_names <- unlist(lapply(file, function(pkc_file) {
    base_pkc_name <- gsub(".pkc", "", trimws(basename(pkc_file)))
    return(base_pkc_name)
  }))
  names(pkc_json_list) <- pkc_names
  pkc_modules <- basename(unlist(lapply(pkc_names, sub, pattern = "_[^_]+$",
                                        replacement = "")))
  names(pkc_modules) <- pkc_names
  header <- list(PKCFileName = sapply(pkc_json_list, function(list) list[["Name"]]),
                 PKCModule = pkc_modules,
                 PKCFileVersion = sapply(pkc_json_list, function(list) list[["Version"]]),
                 PKCFileDate = sapply(pkc_json_list, function(list) list[["Date"]]),
                 AnalyteType = sapply(pkc_json_list, function(list) list[["AnalyteType"]]),
                 MinArea = sapply(pkc_json_list, function(list) list[["MinArea"]]),
                 MinNuclei = sapply(pkc_json_list, function(list) list[["MinNuclei"]]))

  multi_idx <- duplicated(header[["PKCModule"]])
  multi_mods <- unique(header[["PKCModule"]][multi_idx])
  if (length(multi_mods) < 1) {
    if (!is.null(default_pkc_vers)) {
      warning("Only one version found per PKC module. ",
              "No PKCs need to be combined. ", "Therefore, no default PKC versions will be used.")
    }
  }
  else {
    use_pkc_names <- lapply(multi_mods, function(mod) {
      mod_idx <- header[["PKCModule"]] == mod
      max_vers <- as.numeric(as.character(max(as.numeric_version(header[["PKCFileVersion"]][mod_idx]))))
      max_name <- names(header[["PKCFileVersion"]][header[["PKCFileVersion"]] == max_vers])
      return(max_name)
    })
    names(use_pkc_names) <- multi_mods
    if (!is.null(default_pkc_vers)) {
      default_names <- unlist(lapply(default_pkc_vers, function(pkc_file) {
        base_pkc_name <- gsub(".pkc", "", trimws(basename(pkc_file)))
        return(base_pkc_name)
      }))
      default_mods <- unlist(lapply(default_names, sub, pattern = "_[^_]+$",
                                replacement = ""))
      dup_defaults <- default_names[duplicated(default_mods) |
                                      duplicated(default_mods, fromLast = TRUE)]
      if (!all(default_names %in% names(header[["PKCFileName"]]))) {
        removed_pkcs <- default_pkc_vers[!default_names %in%
                                           names(header[["PKCFileName"]])]
        stop("Could not match all default PKC versions with a PKC file name. ",
             "Check default file names match exactly to a PKC file name.\n",
             paste0("Unmatched default PKC versions: ",
                    removed_pkcs))
      }
      else if (length(dup_defaults) > 0) {
        stop("There should only be one default PKC version per module. ",
             "Ensure only one version per module in default PKCs list.\n",
             "Multiple default PKC version conflicts: ",
             paste(dup_defaults, collapse = ", "))
      }
      else {
        use_pkc_names[default_mods] <- default_names
      }
    }
  }
  rtsid_lookup_df <- generate_pkc_lookup(pkc_json_list)
  rtsid_lookup_df$Negative <- grepl("Negative", rtsid_lookup_df$CodeClass)
  rtsid_lookup_df$RTS_ID <- gsub("RNA", "RTS00", rtsid_lookup_df[["RTS_ID"]])
  rtsid_lookup_df <- S4Vectors::DataFrame(rtsid_lookup_df)
  if (length(multi_mods) > 0) {
    for (mod in names(use_pkc_names)) {
      mod_vers <- names(header[["PKCModule"]])[header[["PKCModule"]] ==
                                                 mod]
      mod_lookup <- rtsid_lookup_df[rtsid_lookup_df$Module %in%
                                      mod_vers, ]
      mod_tab <- table(mod_lookup$RTS_ID)
      remove_rts <- names(mod_tab[mod_tab != length(mod_vers)])
      if (length(remove_rts) > 0) {
        warning("The following probes were removed from analysis",
                " as they were not found in all PKC module versions used.\n",
                paste(utils::capture.output(print(subset(rtsid_lookup_df,
                                                  subset = RTS_ID %in% remove_rts))), collapse = "\n"))
        rtsid_lookup_df <- subset(rtsid_lookup_df, subset = !RTS_ID %in%
                                    remove_rts)
      }
      remove_vers <- mod_vers[mod_vers != use_pkc_names[mod]]
      rtsid_lookup_df <- subset(rtsid_lookup_df, subset = !Module %in%
                                  remove_vers)
      warning("The following PKC versions were removed from analysis",
              " as they were overridden by a newer PKC version or",
              " were overridden by a user-defined default PKC version.\n",
              paste(remove_vers, collapse = ", "))
      header <- lapply(header, function(elem) {
        elem[!names(elem) %in% remove_vers]
      })
    }
  }
  S4Vectors::metadata(rtsid_lookup_df) <- header
  return(rtsid_lookup_df)
}

#' readDCC
#'
#' Read a NanoString GeoMx Digital Count Conversion (DCC) file.
#'
#' @param file A character string containing the path to the DCC file.
#'
#' @noRd
readDCC <- function(file)
{
  lines <- trimws(readLines(file))
  trimGalore <- grep("trimGalore", lines)
  if (length(trimGalore) > 0) {
    Raw <- grep("Raw", lines)
    lines <- lines[-c(trimGalore:(Raw - 1))]
  }
  lines <- gsub("SoftwareVersion,\"GeoMx_NGS_Pipeline_ ",
                "SoftwareVersion,", lines)
  lines <- gsub("SoftwareVersion,\"GeoMx_NGS_Pipeline_", "SoftwareVersion,",
                lines)
  lines <- gsub("SoftwareVersion,DRAGEN_GeoMx_", "SoftwareVersion,",
                lines)
  lines[grepl("SoftwareVersion", lines)] <- gsub("\"", "",
                                                 lines[grepl("SoftwareVersion", lines)])
  tags <- c("Header", "Scan_Attributes", "NGS_Processing_Attributes", "Code_Summary")
  output <- sapply(tags, function(tag) {
    bounds <- charmatch(sprintf(c("<%s>", "</%s>"), tag),
                        lines)
    if (anyNA(bounds) || bounds[1L] + 1L >= bounds[2L])
      lines[integer(0)]
    else lines[(bounds[1L] + 1L):(bounds[2L] - 1L)]
  }, simplify = FALSE)
  for (tag in c("Header", "Scan_Attributes", "NGS_Processing_Attributes")) {
    while (length(bad <- grep(",", output[[tag]], invert = TRUE)) >
           0L) {
      bad <- bad[1L]
      if (bad == 1L)
        stop(sprintf("%s section has malformed first line",
                     tag))
      fixed <- output[[tag]]
      fixed[bad - 1L] <- sprintf("%s %s", fixed[bad -
                                                  1L], fixed[bad])
      output[[tag]] <- fixed[-bad]
    }
    output[[tag]] <- strsplit(output[[tag]], split = ",")
    output[[tag]] <- structure(lapply(output[[tag]], function(x) if (length(x) ==
                                                                     1L)
      ""
      else x[2L]), names = lapply(output[[tag]], `[`, 1L),
      class = "data.frame", row.names = basename(file))
  }
  cols <- c("FileVersion", "SoftwareVersion")
  if (!(all(cols %in% colnames(output[["Header"]]))))
    stop("Header section must contain \"FileVersion\" and \"SoftwareVersion\"")
  output[["Header"]][, cols] <- lapply(output[["Header"]][,
                                                          cols], numeric_version)
  fileVersion <- output[["Header"]][1L, "FileVersion"]
  if (!(numeric_version(fileVersion) %in% numeric_version(c("0.01",
                                                            "0.02"))))
    stop("\"FileVersion\" in Header section must be 0.01 or 0.02")
  for (section in c("Header", "Scan_Attributes", "NGS_Processing_Attributes")) {
    valid <- .validDccSchema(output[[section]], fileVersion,
                             section)
    if (!isTRUE(valid))
      stop(valid)
  }
  output[["Header"]][["Date"]] <- as.Date(output[["Header"]][["Date"]],
                                          format = "%Y-%m-%d")
  cols <- c("Raw", "Trimmed", "Stitched", "Aligned", "umiQ30",
            "rtsQ30")
  output[["NGS_Processing_Attributes"]][, cols] <- lapply(output[["NGS_Processing_Attributes"]][,
                                                                                                cols], as.numeric)
  names(output[["Scan_Attributes"]])[names(output[["Scan_Attributes"]]) ==
                                       "ID"] <- "SampleID"
  output[["Code_Summary"]] <- paste0("RTS_ID,Count\n", paste(output[["Code_Summary"]],
                                                             collapse = "\n"))
  output[["Code_Summary"]] <- utils::read.csv(textConnection(output[["Code_Summary"]]),
                                              colClasses = c(RTS_ID = "character", Count = "numeric"))
  output[["Code_Summary"]][["Count"]] <- as.integer(round(output[["Code_Summary"]][["Count"]]))
  rn <- output[["Code_Summary"]][["RTS_ID"]]
  if ((ndups <- anyDuplicated(rn)) > 0L) {
    warning(sprintf("removed %d rows from \"Code_Summary\" due to duplicate rownames",
                    ndups))
    ok <- which(!duplicated(rn, fromLast = FALSE) & !duplicated(rn,
                                                                fromLast = TRUE))
    rn <- rn[ok]
    output[["Code_Summary"]] <- output[["Code_Summary"]][ok,
                                                         , drop = FALSE]
  }
  rownames(output[["Code_Summary"]]) <- rn
  output[["NGS_Processing_Attributes"]][, "DeduplicatedReads"] <- sum(output[["Code_Summary"]][["Count"]])
  return(output)
}


#' importGeoMxSegments
#'
#' Import ROI polygons from the OME.TIFF file
#'
#' @param ome.tiff the OME.TIFF file of the GeoMx Experiment
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#'
#' @noRd
importGeoMxSegments <- function(ome.tiff, summary, imageinfo){

  # check file
  if(file.exists(ome.tiff)){
    if(grepl(".ome.tiff$|.ome.tif$", ome.tiff)){
      if (!requireNamespace('RBioFormats'))
        stop("Please install RBioFormats package to extract xml from the ome.tiff file!")
      if (!requireNamespace('XML'))
        stop("Please install XML package to extract xml from the ome.tiff file!")
      omexml <- RBioFormats::read.omexml(ome.tiff)
    } else if(grepl(".xml$", ome.tiff)){
      omexml <- XML::xmlParse(file = ome.tiff)
    } else {
      stop("Please provide either an ome.tiff or .xml file!")
    }
    omexml <- XML::xmlToList(omexml, simplify = TRUE)
  } else {
    stop("There are no files named ", ome.tiff," in the path")
  }

  # get ROIs
  ROIs <- omexml[which(names(omexml) == "ROI")]

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
#' Rescale GeoMx point (center or polygon corners of ROI) coordinates
#'
#' @param pts coordinates of the cells from the Xenium assays
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#'
#' @noRd
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

#' importGeoMxChannels
#'
#' Rescale GeoMx channels with respect to the scan image
#'
#' @param ome.tiff the OME.TIFF file of the GeoMx Experiment
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#' @param resolution_level the resolution level (1-7) of the image parsed from the OME.TIFF file
#'
#' @param RBioFormats read.image
#' @param EBImage as.Image
#' @param grDevices as.raster
#' @param magick image_read
#'
#' @noRd
importGeoMxChannels <- function(ome.tiff, summary, imageinfo, resolution_level){

  # check file
  if(file.exists(ome.tiff)){
    if(grepl(".ome.tiff$|.ome.tif$", ome.tiff)){
      if (!requireNamespace('RBioFormats'))
        stop("Please install RBioFormats package to extract xml from the ome.tiff file!")
      if (!requireNamespace('XML'))
        stop("Please install XML package to extract xml from the ome.tiff file!")
      omexml <- RBioFormats::read.omexml(ome.tiff)
      omexml <- XML::xmlToList(omexml, simplify = TRUE)
    } else {
      warning("ome.tiff format not found!")
      return(NULL)
    }
  } else {
    stop("There are no files named ", ome.tiff," in the path")
  }

  # get all channels
  ome.tiff <- RBioFormats::read.image(ome.tiff, resolution = resolution_level)

  # get frame information
  nframes <- ome.tiff@metadata$coreMetadata$imageCount
  frames <- EBImage::getFrames(ome.tiff)
  frames <- lapply(frames, function(x){
    img <- magick::image_read(grDevices::as.raster(EBImage::as.Image(x)))
    rescaleGeoMxImage(img, summary, imageinfo, resolution_level = resolution_level)
  })

  # get channel names
  omexml <- omexml$StructuredAnnotations[seq(from = 1, by = 2, length.out = nframes)]
  channel_names <- lapply(omexml, function(x){
    x$Value$ChannelInfo$BiologicalTarget
  })
  names(frames) <- channel_names

  # return frames
  return(frames)
}

#' rescaleGeoMxImage
#'
#' Rescale GeoMx channels with respect to the scan image
#'
#' @param img coordinates of the cells from the Xenium assays
#' @param summary segmentation summary data frame
#' @param imageinfo image information
#' @param resolution_level the resolution level (1-7) of the image parsed from the OME.TIFF file
#'
#' @param magick image_crop
#'
#' @noRd
rescaleGeoMxImage <- function(img, summary, imageinfo, resolution_level){

  # adjust offset and scan size to the resolution level
  offset.x <- summary$Scan.Offset.X[1]/(2*(resolution_level-1))
  offset.y <- summary$Scan.Offset.Y[1]/(2*(resolution_level-1))
  scan.width <- summary$Scan.Width[1]/(2*(resolution_level-1))
  scan.height <- summary$Scan.Height[1]/(2*(resolution_level-1))

  # crop image given adjusted offsets
  new_img <- magick::image_crop(img, geometry = paste0(scan.width, "x", scan.height, "+", offset.x, "+", offset.y))

  # resize image
  new_img <- magick::image_resize(new_img, geometry = paste0(imageinfo$width, "x", imageinfo$height))

  # return
  return(new_img)
}

####
## CosMx ####
####

#' importCosMx
#'
#' Import CosMx data
#'
#' @param tiledbURI the path to the tiledb folder
#' @param assay_name the assay name, default: CosMx
#' @param image the reference morphology image of the CosMx assay
#' @param image_name the image name of the CosMx assay, Default: main
#' @param ome.tiff the OME.TIFF file of the CosMx experiment if exists
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom data.table data.table
#' @importFrom ids random_id
#'
#' @export
#'
importCosMx <- function(tiledbURI, assay_name = "CosMx",
                        image = NULL, image_name = "main", ome.tiff = NULL, import_molecules = FALSE, ...)
{
  # check tiledb and tiledbsc
  if (!requireNamespace("tiledb", quietly = TRUE))
    stop("Please install the tiledb package: \n
         remotes::install_github('TileDB-Inc/TileDB-R', force = TRUE,
                            ref = '0.17.0')")
  if (!requireNamespace("tiledbsc", quietly = TRUE))
    stop("Please install the tiledbsc package: \n
         remotes::install_github('tiledb-inc/tiledbsc', force = TRUE,
                            ref = '8157b7d54398b1f957832f37fff0b173d355530e')")

  # get tiledb
  message("Scanning TileDB array for cell data ...")
  tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, verbose = FALSE)

  # raw counts
  counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE)
  counts <- as.matrix(counts)

  # cell metadata
  metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()

  # coordinates
  coords <- as.matrix(metadata[,c("x_slide_mm", "y_slide_mm")])
  colnames(coords) <- c("x","y")

  # transcripts
  if(import_molecules){
    message("Scanning TileDB array for molecule data ...")
    subcellular <- tiledb::tiledb_array(
      tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
      return_as="data.table")[]
    colnames(subcellular)[colnames(subcellular)=="target"] <- "gene"
  }

  # get slides and construct VoltRon objects for each slides
  slides <- unique(metadata$slide_ID_numeric)

  # for each slide create a VoltRon object with combined layers
  vr_list <- list()
  for(slide in slides){

    # cell assay
    message("Creating cell level assay for slide ", slide, " ...")

    # slide info
    cur_coords <- coords[metadata$slide_ID_numeric == slide,]
    cur_counts <- counts[,rownames(cur_coords)]
    cur_metadata <- metadata[rownames(cur_coords),]

    # create VoltRon object
    cell_object <- formVoltRon(data = cur_counts, metadata = cur_metadata, image = image, coords = cur_coords, main.assay = assay_name, assay.type = "cell", image_name = image_name, ...)
    cell_object$Sample <- paste0("Slide", slide)

    # molecule assay
    if(import_molecules){

      # get slide
      message("Creating molecule level assay for slide ", slide, " ...")
      cur_subcellular <- subset(subcellular, slideID == slide)

      # coordinates
      mol_coords <- as.matrix(cur_subcellular[,c("x_FOV_px", "y_FOV_px")])
      colnames(mol_coords) <- c("x", "y")

      # get subcellular data components
      mol_metadata <- cur_subcellular[,colnames(cur_subcellular)[!colnames(cur_subcellular) %in% c("CellId", "cell_id", "x_FOV_px", "y_FOV_px")], with = FALSE]
      set.seed(nrow(mol_metadata))
      # entity_ID <- paste0(1:nrow(mol_metadata), "_", ids::random_id(bytes = 3, use_openssl = FALSE))
      entity_ID <- paste0(rownames(mol_metadata), "_", ids::random_id(bytes = 3, use_openssl = FALSE))
      mol_metadata <- data.table::data.table(id = entity_ID, assay_id = "Assay1", mol_metadata)

      # coord names
      rownames(mol_coords) <- entity_ID

      # create VoltRon assay for molecules
      mol_assay <- formAssay(coords = mol_coords, image = image, type = "molecule", main_image = image_name)

      # merge assays in one section
      message("Merging assays for slide ", slide, " ...")
      sample.metadata <- SampleMetadata(cell_object)
      cell_object <- addAssay(cell_object,
                              assay = mol_assay,
                              metadata = mol_metadata,
                              assay_name = paste0(assay_name, "_mol"),
                              sample = sample.metadata["Assay1", "Sample"],
                              layer = sample.metadata["Assay1", "Layer"])
    }
    vr_list <- append(vr_list, cell_object)
  }

  # return
  message("Merging slides ...")
  vr <- merge(vr_list[[1]], vr_list[-1])
}

#' generateCosMxImage
#'
#' Generates a low resolution Morphology image of the CosMx experiment
#'
#' @param dir.path CosMx output folder
#' @param increase.contrast increase the contrast of the image before writing
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param ... additional parameters passed to the \link{writeImage} function
#'
#' @importFrom magick image_read image_contrast
#' @importFrom EBImage writeImage
#' @importFrom reshape2 acast
#'
#' @export
#'
generateCosMxImage <- function(dir.path, increase.contrast = TRUE, output.path = NULL, ...) {

  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/CellComposite_lowres.tif")
  output.file <- paste0(output.path, "/CellComposite_lowres.tif")

  # check if the file exists in either Xenium output folder, or the specified location
  if(file.exists(file.path) | file.exists(paste0(output.file))){
    cat("CellComposite_lowres.tif already exists! \n")
    return(NULL)
  }

  # FOV positions of CosMx
  list_of_files <- list.files(dir.path)
  fov_positions_path <- paste0(dir.path, "/", list_of_files[grepl("fov_positions_file.csv$",list_of_files)][1])
  fov_positions <- read.csv(fov_positions_path)

  # manipulate fov positions matrix
  cat("Getting FOV Positions \n")
  relative_fov_positions <- fov_positions
  x_min <- min(relative_fov_positions$x_global_px)
  y_min <- min(relative_fov_positions$y_global_px)
  x_gap <- diff(unique(fov_positions$x_global_px))[1]
  y_gap <- diff(unique(fov_positions$y_global_px))[1]
  relative_fov_positions[,c("x_global_px","y_global_px")] <- t(apply(relative_fov_positions[,c("x_global_px","y_global_px")], 1, function(cell){
    c((cell[1]-x_min)/x_gap,(cell[2]-y_min)/y_gap)
  }))
  relative_fov_positions <- relative_fov_positions[order(relative_fov_positions$y_global_px, decreasing = TRUE),]

  # Combine Images of the FOV grid
  cat("Loading FOV tif files \n")
  image.dir.path <- paste0(dir.path,"/CellComposite/")
  morphology_image_data <- NULL
  for(i in relative_fov_positions$fov){
    image_path <- paste0(image.dir.path, "CellComposite_F", str_pad(as.character(i), 3, pad = 0), ".jpg")
    image_data <- magick::image_read(image_path) %>% magick::image_resize("x500") %>% magick::image_raster()
    if(is.null(morphology_image_data))
      dim_image <- apply(image_data[,1:2], 2, max)
    scale_dim <- relative_fov_positions[i,2:3]*dim_image
    image_data[,1:2] <- image_data[,1:2] +
      rep(1, nrow(image_data)) %o% as.matrix(scale_dim)[1,]
    morphology_image_data <- rbind(morphology_image_data, image_data)
  }
  morphology_image_data_array <- reshape2::acast(morphology_image_data, y ~ x)
  morphology_image <- magick::image_read(morphology_image_data_array) %>% magick::image_resize("x800")

  # pick a resolution level
  morphology_image_info <- image_info(morphology_image)
  cat(paste0("Image Resolution (X:", morphology_image_info$width, " Y:", morphology_image_info$height, ") \n"))

  # increase contrast using EBImage
  if(increase.contrast) {
    cat("Increasing Contrast \n")
    morphology_image <- magick::image_contrast(morphology_image, sharpen = 1)
  }

  # write to the same folder
  cat("Writing Tiff File \n")
  if(is.null(output.path)){
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = file.path, ...)
  } else {
    EBImage::writeImage(magick::as_EBImage(morphology_image), file = output.file, ...)
  }

  return(NULL)
}

####
# Image Data ####
####

#' importImageData
#'
#' import an image as VoltRon object
#'
#' @param image the path to an image file
#' @param tile.size the size of tiles
#' @param stack.id the id of the stack when the magick image composed of multiple layers
#' @param segments Either a list of segments or a GeoJSON file. This will result in a second assay in the VoltRon object to be created
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom data.table data.table
#'
#' @export
#'
importImageData <- function(image, tile.size = 10, stack.id = 1, segments = NULL, ...){

  # get image
  if(!inherits(image, "magick-image")){
    if(!is.character(image)){
      stop("image should either be a magick-image object or a file.path")
    } else{
      if(file.exists(image)){
        image <- magick::image_read(image)
      } else {
        stop(image, " is not found!")
      }
    }
  }

  # get image layer from stacked magick images
  if(length(image) > 1){
    if(stack.id > length(image)){
      stop("The stack.id should be an integer between 1 and ", length(image))
    } else {
      image <- image[stack.id]
    }
  }

  # image info
  imageinfo <- magick::image_info(image)

  # coordinates
  even_odd_correction <- (!tile.size%%2)*(0.5)
  # x_coords <- seq((tile.size/2) + even_odd_correction, length.out = imageinfo$width %/% tile.size)
  # y_coords <- seq(imageinfo$height %/% tile.size, 1)*(tile.size/2)
  x_coords <- seq((tile.size/2) + even_odd_correction, imageinfo$width, tile.size)[1:(imageinfo$width %/% tile.size)]
  y_coords <- seq((tile.size/2) + even_odd_correction, imageinfo$height, tile.size)[1:(imageinfo$height %/% tile.size)]
  y_coords <- rev(y_coords)
  coords <- as.matrix(expand.grid(x_coords, y_coords))
  colnames(coords) <- c("x", "y")
  rownames(coords) <- paste0("tile", 1:nrow(coords))

  # metadata
  metadata <- data.table(id = rownames(coords))

  # create voltron object with tiles
  object <- formVoltRon(data = NULL, metadata = metadata, image = image, coords, main.assay = "ImageData", assay.type = "tile", params = list(tile.size = tile.size), ...)
  
  # check if segments are defined
  if(is.null(segments)){
    return(object)
  } else {
    
    # check if segments are paths
    if(inherits(segments, "character")){
      if(grepl(".geojson$", segments)){
        segments <- generateSegmentsFromGeoJSON(geojson.file = segments)
      } else {
        stop("Only lists or GeoJSON files are accepted as segments input!")
      }
    }
    
    # make coordinates out of segments
    coords <- t(sapply(segments, function(dat){
      apply(dat[,c("x", "y")], 2, mean)
    }))
    rownames(coords) <- names(segments)
    
    # make segment assay
    assay <- formAssay(coords = coords, segments = segments, image = image, type = "ROI") 
    
    # add segments as assay
    sample_metadata <- SampleMetadata(object)
    object <- addAssay(object,
                       assay = assay,
                       metadata = data.frame(row.names = paste0(names(segments), "_Assay2")),
                       assay_name = "ROIAnnotation",
                       sample = sample_metadata$Sample, 
                       layer = sample_metadata$Layer)
    
    # return
    return(object)
  }
}


#' generateSegmentsFromGeoJSON
#' 
#' The function to import segments from a json data
#'
#' @param geojson.file the GeoJSON file, typically generated by QuPath software
#'
#' @importFrom rjson fromJSON
#' @importFrom dplyr tibble
#' 
#' @export
generateSegmentsFromGeoJSON <- function(geojson.file){
  
  # get segments
  if(inherits(geojson.file, "character")){
    if(file.exists(geojson.file)){
      segments <- rjson::fromJSON(file = geojson.file)
    } else {
      stop("geojson.file doesn't exist!")
    }
  } else {
    stop("geojson.file should be the path to the GeoJSON file!")
  }
  
  # parse polygons as segments/ROIs
  segments <- lapply(segments, function(x){
    type <- x$geometry$type
    poly <- x$geometry$coordinates
    if(grepl("Polygon", type)){
      poly <- as.data.frame(matrix(unlist(poly[[1]]), ncol = 2, byrow = TRUE))
    }
    colnames(poly) <- c("x", "y")
    dplyr::tibble(poly)
  })
  
  # attach names to segments
  segments <- mapply(function(x, sgt){
    dplyr::tibble(data.frame(id = x, sgt))
  }, 1:length(segments), segments, SIMPLIFY = FALSE)
  
  # generate ROI names
  names(segments) <- paste0("ROI", 1:length(segments))

  # return
  return(segments)
}

#' generateSegmentsFromGeoJSON
#' 
#' The function to import segments from a json data
#'
#' @param segments the segments, typically from \link{vrSegments}.
#' @param geojson.file the GeoJSON file, typically to be used by QuPath software.
#'
#' @importFrom rjson fromJSON
#' 
#' @export
generateGeoJSONFromSegments <- function(segments, geojson.file){
  
  if(!requireNamespace('geojsonR'))
    stop("Please install geojsonR package for using geojsonR functions")
  
  # get segments
  if(!inherits(geojson.file, "character")){
    stop("geojson.file should be the path to the GeoJSON file!")
  } 
  
  # reshape segments
  segments <- mapply(function(id, sgt){
    poly <- as.list(data.frame(t(as.matrix(sgt[,c("x", "y")]))))
    names(poly) <- NULL
    init <- geojsonR::TO_GeoJson$new()
    geometry <- init$Polygon(list(poly), stringify = FALSE)
    feature <- list(id = id, geometry = geometry, properties = list(objectType = "annotation"))
    feature
  }, names(segments), segments, SIMPLIFY = FALSE, USE.NAMES = FALSE)

  # save as json
  segments <- rjson::toJSON(segments)
  write(segments, file = geojson.file)
  
  # return
  return(segments)
}