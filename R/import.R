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
#' @param verbose verbose
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @importFrom ids random_id
#'
#' @export
#'
importXenium <- function (dir.path, selected_assay = "Gene Expression", assay_name = "Xenium", sample_name = NULL, use_image = TRUE, 
                          morphology_image = "morphology_lowres.tif", resolution_level = 7, overwrite_resolution = TRUE, 
                          image_name = "main", channel_name = "DAPI", import_molecules = FALSE, verbose = TRUE, ...)
{
  # cell assay
  if(verbose)
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
    stop("There are no files named 'cell_feature_matrix.h5' in the path")
  }

  # image
  if(use_image){
    generateXeniumImage(dir.path, file.name = morphology_image, resolution_level = resolution_level, overwrite_resolution = overwrite_resolution,
                        verbose = verbose)
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
    segments <- segments[,c("cell_id", "vertex_x", "vertex_y")]
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
  cell_object <- formVoltRon(rawdata, metadata = NULL, image = image, coords, segments = segments, main.assay = assay_name, 
                             assay.type = "cell", image_name = image_name, main_channel = channel_name, sample_name = sample_name, 
                             feature_name = ifelse(selected_assay == "Gene Expression", "RNA", "main"), ...)

  # molecule assay
  if(!import_molecules){
    return(cell_object)
  } else {
    if(verbose)
      message("Creating molecule level assay ...")
    # transcripts
    transcripts_file <- paste0(dir.path, "/transcripts.parquet")
    if(!file.exists(transcripts_file)){
      stop("There are no file named 'transcripts.csv.gz' in the path")
    } else {
      if (!requireNamespace('arrow'))
        stop("Please install arrow package to extract molecule data!: install.packages('arrow')")
      
      # get subcellur data components
      # subcellular_data <- data.table::fread(transcripts_file)
      subcellular_data <- data.table::as.data.table(arrow::read_parquet(transcripts_file, as_data_frame = FALSE))
      subcellular_data <- subcellular_data[,c("transcript_id", colnames(subcellular_data)[!colnames(subcellular_data) %in% "transcript_id"]), with = FALSE]
      colnames(subcellular_data)[colnames(subcellular_data)=="transcript_id"] <- "id"
      colnames(subcellular_data)[colnames(subcellular_data)=="feature_name"] <- "gene"
      subcellular_data <- subcellular_data[subcellular_data$qv >= 20, ]

      # coordinates
      coords <- as.matrix(subcellular_data[,c("x_location", "y_location", "z_location")])
      colnames(coords) <- c("x", "y", "z")
      if(use_image){
        coords[,c("x", "y")] <- coords[,c("x", "y")]/scaleparam
      }
      coords[,"y"] <- range_coords[2] - coords[,"y"]  + range_coords[1]

      # metadata
      mol_metadata <- subcellular_data[,colnames(subcellular_data)[!colnames(subcellular_data) %in% c("cell_id", "transcript_id", "x_location", "y_location", "z_location")], with = FALSE]
      set.seed(nrow(mol_metadata$id))
      mol_metadata[, postfix:=paste0("_", ids::random_id(bytes = 3, use_openssl = FALSE))]
      mol_metadata[, assay_id:="Assay1"]
      mol_metadata[, id:=do.call(paste0,.SD), .SDcols=c("id", "postfix")]
      
      # coord names
      rownames(coords) <- mol_metadata$id

      # create VoltRon assay for molecules
      mol_assay <- formAssay(coords = coords, image = image, type = "molecule", main_image = image_name, main_channel = channel_name)

      # merge assays in one section
      if(verbose)
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
        connectivity[, cell_assay_id:="_Assay1"]
        connectivity[, cell_id:=do.call(paste0,.SD), .SDcols=c("cell_id", "cell_assay_id")]
        connectivity$cell_assay_id <- NULL
      }

      # add connectivity of spatial points across assays
      object <- addLayerConnectivity(object,
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
#' @param verbose verbose
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
generateXeniumImage <- function(dir.path, increase.contrast = TRUE, resolution_level = 7, overwrite_resolution = FALSE, 
                                output.path = NULL, file.name = "morphology_lowres.tif", verbose = TRUE, ...) {
  
  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/", file.name)
  output.file <- paste0(output.path, "/", file.name)
  
  # check if the file exists in either Xenium output folder, or the specified location
  if((file.exists(file.path) | file.exists(paste0(output.file))) & !overwrite_resolution){
    if(verbose)
      message(file.name, " already exists!")
  } else {
    if (!requireNamespace('RBioFormats'))
      stop("Please install RBioFormats package to extract xml from the ome.tiff file!: BiocManager::install('RBioFormats')")
    if(dir.exists(paste0(dir.path, "/morphology_focus"))){
      if(verbose)
        message("Loading morphology_focus_0000.ome.tif ...")
      morphology_image_lowres <- RBioFormats::read.image(paste0(dir.path, "/morphology_focus/morphology_focus_0000.ome.tif"),
                                                         resolution = resolution_level,
                                                         subset=list(C=1))
    } else if(file.exists(paste0(dir.path, "/morphology_mip.ome.tif"))) {
      if(verbose)
        message("Loading morphology_mip.ome.tif ...")
      morphology_image_lowres <- RBioFormats::read.image(paste0(dir.path, "/morphology_mip.ome.tif"), resolution = resolution_level)
    }
    
    # pick a resolution level
    image_info <- morphology_image_lowres@metadata$coreMetadata
    if(verbose)
      message("  Image Resolution (X:", image_info$sizeX, " Y:", image_info$sizeY, ") ...")
    
    # increase contrast using EBImage
    if(increase.contrast) {
      if(verbose)
        message("  Increasing Contrast ...")
      morphology_image_lowres <- (morphology_image_lowres/max(morphology_image_lowres))
    }
    
    # write to the same folder
    if(verbose)
      message("  Writing Tiff File ...")
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
#' @importFrom magick image_read image_info
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
  if(length(coords_file) > 0){
    if(length(coords_file) > 1) {
      message("There are more than 1 position files in the path, using the first!")
      coords_file <- coords_file[1]
    } 
    coords <- utils::read.csv(file = coords_file, header = FALSE)
    if("barcode" %in% as.vector(coords[1,,drop = TRUE])){
      coords <- utils::read.csv(file = coords_file, header = FALSE, skip = 1)
    } 
    colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    # if(grepl("tissue_positions_list.csv", coords_file)) {
    #   coords <- utils::read.csv(file = coords_file, header = FALSE)
    #   colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    # } else {
    #   coords <- utils::read.csv(file = coords_file, header = TRUE)
    # }
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
    params <- list(
      nearestpost.distance = 200*scales, # distance to nearest spot
      spot.radius = scalefactors$spot_diameter_fullres*scales,
      vis.spot.radius = scalefactors$spot_diameter_fullres*scales*2,
      spot.type = "circle")
    coords <- coords*scales
    coords[,2] <- info$height - coords[,2]
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }
  
  # create VoltRon
  formVoltRon(rawdata, metadata = NULL, image, coords, main.assay = assay_name, params = params, assay.type = "spot", 
              image_name = image_name, main_channel = channel_name, sample_name = sample_name, 
              feature_name = ifelse(selected_assay == "Gene Expression", "RNA", "main"), ...)
}

#' importVisiumHD
#'
#' Importing VisiumHD data
#'
#' @param dir.path path to Visium output folder
#' @param bin.size bin size of the VisiumHD output (Exp: "2", "8" and "16")
#' @param selected_assay selected assay from Visium
#' @param assay_name the assay name
#' @param sample_name the name of the sample
#' @param image_name the image name of the Visium assay, Default: main
#' @param channel_name the channel name of the image of the Visium assay, Default: H&E
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#' @param resolution_level the level of resolution of Visium image: "lowres" (default) or "hires"
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom rjson fromJSON
#' @importFrom utils read.csv
#'
#' @export
#'
importVisiumHD <- function(dir.path, bin.size = "8", selected_assay = "Gene Expression", assay_name = "VisiumHD", sample_name = NULL, image_name = "main", channel_name = "H&E", inTissue = TRUE, resolution_level = "lowres", ...){
  
  # raw counts
  bin.size <- formatC(as.numeric(bin.size), width = 3, format = "d", flag = "0")
  dir.path <- paste0(dir.path, "binned_outputs/square_", bin.size, "um/")
  listoffiles <- list.files(dir.path)
  datafile <- listoffiles[grepl("filtered_feature_bc_matrix.h5", listoffiles)][1]
  datafile <- paste0(dir.path, "/", datafile)
  if(file.exists(datafile)){
    rawdata <- import10Xh5(filename = datafile)
    if(any(names(rawdata) %in% selected_assay)){
      rawdata <- as(rawdata[[selected_assay]], "CsparseMatrix")
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
    info <- magick::image_info(image)
  } else {
    stop("There are no spatial image files in the path")
  }
  
  # coordinates
  if(!requireNamespace("arrow"))
    stop("Please install arrow package to import VisiumHD!: install.packages('arrow')")
  coords_file <- list.files(paste0(dir.path, "/spatial/"), full.names = TRUE)
  coords_file <- coords_file[grepl("tissue_positions",coords_file)]
  if(length(coords_file) == 1){
    coords <- data.table::as.data.table(arrow::read_parquet(coords_file, as_data_frame = FALSE))
  } else if(length(coords_file) > 1) {
    stop("There are more than 1 position files in the path")
  } else {
    stop("There are no files named 'tissue_positions.parquet' in the path")
  }
  if(inTissue){
    coords <- subset(coords, in_tissue == 1)
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
    params <- list(
      nearestpost.distance = scalefactors$spot_diameter_fullres*scales*(3/sqrt(2)),
      spot.radius = scalefactors$spot_diameter_fullres*scales,
      vis.spot.radius = scalefactors$spot_diameter_fullres*scales,
      spot.type = "square")
    coords <- coords*scales
    coords[,2] <- info$height - coords[,2]
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }
  
  # create VoltRon
  formVoltRon(rawdata, metadata = NULL, image, coords, main.assay = assay_name, params = params, assay.type = "spot", 
              image_name = image_name, main_channel = channel_name, sample_name = sample_name, 
              feature_name = ifelse(selected_assay == "Gene Expression", "RNA", "main"), ...)
}

####
## Auxiliary ####
####

#' import10Xh5
#'
#' import the sparse matrix from the H5 file
#'
#' @param filename the path to h5 file
#'
#' @importFrom Matrix sparseMatrix
#'
#' @export
import10Xh5 <- function(filename){

  # check package
  if(!requireNamespace("rhdf5")){
    stop("You have to install the rhdf5 package!: BiocManager::install('rhdf5')")
  }
  
  # check file
  if(file.exists(filename)){
    matrix.10X <- rhdf5::h5read(file = filename, name = "matrix")
  } else {
    stop("There are no files named ", filename," in the path")
  }

  # get data, barcodes and feature
  features <- matrix.10X[["features"]][["name"]]
  feature_type <- matrix.10X[["features"]][["feature_type"]]
  cells <- matrix.10X[["barcodes"]]
  mat <- matrix.10X[["data"]]
  indices <- matrix.10X[["indices"]]
  indptr <- matrix.10X[["indptr"]]
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
#' @param image_name the image name of the Visium assay, Default: main
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom dplyr %>% full_join
#' @importFrom utils read.csv
#' @importFrom magick image_info image_read
#'
#' @export
#'
importGeoMx <- function(dcc.path, pkc.file, summarySegment, summarySegmentSheetName, assay_name = "GeoMx",
                        image = NULL, segment_polygons = FALSE, ome.tiff = NULL, resolution_level = 3, image_name = "main", ...)
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
  for(i in seq_len(length(dccData))){
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
        stop("Please install xlsx package for using the read.xlsx function!: install.packages('xlsx')")
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
  object <- formVoltRon(rawdata, metadata = segmentsummary, image, coords, segments, main.assay = assay_name, assay.type = "ROI", 
                        feature_name = "RNA", image_name = image_name, ...)
  
  # add negative probe assay as new feature set
  object <- addFeature(object, assay = assay_name, data = rawdata_neg, feature_name = "NegProbe")

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
#' @importFrom rjson fromJSON
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
  # rtsid_lookup_df <- S4Vectors::DataFrame(rtsid_lookup_df)
  rtsid_lookup_df <- data.table::data.table(rtsid_lookup_df)
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
                )
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
  # S4Vectors::metadata(rtsid_lookup_df) <- header
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
        stop("Please install RBioFormats package to extract xml from the ome.tiff file!: BiocManager::install('RBioFormats')")
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
  for(i in seq_len(length(ROIs))){
    cur_ROI <- ROIs[[i]]

    # if the shape is a polygon
    if("Polygon" %in% names(cur_ROI$Union)){
      coords <- strsplit(cur_ROI$Union$Polygon, split = "\\n")
      coords <- strsplit(coords$Points, split = " ")[[1]]
      coords <- sapply(coords, function(x) as.numeric(strsplit(x, split = ",")[[1]]), USE.NAMES = FALSE)
      coords <- as.data.frame(t(coords))
      colnames(coords) <- c("x", "y")
      coords <- rescaleGeoMxPoints(coords, summary, imageinfo)
      coords <- data.frame(id = cur_ROI$Union$Label[["Text"]], coords)
      mask_lists[[cur_ROI$Union$Label[["Text"]]]] <- data.frame(coords)

    # if the shape is an ellipse
    } else if("Ellipse" %in% names(cur_ROI$Union)){
      coords <- as.numeric(cur_ROI$Union$Ellipse[c("X","Y", "RadiusX", "RadiusY")])
      coords <-  as.data.frame(matrix(coords, nrow = 1))
      colnames(coords) <- c("x", "y", "rx", "ry")
      coords[,c("x", "y")] <- rescaleGeoMxPoints(coords[,c("x", "y")], summary, imageinfo)
      coords$rx <- coords$rx * imageinfo$width/summary$Scan.Width[1]
      coords$ry <- coords$ry * imageinfo$height/summary$Scan.Height[1]
      coords <- data.frame(id = cur_ROI$Union$Label[["Text"]], coords)
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
        stop("Please install RBioFormats package to extract xml from the ome.tiff file!: BiocManager::install('RBioFormats')")
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
#' @param path the path to the tiledb folder
#' @param assay_name the assay name, default: CosMx
#' @param image the reference morphology image of the CosMx assay
#' @param image_name the image name of the CosMx assay, Default: main
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param verbose verbose
#' @param method the approach for importing the CosMx assay either by the folder of CSVs or with TileDB array.
#' @param feature_name the name/key of the feature set.
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @export
importCosMx <- function(path, 
                        assay_name = "CosMx",
                        image = NULL, 
                        image_name = "main", 
                        import_molecules = FALSE, 
                        verbose = TRUE, 
                        method = "CSV", 
                        feature_name = NULL, ...)
{
  if(method == "CSV"){
    vr <- importCosMxCSV(path = path, assay_name = assay_name,
                         image = image, image_name = image_name, 
                         import_molecules = import_molecules, 
                         feature_name = feature_name, verbose = verbose, ...)
  } else if(method == "TileDB"){
    stop("TileDB importer is currently deprecated!")
    # vr <- importCosMxTileDB(tiledbURI = path, assay_name = assay_name,
    #                         image = image, image_name = image_name, 
    #                         import_molecules = import_molecules, 
    #                         feature_name = feature_name, verbose = verbose, ...)
  } else {
    stop("method should be either 'CSV' or 'TileDB'!")
  }
  vr
}

#' importCosMx
#'
#' Import CosMx data
#'
#' @param path the path to the tiledb folder
#' @param assay_name the assay name, default: CosMx
#' @param image the reference morphology image of the CosMx assay
#' @param image_name the image name of the CosMx assay, Default: main
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param feature_name the name/key of the feature set
#' @param verbose verbose
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom data.table data.table
#' @importFrom ids random_id
#'
#' @noRd
importCosMxCSV <- function(path, 
                           assay_name = "CosMx",
                           image = NULL, 
                           image_name = "main", 
                           import_molecules = FALSE, 
                           feature_name = NULL,
                           verbose = TRUE, ...)
{
  list_of_files <- list.files(path, full.names = TRUE)

  # metadata
  metadata_file <- list_of_files[grepl("metadata_file", list_of_files)]
  if(file.exists(metadata_file)){
    metadata <- data.table::fread(metadata_file, data.table = FALSE)
    rownames(metadata) <- metadata$cell
  } else {
    stop("There are no file with mattern 'metadata_file' in the path")
  }

  # raw counts
  if(verbose)
    message("Reading cell data from CosMx folder ...")
  datafile <- list_of_files[grepl("exprMat_file", list_of_files)]
  if(file.exists(datafile)){
    tmp <- data.table::fread(datafile)
    counts <- t(as(as.matrix(tmp[,-(1:2)]), "dgCMatrix"))
    colnames(counts) <- metadata$cell
  } else {
    stop("There are no files with pattern 'exprMat_file' in the path")
  }
  
  # coordinates
  coords <- as.matrix(metadata[,c("CenterX_global_px", "CenterY_global_px")])
  colnames(coords) <- c("x","y")
  rownames(coords) <- metadata$cell

  # segments
  segments_file <- list_of_files[grepl("polygons", list_of_files)]
  if(file.exists(segments_file)){
    segments <- as.data.frame(data.table::fread(segments_file))
    segments <- segments[,c("cell", "x_global_px", "y_global_px")]
    colnames(segments) <- c("cell_id", "x", "y")
    segments <- segments %>% dplyr::group_split(cell_id)
    segments <- as.list(segments)
    sgt <- do.call(rbind, segments)
    names(segments) <- unique(sgt$cell_id)
  } else {
    stop("There are no file with mattern 'polygons' in the path")
  }
  
  # transcripts
  if(import_molecules){
    if(verbose)
      message("Reading molecule data from CosMx folder ...")
    molecule_file <- list_of_files[grepl("tx_file", list_of_files)]
    if(file.exists(molecule_file)){
      molecules <- data.table::fread(molecule_file)
      colnames(molecules)[colnames(molecules)=="target"] <- "gene"
    } else {
      stop("There are no file with mattern 'tx_file' in the path")
    }
  }
  
  # get slides and construct VoltRon objects for each slides
  slides <- unique(metadata$slide_ID)
  
  # for each slide create a VoltRon object with combined layers
  vr_list <- list()
  for(slide in slides){
    
    # cell assay
    if(verbose)
      message("Creating cell level assay for slide ", slide, " ...")
    
    # slide info
    ind <- metadata$slide_ID == slide
    cur_metadata <- metadata[ind,]
    cur_coords <- coords[ind,]
    cur_counts <- counts[,ind]
    cur_segments <- segments[rownames(cur_coords)]
    
    # create VoltRon object
    cell_object <- formVoltRon(data = cur_counts, metadata = cur_metadata, image = image, coords = cur_coords, segments = cur_segments,
                               main.assay = assay_name, assay.type = "cell", image_name = image_name, feature_name = feature_name, ...)
    cell_object$Sample <- paste0("Slide", slide)
    
    # molecule assay
    if(import_molecules){
      
      # get slide
      if(verbose)
        message("Creating molecule level assay for slide ", slide, " ...")
      if(!"slideID" %in% colnames(molecules)){
        cur_molecules <- molecules
      } else {
        cur_molecules <- subset(molecules, slideID == slide)
      }
      
      # coordinates
      mol_coords <- as.matrix(cur_molecules[,c("x_global_px", "y_global_px")])
      colnames(mol_coords) <- c("x", "y")
      
      # get molecules data components
      mol_metadata <- cur_molecules[,colnames(cur_molecules)[!colnames(cur_molecules) %in% c("CellId", "x_global_px", "y_global_px")], with = FALSE]
      set.seed(nrow(mol_metadata))
      mol_metadata[, id:=1:.N]
      mol_metadata[, assay_id:="Assay1"]
      mol_metadata[, postfix:=paste0("_", ids::random_id(bytes = 3, use_openssl = FALSE))]
      mol_metadata[, id:=do.call(paste0,.SD), .SDcols=c("id", "postfix")]
      
      # coord names
      rownames(mol_coords) <- mol_metadata$id
      
      # create VoltRon assay for molecules
      mol_assay <- formAssay(coords = mol_coords, image = image, type = "molecule", main_image = image_name)
      
      # merge assays in one section
      if(verbose)
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
  if(verbose)
    message("Merging slides ...")
  if(length(vr_list) > 1){
    vr <- merge(vr_list[[1]], vr_list[-1])
  } else {
    vr <- vr_list[[1]]
  }
  vr
}

# importCosMxTileDB <- function(tiledbURI,
#                               assay_name = "CosMx",
#                               image = NULL,
#                               image_name = "main",
#                               import_molecules = FALSE,
#                               feature_name = NULL,
#                               verbose = TRUE, ...)
# {
#   # check tiledb and tiledbsc
#   if (!requireNamespace("tiledb", quietly = TRUE))
#     stop("Please install the tiledb package: \n
#          remotes::install_github('TileDB-Inc/TileDB-R', force = TRUE, ref = '0.17.0')")
#   if (!requireNamespace("tiledbsc", quietly = TRUE))
#     stop("Please install the tiledbsc package: \n
#          remotes::install_github('tiledb-inc/tiledbsc', force = TRUE, ref = '8157b7d54398b1f957832f37fff0b173d355530e')")
# 
#   # get tiledb
#   if(verbose)
#     message("Scanning TileDB array for cell data ...")
#   tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, verbose = FALSE)
# 
#   # raw counts
#   counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE)
#   counts <- as.matrix(counts)
# 
#   # cell metadata
#   metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
# 
#   # coordinates
#   coords <- as.matrix(metadata[,c("x_slide_mm", "y_slide_mm")])
#   colnames(coords) <- c("x","y")
# 
#   # transcripts
#   if(import_molecules){
#     if(verbose)
#       message("Scanning TileDB array for molecule data ...")
#     subcellular <- tiledb::tiledb_array(
#       tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
#       return_as="data.table")[]
#     colnames(subcellular)[colnames(subcellular)=="target"] <- "gene"
#   }
# 
#   # get slides and construct VoltRon objects for each slides
#   slides <- unique(metadata$slide_ID_numeric)
# 
#   # for each slide create a VoltRon object with combined layers
#   vr_list <- list()
#   for(slide in slides){
# 
#     # cell assay
#     if(verbose)
#       message("Creating cell level assay for slide ", slide, " ...")
# 
#     # slide info
#     cur_coords <- coords[metadata$slide_ID_numeric == slide,]
#     cur_counts <- counts[,rownames(cur_coords)]
#     cur_metadata <- metadata[rownames(cur_coords),]
# 
#     # create VoltRon object
#     cell_object <- formVoltRon(data = cur_counts, metadata = cur_metadata, image = image, coords = cur_coords,
#                                main.assay = assay_name, assay.type = "cell", image_name = image_name, feature_name = feature_name, ...)
#     cell_object$Sample <- paste0("Slide", slide)
# 
#     # molecule assay
#     if(import_molecules){
# 
#       # get slide
#       if(verbose)
#         message("Creating molecule level assay for slide ", slide, " ...")
#       if("slideID" %in% colnames(subcellular)){
#         cur_subcellular <- subcellular
#       } else {
#         cur_subcellular <- subset(subcellular, slideID == slide)
#       }
# 
#       # coordinates
#       mol_coords <- as.matrix(cur_subcellular[,c("x_global_px", "y_global_px")])
#       colnames(mol_coords) <- c("x", "y")
# 
#       # get subcellular data components
#       mol_metadata <- cur_subcellular[,colnames(cur_subcellular)[!colnames(cur_subcellular) %in% c("CellId", "cell_id", "x_global_px", "y_global_px")], with = FALSE]
#       set.seed(nrow(mol_metadata))
#       mol_metadata[, id:=1:.N]
#       mol_metadata[, assay_id:="Assay1"]
#       mol_metadata[, postfix:=paste0("_", ids::random_id(bytes = 3, use_openssl = FALSE))]
#       mol_metadata[, id:=do.call(paste0,.SD), .SDcols=c("id", "postfix")]
# 
#       # coord names
#       rownames(mol_coords) <- mol_metadata$id
# 
#       # create VoltRon assay for molecules
#       mol_assay <- formAssay(coords = mol_coords, image = image, type = "molecule", main_image = image_name)
# 
#       # merge assays in one section
#       if(verbose)
#         message("Merging assays for slide ", slide, " ...")
#       sample.metadata <- SampleMetadata(cell_object)
#       cell_object <- addAssay(cell_object,
#                               assay = mol_assay,
#                               metadata = mol_metadata,
#                               assay_name = paste0(assay_name, "_mol"),
#                               sample = sample.metadata["Assay1", "Sample"],
#                               layer = sample.metadata["Assay1", "Layer"])
#     }
#     vr_list <- append(vr_list, cell_object)
#   }
# 
#   # return
#   if(verbose)
#     message("Merging slides ...")
#   if(length(vr_list) > 1){
#     vr <- merge(vr_list[[1]], vr_list[-1])
#   } else {
#     vr <- vr_list[[1]]
#   }
#   vr
# }

#' generateCosMxImage
#'
#' Generates a low resolution Morphology image of the CosMx experiment
#'
#' @param dir.path CosMx folder of images
#' @param fov.position.file the file providing FOV positions
#' @param increase.contrast increase the contrast of the image before writing
#' @param output.path The path to the new morphology image created if the image should be saved to a location other than Xenium output folder.
#' @param verbose verbose
#' @param ... additional parameters passed to the \link{writeImage} function
#'
#' @importFrom magick image_read image_contrast image_composite image_blank
#' @importFrom EBImage writeImage
#' @importFrom stringr str_pad
#'
#' @export
generateCosMxImage <- function(dir.path, fov.position.file, increase.contrast = FALSE, output.path = NULL, verbose = TRUE, ...) {
  
  # check package
  if(!requireNamespace("reshape2")){
    stop("You have to install the reshape2 package!: install.packages('reshape2')")
  }
  
  # file path to either Xenium output folder or specified folder
  file.path <- paste0(dir.path, "/CellComposite_lowres.tif")
  output.file <- paste0(output.path, "/CellComposite_lowres.tif")
  
  # check if the file exists in either Xenium output folder, or the specified location
  if(file.exists(file.path) | file.exists(paste0(output.file))){
    if(verbose)
      message("CellComposite_lowres.tif already exists! \n")
    return(NULL)
  }
  
  # Combine Images of the FOV grid
  if(verbose)
    message("Parsing tif files ...")
  image.dir.path <- paste0(dir.path,"/CellComposite/")
  image.files <- list.files(image.dir.path)
  image.files <- image.files[grepl("^CellComposite_F", image.files)]
  image.files_FOV_ID <- vapply(image.files, function(x) strsplit(x, split = "CellComposite_F[0]+")[[1]][2], character(1))
  image.files_FOV_ID <- as.numeric(gsub(".jpg", "", image.files_FOV_ID))
  
  # FOV positions of CosMx
  if(verbose)
    message("Getting FOV Positions ...")
  fov_positions_path <- fov.position.file
  fov_positions <- read.csv(fov_positions_path)
  fov_positions[,3] <- max(fov_positions[,3]) - fov_positions[,3]
  fov_positions[,c(2,3)] <- floor(fov_positions[,c(2,3)]/10)
  fov_info <- colnames(fov_positions)[grepl("fov|FOV", colnames(fov_positions))]
  if(length(fov_info) > 0){
    fov_info <- fov_positions[[fov_info]]
    fov_info <- fov_info[fov_info %in% image.files_FOV_ID]
  } else {
    stop("no FOV info is found in fov_positions_file.csv")
  }
  
  # Combine images 
  if(verbose)
    message("Combining Images ...")
  extent <- apply(fov_positions[,c(2,3)], 2, max)
  extent <- extent + 600
  morphology_image <- magick::image_blank(extent[1], 
                                          extent[2], 
                                          color = "black")
  n <- length(fov_info)
  for(i in 1:n){
    
    # progress 
    cat(paste0(round(i / n * 100), '% completed \n'))
    Sys.sleep(.05)
    if (i == n) cat('Done\n') else cat("\033[2K")
    
    # add image
    image_path <- image.files[grepl(paste0("F[0]+", fov_info[i],".jpg$"), image.files)][1]
    imagedata <- magick::image_read(file.path(image.dir.path, image_path)) %>% 
      magick::image_resize(magick::geometry_size_percent(10))
    off_set <- fov_positions[i,c(2,3)]
    morphology_image <- morphology_image %>% 
      image_composite(imagedata, offset = paste0("+", off_set[1], "+", off_set[2]))
    rm(imagedata)
  }
  
  # pick a resolution level
  morphology_image_info <- image_info(morphology_image)
  if(verbose)
    message("Image Resolution ",
            "(X:", morphology_image_info$width, 
            " Y:", morphology_image_info$height, ") \n")
  
  # increase contrast
  if(increase.contrast) {
    if(verbose)
      message("Increasing Contrast \n")
    morphology_image <- magick::image_contrast(morphology_image, sharpen = 1)
  }
  
  # write to the same folder
  if(verbose)
    message("Writing tif file ...")
  if(is.null(output.path)){
    EBImage::writeImage(magick::as_EBImage(morphology_image), 
                        file = file.path, ...)
  } else {
    EBImage::writeImage(magick::as_EBImage(morphology_image), 
                        file = output.file, ...)
  }
  
  # return
  return(morphology_image)
}

####
# Spatial Genomics ####
####

####
## GenePs ####
####

#' importGenePS
#'
#' Importing GenePS data
#' 
#' @param dir.path path to Xenium output folder
#' @param assay_name the assay name of the SR object
#' @param sample_name the name of the sample
#' @param use_image if TRUE, the DAPI image will be used.
#' @param resolution_level the level of resolution within TIFF image. Default: 7 (971x638)
#' @param image_name the image name of the Xenium assay, Default: main
#' @param channel_name the channel name of the image of the Xenium assay, Default: DAPI
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param verbose verbose
#' @param ... additional parameters passed to \link{formVoltRon}
#' 
#' @importFrom magick image_read image_info
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @importFrom ids random_id
#' @importFrom grDevices as.raster
#' @importFrom EBImage as.Image
#'
#' @export
#'
importGenePS <- function (dir.path, assay_name = "GenePS", sample_name = NULL, use_image = TRUE, 
                          resolution_level = 7, image_name = "main", channel_name = "DAPI", import_molecules = FALSE, 
                          verbose = TRUE, ...)
{
  # cell assay
  if(verbose)
     message("Creating cell level assay ...")
  
  # raw counts
  DataOutputfiles <- list.files(paste0(dir.path, "DataOutput/"))
  datafile <- DataOutputfiles[grepl("_cellxgene.csv", DataOutputfiles)]
  if(length(datafile) > 0){
    rawdata <- utils::read.csv(file = paste0(dir.path, "DataOutput/", datafile), row.names = NULL)
    rawdata <- rawdata[,colnames(rawdata)[!colnames(rawdata) %in% "X"]]
    rawdata <- t(rawdata)
  } else {
    stop("There are no files ending with '_cellxgene.csv' in the path")
  }
  
  # image
  if(use_image){
    
    # check RBioFormats
    if (!requireNamespace('RBioFormats'))
      stop("Please install RBioFormats package to extract xml from the ome.tiff file!: BiocManager::install('RBioFormats')")
    
    # check image
    image_file <- paste0(dir.path, "/images/DAPI.tiff")
    if(file.exists(image_file)){
      image <- RBioFormats::read.image(image_file, resolution = resolution_level)
      image <- EBImage::as.Image(image)
      image <- grDevices::as.raster(image)
      image <- magick::image_read(image)
    } else {
      stop("There are no image files in the path")
    }
    # scale the xenium image instructed by 10x Genomics help page
    scaleparam <- ceiling(2^(resolution_level-1))
  } else {
    image <- NULL
  }
  
  # coordinates
  coordsfile <- DataOutputfiles[grepl("_cellcoordinate.csv", DataOutputfiles)]
  if(length(coordsfile) > 0){
    coords <- utils::read.csv(file = paste0(dir.path, "DataOutput/", coordsfile))
    coords <- as.matrix(coords[,c("center_x", "center_y")])
    colnames(coords) <- c("x", "y")
    if(use_image) {
      coords <- coords/scaleparam
      imageinfo <- unlist(magick::image_info(image)[c("height")])
      range_coords <- c(0,imageinfo)
    } else {
      range_coords <- range(coords[,2])
    }
    coords[,2] <- range_coords[2] - coords[,2] + range_coords[1]
  } else {
    stop("There are no files ending with '_cellcoordinate.csv' in the path")
  }
  
  # create VoltRon object for cells
  cell_object <- formVoltRon(rawdata, metadata = NULL, image = image, coords, main.assay = assay_name, 
                             assay.type = "cell", image_name = image_name, main_channel = channel_name, sample_name = sample_name, 
                             feature_name = "RNA", ...)
  
  # add molecules
  if(!import_molecules){
    return(cell_object)
  } else {
    if(verbose)
      message("Creating molecule level assay ...")
    transcripts_file <- DataOutputfiles[grepl("_transcriptlocation.csv", DataOutputfiles)]
    if(length(transcripts_file) == 0){
      stop("There are no files ending with '_transcriptlocation.csv' in the path")
    } else {
      
      # get subcellur data components
      subcellular_data <- data.table::fread(paste0(dir.path, "DataOutput/", transcripts_file))
      subcellular_data$id <- seq_len(nrow(subcellular_data))
      subcellular_data <- subcellular_data[,c("id", colnames(subcellular_data)[!colnames(subcellular_data) %in% "id"]), with = FALSE]
      colnames(subcellular_data)[colnames(subcellular_data)=="name"] <- "gene"
      
      # coordinates
      coords <- as.matrix(subcellular_data[,c("x", "y")])
      colnames(coords) <- c("x", "y")
      if(use_image){
        coords <- coords/scaleparam
      }
      coords[,"y"] <- range_coords[2] - coords[,"y"]  + range_coords[1]
      
      # metadata
      mol_metadata <- subcellular_data[,colnames(subcellular_data)[!colnames(subcellular_data) %in% c("cell", "x", "y")], with = FALSE]
      set.seed(nrow(mol_metadata$id))
      mol_metadata$postfix <- paste0("_", ids::random_id(bytes = 3, use_openssl = FALSE))
      mol_metadata$assay_id <- "Assay1"
      mol_metadata[, id:=do.call(paste0,.SD), .SDcols=c("id", "postfix")]
      
      # coord names
      rownames(coords) <- mol_metadata$id

      # create VoltRon assay for molecules
      mol_assay <- formAssay(coords = coords, image = image, type = "molecule", main_image = image_name, main_channel = channel_name)
      
      # merge assays in one section
      if(verbose)
        message("Merging assays ...")
      sample.metadata <- SampleMetadata(cell_object)
      object <- addAssay(cell_object,
                         assay = mol_assay,
                         metadata = mol_metadata,
                         assay_name = paste0(assay_name, "_mol"),
                         sample = sample.metadata["Assay1", "Sample"],
                         layer = sample.metadata["Assay1", "Layer"])
      
      # connectivity
      connectivity <- data.table::data.table(transcript_id = mol_metadata$id, cell_id = subcellular_data[["cell"]])
      connectivity <- subset(connectivity, cell_id != 0)
      connectivity[["cell_id"]] <- vrSpatialPoints(cell_object)[connectivity[["cell_id"]]]
      
      # add connectivity of spatial points across assays
      object <- addLayerConnectivity(object,
                                connectivity = connectivity,
                                sample = sample.metadata["Assay1", "Sample"],
                                layer = sample.metadata["Assay1", "Layer"])
      
      # return
      return(object)
    }
  }
}
  
####
# BGI Genomics ####
####

####
## STOmics ####
####

#' importSTOmics
#'
#' Importing STOmics (Stereo-Seq) data
#'
#' @param h5ad.path path to h5ad file of STOmics output
#' @param assay_name the assay name
#' @param sample_name the name of the sample
#' @param image_name the image name of the Visium assay, Default: main
#' @param channel_name the channel name of the image of the Visium assay, Default: H&E
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom methods as
#' @export
#'
importSTOmics <- function(h5ad.path, assay_name = "STOmics", sample_name = NULL, image_name = "main", channel_name = "H&E", ...)
{
  # check package
  if(!requireNamespace('anndataR'))
    stop("Please install anndataR package!: devtools::install_github('scverse/anndataR')")
  
  # get h5ad data
  stdata <- anndataR::read_h5ad(h5ad.path, to = "HDF5AnnData")
  
  # observation and feature names
  obs_names <- stdata$obs_names
  var_names <- stdata$var_names
  
  # raw counts
  rawdata <- Matrix::t(stdata$X)
  rownames(rawdata) <- var_names
  colnames(rawdata) <- obs_names
  rawdata <- methods::as(rawdata, 'CsparseMatrix')
  
  # metadata
  metadata <- stdata$obs
  rownames(metadata) <- obs_names
  
  # coordinates
  coords <- stdata$obsm$spatial
  rownames(coords) <- obs_names
  
  # scale coordinates
  binsize <- stdata$uns$bin_size
  params <- list(
    spot.radius = 0.5 + (binsize-1),
    vis.spot.radius = 0.5 + (binsize-1))
  
  # create VoltRon
  formVoltRon(rawdata, metadata = metadata, coords = coords, main.assay = assay_name, params = params, assay.type = "spot", 
              image_name = image_name, main_channel = channel_name, sample_name = sample_name, feature_name = "RNA", ...)
}

####
# Akoya ####
####

####
## PhenoCycler ####
####

#' importPhenoCycler
#'
#' Importing PhenoCycler data
#' 
#' @param dir.path path to PhenoCycler output folder
#' @param assay_name the assay name of the SR object
#' @param sample_name the name of the sample
#' @param image_name the image name of the Xenium assay, Default: main
#' @param type Specify which type matrix is being provided.
#' \itemize{
#'  \item \dQuote{\code{processor}}: matrix generated by CODEX Processor
#'  \item \dQuote{\code{inform}}: matrix generated by inForm
#'  \item \dQuote{\code{qupath}}: matrix generated by QuPath
#' }
#' @param filter A pattern to filter features by; pass \code{NA} to skip feature filtering
#' @param inform.quant When \code{type} is \dQuote{\code{inform}}, the quantification level to read in 
#' @param verbose verbose
#' @param ... additional parameters passed to \link{formVoltRon} 
#' 
#' @importFrom magick image_info image_read
#'
#' @export
importPhenoCycler <- function(dir.path, assay_name = "PhenoCycler", sample_name = NULL, image_name = "main", 
                              type = c('inform', 'processor', 'qupath'), filter = 'DAPI|Blank|Empty', 
                              inform.quant = c('mean', 'total', 'min', 'max', 'std'), verbose = TRUE, ...){
  
  # raw counts, metadata and coordinates 
  listoffiles <- list.files(paste0(dir.path, "/processed/segm/segm-1/fcs/compensated/"), full.names = TRUE)
  datafile <- listoffiles[grepl("_compensated.csv$", listoffiles)][1]
  if(file.exists(datafile)){
    rawdata <- readPhenoCyclerMat(filename = datafile, type = type, filter = filter, inform.quant = inform.quant)
  } else {
    stop("There are no files named ending with '_compensated.csv' in the processed/segm/segm-1/fcs/compensated/ subfolder")
  }
  
  # cell id 
  cellid <- paste0("cell", colnames(rawdata$matrix))
  
  # coordinates
  coords <- rawdata$centroids[,c("x", "y")]
  rownames(coords) <- cellid
  coords <- as.matrix(coords)
  
  # metadata
  metadata <- rawdata$metadata
  rownames(metadata) <- cellid
  
  # data
  rawdata <- rawdata$matrix
  colnames(rawdata) <- cellid
  rownames(rawdata) <- gsub("(\t|\r|\n)", "", rownames(rawdata))
  
  # images 
  image_dir <- paste0(dir.path, "/processed/stitched/reg001/")
  list_files <- list.files(image_dir)
  if(!dir.exists(image_dir)){
    if(verbose)
      message("There are no images of channels!")
    image_list <- NULL
  } else {
    if(!any(grepl(".tif$", list_files))){
      stop("The folder doesnt have any images associated with channels!")
    } else{
      image_channel_names <- sapply(list_files, function(x) {
        name <- strsplit(x, split = "_")[[1]]
        name <- gsub(".tif$", "", name[length(name)])
        name <- make.unique(name)
        return(name)
      })
      image_channel_names <- c(image_channel_names[grepl("DAPI", image_channel_names)][1], 
                               image_channel_names[image_channel_names %in% rownames(rawdata)])
      list_files <- names(image_channel_names)
      image_list <- lapply(list_files, function(x){
        magick::image_read(paste0(image_dir, "/", x))
      })
      names(image_list) <- image_channel_names   
      coords[,2] <- magick::image_info(image_list[[1]])$height - coords[,2]
    }
  }
  
  # voltron object
  object <- formVoltRon(data = rawdata, metadata = metadata, image = image_list, coords = coords, assay.type = "cell", 
                        sample_name = sample_name, main.assay = assay_name, image_name = image_name, feature_name = "RNA", ...)
  
  # return
  object
}

#' readPhenoCyclerMat
#' 
#' Read and Load Akoya CODEX data, adapted from the \code{ReadAkoya} function \code{Seurat} package
#'
#' @param filename Path to matrix generated by upstream processing.
#' @param type Specify which type matrix is being provided.
#' \itemize{
#'  \item \dQuote{\code{processor}}: matrix generated by CODEX Processor
#'  \item \dQuote{\code{inform}}: matrix generated by inForm
#'  \item \dQuote{\code{qupath}}: matrix generated by QuPath
#' }
#' @param filter A pattern to filter features by; pass \code{NA} to skip feature filtering
#' @param inform.quant When \code{type} is \dQuote{\code{inform}}, the quantification level to read in
#' 
#' @importFrom methods as
#' 
#' @noRd
readPhenoCyclerMat <- function(
    filename,
    type = c('inform', 'processor', 'qupath'),
    filter = 'DAPI|Blank|Empty',
    inform.quant = c('mean', 'total', 'min', 'max', 'std')
) {
  # Check arguments
  if (!file.exists(filename)) {
    stop("Can't find file: ", filename)
  }
  type <- tolower(x = type[1L])
  type <- match.arg(arg = type)
  ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
  # Preload matrix
  sep <- switch(EXPR = type, 'inform' = '\t', ',')
  mtx <- data.table::fread(
    file = filename,
    sep = sep,
    data.table = FALSE,
    verbose = FALSE
  )
  # Assemble outputs
  outs <- switch(
    EXPR = type,
    'processor' = {
      # Create centroids data frame
      centroids <- data.frame(
        x = mtx[['x:x']],
        y = mtx[['y:y']],
        cell = as.character(x = mtx[['cell_id:cell_id']]),
        stringsAsFactors = FALSE
      )
      rownames(x = mtx) <- as.character(x = mtx[['cell_id:cell_id']])
      # Create metadata data frame
      md <- mtx[, !grepl(pattern = '^cyc', x = colnames(x = mtx)), drop = FALSE]
      colnames(x = md) <- vapply(
        X = strsplit(x = colnames(x = md), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        2L
      )
      # Create expression matrix
      mtx <- mtx[, grepl(pattern = '^cyc', x = colnames(x = mtx)), drop = FALSE]
      colnames(x = mtx) <- vapply(
        X = strsplit(x = colnames(x = mtx), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        2L
      )
      if (!is.na(x = filter)) {
        mtx <- mtx[, !grepl(pattern = filter, x = colnames(x = mtx)), drop = FALSE]
      }
      mtx <- t(x = mtx)
      if ((sum(mtx == 0) / length(x = mtx)) > ratio) {
        # mtx <- as.sparse(x = mtx)
        mtx <- methods::as(mtx, 'CsparseMatrix')
      }
      list(matrix = mtx, centroids = centroids, metadata = md)
    },
    'inform' = {
      inform.quant <- tolower(x = inform.quant[1L])
      inform.quant <- match.arg(arg = inform.quant)
      expr.key <- c(
        mean = 'Mean',
        total = 'Total',
        min = 'Min',
        max = 'Max',
        std = 'Std Dev'
      )[inform.quant]
      expr.pattern <- '\\(Normalized Counts, Total Weighting\\)'
      rownames(x = mtx) <- mtx[['Cell ID']]
      mtx <- mtx[, setdiff(x = colnames(x = mtx), y = 'Cell ID'), drop = FALSE]
      # Create centroids
      centroids <- data.frame(
        x = mtx[['Cell X Position']],
        y = mtx[['Cell Y Position']],
        cell  = rownames(x = mtx),
        stringsAsFactors = FALSE
      )
      # Create metadata
      cols <- setdiff(
        x = grep(
          pattern = expr.pattern,
          x = colnames(x = mtx),
          value = TRUE,
          invert = TRUE
        ),
        y = paste('Cell', c('X', 'Y'), 'Position')
      )
      md <- mtx[, cols, drop = FALSE]
      # Create expression matrices
      exprs <- data.frame(
        cols = grep(
          pattern = paste(expr.key, expr.pattern),
          x = colnames(x = mtx),
          value = TRUE
        )
      )
      exprs$feature <- vapply(
        X = trimws(x = gsub(
          pattern = paste(expr.key, expr.pattern),
          replacement = '',
          x = exprs$cols
        )),
        FUN = function(x) {
          x <- unlist(x = strsplit(x = x, split = ' '))
          x <- x[length(x = x)]
          return(gsub(pattern = '\\(|\\)', replacement = '', x = x))
        },
        FUN.VALUE = character(length = 1L)
      )
      exprs$class <- tolower(x = vapply(
        X = strsplit(x = exprs$cols, split = ' '),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      ))
      classes <- unique(x = exprs$class)
      outs <- vector(
        mode = 'list',
        length = length(x = classes) + 2L
      )
      names(x = outs) <- c(
        'matrix',
        'centroids',
        'metadata',
        setdiff(x = classes, y = 'entire')
      )
      outs$centroids <- centroids
      outs$metadata <- md
      for (i in classes) {
        df <- exprs[exprs$class == i, , drop = FALSE]
        expr <- mtx[, df$cols]
        colnames(x = expr) <- df$feature
        if (!is.na(x = filter)) {
          expr <- expr[, !grepl(pattern = filter, x = colnames(x = expr)), drop = FALSE]
        }
        expr <- t(x = expr)
        if ((sum(expr == 0, na.rm = TRUE) / length(x = expr)) > ratio) {
          # expr <- as.sparse(x = expr)
          expr <- methods::as(expr, 'CsparseMatrix')
        }
        outs[[switch(EXPR = i, 'entire' = 'matrix', i)]] <- expr
      }
      outs
    },
    'qupath' = {
      rownames(x = mtx) <- as.character(x = seq_len(length.out = nrow(x = mtx)))
      # Create centroids
      xpos <- sort(
        x = grep(pattern = 'Centroid X', x = colnames(x = mtx), value = TRUE),
        decreasing = TRUE
      )[1L]
      ypos <- sort(
        x = grep(pattern = 'Centroid Y', x = colnames(x = mtx), value = TRUE),
        decreasing = TRUE
      )[1L]
      centroids <- data.frame(
        x = mtx[[xpos]],
        y = mtx[[ypos]],
        cell = rownames(x = mtx),
        stringsAsFactors = FALSE
      )
      # Create metadata
      cols <- setdiff(
        x = grep(
          pattern = 'Cell: Mean',
          x = colnames(x = mtx),
          ignore.case = TRUE,
          value = TRUE,
          invert = TRUE
        ),
        y = c(xpos, ypos)
      )
      md <- mtx[, cols, drop = FALSE]
      # Create expression matrix
      idx <- which(x = grepl(
        pattern = 'Cell: Mean',
        x = colnames(x = mtx),
        ignore.case = TRUE
      ))
      mtx <- mtx[, idx, drop = FALSE]
      colnames(x = mtx) <- vapply(
        X = strsplit(x = colnames(x = mtx), split = ':'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      )
      if (!is.na(x = filter)) {
        mtx <- mtx[, !grepl(pattern = filter, x = colnames(x = mtx)), drop = FALSE]
      }
      mtx <- t(x = mtx)
      if ((sum(mtx == 0) / length(x = mtx)) > ratio) {
        # mtx <- as.sparse(x = mtx)
        mtx <- methods::as(mtx, 'CsparseMatrix')
      }
      list(matrix = mtx, centroids = centroids, metadata = md)
    },
    stop("Unknown matrix type: ", type)
  )
  return(outs)
}

####
# Non-Commercial ####
####

####
## OpenST ####
####

#' importOpenST
#'
#' Importing OpenST data
#'
#' @param h5ad.path path to h5ad file of STOmics output
#' @param assay_name the assay name
#' @param sample_name the name of the sample
#' @param image_name the image name of the Visium assay, Default: main
#' @param channel_name the channel name of the image of the Visium assay, Default: H&E
#' @param verbose verbose
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom methods as
#' @importFrom Matrix t
#' 
#' @export
importOpenST <- function(h5ad.path, assay_name = "OpenST", sample_name = NULL, image_name = "main", channel_name = "H&E", verbose = TRUE, ...)
{
  # check package
  if(!requireNamespace('anndataR'))
    stop("Please install anndataR package!: devtools::install_github('scverse/anndataR')")
  
  # get h5ad data
  stdata <- anndataR::read_h5ad(h5ad.path, to = "HDF5AnnData")
  
  # observation and feature names
  obs_names <- stdata$obs_names
  var_names <- stdata$var_names
  
  # raw counts
  rawdata <- Matrix::t(stdata$X)
  rownames(rawdata) <- var_names
  colnames(rawdata) <- obs_names
  rawdata <- methods::as(rawdata, 'CsparseMatrix')
  
  # metadata
  metadata <- stdata$obs
  rownames(metadata) <- obs_names
  
  # coordinates
  obsm <- stdata$obsm
  if("spatial_3d_aligned" %in% names(obsm)) {
    coords <- stdata$obsm$spatial_3d_aligned
    rownames(coords) <- obs_names
    zlocation <- unique(coords[,3])
    sections <- unique(metadata$n_section)
    zlocation <- zlocation[order(sections)]
    sections <- sections[order(sections)]
  } else {
    coords <- stdata$obsm$spatial
    rownames(coords) <- obs_names
    zlocation <- 0
    metadata$n_sections <- sections <- 1
  }
  
  # get individual sections as voltron data if there are any
  vr_data_list <- list()
  if(verbose)
    message("Creating Layers ...")
  for(i in seq_len(length(sections))){
    ind <- metadata$n_section == sections[i]
    spatialpoints <- rownames(metadata[metadata$n_section == sections[i],])
    cur_data <- rawdata[,spatialpoints]
    cur_metadata <- metadata[spatialpoints,]
    cur_coords <- coords[ind,c(1,2)]
    rownames(cur_coords) <- spatialpoints
    vr_data_list[[i]] <- 
      formVoltRon(data = cur_data, 
                  metadata = cur_metadata, 
                  coords = cur_coords, 
                  main.assay = assay_name, 
                  sample_name = paste0("Section", sections[i]),
                  image_name = image_name, 
                  main_channel = channel_name, 
                  feature_name = "RNA", ...)
  }
  
  # create VoltRon
  sample_name <- ifelse(is.null(sample_name), "Sample", sample_name)
  
  # set zlocations and adjacency of layer in the vrBlock
  if(length(vr_data_list) > 1){
    vr_data <- mergeVoltRon(vr_data_list[[1]], 
                            vr_data_list[-1], 
                            samples = sample_name)
    connectivity <- data.frame(Var1 = rep(seq_len(length(sections)), 
                                          length(sections)), 
                               Var2 = rep(seq_len(length(sections)), 
                                          each = length(sections)))
    vr_data <- addBlockConnectivity(vr_data, 
                                    connectivity = connectivity, 
                                    zlocation = zlocation, 
                                    sample = sample_name)
  } else {
    vr_data <- vr_data_list[[1]]
  }
  
  # return
  vr_data
}

####
## DBIT-Seq ####
####
 
#' importDBITSeq
#'
#' Importing DBIT-Seq data
#'
#' @param path.rna path to rna count matrix
#' @param path.prot path to protein count matrix
#' @param size the size of the in situ pixel (defualt is 10 (micron))
#' @param assay_name the assay name
#' @param sample_name the name of the sample
#' @param image_name the image name of the Visium assay, Default: main
#' @param channel_name the channel name of the image of the Visium assay, Default: H&E
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom utils read.table
#'
#' @export
importDBITSeq <- function(path.rna, path.prot = NULL, size = 10, assay_name = "DBIT-Seq", sample_name = NULL, image_name = "main", channel_name = "H&E", ...)
{
  # count matrix RNA
  rnadata <- utils::read.table(path.rna, header = TRUE, sep = "\t", row.names = 1)
  rnadata <- t(as.matrix(rnadata))
  
  # count matrix Protein
  protdata <- utils::read.table(path.prot, header = TRUE, sep = "\t", row.names = 1)
  protdata <- t(as.matrix(protdata))
  
  # coords 
  coords <- sapply(colnames(rnadata), function(x) as.numeric(strsplit(x, split = "x")[[1]]))
  coords <- t(coords)
  colnames(coords) <- c("x", "y")
  
  # get DBIT-Seq parameters
  params <- list(
    spot.radius = size/2,
    vis.spot.radius = size,
    nearestpost.distance = size*2*sqrt(2))
  coords <- coords*size*3/2
  
  # make voltron object
  object <- formVoltRon(data = rnadata, coords = coords, image = NULL, assay.type = "spot", params = params, image_name = "main", 
                        main.assay = assay_name, sample_name = sample_name, feature_name = "RNA", ...)
  
  # add protein assay
  if(!is.null(path.prot)){
   
    # # create protein assay
    # new_assay <- formAssay(data = protdata,
    #                        coords = coords,
    #                        type = "spot")
    # new_assay@image <- object[["Assay1"]]@image
    # sample.metadata <- SampleMetadata(object)
    
    # # add new assay
    # object <- addAssay(object,
    #                    assay = new_assay,
    #                    metadata = Metadata(object),
    #                    assay_name = paste(assay_name, "Prot", sep = "-"),
    #                    sample = sample.metadata["Assay1", "Sample"],
    #                    layer = sample.metadata["Assay1", "Layer"])
    # 
    # # add connectivity of spatial points across assays
    # connectivity <- cbind(vrSpatialPoints(object, assay = "Assay1"),
    #                       vrSpatialPoints(object, assay = "Assay2"))
    # object <- addConnectivity(object,
    #                           connectivity = connectivity,
    #                           sample = sample.metadata["Assay1", "Sample"],
    #                           layer = sample.metadata["Assay1", "Layer"])
    
    # add negative probe assay as new feature set
    object <- addFeature(object, assay = assay_name, data = protdata, feature_name = "Protein")

  }
  
  # return 
  return(object)
}

####
# Image Data ####
####

#' importImageData
#'
#' import an image as VoltRon object
#'
#' @param image a single or a list of image paths or magick-image objects
#' @param tile.size the size of tiles
#' @param segments Either a list of segments or a GeoJSON file. This will result in a second assay in the VoltRon object to be created
#' @param image_name the image name of the Image assay, Default: main
#' @param channel_names the channel names of the images if multiple images are provided
#' @param series the series IDs of the pyramidal image, 
#' typical an integer starting from 1
#' @param resolution the resolution IDs of the 
#' pyramidal image, typical an integer starting from 1
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom data.table data.table
#'
#' @examples
#' # single image
#' imgfile <- system.file("extdata", "DAPI.tif", package = "VoltRon")
#' vrdata <- importImageData(imgfile, image_name = "main")
#' 
#' # multiple images
#' imgfile <- c(system.file("extdata", "DAPI.tif", package = "VoltRon"), 
#'              system.file("extdata", "DAPI.tif", package = "VoltRon"))
#' vrdata <- importImageData(imgfile, image_name = "main", channel_name = c("DAPI", "DAPI2"))
#' 
#' @export
importImageData <- function(image, tile.size = 10, segments = NULL, 
                            image_name = "main", channel_names = NULL, 
                            series = 1, resolution = NULL, ...){
  
  # # images and channel names
  # if(!is.null(channel_names)){
  #   if(length(image) != length(channel_names))
  #     stop("Provided channel names should of the same length as the images!")
  #   if(any(!is.character(channel_names)))
  #     stop("Invalid channel names!")  
  # }
  # 
  # # get image
  # if(!is.list(image)){}
  #   image <- as.list(image)
  # image <- sapply(image, function(img){
  #   if(!inherits(img, "magick-image")){
  #     if(!is.character(img)){
  #       stop("image should either be a magick-image object or a file.path")
  #     } else{
  #       if(file.exists(img)){
  #         img <- magick::image_read(img)
  #       } else {
  #         stop(img, " is not found!")
  #       }
  #     }
  #   }
  #   img
  # }, USE.NAMES = TRUE, simplify = FALSE)
  # 
  # # channel names
  # if(!is.null(channel_names)){
  #   names(image) <- channel_names
  # }
  image <- importImage(image, channel_names = channel_names, 
                       series = series, resolution = resolution)
  
  # check image size
  imageinfo <- vapply(image, function(img) {
    info <- magick::image_info(img)
    c(info$width, info$height)
  }, numeric(2))
  unique_width <- unique(imageinfo[1,])
  unique_height <- unique(imageinfo[2,])
  if(length(unique_width) == 1 && length(unique_height) == 1){
    imageinfo <- list(width = imageinfo[1,1], height = imageinfo[2,1])
  }

  # coordinates
  row_tile_size <- (imageinfo$height %/% tile.size)
  col_tile_size <- (imageinfo$width %/% tile.size)
  x_coords <- tile.size/2 + (seq_len(col_tile_size) - 1) * tile.size
  y_coords <- tile.size/2 + (seq_len(row_tile_size) - 1) * tile.size
  y_coords <- imageinfo$height - y_coords
  coords <- as.matrix(expand.grid(x_coords, y_coords))
  colnames(coords) <- c("x", "y")
  rownames(coords) <- paste0("tile", seq_len(nrow(coords)))

  # metadata
  metadata <- data.table::data.table(id = rownames(coords))

  # create voltron object with tiles
  object <- formVoltRon(data = NULL, metadata = metadata, image = image, coords, main.assay = "ImageData", assay.type = "tile", params = list(tile.size = tile.size), image_name = image_name, ...)
  
  # check if segments are defined
  if(is.null(segments)){
    return(object)
  } else {
    
    # check if segments are paths
    if(inherits(segments, "character")){
      if(grepl(".geojson$", segments)){
        segments <- generateSegments(geojson.file = segments)
      } else {
        stop("Only lists or GeoJSON files are accepted as segments input!")
      }
    }
    
    # make coordinates out of segments
    coords <- t(vapply(segments, function(dat){
      apply(dat[,c("x", "y")], 2, mean)
    }, numeric(2)))
    rownames(coords) <- names(segments)
    
    # make segment assay
    assay <- formAssay(coords = coords, segments = segments, image = image, type = "ROI", main_image = image_name) 
    
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

#' importImage
#'
#' import an image to be used in \code{importImageData}
#'
#' @param image a single or a list of image paths or magick-image objects
#' @param channel_names the channel names of the images if multiple images are provided
#' @param series the series IDs of the pyramidal image, 
#' typical an integer starting from 1
#' @param resolution the resolution IDs of the 
#' pyramidal image, typical an integer starting from 1
#'
#' @importFrom magick image_read image_info
#' @importFrom data.table data.table
#'
#' @noRd
importImage <- function(image, channel_names = NULL, 
                        series = 1, resolution = NULL){
  
  # images and channel names
  if(!is.null(channel_names)){
    if(length(image) != length(channel_names))
      stop("Provided channel names should of the same length as the images!")
    if(any(!is.character(channel_names)))
      stop("Invalid channel names!")
  }

  # check if image is ome.tiff
  if(is.character(image)){
    if(any(grepl(".ome.tiff$|.ome.tif$", image))){
      if (!requireNamespace('RBioFormats'))
        stop("Please install RBioFormats package to images from the ome.tiff file!: BiocManager::install('RBioFormats')")
      if(is.null(resolution))
        stop("For importing images from ome.tiff files, please specify resolution. ", 
             "See help(read.metadata) from RBioFormats package.")
      if(length(image) > 1)
        stop("Only a single ome.tiff file")
    }
  } else {
    if(!inherits(image, "magick-image"))
      stop("image should either be a magick-image object or a file.path")
  }
    
  # get image
  if(!is.list(image))
    image <- as.list(image)
  image <- sapply(image, function(img){
    # check if image exists
    if(is.character(img)){
      if(file.exists(img)){
        # ome.tiff images
        if(grepl(".ome.tiff$|.ome.tif$", img)){
          omexml <- RBioFormats::read.omexml(img)
          omexml <- XML::xmlToList(omexml, simplify = TRUE)
          meta <- RBioFormats::read.metadata(img)
          img <- RBioFormats::read.image(img, 
                                         series = series, 
                                         resolution = resolution, 
                                         normalize = TRUE)
          img <- EBImage::as.Image(img)
          # check if there are more than 3 or 1 channels
          if(length(d <- dim(img)) > 2){
            if(d[3] != 3){
              img <- sapply(seq_len(d[3]), function(i){
                tmp <- img[,,i]
                magick::image_read(grDevices::as.raster(tmp))
              })
            } else {
              img <- magick::image_read(grDevices::as.raster(img))
            }
          } else {
            img <- magick::image_read(grDevices::as.raster(img))
          }
        # regular tiff images
        } else {
          img <- magick::image_read(img)
        }
      } else {
        stop(img, " is not found!")
      }
    }
    img
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  # flatten nested image lists
  image <- unlist(image, recursive = FALSE)
  
  # channel names
  if(!is.null(channel_names)){
    if(length(image) == length(channel_names)){
      names(image) <- channel_names
    } else {
      stop("Provided channel_names should be the same length as the images!")
    }
  }
  
  # return
  return(image)
}


#' generateSegments
#' 
#' The function to import segments from a geojson file
#'
#' @param geojson.file the GeoJSON file, typically generated by QuPath software
#'
#' @importFrom rjson fromJSON
#' @importFrom dplyr tibble
#' 
#' @export
generateSegments <- function(geojson.file){
  
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
  }, seq_len(length(segments)), segments, SIMPLIFY = FALSE)
  
  # generate ROI names
  names(segments) <- paste0("ROI", seq_len(length(segments)))

  # return
  return(segments)
}

#' generateGeoJSON
#' 
#' generating geojson files from segments
#'
#' @param segments the segments, typically from \link{vrSegments}.
#' @param file the GeoJSON file, typically to be used by QuPath software.
#'
#' @importFrom rjson fromJSON
#' 
#' @export
generateGeoJSON <- function(segments, file){
  
  if(!requireNamespace('geojsonR'))
    stop("Please install geojsonR package for using geojsonR functions")
  
  # get segments
  if(!inherits(file, "character")){
    stop("file should be the path to the GeoJSON file!")
  } 
  
  # reshape segments
  segments <- mapply(function(id, sgt){
    poly <- na.omit(as.matrix(sgt[,c("x", "y")]))
    poly <- rbind(poly, poly[1,,drop = FALSE])
    poly <- as.list(data.frame(t(poly)))
    names(poly) <- NULL
    init <- geojsonR::TO_GeoJson$new()
    geometry <- init$Polygon(list(poly), stringify = TRUE)
    feature <- list(type = "Feature", 
                    id = id, 
                    geometry = geometry[!names(geometry) %in% "json_dump"],
                    properties = list(objectType = "annotation"))
    feature
  }, names(segments), segments, SIMPLIFY = FALSE, USE.NAMES = FALSE)

  # save as json
  segments <- rjson::toJSON(segments)
  write(segments, file = file)
}