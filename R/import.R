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
#' @param use_image if TRUE, the DAPI image will be used.
#' @param morphology_image the name of the lowred morphology image. Default: morphology_lowres.tif
#' @param resolution_level the level of resolution within Xenium OME-TIFF image, see \code{generateXeniumImage}. Default: 7 (553x402)
#' @param overwrite_resolution if TRUE, the image "file.name" will be generated again although it exists at "dir.path"
#' @param image_name the image name of the Xenium assay, Default: DAPI
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param ... additional parameters passed to \code{formVoltRon}
#'
#' @importFrom magick image_read image_info
#' @importFrom utils read.csv
#' @importFrom data.table fread
#' @importFrom ids random_id
#'
#' @export
#'
importXenium <- function (dir.path, selected_assay = "Gene Expression", assay_name = "Xenium", use_image = TRUE, morphology_image = "morphology_lowres.tif", resolution_level = 7, overwrite_resolution = FALSE, image_name = "DAPI", import_molecules = FALSE, ...)
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
  cell_object <- formVoltRon(rawdata, metadata = NULL, image = image, coords, segments = segments, main.assay = assay_name, assay.type = "cell", image_name = image_name, ...)

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
      entity_ID <- paste0(mol_metadata$id, "_", ids::random_id(bytes = 3, use_openssl = FALSE))
      mol_metadata <- data.table::data.table(id = entity_ID, assay_id = "Assay1", mol_metadata)

      # coord names
      rownames(coords) <- entity_ID
    }

    # create VoltRon assay for molecules
    mol_assay <- formAssay(coords = coords, image = image, type = "molecule", main_image = image_name)

    # merge assays in one section
    message("Merging assays ...")
    sample.metadata <- SampleMetadata(cell_object)
    object <- addAssay(cell_object,
                       assay = mol_assay,
                       metadata = mol_metadata,
                       assay_name = paste0(assay_name, "_mol"),
                       sample = sample.metadata["Assay1", "Sample"],
                       layer = sample.metadata["Assay1", "Layer"])
    return(object)
  }
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
#' @param image_name the image name of the Visium assay, Default: H&E
#' @param inTissue if TRUE, only barcodes that are in the tissue will be kept (default: TRUE)
#' @param ... additional parameters passed to \code{formVoltRon}
#'
#' @importFrom magick image_read
#' @importFrom rjson fromJSON
#' @importFrom utils read.csv
#'
#' @export
#'
importVisium <- function(dir.path, selected_assay = "Gene Expression", assay_name = "Visium", image_name = "H&E", inTissue = TRUE, ...)
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
  # coords$pxl_row_in_fullres <- max(coords$pxl_row_in_fullres) - coords$pxl_row_in_fullres + min(coords$pxl_row_in_fullres)
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
    scales <- scalefactors$tissue_lowres_scalef
    # spot.radius is the half of the diameter, but we visualize by a factor of 1.5 larger
    # params <- list(spot.radius = scalefactors$spot_diameter_fullres*scalefactors$tissue_lowres_scalef*1.5)
    params <- list(spot.radius = scalefactors$spot_diameter_fullres*scalefactors$tissue_lowres_scalef*2)
    coords <- coords*scales
    coords[,2] <- info$height - coords[,2]
  } else {
    stop("There are no files named 'scalefactors_json.json' in the path")
  }

  # create VoltRon
  formVoltRon(rawdata, metadata = NULL, image, coords, main.assay = assay_name, params = params, assay.type = "spot", image_name = image_name, ...)
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
#' @param ... additional parameters passed to \code{formVoltRon}
#'
#' @importFrom dplyr %>% full_join
#' @importFrom utils read.csv
#' @importFrom magick image_info image_read
#'
#' @export
#'
importGeoMx <- function(dcc.path, pkc.file, summarySegment, summarySegmentSheetName, assay_name = "GeoMx",
                        image = NULL, segment_polygons = FALSE, ome.tiff = NULL, ...)
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

  # get genes
  NegProbes <- pkcdata$RTS_ID[pkcdata$Target == "NegProbe-WTX"]
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
    segments <- importGeoMxSegments(ome.tiff, segmentsummary, geomx_image_info)
    segments <- segments[rownames(coords)]
  }

  # create VoltRon
  formVoltRon(rawdata, metadata = segmentsummary, image, coords, segments, main.assay = assay_name, assay.type = "ROI", ...)
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
#' @importFrom XML xmlToList
#'
importGeoMxSegments <- function(ome.tiff, summary, imageinfo){

  # check file
  if(file.exists(ome.tiff)){
    options(java.parameters = "-Xmx4g")
    if(grepl(".ome.tiff$|.ome.tif$", ome.tiff)){
      if (!requireNamespace('RBioFormats'))
        stop("Please install RBioFormats package extract xml from the ome.tiff file!")
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
#' Rescale GeoMx point (center or polygon corners of ROI) coordinates for image registration
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
#' @param image_name the image name of the CosMx assay, Default: DAPI
#' @param ome.tiff the OME.TIFF file of the CosMx experiment if exists
#' @param import_molecules if TRUE, molecule assay will be created along with cell assay.
#' @param ... additional parameters passed to \code{formVoltRon}
#'
#' @importFrom data.table data.table
#' @importFrom ids random_id
#'
#' @export
#'
importCosMx <- function(tiledbURI, assay_name = "CosMx",
                        image = NULL, image_name = "DAPI", ome.tiff = NULL, import_molecules = FALSE, ...)
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
      entity_ID <- paste0(1:nrow(mol_metadata), "_", ids::random_id(bytes = 3, use_openssl = FALSE))
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
  vr <- merge(vr_list[[1]], vr_list[-1])
}

####
# Image Data ####
####

#' importImageData
#'
#' import an image as VoltRon object
#'
#' @param img an image
#'
#' @importFrom magick image_read image_raster
#' @importFrom data.table data.table
#'
importImageData <- function(image, ...){

  # get image
  if(file.exists(image)){
    image <- magick::image_read(image)
  } else {
    stop(image, " is not found!")
  }

  # coordinates
  image_data <- magick::image_raster(image)

  # metadata
  metadata <- data.table(id = rownames(image_data))

  # coordinates
  coords <- as.matrix(image_data[,c("x","y")])

  # create voltron object
  formVoltRon(data = NULL, metadata = metadata, image = image, coords, main.assay = "ImageData", assay.type = "tile", ...)
}
