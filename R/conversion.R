####
# Seurat ####
####

#' @param type the spatial data type of Seurat object: "image" or "spatial"
#' @param assay_name the assay name
#' @param ... Additional parameter passed to \link{formVoltRon}
#'
#' @rdname as.VoltRon
#' @method as.VoltRon Seurat
#'
#' @importFrom stringr str_replace str_extract
#' @export
#'
as.VoltRon.Seurat <- function(object, type = c("image", "spatial"), assay_name = NULL, ...){

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects!: install.packages('Seurat')")

  # raw counts
  rawdata <- SeuratObject::LayerData(object, assay = Seurat::DefaultAssay(object), layer = "counts")

  # metadata
  metadata <- object@meta.data

  # embeddings
  if(length(object@reductions) > 0){
    embeddings_flag <- TRUE
    embedding_list <- sapply(object@reductions, Seurat::Embeddings, USE.NAMES = TRUE)
  } else {
    embeddings_flag <- FALSE
  }

  # image
  voltron_list <- list()
  spatialobjectlist <- object@images
  fov_names <- names(spatialobjectlist)
  if(length(spatialobjectlist) > 0){
    for(fn in fov_names){

      # message
      message("Converting FOV: ", fn, " ...")

      # image object
      spatialobject <- spatialobjectlist[[fn]]

      # cells
      cells <- Seurat::Cells(spatialobject)
      cells_nopostfix <- gsub("_Assay[0-9]+$", "", cells)

      # count
      cur_rawdata <- as.matrix(rawdata[,cells])
      colnames(cur_rawdata) <- cells_nopostfix

      # metadata
      cur_metadata <- metadata[cells,]
      rownames(cur_metadata) <- cells_nopostfix

      # coords
      coords <- as.matrix(Seurat::GetTissueCoordinates(spatialobject))[,1:2]
      coords <- apply(coords, 2, as.numeric)
      colnames(coords) <- c("x", "y")
      rownames(coords) <- cells_nopostfix

      # form voltron
      params <- list()
      assay.type <- "cell"
      assay_name <- "FOV"
      voltron_list[[fn]] <- formVoltRon(data = cur_rawdata, metadata = cur_metadata, coords = coords, main.assay = assay_name, params = params, assay.type = assay.type, sample_name = fn, ...)

      # embeddings
      spatialpoints <- vrSpatialPoints(voltron_list[[fn]])
      spatialpoints_nopostfix <- stringr::str_replace(spatialpoints, "_Assay[0-9]+$", "")
      spatialpoints_assay <- stringr::str_extract(spatialpoints, "Assay[0-9]+$")
      if(embeddings_flag){
        for(embed_name in names(embedding_list)){
          cur_embedding <- embedding_list[[embed_name]][cells,]
          rownames(cur_embedding) <- spatialpoints
          vrEmbeddings(voltron_list[[fn]], type = embed_name) <- cur_embedding
        }
      }
    }

    # merge object
    message("Merging object ...")
    vrobject <- merge(voltron_list[[1]], voltron_list[-1])
  } else{
    image <- NULL
    stop("There are no spatial objects available in this Seurat object")
  }

  return(vrobject)
}

#' as.Seurat
#'
#' Converting a VoltRon object into a Seurat object
#'
#' @param object a VoltRon object
#' @param cell.assay the name(type) of the cell assay to be converted
#' @param molecule.assay the name(type) of the molecule assay to be added to the cell assay in Seurat object
#' @param image_key the name (or prefix) of the image(s)
#' @param type the spatial data type of Seurat object: "image" or "spatial"
#' @param reg if TRUE, registered coordinates will be used
#'
#' @rdname as.Seurat
#'
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_replace
#'
#' @export
as.Seurat <- function(object, cell.assay = NULL, molecule.assay = NULL, image_key = "fov", type = c("image", "spatial"), reg = FALSE){
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")
  
  # check the number of assays
  if(is.null(cell.assay)){
    if(length(unique(sample_metadata[["Assay"]])) > 1){
      stop("You can only convert a single VoltRon assay into a Seurat object!")
    } else {
      cell.assay <- sample_metadata[["Assay"]]
    }
  } else {
    vrMainAssay(object) <- cell.assay
  }
  
  # check the number of assays
  if(unique(vrAssayTypes(object, assay = cell.assay)) %in% c("spot","ROI")) {
    stop("Conversion of Spot or ROI assays into Seurat is not yet permitted!")
  }
  
  # data
  data <- vrData(object, assay = cell.assay, norm = FALSE)
  
  # metadata
  metadata <- Metadata(object, assay = cell.assay)
  
  # Seurat object
  seu <- Seurat::CreateSeuratObject(counts = data, meta.data = metadata, assay = cell.assay)
  
  # add embeddings
  if(length(vrEmbeddingNames(object)) > 0){
    for(embd in vrEmbeddingNames(object)){
      embd_data <- vrEmbeddings(object, type = embd)
      colnames(embd_data) <- paste0(embd, 1:ncol(embd_data))
      seu[[embd]] <- Seurat::CreateDimReducObject(embd_data, key = paste0(embd, "_"), assay = Seurat::DefaultAssay(seu))
    }
  }
  
  # get image objects for each assay
  for(assy in vrAssayNames(object)){
    assay_object <- object[[assy]]
    if(type == "image"){
      coords <- vrCoordinates(assay_object, reg = reg)
      image.data <- list(centroids = SeuratObject::CreateCentroids(coords))
      if(!is.null(molecule.assay)){
        assay_metadata <- sample_metadata[assy,]
        molecule.assay.id <- rownames(sample_metadata)[sample_metadata$Assay == molecule.assay & (assay_metadata$Layer == sample_metadata$Layer & assay_metadata$Sample == sample_metadata$Sample)]
        if(length(molecule.assay.id) > 0){
          molecules_metadata <- Metadata(object, assay = molecule.assay.id)
          molecules_coords <- vrCoordinates(object, assay = molecule.assay.id, reg = reg)
          molecules <- dplyr::bind_cols(molecules_metadata, molecules_coords)
          rownames(molecules) <- stringr::str_replace(rownames(molecules), pattern = molecule.assay.id, replacement = assy)
          colnames(molecules)[colnames(molecules) %in% "feature_name"] <- "gene"
        }
      } else {
        molecules <- NULL
      }
      image.data <- SeuratObject::CreateFOV(coords = image.data, type = c("centroids"), molecules = molecules, assay = cell.assay)
      image <- paste0(image_key, assy)
      seu[[image]] <- image.data
    } else if(type == "spatial"){
      stop("Currently VoltRon does not support converting into Spatial-type (e.g. VisiumV1) Spatial objects!")
    }
  }
  
  
  # return
  seu
}

####
# AnnData ####
####

#' convertAnnDataToVoltRon
#'
#' converting AnnData h5ad files to VoltRon objects
#'
#' @param file h5ad file
#' @param AssayID the ID assays in the h5ad file
#' @param ... additional parameters passed to \link{formVoltRon}
#'
#' @export
#'
convertAnnDataToVoltRon <- function(file, AssayID = NULL, ...){
  
  # check Seurat package
  if(!requireNamespace('anndata'))
    stop("Please install anndata package")
  
  # read anndata
  adata <- anndata::read_h5ad(file)
  
  # raw counts
  rawdata <- as.matrix(t(adata$X))
  
  # metadata
  metadata <- adata$obs
  
  # coordinates and subcellular
  if(is.null(AssayID)){
    coords <- data.frame(adata$obsm, row.names = colnames(rawdata))
    coords <- apply(coords, 2, as.numeric)
    colnames(coords) <- c("x", "y")
    
    # scale coordinates and assay.type
    params <- list()
    assay.type <- "cell"
    assay_name <- "Xenium"
    
    # create VoltRon
    object <- formVoltRon(rawdata, metadata, image = NULL, coords, main.assay = assay_name, params = params, assay.type = assay.type, ...)
    
    # return
    return(object)
  } else {
  }
}

#' as.AnnData
#'
#' Converting a VoltRon object into a AnnData (.h5ad) object
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param file the name of the h5ad file.
#' @param flip_coordinates if TRUE, the spatial coordinates (including segments) will be flipped.
#' @param method the package to use for conversion: "anndataR" or "anndata".
#' @param create.ometiff should an ometiff file be generated of default image of the object
#' @param ... additional parameters passed to \link{vrImages}.
#' 
#' @details
#' This function converts a VoltRon object into an AnnData object (.h5ad file). It extracts assay data,
#' spatial coordinates, and optionally flips coordinates. Images associated with the assay can be included in the 
#' resulting AnnData file, with additional customization parameters like channel, scale.perc.
#' 
#' @rdname as.AnnData
#'
#' @importFrom stringr str_extract
#' @importFrom magick image_data
#' 
#' @export
as.AnnData <- function(object, 
                       file, 
                       assay = NULL, 
                       flip_coordinates = FALSE, 
                       method = "anndata", 
                       create.ometiff = FALSE, 
                       ...) {
  
  # Check the number of assays
  if (is.null(assay)) {
    if (length(unique(SampleMetadata(object)[["Assay"]])) > 1) {
      stop("You can only convert a single VoltRon assay into a Anndata object!")
    } else {
      assay <- SampleMetadata(object)[["Assay"]]
    }
  }
  assay <- vrAssayNames(object, assay = assay)
  
  # Check the number of assays
  if (unique(vrAssayTypes(object, assay = assay)) %in% c("ROI", "tile")) {
    stop("Conversion of tile or ROI assays into Anndata is not permitted!")
  }
  
  # Data
  data <- vrData(object, assay = assay, norm = FALSE)
  
  # Metadata
  metadata <- Metadata(object, assay = assay)
  metadata[["library_id"]] <- stringr::str_extract(rownames(metadata), "_Assay[0-9]+$")
  metadata[["library_id"]] <- gsub("^_", "", metadata[["library_id"]])
  
  # Embeddings
  obsm <- list()
  if (length(vrEmbeddingNames(object, assay = assay)) > 0) {
    for (embed_name in vrEmbeddingNames(object, assay = assay)) {
      obsm[[embed_name]] <- vrEmbeddings(object, assay = assay, type = embed_name)
    }
  }
  
  # Flip coordinates
  if (flip_coordinates) {
    object <- flipCoordinates(object, assay = assay)
  }
  
  # Coordinates
  coords <- vrCoordinates(object, assay = assay)
  
  # Segments
  segments <- vrSegments(object, assay = assay)
  if(length(segments) > 0){
    max_vertices <- max(sapply(segments, nrow))
    num_cells <- length(segments)
    segmentations_array <- array(NA, dim = c(num_cells, max_vertices, 2))
    cell_ids <- names(segments)
    for (i in seq_along(cell_ids)) {
      seg <- segments[[i]]
      seg_matrix <- as.matrix(seg[, c("x", "y")])
      nrow_diff <- max_vertices - nrow(seg_matrix)
      if (nrow_diff > 0) {
        seg_matrix <- rbind(seg_matrix, matrix(NA, nrow = nrow_diff, ncol = 2))
      }
      segmentations_array[i, , ] <- seg_matrix
    }
    for (k in 1:2) {
      segmentations_array[,,k] <- t(apply(segmentations_array[,,k], 1, fill_na_with_preceding))
    } 
  } else {
    segmentations_array <- array(dim = nrow(coords))
  }

  # Images
  images_mgk <- vrImages(object, assay = assay, ...)
  if(!is.list(images_mgk)){
    images_mgk <- list(images_mgk)
    names(images_mgk) <- vrAssayNames(object, assay = assay)  
  }
  image_list <- lapply(images_mgk, function(img) {
    list(images = list(hires = as.numeric(magick::image_data(img, channels = "rgb"))),
         scalefactors = list(tissue_hires_scalef = 1, spot_diameter_fullres = 0.5))
  })
  
  # obsm
  # TODO: currently embedding and spatial dimensions should of the same size, but its not 
  # always the case in VoltRon objects
  obsm <- c(obsm, list(spatial = coords, 
                       spatial_AssayID = coords, 
                       segmentation = segmentations_array))
  
  # save as zarr
  if(grepl(".zarr[/]?$", file)){
    
    # check packages
    if(!requireNamespace('DelayedArray'))
      stop("Please install DelayedArray package for using DelayedArray functions")
    
    # run basilisk to call zarr methods
    proc <- basilisk::basiliskStart(py_env)
    on.exit(basilisk::basiliskStop(proc))
    success <- basilisk::basiliskRun(proc, function(data, metadata, obsm, coords, segments, image_list, file) {
      zarr <- reticulate::import("zarr")
      anndata <- reticulate::import("anndata")
      make_numpy_friendly <- function(x) {
        if (DelayedArray::is_sparse(x)) {
          methods::as(x, "dgCMatrix")
        }
        else {
          as.matrix(x)
        }
      }
      X <- make_numpy_friendly(t(data))
      obsm <- list(spatial = coords, 
                   spatial_AssayID = coords, 
                   segmentation = segmentations_array)
      adata <- anndata$AnnData(X = X, 
                               obs = metadata, 
                               obsm = obsm, 
                               uns = list(spatial = image_list))
      adata <- reticulate::r_to_py(adata)
      adata$write_zarr(file)
      return(TRUE)
    }, data = data, metadata = metadata, obsm = obsm, coords = coords, segments = segmentations_array, image_list = image_list, file = file)
    if(create.ometiff){
      success2 <- as.OmeTiff(images_mgk[[1]], out_path = gsub("zarr[/]?$", "ome.tiff", file)) 
      success <- success & success2
    } 
    
    return(success)
    
  # save as h5ad
  } else if(grepl(".h5ad$", file)) {
    
    # Check and use a package for saving h5ad
    if (method == "anndataR") {
      if (!requireNamespace('anndataR', quietly = TRUE)) {
        stop("The anndataR package is not installed. Please install it or choose the 'anndata' method.")
      }
      
      # Create anndata using anndataR
      adata <- anndataR::AnnData(obs_names = rownames(metadata), 
                                 var_names = rownames(data), 
                                 X = t(data), 
                                 obs = metadata, 
                                 obsm = list(spatial = coords, 
                                             spatial_AssayID = coords, 
                                             segmentation = segmentations_array),
                                 uns = list(spatial = image_list))
      
      # Write to h5ad file using anndataR
      anndataR::write_h5ad(adata, path = file)
      
    } else if (method == "anndata") {
      if (!requireNamespace('anndata', quietly = TRUE)) {
        stop("The anndata package is not installed. Please install it or choose the 'anndataR' method.")
      }
      
      # Create anndata using anndata
      adata <- anndata::AnnData(X = t(data), 
                                obs = metadata, 
                                obsm = list(spatial = coords,
                                            spatial_AssayID = coords,
                                            segmentation = segmentations_array),
                                uns = list(spatial = image_list))
      
      
      # Write to h5ad file using anndata
      anndata::write_h5ad(adata, filename = file)
      
    } else {
      stop("Invalid method selected. Please choose either 'anndataR' or 'anndata'.")
    }
    
  } else {
    stop("the 'file' should have an .h5ad, .zarr or .zarr/ extension")
  }
}  

####
# AnnData (Zarr) ####
####

#' @rdname as.Zarr
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
as.Zarr.VoltRon <- function (object, out_path, image_id = "image_1")
{
  
  # check packages
  if(!requireNamespace('DelayedArray'))
    stop("Please install DelayedArray package for using DelayedArray functions")
  
  # object data
  datax <- vrData(object, norm = FALSE)
  metadata <- Metadata(object)

  # obsm
  obsm <- list()
  coords <- vrCoordinates(object)
  obsm[["spatial"]] <- t(as.matrix(coords))
  if (length(vrEmbeddingNames(object)) > 0) {
    for (embed_name in vrEmbeddingNames(object)) {
      embedding <- vrEmbeddings(object, type = embed_name)
      obsm[[embed_name]] <- t(as.matrix(embedding))
    }
  }
  
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(datax, metadata, obsm, out_path) {
    anndata <- reticulate::import("anndata")
    zarr <- reticulate::import("zarr")
    make_numpy_friendly <- function(x, transpose = TRUE) {
      if (transpose) {
        x <- Matrix::t(x)
      }
      if (DelayedArray::is_sparse(x)) {
        methods::as(x, "dgCMatrix")
      }
      else {
        as.matrix(x)
      }
    }
    X <- make_numpy_friendly(datax)
    adata <- anndata$AnnData(X = X, obs = metadata)
    if (length(obsm) > 0) {
      obsm <- lapply(obsm, make_numpy_friendly)
      adata$obsm <- obsm
    }
    adata$write_zarr(out_path)
    return(TRUE)
  }, datax = datax, metadata = metadata, obsm = obsm, out_path = out_path)
  return(success)
}

#' @rdname as.Zarr
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
"as.Zarr.magick-image" <- function (object, out_path, image_id = "image_1")
{
  img_arr <- apply(as.matrix(magick::image_raster(object, tidy = FALSE)), c(1, 2), col2rgb)
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(img_arr, image_id, out_path) {
    zarr <- reticulate::import("zarr")
    ome_zarr <- reticulate::import("ome_zarr")
    z_root <- zarr$open_group(out_path, mode = "w")
    obj_list <- function(...) {
      retval <- stats::setNames(list(), character(0))
      param_list <- list(...)
      for (key in names(param_list)) {
        retval[[key]] = param_list[[key]]
      }
      retval
    }
    default_window <- obj_list(start = 0, min = 0, max = 255, end = 255)
    ome_zarr$writer$write_image(image = img_arr,
                                group = z_root,
                                axes = "cyx",
                                omero = obj_list(name = image_id, version = "0.3",
                                                 rdefs = obj_list(),
                                                 channels = list(obj_list(label = "r", color = "FF0000", window = default_window),
                                                                 obj_list(label = "g", color = "00FF00", window = default_window),
                                                                 obj_list(label = "b", color = "0000FF", window = default_window))))
    return(TRUE)
  }, img_arr = img_arr, image_id = image_id, out_path = out_path)
  return(success)
}

####
# OME ####
####

#' as.OmeTiff
#'
#' Converting VoltRon (magick) images to ome.tiff
#'
#' @param object a magick-image object
#' @param out_path output path to ome.tiff file
#' @param image_id image name
#' 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
as.OmeTiff <- function (object, out_path, image_id = "image_1"){
  
  # get image and transpose the array
  img_arr <- apply(as.matrix(magick::image_raster(object, tidy = FALSE)), c(1, 2), col2rgb)
  img_arr <- aperm(img_arr, c(2,3,1))
  
  # run basilisk
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(img_arr, image_id, out_path, e) {
    
    # set up environment
    e <- new.env()
    options("reticulate.engine.environment" = e)
    
    # get image data to python environment
    img_arr <- reticulate::r_to_py(img_arr)
    assign("img_arr_py", img_arr, envir = e)

    # save ome.tiff
    reticulate::py_run_string(
      paste0("import numpy as np
import tifffile
tifimage = r.img_arr_py.astype('uint8')
# tifimage = np.random.randint(0, 255, (32, 32, 3), 'uint8')
np.savetxt('data.csv', tifimage[:,:,0], delimiter=',')
with tifffile.TiffWriter('", out_path, "') as tif: tif.write(tifimage, photometric='rgb')"
      ))
    
    return(TRUE)
  }, img_arr = img_arr, image_id = image_id, out_path = out_path)
  return(success)
}

#' as.OmeZarr
#'
#' Converting VoltRon (magick) images to ome.tiff
#'
#' @param object a magick-image object
#' @param out_path output path to ome.tiff file
#' @param image_id image name
#' 
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
as.OmeZarr <- function (object, out_path, image_id = "image_1"){
  
  # get image and transpose the array
  img_arr <- apply(as.matrix(magick::image_raster(object, tidy = FALSE)), c(1, 2), col2rgb)

  # run basilisk
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(img_arr, image_id, out_path) {
    zarr <- reticulate::import("zarr")
    ome_zarr <- reticulate::import("ome_zarr")
    z_root <- zarr$open_group(out_path, mode = "w")
    obj_list <- function(...) {
      retval <- stats::setNames(list(), character(0))
      param_list <- list(...)
      for (key in names(param_list)) {
        retval[[key]] = param_list[[key]]
      }
      retval
    }
    default_window <- obj_list(start = 0, min = 0, max = 255, end = 255)
    ome_zarr$writer$write_image(image = img_arr,
                                group = z_root,
                                axes = "cyx",
                                omero = obj_list(name = image_id, version = "0.3",
                                                 rdefs = obj_list(),
                                                 channels = list(obj_list(label = "r", color = "FF0000", window = default_window),
                                                                 obj_list(label = "g", color = "00FF00", window = default_window),
                                                                 obj_list(label = "b", color = "0000FF", window = default_window))))
    return(TRUE)
  }, img_arr = img_arr, image_id = image_id, out_path = out_path)
  return(success)
}

####
# Giotto ####
####

#' as.Giotto
#'
#' Converting a VoltRon object into a Giotto object
#'
#' @param object a VoltRon object
#' @param assay the name of the assay to be converted
#' @param reg if TRUE, registered coordinates will be used
#'
#' @rdname as.Giotto
#'
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_replace str_extract
#' @importFrom magick image_write
#'
#' @export
as.Giotto <- function(object, assay = NULL, reg = FALSE){
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # check Seurat package
  if(!requireNamespace('Giotto'))
    stop("Please install Giotto package!")
  
  # check the number of assays
  if(is.null(assay)){
    if(length(unique(sample_metadata[["Assay"]])) > 1){
      stop("You can only convert a single VoltRon assay into a Seurat object!")
    } else {
      assay <- sample_metadata[["Assay"]]
    }
  } else {
    vrMainAssay(object) <- assay
  }
  
  # check the number of assays
  if(!unique(vrAssayTypes(object, assay = assay)) %in% c("cell")) {
    stop("Conversion of assay types other than cells into Giotto is not yet permitted!")
  }
  
  # data
  rowdata <- vrData(object, assay = assay, norm = FALSE)
  
  # metadata
  metadata <- Metadata(object, assay = assay)
  metadata$cell_ID <- rownames(metadata)
  assays <- stringr::str_extract(rownames(metadata), pattern = "_Assay[0-9]+$")
  assays <- gsub("^_", "", assays)
  
  # coordinates
  coords <- vrCoordinates(object, assay = assay, reg = reg)
  
  # Seurat object
  gio <- Giotto::createGiottoObject(expression = rowdata, 
                                    spatial_locs = coords, 
                                    cell_metadata = metadata)
  
  # get image objects for each assay
  for(assy in vrAssayNames(object)){
    assay_object <- object[[assy]]
    if(vrAssayTypes(assay_object) == "cell"){
      img <- vrImages(assay_object)
      gio_img <- Giotto::createGiottoImage(gio, 
                                   spat_unit = "cell", 
                                   mg_object = img)
      gio <- Giotto::addGiottoImage(gio, images = list(gio_img), spat_loc_name = "cell")
    } else {
      stop("Currently VoltRon does only support converting cell type spatial data sets into SpatialExperiment objects!")
    }
  }
  
  # return
  gio
}

####
# SpatialExperiment ####
####

#' @param type the spatial data type of Seurat object: "image" or "spatial"
#' @param assay_type one of two types, 'cell' or 'spot' etc.
#' @param assay_name the assay name of the voltron assays (e.g. Visium, Xenium etc.)
#' @param image_id select image_id names if needed.
#' @param ... Additional parameter passed to \link{formVoltRon}
#'
#' @rdname as.VoltRon
#' @method as.VoltRon SpatialExperiment
#' 
#' @importFrom magick image_read
#'
#' @export
as.VoltRon.SpatialExperiment <- function(object, assay_type = "cell", assay_name = NULL, image_id = NULL, ...){
  
  # check SpatialExperiment package
  if(!requireNamespace('SpatialExperiment'))
    stop("Please install SpatialExperiment package for using SpatialExperiment objects!: BiocManager::install('SpatialExperiment')")
  
  # raw counts
  data <- SummarizedExperiment::assay(object, i = "counts")

  # metadata
  metadata <- as.data.frame(SummarizedExperiment::colData(object))
  
  # embeddings
  dim_names <- SingleCellExperiment::reducedDimNames(object)
  if(length(dim_names) > 0){
    embeddings_flag <- TRUE
    embedding_list <- sapply(dim_names, function(x) {
      SingleCellExperiment::reducedDim(object, type = x)
    }, USE.NAMES = TRUE)
  } else {
    embeddings_flag <- FALSE
  }
  
  # coords
  coords <- SpatialExperiment::spatialCoords(object)
  colnames(coords) <- c("x", "y")
  
  # img data
  imgdata <- SpatialExperiment::imgData(object)
  
  # image
  voltron_list <- list()
  sample_names <- unique(metadata$sample_id)
  for(samp in sample_names){
    
    # spatial points
    sppoints <- rownames(metadata)[metadata$sample_id == samp]
    sppoints_nopostfix <- gsub("_Assay[0-9]+$", "", sppoints)
    
    # metadata 
    cur_metadata <- metadata[sppoints,]
    
    # data
    cur_data <- data[,sppoints]
    
    # coords
    cur_coords <- coords[sppoints,]
    
    # image
    if(nrow(imgdata) > 0){
      if(is.null(image_id)){
        image_names <- imgdata$image_id[imgdata$sample_id == samp]
      } else {
        image_names <- image_id
      }
      img_list <- sapply(image_names, function(img){ 
        imgraster <- SpatialExperiment::imgRaster(object, 
                                                  sample_id = samp, 
                                                  image_id = img)
        magick::image_read(imgraster)
      }, USE.NAMES = TRUE)
    } else {
      img_list <- NULL
    }
    
    # get params
    if(assay_type == "spot"){
      vis.spot.radius <- 1 
      spot.radius <- 1
    } else {
      params <- list()
    }
    
    # form voltron
    assay_name <- assay_name
    assay_type <- assay_type
    voltron_list[[samp]] <- formVoltRon(data = cur_data, metadata = cur_metadata, coords = cur_coords, 
                                      main.assay = assay_name, image = img_list, params = params, 
                                      assay.type = assay_type, sample_name = samp, ...)
    
    # add embeddings
    spatialpoints <- vrSpatialPoints(voltron_list[[samp]])
    spatialpoints_nopostfix <- stringr::str_replace(spatialpoints, "_Assay[0-9]+$", "")
    spatialpoints_assay <- stringr::str_extract(spatialpoints, "Assay[0-9]+$")
    if(embeddings_flag){
      for(embed_name in names(embedding_list)){
        cur_embedding <- embedding_list[[embed_name]][cells,]
        rownames(cur_embedding) <- spatialpoints
        vrEmbeddings(voltron_list[[samp]], type = embed_name) <- cur_embedding
      }
    }
  }
  
  # merge object
  message("Merging object ...")
  if(length(voltron_list) > 1){
    vrobject <- merge(voltron_list[[1]], voltron_list[-1])
  } else {
    vrobject <- voltron_list[[1]]
  }

  return(vrobject)
}

#' as.SpatialExperiment
#'
#' Converting a VoltRon object into a SpatialExperiment object
#'
#' @param object a VoltRon object
#' @param assay the name of the assay to be converted
#' @param reg if TRUE, registered coordinates will be used
#'
#' @rdname as.SpatialExperiment
#'
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_replace str_extract
#' @importFrom magick image_write
#'
#' @export
as.SpatialExperiment <- function(object, assay = NULL, reg = FALSE){
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # check Seurat package
  if(!requireNamespace('SpatialExperiment'))
    stop("Please install SpatialExperiment package!: BiocManager::install('SpatialExperiment'')")
  
  # check the number of assays
  if(is.null(assay)){
    if(length(unique(sample_metadata[["Assay"]])) > 1){
      stop("You can only convert a single VoltRon assay into a Seurat object!")
    } else {
      assay <- sample_metadata[["Assay"]]
    }
  } else {
    vrMainAssay(object) <- assay
  }
  
  # check the number of assays
  if(unique(vrAssayTypes(object, assay = assay)) %in% c("ROI", "molecule", "tile")) {
    stop("Conversion of ROI, molecule and tile assays into SpatialExperiment is not yet permitted!")
  }
  
  # data
  rawdata <- as.matrix(vrData(object, assay = assay, norm = FALSE))
  
  # metadata
  metadata <- Metadata(object, assay = assay)
  if(is.null(rownames(metadata)))
    rownames(metadata) <- metadata$id
  metadata <- metadata[colnames(rawdata),]
  assays <- stringr::str_extract(rownames(metadata), pattern = "_Assay[0-9]+$")
  assays <- gsub("^_", "", assays)
  
  # Embeddings
  reduceddims <- list()
  if (length(vrEmbeddingNames(object, assay = assay)) > 0) {
    for (embed_name in vrEmbeddingNames(object, assay = assay)) {
      reduceddims[[embed_name]] <- vrEmbeddings(object, assay = assay, type = embed_name)
    }
  }
  
  # coordinates
  coords <- as.matrix(vrCoordinates(flipCoordinates(object, assay = assay), assay = assay, reg = reg))
  coords <- coords[colnames(rawdata),]
  
  # Seurat object
  spe <- SpatialExperiment::SpatialExperiment(assay=list(counts = rawdata),
                                              colData=metadata, 
                                              reducedDims = reduceddims,
                                              sample_id=assays,
                                              spatialCoords=coords)
  
  # get image objects for each assay
  for(assy in vrAssayNames(object)){
    assay_object <- object[[assy]]
    img <- vrImages(assay_object)
    imgfile <- tempfile(fileext='.png')
    magick::image_write(image = img, path = imgfile, format = 'png')
    spe <- SpatialExperiment::addImg(spe,
                                     sample_id = vrAssayNames(assay_object),
                                     image_id = "main",
                                     imageSource = imgfile,
                                     scaleFactor = 1,
                                     load = TRUE)
    file.remove(imgfile)
  }
  
  # return
  spe
}