####
# Conversion into a VoltRon class ####
####

#' @param object A Seurat object
#' @param type the spatial data type of Seurat object: "image" or "spatial"
#' @param assay_name the assay name
#' @param ... Additional parameter passed to \code{formVoltRon}
#'
#' @rdname as.VoltRon
#' @method as.VoltRon Seurat
#'
#' @export
#'
as.VoltRon.Seurat <- function(object, type = c("image", "spatial"), assay_name = NULL, ...){

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # raw counts
  rawdata <- as.matrix(object[[Seurat::DefaultAssay(object)]]@counts)

  # metadata
  metadata <- object@meta.data

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

      # count
      cur_rawdata <- rawdata[,cells]

      # metadata
      cur_metadata <- metadata[cells,]

      # coords
      coords <- as.matrix(Seurat::GetTissueCoordinates(spatialobject))[,1:2]
      coords <- apply(coords, 2, as.numeric)
      colnames(coords) <- c("x", "y")

      # from voltron
      params <- list()
      assay.type <- "cell"
      assay_name <- "FOV"
      voltron_list[[fn]] <- formVoltRon(data = cur_rawdata, metadata = cur_metadata, coords = coords, main.assay = assay_name, params = params, assay.type = assay.type, sample_name = fn, ...)
    }

    # merge object
    message("Merging object ...")
    vrobject <- merge(voltron_list[[1]], voltron_list[-1])
  } else{
    image <- NULL
    warning("There are no spatial objects available in this Seurat object")
  }

  return(vrobject)
}

#' convertAnnDataToVoltRon
#'
#' converting AnnData h5ad files to VoltRon objects
#'
#' @param file h5ad file
#' @param AssayID the ID assays in the h5ad file
#' @param ... additional parameters passed to \code{formVoltRon}
#'
#' @importFrom anndata AnnData read_h5ad
#'
#' @export
#'
convertAnnDataToVoltRon <- function(file, AssayID = NULL, ...){

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

####
# Other Packages ####
####

#' @param cell.assay the name(type) of the cell assay to be converted
#' @param molecule.assay the name(type) of the molecule assay to be added to the cell assay in Seurat object
#' @param image_key the name (or prefix) of the image(s)
#' @param type the spatial data type of Seurat object: "image" or "spatial"
#' @param reg if TRUE, registered coordinates will be used
#'
#' @rdname as.Seurat
#' @method as.Seurat VoltRon
#'
#' @importFrom dplyr bind_cols
#' @importFrom stringr str_replace
#'
#' @export
#'
as.Seurat.VoltRon <- function(object, cell.assay = NULL, molecule.assay = NULL, image_key = "fov", type = c("image", "spatial"), reg = FALSE){

  # sample metadata
  sample_metadata <- SampleMetadata(object)

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # check the number of assays
  if(is.null(cell.assay)){
    if(length(unique(sample.metadata[["Assay"]])) > 1){
      stop("You can only convert a single VoltRon assay into a Seurat object!")
    } else {
      cell.assay <- sample.metadata[["Assay"]]
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

#' @param assay the name(type) of the assay to be converted
#' @param file the name of the h5ad file
#' @param image_key the name (or prefix) of the image(s)
#' @param type the spatial data type of Seurat object: "image" or "spatial"
#'
#' @rdname as.AnnData
#' @method as.AnnData VoltRon
#'
#' @importFrom anndata AnnData write_h5ad
#' @importFrom stringr str_extract
#'
#' @export
#'
as.AnnData.VoltRon <- function(object, file, assay = NULL, image_key = "fov", type = c("image", "spatial")){

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # check the number of assays
  if(is.null(assay)){
    if(length(unique(SampleMetadata(object)[["Assay"]])) > 1){
      stop("You can only convert a single VoltRon assay into a Seurat object!")
    } else {
      assay <- SampleMetadata(object)[["Assay"]]
    }
  } else {
    vrMainAssay(object) <- assay
  }

  # check the number of assays
  if(unique(vrAssayTypes(object, assay = assay)) %in% c("spot","ROI")) {
    stop("Conversion of Spot or ROI assays into Seurat is not permitted!")
  }

  # data
  data <- vrData(object, assay = assay, norm = FALSE)

  # metadata
  metadata <- Metadata(object, assay = assay)
  metadata$AssayID <- stringr::str_extract(rownames(metadata), "_Assay[0-9]+$")

  # coordinates
  coords <- vrCoordinates(object, assay = assay, reg = TRUE)

  # create anndata
  adata <- anndata::AnnData(X = t(data), obs = metadata, obsm = list(spatial = coords, spatial_AssayID = coords))

  # create anndata file
  anndata::write_h5ad(adata, filename = file)

  # return
  NULL
}

#' @param object VoltRon image
#' @param out_path output path to ome.zarr
#' @param image_id image name
#'
#' @rdname as.Zarr
#' @method as.Zarr VoltRon
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
#'
as.Zarr.VoltRon <- function (object, out_path, image_id = "main_image")
{
  # object data
  datax <- vrData(object, norm = FALSE)
  metadata <- Metadata(object)
  # feature.metadata <- vrFeatureData(object)

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

#' @param object a VoltRon image
#' @param out_path output path to ome.zarr
#' @param image_id image name
#'
#' @rdname as.Zarr
#' @method as.Zarr magick-image
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
#'
"as.Zarr.magick-image" <- function (object, out_path, image_id = "main_image")
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
