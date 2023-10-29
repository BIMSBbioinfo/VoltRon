####
# Conversion into a VoltRon class ####
####

#' @param object A Seurat object
#' @param ... Additional parameter passed to \code{formVoltRon}
#'
#' @rdname as.VoltRon
#' @method as.VoltRon Seurat
#'
#' @export
#'
as.VoltRon.Seurat <- function(object, ...){

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # raw counts
  rawdata <- as.matrix(object[[Seurat::DefaultAssay(object)]]@counts)

  # metadata
  metadata <- object@meta.data

  # coordinates and subcellular
  if(grepl("Visium", class(object@images[[1]]))){
    coords <- as.matrix(Seurat::GetTissueCoordinates(object))[,2:1]
    colnames(coords) <- c("x", "y")
  } else{
    coords <- as.matrix(Seurat::GetTissueCoordinates(object))[,1:2]
    coords <- apply(coords, 2, as.numeric)
    colnames(coords) <- c("x", "y")
  }

  # image
  spatialobjectlist <- object@images
  if(length(spatialobjectlist) > 0){
    spatialobject <- spatialobjectlist[[1]]
    if("image" %in% slotNames(spatialobject)){
      image <-  magick::image_read(spatialobject@image)
      info <- image_info(image)
      coords[,2] <- info$height - coords[,2]
    } else {
      image <- NULL
      warning("There are no images available in this Seurat object")
    }
  } else{
    image <- NULL
    warning("There are no images available in this Seurat object")
  }

  # scale coordinates and assay.type
  if(grepl("Visium", class(object@images[[1]]))){
    params <- list(spot.radius = Seurat::Radius(object@images[[1]])*max(info$width, info$height))
    assay.type <- "spot"
    assay_name <- "Visium"
  } else{
    params <- list()
    assay.type <- "cell"
    assay_name <- "Xenium"
  }

  # create VoltRon
  formVoltRon(rawdata, metadata, image, coords, main.assay = assay_name, params = params, assay.type = assay.type, ...)
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
convertAnnDataToVoltRon <- function(file, AssayID = NULL, Sample = NULL, ...){

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

#' @param assay the name(type) of the assay to be converted
#' @param image_key the name (or prefix) of the image(s)
#' @param type the spatial data type of Seurat object: "image" or "spatial"
#'
#' @rdname as.Seurat
#' @method as.Seurat VoltRon
#'
#' @export
#'
as.Seurat.VoltRon <- function(object, assay = NULL, image_key = "fov", type = c("image", "spatial")){

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
  # colnames(data) <- gsub("_", "-", colnames(data))

  # metadata
  metadata <- Metadata(object, assay = assay)
  # rownames(metadata) <- gsub("_", "-", rownames(metadata))

  # Seurat object
  seu <- Seurat::CreateSeuratObject(counts = data, meta.data = metadata, assay = assay)

  # get image objects for each assay
  for(assy in vrAssayNames(object)){
    assay_object <- object[[assy]]
    if(type == "image"){
      coords <- vrCoordinates(assay_object, reg = TRUE)
      # rownames(coords) <- gsub("_", "-", rownames(coords))
      image.data <- list(centroids = SeuratObject::CreateCentroids(coords))
      subcellular <- vrSubcellular(assay_object, reg = TRUE)
      if(nrow(subcellular) > 0){
        colnames(subcellular)[colnames(subcellular) %in% "feature_name"] <- "gene"
      } else {
        subcellular <- NULL
      }
      image.data <- SeuratObject::CreateFOV(coords = image.data, type = c("centroids"), molecules = subcellular, assay = assay)
      image <- paste0(image_key, "_", assy)
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
