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
