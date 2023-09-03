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

  # coordinates
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

#' @param object A Giotto object
#'
#' @rdname as.VoltRon
#' @method as.VoltRon Giotto
#'
#' @export
#'
as.VoltRon.Giotto <- function(object){

}

####
# Other Packages ####
####

#' @param object A VoltRon object
#'
#' @rdname as.Seurat
#' @method as.Seurat VoltRon
#'
#' @export
#'
as.Seurat.VoltRon <- function(object, image = "fov"){

  # check Seurat package
  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # check the number of assays
  if(nrow(SampleMetadata(object)) > 1)
    stop("You can only convert a single VoltRon assay into a Seurat object!")

  # check the number of assays
  if(vrAssayTypes(object) == "spot") {
    stop("Conversion of Spot assays into Seurat is not permitted!")
  } else {
    assay = "Xenium"
  }

  # data
  data <- vrData(object, norm = FALSE)
  colnames(data) <- gsub("_Assay[0-9]+", "", colnames(data))

  # metadata
  metadata <- Metadata(object)
  rownames(metadata) <- gsub("_Assay[0-9]+", "", rownames(metadata))

  # Seurat object
  seu <- CreateSeuratObject(counts = data, meta.data = metadata, assay = assay)

  # get coordinates
  coords <- vrCoordinates(object, reg = TRUE)
  rownames(coords) <- gsub("_Assay[0-9]+", "", rownames(coords))

  # define image
  image.data <- list(centroids = CreateCentroids(coords))
  image.data <- CreateFOV(coords = image.data, type = c("centroids"), assay = assay)
  seu[[image]] <- image.data

  # return
  seu
}

#' @param object A VoltRon object
#'
#' @rdname as.Giotto
#' @method as.Giotto VoltRon
#'
#' @export
#'
as.Giotto.VoltRon <- function(object){

}
