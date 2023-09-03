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

  if(!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # raw counts
  rawdata <- as.matrix(object[[Seurat::DefaultAssay(object)]]@counts)

  # metadata
  metadata <- object@meta.data

  # image
  spatialobjectlist <- object@images
  if(length(spatialobjectlist) > 0){
    spatialobject <- spatialobjectlist[[1]]
    if("image" %in% slotNames(spatialobject)){
      image <-  magick::image_read(spatialobject@image)
      info <- image_info(image)
    } else {
      stop("There are no images available in this Seurat object")
    }
  } else{
    stop("There are no images available in this Seurat object")
  }

  # coordinates
  coords <- as.matrix(Seurat::GetTissueCoordinates(object))[,2:1]
  coords[,2] <- info$height - coords[,2]
  colnames(coords) <- c("x", "y")

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
as.Seurat.VoltRon <- function(object){

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
