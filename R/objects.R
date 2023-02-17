# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))

#' The FOVImage class
#'
#' The FOVImage stores the morphology image from the Xenium assay
#'
#' @slot image A magick-image class object from magick package
#'
#' @name FOVImage-class
#' @rdname FOVImage-class
#' @exportClass FOVImage
#'
FOVImage <- setClass(
  Class = 'FOVImage',
  slots = list(
    'image' = 'magick-image'
  )
)

# set method for FOVImage class
setMethod(
  f = 'show',
  signature = 'FOVImage',
  definition = function(object) {
    cat("FOV Morphology Image \n")
    cat(paste(format(image_info(object@image)), collapse = " \n "))
    return(invisible(x = NULL))
  }
)


#' LoadNanostring
#'
#' Load Nanostring CosMx data
#'
#' @param data.dir Path to folder containing Nanostring SMI outputs
#' @param fov FOV name
#' @param assay assay
#' @param type the type of boundaries used for cells
#'
#' @export
#'
LoadNanostring <- function (data.dir, fov, assay = "Nanostring", type = c("centroids", "segmentations"))
{
  data <- ReadNanostring(data.dir = data.dir, type = type)
  segmentations.data <- list()
  if("segmentations" %in% type){
    segmentations.data[["segmentation"]] <- CreateSegmentation(data$segmentations)
  }
  if("centroids" %in% type){
    segmentations.data[["centroids"]] <- CreateCentroids(data$centroids)
  }
  coords <- CreateFOV(coords = segmentations.data, type = type, molecules = data$pixels, assay = assay)
  obj <- CreateSeuratObject(counts = data$matrix, assay = assay)
  if(all(c("centroids", "segmentations") %in% type)){
    cells <- intersect(Cells(x = coords, boundary = "segmentation"),
                       Cells(x = coords, boundary = "centroids"))
  } else {
    cells <- Cells(x = coords, boundary = "centroids")
  }
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[[fov]] <- coords
  return(obj)
}
