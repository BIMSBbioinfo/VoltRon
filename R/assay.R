####
# Objects and Classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))

## vrAssay ####

#' The vrAssay (VoltRon Assay) Class
#'
#' @slot rawdata raw count table
#' @slot normdata normalized count table
#' @slot coord spatial coordinates of the assay
#' @slot coord_reg spatial coordinates of the registered assay
#' @slot segments spatial coordinates of the segments
#' @slot segments_reg spatial coordinates of the registered segments
#' @slot image image of the spatial assay
#' @slot params additional parameters used by different assay types
#' @slot type the type of the assay (cell, spot, ROI)
#'
#' @name vrAssay-class
#' @rdname vrAssay-class
#' @exportClass vrAssay
#'
vrAssay <- setClass(
  Class = 'vrAssay',
  slots = c(
    rawdata = 'matrix',
    normdata = 'matrix',
    featuredata = 'data.frame',
    embeddings = "list",
    coords = 'matrix',
    coords_reg = 'matrix',
    segments = 'list',
    segments_reg = 'list',
    image = 'magick-image',
    params = "list",
    type = "character"
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'vrAssay',
  definition = function(object) {
    cat("vrAssay (VoltRon Assay) of", ncol(object@rawdata), "spatial entities and", nrow(object@rawdata), "features. \n")
    return(invisible(x = NULL))
  }
)

####
# Methods ####
####

### Subset vrAssay objects ####

#' @method subset vrAssay
#'
#' @importFrom rlang enquo
#'
#' @export
#'
subset.vrAssay <- function(object, subset, entities = NULL, features = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(entities)){

    # data
    object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% entities]
    object@normdata  <- object@normdata[,colnames(object@normdata) %in% entities]

    # coordinates
    object@coords  <- object@coords[rownames(object@coords) %in% entities,]
    if(nrow(object@coords_reg) > 0)
      object@coords_reg  <- object@coords_reg[rownames(object@coords_reg) %in% entities,]

    # segments
    if(length(object@segments) > 0)
      object@segments  <- object@segments[names(object@segments) %in% entities]
    if(length(object@segments_reg) > 0)
      object@segments_reg  <- object@segments_reg[names(object@segments_reg) %in% entities]

  } else if(!is.null(features)){

    # select features
    object@rawdata <- object@rawdata[rownames(object@rawdata) %in% features,]
    object@normdata <- object@normdata[rownames(object@normdata) %in% features,]

  } else if(!is.null(image)) {

    # coordinates
    cropped_coords <- subsetCoordinates(object@coords, object@image, image)
    object@coords <- cropped_coords

    # segments
    cropped_segments <- object@segments[rownames(cropped_coords)]
    if(length(object@segments) > 0){
      object@segments[rownames(cropped_coords)] <- subsetSegments(cropped_segments, object@image, image)
    }

    # image
    object <- subset.vrAssay(object, entities = rownames(cropped_coords))
    object@image <- image_crop(object@image, image)

  }

  # set VoltRon class
  return(object)
}

#' subsetCoordinates
#'
#' subsetting coordinates given cropping parameters of a magick image objects
#'
#' @param coords coordinates
#' @param image the magick image associated with the coordinates
#' @param crop_info image crop parameters in character, see \code{image_crop}
#'
subsetCoordinates <- function(coords, image, crop_info){

  # image
  imageinfo <- image_info(image)

  # get crop information
  crop_info <- strsplit(crop_info, split = "\\+")[[1]]
  crop_info <- unlist(lapply(crop_info, function(x) strsplit(x, "x")))
  crop_info <- as.numeric(crop_info)

  # get uncropped entities
  xlim <- c(crop_info[3], crop_info[3]+crop_info[1])
  ylim <- c(crop_info[4], crop_info[4]+crop_info[2])
  ylim <- rev(imageinfo$height - ylim)

  # adjust for maximum res
  if(ylim[2] < 0){
    ylim[2] <- 0
    ylim[1] <- ylim[2] - imageinfo$height + crop_info[2]
  }
  if(xlim[2] > imageinfo$width){
    xlim[2] <- imageinfo$width
    xlim[1] <- xlim[2] - crop_info[1]
  }

  # get inside coords
  inside <- (coords[,1] > xlim[1] & coords[,1] < xlim[2]) & (coords[,2] > ylim[1] & coords[,2] < ylim[2])
  coords <- coords[inside,]

  # adjust coordinates
  coords[,1] <- coords[,1] - xlim[1]
  coords[,2] <- coords[,2] - ylim[1]

  # return new coords
  return(coords)
}

#' subsetSegments
#'
#' subsetting segments given cropping parameters of a magick image objects
#'
#' @param coords coordinates
#' @param image the magick image associated with the coordinates
#' @param crop_info image crop parameters in character, see \code{image_crop}
#'
subsetSegments <- function(segments, image, crop_info){
  for(i in 1:length(segments)){
    if(nrow(segments[[i]]) > 1){
      segments[[i]] <- subsetCoordinates(segments[[i]], image, crop_info)
    } else {
      segments[[i]][,c("x","y")] <- subsetCoordinates(segments[[i]][,c("x","y")], image, crop_info)
    }
  }
  segments
}

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrAssay
#'
#' @export
#'
vrSpatialPoints.vrAssay <- function(object, ...) {
  colnames(object@rawdata)
}

#' @rdname vrFeatures
#' @method vrFeatures vrAssay
#'
#' @export
#'
vrFeatures.vrAssay <- function(object, ...) {
  return(rownames(object@rawdata))
}

#' @rdname FeatureData
#' @method FeatureData vrAssay
#'
#' @export
#'
FeatureData.vrAssay <- function(object, ...) {
  return(object@featuredata)
}

#' @rdname FeatureData
#' @method FeatureData<- vrAssay
#'
#' @export
#'
"FeatureData<-.vrAssay" <- function(object, ..., value) {
  object@featuredata <- value
  return(object)
}

#' @rdname AssayNames
#' @method AssayNames vrAssay
#'
#' @export
#'
AssayNames.vrAssay <- function(object, ...) {
  assay_ids <- stringr::str_extract(vrSpatialPoints(object), "Assay[0-9]+")
  assay_id <- unique(assay_ids)
  return(assay_id)
}

#' @rdname AssayNames
#' @method AssayNames<- vrAssay
#'
#' @export
#'
"AssayNames<-.vrAssay" <- function(object, ..., value){

  # change assay names
  colnames(object@rawdata) <- gsub("Assay[0-9]+$", value, colnames(object@rawdata))
  colnames(object@normdata) <- gsub("Assay[0-9]+$", value, colnames(object@normdata))

  # coordinates
  rownames(object@coords)  <- gsub("Assay[0-9]+$", value, rownames(object@coords))
  if(nrow(object@coords_reg) > 0)
    rownames(object@coords_reg) <- gsub("Assay[0-9]+$", value, rownames(object@coords_reg))

  # segments
  if(length(object@segments) > 0)
    names(object@segments) <- gsub("Assay[0-9]+$", value, names(object@segments))
  if(length(object@segments_reg) > 0)
    names(object@segments_reg) <- gsub("Assay[0-9]+$", value, names(object@segments_reg))

  # return
  return(object)
}

#' @rdname AssayTypes
#' @method AssayTypes vrAssay
#'
#' @export
#'
AssayTypes.vrAssay <- function(object, ...) {
  return(object@type)
}

#' @rdname vrData
#' @method vrData vrAssay
#'
#' @export
#'
vrData.vrAssay <- function(object, norm = FALSE) {
  if(norm){
    return(object@normdata)
  } else {
    return(object@rawdata)
  }
}

#' @rdname Coordinates
#' @method Coordinates vrAssay
#'
#' @export
#'
Coordinates.vrAssay <- function(object, reg = FALSE, ...) {
  if(reg){
    if(nrow(object@coords_reg) < 1) {
      return(object@coords)
    } else {
      return(object@coords_reg)
    }
  } else {
    return(object@coords)
  }
}

#' @rdname Coordinates
#' @method Coordinates<- vrAssay
#'
#' @export
#'
"Coordinates<-.vrAssay" <- function(object, reg = FALSE, ..., value) {

  # get coordinates
  coords <- Coordinates(object, ...)

  # stop if the rownames are not matching
  if(any(sapply(rownames(values),is.null)))
    stop("Provided coordinates data does not have cell/spot/ROI names")

  if(!all(rownames(values) %in% rownames(coords)))
    stop("Cant overwrite coordinates, non-existing cells/spots/ROIs!")

  # stop if the colnames are not matching
  if(!all(colnames(values) %in% colnames(coords)))
    stop("Cant overwrite coordinates, only x or y coordinates should be provided!")

  if(reg){
    slot(object = object, name = 'coords_reg') <- value
  } else{
    slot(object = object, name = 'coords') <- value
  }
  return(object)
}

#' @rdname Segments
#' @method Segments vrAssay
#'
#' @export
#'
Segments.vrAssay <- function(object, reg = FALSE, ...) {
  if(reg){
    if(length(object@segments_reg) < 1) {
      return(object@segments)
    } else {
      return(object@segments_reg)
    }
  } else {
    return(object@segments)
  }
}

#' @rdname Segments
#' @method Segments<- vrAssay
#'
#' @export
#'
"Segments<-.vrAssay" <- function(object, reg = FALSE, ..., value) {

  # get coordinates
  segts <- Segments(object, ...)

  # stop if the names are not matching
  if(any(sapply(names(values),is.null)))
    stop("Provided coordinates data does not have cell/spot/ROI names")

  if(!all(names(values) %in% names(coords)))
    stop("Cant overwrite coordinates, non-existing cells/spots/ROIs!")

  if(reg){
    slot(object = object, name = 'segments_reg') <- value
  } else{
    slot(object = object, name = 'segments') <- value
  }
  return(object)
}

#' @rdname Distances
#' @method Distances vrAssay
#'
#' @export
#'
Distances.vrAssay <- function(object, reg = FALSE, method = "euclidean", ...) {
  coords <- Coordinates(object, reg = reg, ...)
  return(as.matrix(dist(coords, method = method)))
}

#' @rdname Embeddings
#' @method Embeddings vrAssay
#'
#' @export
#'
Embeddings.vrAssay <- function(object, type = "pca") {
  return(object@embeddings[[type]])
}

#' @rdname Embeddings
#' @method Embeddings<- vrAssay
#'
#' @export
#'
"Embeddings<-.vrAssay" <- function(object, type = "pca", ..., value) {
  object@embeddings[[type]] <- value
  return(object)
}
