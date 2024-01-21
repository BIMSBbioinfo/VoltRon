#' @include zzz.R
NULL

####
# Objects and Classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('magick-image'))
setOldClass(Classes = c('raster'))
setOldClass(Classes = c('bitmap'))

## vrAssay ####

#' The vrAssay (VoltRon Assay) Class
#'
#' @slot rawdata raw count table
#' @slot normdata normalized count table
#' @slot featuredata feature metadata
#' @slot embeddings list of embeddings
#' @slot image a list of vrImage objects
#' @slot params additional parameters used by different assay types
#' @slot type the type of the assay (cell, spot, ROI)
#' @slot main_image the key of the main image
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
    image = "list",
    params = "list",
    type = "character",
    name = "character",
    main_image = "character"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrAssay',
  definition = function(object) {
    cat("vrAssay (VoltRon Assay) of", nrow(vrCoordinates(object)), "spatial points and", nrow(object@rawdata), "features. \n")
    return(invisible(x = NULL))
  }
)

####
# Methods ####
####

### Create vrAssay Object ####

#' formAssay
#'
#' Create a vrAssay (VoltRon assay) object
#'
#' @param data the count table
#' @param coords the coordinates of the spatial points
#' @param segments the segments of the spatial points (optional)
#' @param image the image of the data
#' @param params additional parameters of the object
#' @param type the type of the assay (cells, spots, ROIs)
#' @param name the name of the assay
#' @param main_image the name of the main_image
#'
#' @importFrom methods new
#'
#' @export
#'
formAssay <- function(data = NULL, coords, segments = list(), image = NULL, params = list(), type = "ROI", name = "Assay1", main_image = "image_1"){

  # get data
  if(is.null(data)){
    data <- matrix(nrow = 0, ncol = nrow(coords))
    colnames(data) <- rownames(coords)
  }

  # get image object
  image <- formImage(coords = coords, segments = segments, image = image)
  image <- list(image)
  names(image) <- main_image

  # make vrAssay object
  methods::new("vrAssay", rawdata = data, normdata = data,
               image = image, params = params, type = type, name = name, main_image = main_image)
}

### Subset vrAssay objects ####

#' Subsetting vrAssay objects
#'
#' Given a vrAssay object, subset the object given one of the attributes
#'
#' @param object A vrAssay object
#' @param subset Logical statement for subsetting
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#' @param image the subseting string passed to \code{magick::image_crop}
#'
#' @method subset vrAssay
#'
#' @importFrom rlang enquo
# #' @importFrom magick image_crop
#'
#' @export
#'
subset.vrAssay <- function(object, subset, spatialpoints = NULL, features = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- rlang::enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(features)){

    # select features
    nonmatching_features <- setdiff(features, vrFeatures(object))
    features <- intersect(vrFeatures(object), features)

    if(length(features) > 0){
      object@rawdata <- object@rawdata[rownames(object@rawdata) %in% features,, drop = FALSE]
      object@normdata <- object@normdata[rownames(object@normdata) %in% features,, drop = FALSE]
    } else {
      stop("none of the provided features are found in the assay")
    }

    if(length(nonmatching_features))
      message("the following features are not found in the assay: ", paste(nonmatching_features, collapse = ", "))

  } else {

    if(!is.null(spatialpoints)){

      # check if spatial points are here
      if(length(intersect(spatialpoints, vrSpatialPoints(object))) == 0)
        return(NULL)

      # data
      object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% spatialpoints, drop = FALSE]
      object@normdata  <- object@normdata[,colnames(object@normdata) %in% spatialpoints, drop = FALSE]

      # embeddings
      for(embed in vrEmbeddingNames(object)){
        embedding <- vrEmbeddings(object, type = embed)
        vrEmbeddings(object, type = embed) <- embedding[rownames(embedding) %in% spatialpoints,]
      }

      # images
      for(img in vrImageNames(object))
        object@image[[img]] <- subset.vrImage(object@image[[img]], spatialpoints = spatialpoints)

    } else if(!is.null(image)) {

      # images
      for(img in vrImageNames(object))
        object@image[[img]] <- subset.vrImage(object@image[[img]], image = image)
      spatialpoints <- rownames(vrCoordinates(object@image[[img]]))

      # data
      object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% spatialpoints, drop = FALSE]
      object@normdata  <- object@normdata[,colnames(object@normdata) %in% spatialpoints, drop = FALSE]

      # embeddings
      for(embed in vrEmbeddingNames(object)){
        embedding <- vrEmbeddings(object, type = embed)
        vrEmbeddings(object, type = embed) <- embedding[rownames(embedding) %in% spatialpoints,]
      }
    }
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
#' @param crop_info the subseting string passed to \code{magick::image_crop}
#'
subsetCoordinates <- function(coords, image, crop_info){

  # image
  imageinfo <- image_info(image)

  # get crop information
  crop_info <- strsplit(crop_info, split = "\\+")[[1]]
  crop_info <- unlist(lapply(crop_info, function(x) strsplit(x, "x")))
  crop_info <- as.numeric(crop_info)

  # get uncropped spatial points
  xlim <- c(crop_info[3], crop_info[3]+crop_info[1])
  ylim <- c(crop_info[4], crop_info[4]+crop_info[2])
  ylim <- rev(imageinfo$height - ylim)

  # adjust for maximum res
  if(ylim[2] < 0){
    ylim[2] <- 0
    # ylim[1] <- ylim[2] - imageinfo$height + crop_info[2] # CHANGE THIS LATER ?
  }
  if(xlim[2] > imageinfo$width){
    xlim[2] <- imageinfo$width
    # xlim[1] <- xlim[2] - crop_info[1] # CHANGE THIS LATER ?
  }

  # get inside coords
  inside <- (coords[,1] > xlim[1] & coords[,1] < xlim[2]) & (coords[,2] > ylim[1] & coords[,2] < ylim[2])
  coords <- coords[inside,]

  if(nrow(coords) > 0){
    # adjust coordinates
    coords[,1] <- coords[,1] - xlim[1]
    coords[,2] <- coords[,2] - ylim[1]

    # return new coords
    return(coords)
  } else {
    stop("No spatial points remain after cropping!")
  }
}

#' subsetSegments
#'
#' subsetting segments given cropping parameters of a magick image objects
#'
#' @param segments a list of coordinates for each segment
#' @param image the magick image associated with the coordinates
#' @param crop_info the subseting string passed to \code{magick::image_crop}
#'
subsetSegments <- function(segments, image, crop_info){

  # change strategy based on the length of segments
  if(length(segments) < 200) {
    for(i in 1:length(segments)){
      segments[[i]][,c("x","y")] <- subsetCoordinates(segments[[i]][,c("x","y")], image, crop_info)
      # segments[[i]] <- subsetCoordinates(segments[[i]][,c("x","y")], image, crop_info)
    }
  } else {
    segment_names <- rep(names(segments), sapply(segments, nrow, simplify = TRUE))
    segments <- do.call(rbind,segments)
    rownames(segments) <- 1:nrow(segments)
    segments <- data.frame(segments, row_id = rownames(segments))
    cropped_segments <- subsetCoordinates(segments[,c("x","y")], image, crop_info)
    cropped_segments <- data.frame(cropped_segments, cell_id = segments[rownames(cropped_segments),]$cell_id, row_id = rownames(cropped_segments))
    cropped_segments <- cropped_segments %>% right_join(segments[,c("cell_id", "row_id")], by = c("row_id" = "row_id"))
    cropped_segments <- cropped_segments[,c("cell_id.y", "x", "y")]
    colnames(cropped_segments) <- c("cell_id", "x", "y")
    segments <- cropped_segments %>% dplyr::group_split(cell_id)
  }

  segments
}

### Other Methods ####

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrAssay
#'
#' @export
#'
vrSpatialPoints.vrAssay <- function(object, ...) {
  if(ncol(object@rawdata) > 0){
    return(colnames(object@rawdata))
  } else {
    return(rownames(vrCoordinates(object)))
  }
}

#' @param value new spatial points
#'
#' @rdname vrSpatialPoints
#' @method vrSpatialPoints<- vrAssay
#'
#' @export
#'
"vrSpatialPoints<-.vrAssay" <- function(object, ..., value) {

  # data
  if(length(vrSpatialPoints(object)) != length(value)){
    stop("The number of spatial points is not matching with the input")
  } else {
    if(nrow(object@rawdata) > 0){
      colnames(object@rawdata) <- value
      colnames(object@normdata) <- value
    }
  }

  # images
  for(img_name in vrImageNames(object)){
    vrSpatialPoints(object@image[[img_name]]) <- value
  }

  # embeddings
  embeddings <- object@embeddings
  embed_names <- names(embeddings)
  if(length(embed_names) > 0){
    for(type in embed_names){
      if(nrow(embeddings[[type]]) > 0){
        if(nrow(embeddings[[type]]) != length(value)){
          stop("The number of spatial points is not matching with the input")
        } else {
          rownames(embeddings[[type]]) <- value
        }
        object@embeddings[[type]] <- embeddings[[type]]
      }
    }
  }


  # return
  return(object)
}

#' @rdname vrFeatures
#' @method vrFeatures vrAssay
#'
#' @export
#'
vrFeatures.vrAssay <- function(object, ...) {
  return(rownames(object@rawdata))
}

#' @rdname vrFeatureData
#' @method vrFeatureData vrAssay
#'
#' @export
#'
vrFeatureData.vrAssay <- function(object, ...) {
  return(object@featuredata)
}

#' @param value new feature data
#'
#' @rdname vrFeatureData
#' @method vrFeatureData<- vrAssay
#'
#' @export
#'
"vrFeatureData<-.vrAssay" <- function(object, ..., value) {
  object@featuredata <- value
  return(object)
}

#' @rdname vrAssayNames
#' @method vrAssayNames vrAssay
#'
#' @export
#'
vrAssayNames.vrAssay <- function(object, ...) {

  if(.hasSlot(object, name = "name")){
    if(grep("Assay", object@name)){
      return(object@name)
    } else {
      assay_ids <- stringr::str_extract(vrSpatialPoints(object), "Assay[0-9]+$")
      assay_id <- unique(assay_ids)
      return(assay_id)
    }
  } else {
    assay_ids <- stringr::str_extract(vrSpatialPoints(object), "Assay[0-9]+$")
    assay_id <- unique(assay_ids)
    return(assay_id)
  }
}

#' @param value the new postfix name of the assay
#'
#' @rdname vrAssayNames
#' @method vrAssayNames<- vrAssay
#'
#' @importFrom stringr str_replace
#'
"vrAssayNames<-.vrAssay" <- function(object, ..., value){

  # get original assay name
  assayname <- vrAssayNames(object)

  # change assay names
  vrSpatialPoints(object)  <- stringr::str_replace(vrSpatialPoints(object), assayname, value)
  object@name <- value

  # return
  return(object)
}

#' @rdname vrAssayTypes
#' @method vrAssayTypes vrAssay
#'
#' @export
#'
vrAssayTypes.vrAssay <- function(object) {
  return(object@type)
}

#' @rdname vrAssayParams
#' @method vrAssayParams vrAssay
#'
#' @export
#'
vrAssayParams.vrAssay <- function(object, param = NULL) {
  if(!is.null(param)){
    if(param %in% names(object@params)){
      return(object@params[[param]])
    } else {
      warning(param, " not found in the param list")
      return(NULL)
    }
  } else {
    return(object@params)
  }
}

#' @param features features
#' @param norm TRUE if normalized data should be returned
#' @param ... additonal parameters passed to \code{vrImage}
#'
#' @rdname vrData
#' @method vrData vrAssay
#'
#' @importFrom magick image_raster
#'
#' @export
#'
vrData.vrAssay <- function(object, features = NULL, norm = FALSE, ...) {

  # get assay types
  assay.type <- vrAssayTypes(object)

  # for ROIs, cells and spots
  if(assay.type %in% c("ROI", "cell", "spot")){

    # check if there are features
    if(!is.null(features)){
      if(!all(features %in% vrFeatures(object))){
        stop("Some features are not available in the assay!")
      }
      if(norm){
        return(object@normdata[features,,drop = FALSE])
      } else {
        return(object@rawdata[features,,drop = FALSE])
      }

    # if there are no features requested, return the data
    } else {
      if(norm){
        return(object@normdata)
      } else {
        return(object@rawdata)
      }
    }

  # for tiles and molecules
  } else {

    # check if features are requested
    if(!is.null(features)){
      stop("No features are available for tile and molecule assays!")
    } else{

      # for tile only
      if(assay.type == "tile") {
        image_data <- as.numeric(vrImages(object, as.raster = TRUE, ...))
        # image_data <- image_data[,,1]
        image_data <- (0.299 * image_data[,,1] + 0.587 * image_data[,,2] + 0.114 * image_data[,,3])
        image_data <- split_into_tiles(image_data, tile_size = vrAssayParams(object, param = "tile.size"))
        image_data <- sapply(image_data, function(x) return(as.vector(x)))
        image_data <- image_data*255
        rownames(image_data) <- paste0("pixel", 1:nrow(image_data))
        colnames(image_data) <- vrSpatialPoints(object)
        return(image_data)
      # for molecules only
      } else if(assay.type == "molecule"){
        return(matrix(nrow = 0, ncol = 0))
      }
    }
  }
}

#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrCoordinates
#' @method vrCoordinates vrAssay
#'
#' @export
#'
vrCoordinates.vrAssay <- function(object, image_name = NULL, reg = FALSE) {

  # check main image
  if(is.null(image_name)){
    image_name <- vrMainImage(object)
  }

  # check registered coordinates
  if(reg){
    if(!paste0(image_name, "_reg") %in% vrImageNames(object)){
      warning("There are no registered images with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }

  # check coordinates
  if(!image_name %in% vrImageNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  # return coordinates
  return(vrCoordinates(object@image[[image_name]]))
}

#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D coordinates
#'
#' @rdname vrCoordinates
#' @method vrCoordinates<- vrAssay
#'
#' @importFrom methods slot
#'
#' @export
#'
"vrCoordinates<-.vrAssay" <- function(object, image_name = NULL, reg = FALSE, ..., value) {

  # check main image
  if(is.null(image_name)){
    image_name <- vrMainImage(object)
  }

  # check registered coordinates
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }

  # check coordinates
  if(!image_name %in% vrImageNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  vrCoordinates(object@image[[image_name]]) <- value
  return(object)
}

#' @param image_name the key of the image
#' @param ... additional parameters passed to \code{vrCoordinates}
#'
#' @rdname flipCoordinates
#' @method flipCoordinates vrAssay
#'
#' @importFrom magick image_info
#'
#' @export
#'
flipCoordinates.vrAssay <- function(object, image_name = NULL, ...) {
  imageinfo <- magick::image_info(vrImages(object, name = image_name))
  coords <- vrCoordinates(object, image_name = image_name, ...)
  coords[,"y"] <- imageinfo$height - coords[,"y"]
  vrCoordinates(object, image_name = image_name, ...) <- coords
  return(object)
}

#' @param image_name the key of the image associated with the segments
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrSegments
#' @method vrSegments vrAssay
#'
#' @export
#'
vrSegments.vrAssay <- function(object, image_name = NULL, reg = FALSE) {

  # check main image
  if(is.null(image_name)){
    image_name <- vrMainImage(object)
  }

  # check registered segments
  if(reg){
    if(!paste0(image_name, "_reg") %in% vrImageNames(object)){
      warning("There are no registered images with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }

  # check coordinates
  if(!image_name %in% vrImageNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  # return coordinates
  return(vrSegments(object@image[[image_name]]))
}

#' @param image_name the key of the image associated with the segments
#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D segments for each spatial point
#'
#' @rdname vrSegments
#' @method vrSegments<- vrAssay
#'
#' @importFrom methods slot
#' @export
#'
"vrSegments<-.vrAssay" <- function(object, image_name = NULL, reg = FALSE, ..., value) {

  # check main image
  if(is.null(image_name)){
    image_name <- vrMainImage(object)
  }

  # check registered segments
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }

  # check coordinates
  if(!image_name %in% vrImageNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  vrSegments(object@image[[image_name]]) <- value
  return(object)
}

#' @param reg TRUE if registered segments are being updated
#' @param method the method argument passed to \code{base::dist}
#'
#' @rdname vrDistances
#' @method vrDistances vrAssay
#'
#' @importFrom stats dist
#' @export
#'
vrDistances.vrAssay <- function(object, reg = FALSE, method = "euclidean", ...) {
  coords <- vrCoordinates(object, reg = reg, ...)
  return(as.matrix(stats::dist(coords, method = method)))
}

#' @param type the key name for the embedding
#' @param dims the set of dimensions of the embedding data
#'
#' @rdname vrEmbeddings
#' @method vrEmbeddings vrAssay
#'
#' @export
#'
vrEmbeddings.vrAssay <- function(object, type = "pca", dims = 1:30) {

  # embeddings
  embeddings <- object@embeddings
  embedding_names <- names(embeddings)

  # check embeddings and return
  if(!type %in% embedding_names){
    stop("Embedding type ", type, " is not found!")
  } else{
    embedding <- object@embeddings[[type]]
    if(max(dims) > ncol(embedding)){
      dims <- 1:ncol(embedding)
    }
    return(embedding[,dims])
  }
}

#' @param type the key name for the embedding
#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @method vrEmbeddings<- vrAssay
#'
#' @export
#'
"vrEmbeddings<-.vrAssay" <- function(object, type = "pca", ..., value) {
  object@embeddings[[type]] <- value
  return(object)
}

#' @rdname vrEmbeddingNames
#' @method vrEmbeddingNames vrAssay
#'
#' @export
#'
vrEmbeddingNames.vrAssay <- function(object){
  return(names(object@embeddings))
}
