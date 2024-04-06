#' @include zzz.R
#' @importClassesFrom Matrix dgCMatrix
NULL

####
# Objects and Classes ####
####

## Auxiliary ####

# Set old classes
setClassUnion("data_matrix", members = c("matrix", "dgCMatrix"))

## vrAssay ####

#' The vrAssay (VoltRon Assay) Class
#'
#' @slot rawdata raw count table
#' @slot normdata normalized count table
#' @slot featuredata feature metadata
#' @slot embeddings list of embeddings
#' @slot image a list of vrImage objects
#' @slot params additional parameters used by different assay types
#' @slot type the type of the assay (tile, molecule, cell, spot, ROI)
#' @slot name the assay name
#' @slot main_image the key of the main image
#'
#' @name vrAssay-class
#' @rdname vrAssay-class
#' @exportClass vrAssay
#'
vrAssay <- setClass(
  Class = 'vrAssay',
  slots = c(
    rawdata = 'data_matrix',
    normdata = 'data_matrix',
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
#' @param data the feature matrix of spatialpoints
#' @param coords the coordinates of the spatial points
#' @param segments the list of segments each associated with a spatial point (optional)
#' @param image a singelton or list of images as magick-image objects
#' @param params additional parameters of the object
#' @param type the type of the assay (tile, molecule, cell, spot or ROI)
#' @param name the name of the assay
#' @param main_image the name of the main_image
#' @param ... additional arguements passed to \link{formImage}
#'
#' @importFrom methods new
#'
#' @export
#'
formAssay <- function(data = NULL, coords, segments = list(), image = NULL, params = list(), type = "ROI", name = "Assay1", main_image = "image_1", ...){

  # get data
  if(is.null(data)){
    data <- matrix(nrow = 0, ncol = nrow(coords))
    colnames(data) <- rownames(coords)
  }

  # get image object
  image <- formImage(coords = coords, segments = segments, image = image, ...)
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
#' @param object a vrAssay object
#' @param subset Logical statement for subsetting
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrAssay
#' @order 4
#'
#' @importFrom rlang enquo
#'
#' @export
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
      spatialpoints <- intersect(spatialpoints, vrSpatialPoints(object))
      if(length(spatialpoints) == 0){
        return(NULL)
      }

      # data
      object@rawdata  <- object@rawdata[,spatialpoints, drop = FALSE]
      object@normdata  <- object@normdata[,spatialpoints, drop = FALSE]

      # embeddings
      for(embed in vrEmbeddingNames(object)){
        embedding <- vrEmbeddings(object, type = embed)
        vrEmbeddings(object, type = embed) <- embedding[spatialpoints,, drop = FALSE]
      }

      # image
      # for(img in vrImageNames(object))
      for(img in vrSpatialNames(object))
        object@image[[img]] <- subset.vrImage(object@image[[img]], spatialpoints = spatialpoints)

    } else if(!is.null(image)) {

      # images
      # for(img in vrImageNames(object))
      for(img in vrSpatialNames(object))
        object@image[[img]] <- subset.vrImage(object@image[[img]], image = image)
      spatialpoints <- rownames(vrCoordinates(object@image[[img]]))

      # data
      object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% spatialpoints, drop = FALSE]
      object@normdata  <- object@normdata[,colnames(object@normdata) %in% spatialpoints, drop = FALSE]

      # embeddings
      for(embed in vrEmbeddingNames(object)){
        embedding <- vrEmbeddings(object, type = embed)
        vrEmbeddings(object, type = embed) <- embedding[rownames(embedding) %in% spatialpoints,, drop = FALSE]
      }
    } else {

      # else return empty
      return(NULL)
    }
  }

  # set VoltRon class
  return(object)
}

#' subsetCoordinates
#'
#' subsetting coordinates given cropping parameters of a magick image objects
#'
#' @param coords the coordinates of the spatial points
#' @param image the magick image associated with the coordinates
#' @param crop_info the subseting string passed to \link{image_crop}
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
#' @param segments the list of segments each associated with a spatial point
#' @param image the magick image associated with the coordinates
#' @param crop_info the subseting string passed to \link{image_crop}
#'
subsetSegments <- function(segments, image, crop_info){

  # change strategy based on the length of segments
  if(length(segments) < 200) {
    for(i in 1:length(segments)){
      segments[[i]][,c("x","y")] <- subsetCoordinates(segments[[i]][,c("x","y")], image, crop_info)
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
#' @order 4
#' 
#' @export
vrSpatialPoints.vrAssay <- function(object) {
  if(ncol(object@rawdata) > 0){
    return(colnames(object@rawdata))
  } else {
    return(rownames(vrCoordinates(object)))
  }
}

#' @rdname vrSpatialPoints
#' @order 8
#' @export
"vrSpatialPoints<-.vrAssay" <- function(object, value) {

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
  # for(img in vrImageNames(object))
  for(img in vrSpatialNames(object)){
    vrSpatialPoints(object@image[[img]]) <- value
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
#' @order 3
#' @export
vrFeatures.vrAssay <- function(object) {
  return(rownames(object@rawdata))
}

#' @rdname vrFeatureData
#' @order 3
#' @export
vrFeatureData.vrAssay <- function(object) {
  return(object@featuredata)
}

#' @rdname vrFeatureData
#' @order 5
#' @export
"vrFeatureData<-.vrAssay" <- function(object, value) {
  object@featuredata <- value
  return(object)
}

#' @rdname vrAssayNames
#' @order 4
#' @export
vrAssayNames.vrAssay <- function(object) {

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

#' @param value assay name
#' 
#' @rdname vrAssayNames
#' @order 5
#' @importFrom stringr str_replace
"vrAssayNames<-.vrAssay" <- function(object, value){

  # get original assay name
  assayname <- vrAssayNames(object)

  # change assay names
  spatialpoints <- stringr::str_replace(vrSpatialPoints(object), assayname, value)

  # add assay name if missing
  if(vrAssayTypes(object) %in% c("ROI", "cell", "spot")){
    ind <- !grepl("Assay[0-9]+$", spatialpoints)
    spatialpoints[ind] <- stringr::str_replace(spatialpoints[ind], "$", paste0("_", value))
  }

  # replace spatial point names
  vrSpatialPoints(object) <- spatialpoints
  object@name <- value

  # return
  return(object)
}

#' @rdname vrAssayTypes
#' @order 3
#' @export
vrAssayTypes.vrAssay <- function(object) {
  return(object@type)
}

#' Get assay parameters
#'
#' Given a vrAssay object, if there are any, get a list of parameters of the assay(s)
#'
#' @param object a vrAssay object
#' @param param the parameter value to return
#'
#' @rdname vrAssayParams
#'
#' @export
vrAssayParams <- function(object, param = NULL) {
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

#' @param features the set of features
#' @param norm TRUE if normalized data should be returned
#' 
#' @rdname vrData
#' @order 3
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

#' @rdname vrCoordinates
#' @order 3
#' @export
#'
vrCoordinates.vrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE) {

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    # image_name <- vrMainImage(object)
    image_name <- vrMainSpatial(object)
  }

  # check registered coordinates
  if(reg){
    # if(!paste0(image_name, "_reg") %in% vrImageNames(object)){
    if(!paste0(image_name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }

  # check coordinates
  # if(!image_name %in% vrImageNames(object)){
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  # return coordinates
  return(vrCoordinates(object@image[[image_name]]))
}

#' @rdname vrCoordinates
#' @order 5
#' @importFrom methods slot
#'
#' @export
"vrCoordinates<-.vrAssay" <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    # image_name <- vrMainImage(object)
    image_name <- vrMainSpatial(object)
  }

  # check registered coordinates
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }

  # check coordinates
  # if(!image_name %in% vrImageNames(object)){
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  vrCoordinates(object@image[[image_name]]) <- value
  return(object)
}

#' @rdname flipCoordinates
#' @order 3
#'
#' @importFrom magick image_info
#'
#' @export
#'
flipCoordinates.vrAssay <- function(object, image_name = NULL, spatial_name = NULL, ...) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # get image info
  imageinfo <- magick::image_info(vrImages(object, name = image_name))
  
  # flip coordinates
  coords <- vrCoordinates(object, image_name = image_name, ...)
  coords[,"y"] <- imageinfo$height - coords[,"y"]
  vrCoordinates(object, image_name = image_name, ...) <- coords
  
  # flip segments
  segments <- vrSegments(object, image_name = image_name, ...)
  if(length(segments) > 0){
    name_segments <- names(segments)
    segments <- do.call("rbind", segments)
    segments[,"y"] <- imageinfo$height - segments[,"y"]
    segments <- split(segments, segments[,1])
    names(segments) <- name_segments
    vrSegments(object, image_name = image_name, ...) <- segments
  }
  
  # return
  return(object)
}

#' @rdname vrSegments
#' @order 3
#' @export
vrSegments.vrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE) {

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    # image_name <- vrMainImage(object)
    image_name <- vrMainSpatial(object)
  }

  # check registered segments
  if(reg){
    # if(!paste0(image_name, "_reg") %in% vrImageNames(object)){
    if(!paste0(image_name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }

  # check coordinates
  # if(!image_name %in% vrImageNames(object)){
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  # return coordinates
  return(vrSegments(object@image[[image_name]]))
}

#' @rdname vrSegments
#' @order 6
#' @importFrom methods slot
#' @export
"vrSegments<-.vrAssay" <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    # image_name <- vrMainImage(object)
    image_name <- vrMainSpatial(object)
  }

  # check registered segments
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }

  # check coordinates
  # if(!image_name %in% vrImageNames(object)){
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any image in this vrAssay object")
  }

  vrSegments(object@image[[image_name]]) <- value
  return(object)
}

#' @param type the key name for the embedding, i.e. "pca" or "umap"
#' @param dims the set of dimensions of the embedding data
#'
#' @rdname vrEmbeddings
#' @order 3
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
    return(embedding[,dims, drop = FALSE])
  }
}

#' @rdname vrEmbeddings
#' @order 4
#' @export
"vrEmbeddings<-.vrAssay" <- function(object, type = "pca", value) {
  object@embeddings[[type]] <- value
  return(object)
}

#' @rdname vrEmbeddingNames
#' @order 3
#'
#' @export
vrEmbeddingNames.vrAssay <- function(object){
  return(names(object@embeddings))
}
