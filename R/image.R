####
# Objects and Classes ####
####

# Set class union
suppressWarnings({
  setClassUnion("image_matrix", 
                members = c("matrix", 
                            "data.frame",
                            "dgRMatrix", 
                            "dgeMatrix", 
                            "Array", 
                            if (requireNamespace("BPCells", quietly = TRUE)) "IterableMatrix" else NULL))
})

## vrImage ####

#' The vrImage (VoltRon Image) Class
#'
#' @slot coords spatial coordinates of the assay
#' @slot segments spatial coordinates of the segments, if available
#' @slot image image of the spatial assay, bitmap class
#' @slot main_channel the key of the main channel of vrImage object
#'
#' @name vrImage-class
#' @rdname vrImage-class
#' @exportClass vrImage
#'
vrImage <- setClass(
  Class = 'vrImage',
  slots = c(
    coords = 'image_matrix',
    segments = 'list',
    image = "list",
    main_channel = "character"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrImage',
  definition = function(object) {

    # separate names
    image_names <- names(object@image)
    image_id <- seq_along(image_names)
    image_names_split <- split(image_names, ceiling(image_id/10))

    cat("vrImage (VoltRon Image) Object \n")
    text <- "Channels:"
    for(img in image_names_split){
      cat(text, paste(img, collapse = ", "), "\n")
      text <- "         "
    }
    return(invisible(x = NULL))
  }
)

## vrSpatial ####

#' The vrSpatial (VoltRon Spatial) Class
#'
#' @slot coords spatial coordinates of the assay
#' @slot segments spatial coordinates of the segments, if available
#' @slot image image of the spatial assay, bitmap class
#' @slot main_channel the key of the main channel of vrImage object
#'
#' @name vrSpatial-class
#' @rdname vrSpatial-class
#' @exportClass vrSpatial
#'
vrSpatial <- setClass(
  Class = 'vrSpatial',
  slots = c(
    coords = 'image_matrix',
    segments = 'list',
    image = "list",
    main_channel = "character"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrSpatial',
  definition = function(object) {
    
    # separate names
    image_names <- names(object@image)
    image_id <- seq_along(image_names)
    image_names_split <- split(image_names, ceiling(image_id/10))
    
    cat("vrSpatial (VoltRon Spatial) Object \n")
    text <- "Channels:"
    for(img in image_names_split){
      cat(text, paste(img, collapse = ", "), "\n")
      text <- "         "
    }
    return(invisible(x = NULL))
  }
)

####
# Create vrImage Object ####
####

#' formImage
#'
#' Create a vrImage (VoltRon image) object
#'
#' @param coords the coordinates of the spatial points
#' @param segments the list of segments each associated with a spatial point
#' @param image a singelton or list of images as magick-image objects
#' @param main_channel the key of the main channel of vrImage object
#'
#' @importFrom magick image_data image_read image_info
#' @importFrom methods new
#'
#' @export
#'
formImage <- function(coords, segments = list(), image = NULL, main_channel = NULL){

  # get coordinates
  if(inherits(coords, "data.frame")){
    coords <- as.matrix(coords)
  }
  if(!inherits(coords, c("matrix", "dgCMatrix", "Matrix", "IterableMatrix"))){
    stop("Coordinates table should either of a matrix or data.frame class!")
  }
  if(ncol(coords) == 2){
    coords <- cbind(coords,0)
    colnames(coords) <- c("x", "y", "z")
  }
  if(!ncol(coords) %in% c(2,3)){
    stop("The length of colnames of the coordinates matrix should be either two or three!")
  } 

  # get segments
  if(length(segments) > 0){
    if(length(segments) == length(rownames(coords))){
      names(segments) <- rownames(coords)
    } else {
      stop("Number of segments doesnt match the number of points!")
    }
  }

  # check if the image input is a list
  if(!is.null(image)){
    if(is.list(image)){

      # enter names if there are no names
      if(is.null(names(image)))
        names(image) <- paste("channel_", 1:length(image), sep = "")

      # get image information
      imageinfo <- sapply(image, function(x) magick::image_info(x)[,c("width", "height")], USE.NAMES = TRUE)
      flag <- all(apply(imageinfo, 1, function(x) length(unique(x)) == 1))

      #
      if(!flag){
        stop("When providing multiple images as channels, make sure that all images have the same dimensionality!")
      } else {
        image <- lapply(image, magick::image_data)
        names(image) <- colnames(imageinfo)
        if(is.null(main_channel))
          main_channel <- names(image)[1]
      }
    } else {
      image <- list(magick::image_data(image))
      if(is.null(main_channel))
        main_channel <- "channel_1"
      names(image) <- main_channel
    }
  } else {
    image <- list()
    main_channel <- ""
  }

  # make vrimage object
  methods::new("vrSpatial", coords = coords, segments = segments, image = image, main_channel = main_channel)
}

### Subset vrImage objects ####

#' Subsetting vrImage objects
#'
#' Given a vrImage object, subset the object given one of the attributes.
#'
#' @param object A vrImage object
#' @param subset Logical statement for subsetting
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrImage
#' @order 5
#'
#' @importFrom rlang enquo
#' @importFrom magick image_crop
#'
#' @export
#'
subset.vrImage <- function(object, subset, spatialpoints = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- rlang::enquo(arg = subset)
  }

  # coords and segments
  coords <- vrCoordinates(object)
  segments <- vrSegments(object)

  if(!is.null(spatialpoints)){

    # check if spatial points are here
    spatialpoints <- intersect(spatialpoints, rownames(coords))
    if(length(spatialpoints) == 0){
      return(NULL)
    }

    # coordinates
    vrCoordinates(object) <- coords[spatialpoints,, drop = FALSE]

    # segments
    if(length(segments) > 0)
      vrSegments(object) <- segments[spatialpoints]

  } else if(!is.null(image)) {

    # get one image
    vrimage <- vrImages(object)

    # coordinates
    cropped_coords <- subsetCoordinates(coords, vrimage, image)
    vrCoordinates(object) <- cropped_coords

    # segments
    cropped_segments <- segments[rownames(cropped_coords)]
    if(length(segments) > 0){
      segments[rownames(cropped_coords)] <- subsetSegments(cropped_segments, vrimage, image)
      vrSegments(object) <- segments
    }

    # spatial points
    object <- subset.vrImage(object, spatialpoints = rownames(cropped_coords))

    # image
    for(img in vrImageChannelNames(object)){
      
      # check if the image is either ondisk or inmemory
      img_data <- object@image[[img]]
      if(inherits(img_data, "Image_Array")){
        crop_info_int <- as.integer(strsplit(image, split = "[x|+]")[[1]])
        img_data <- ImageArray::crop(img_data, ind = list(crop_info_int[3]:(crop_info_int[3]+crop_info_int[1]), crop_info_int[4]:(crop_info_int[4]+crop_info_int[2])))
        object@image[[img]] <- img_data
      } else {
        img_data <- magick::image_read(img_data)
        img_data <- magick::image_crop(img_data, image)
        object@image[[img]] <- magick::image_data(img_data) 
      }
    }
  }

  # set VoltRon class
  return(object)
}

#' Subsetting vrSpatial objects
#'
#' Given a vrSpatial object, subset the object given one of the attributes.
#'
#' @param object A vrSpatial object
#' @param subset Logical statement for subsetting
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrSpatial
#' @order 5
#'
#' @importFrom rlang enquo
#' @importFrom magick image_crop
#'
#' @export
#'
subset.vrSpatial <- function(object, subset, spatialpoints = NULL, image = NULL){
  subset.vrImage(object, subset = subset, spatialpoints = spatialpoints, image = image)
}

####
# Methods ####
####

vrImagesVoltRon <- function(object, assay = NULL, name = NULL, reg = FALSE, channel = NULL, as.raster = FALSE, scale.perc = 100){
  
  # get assay names
  if(is.null(assay)){
    assay_names <- vrAssayNames(object, assay = "all")
  } else {
    assay_names <- vrAssayNames(object, assay = assay)
  }
  
  # get images
  images <- sapply(assay_names, function(assy) vrImages(object[[assy]], 
                                                        name = name, 
                                                        reg = reg, 
                                                        channel = channel,
                                                        as.raster = as.raster, 
                                                        scale.perc = scale.perc), USE.NAMES = TRUE)
  if(length(images) == 1){
    return(images[[1]])
  } else {
    return(images)
  }
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param name the name of the main spatial system
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainSpatial}) is requested
#' @param channel the name of the channel associated with the image
#' @param as.raster return as raster
#' @param scale.perc scale percentage if lower resolution image needed
#'
#' @rdname vrImages
#' @order 2
#' @export
setMethod("vrImages", "VoltRon", vrImagesVoltRon)

vrImagesvrAssay <- function(object, name = NULL, reg = FALSE, channel = NULL, as.raster = FALSE, scale.perc = 100){
  
  # check image name
  if(is.null(name)) {
    name <- object@main_image
  }
  
  # get registered image
  if(reg){
    if(!paste0(name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", name, "!")
    } else {
      name <- paste0(name, "_reg")
    }
  }
  
  # check main image
  if(!name %in% vrSpatialNames(object)){
    stop(name, " is not among any image in this vrAssay object")
  }
  
  return(vrImages(object@image[[name]], channel = channel, as.raster = as.raster, scale.perc = scale.perc))
}

#' @rdname vrImages
#' @order 3
#' @export
setMethod("vrImages", "vrAssay", vrImagesvrAssay)

#' @rdname vrImages
#' @order 3
#' @export
setMethod("vrImages", "vrAssayV2", vrImagesvrAssay)

vrImagesReplacevrAssay <- function(object, name = NULL, channel = NULL, reg = FALSE, value) {
  if(is.null(name)) {
    name <- object@main_image
  }
  
  if(reg){
    name <- paste0(name, "_reg")
  }
  
  if(inherits(value, "vrImage") | inherits(value, "vrSpatial")){
    object@image[[name]] <- value
  } else {
    if(!is.null(channel)){
      vrImages(object@image[[name]], channel = channel) <- value
    }
  }
  return(object)
}

#' @param value new image
#' 
#' @rdname vrImages
#'
#' @importFrom magick image_data
#' @order 5
#' @export
setMethod("vrImages<-", "vrAssay", vrImagesReplacevrAssay)

#' @param value new image
#' 
#' @rdname vrImages
#'
#' @importFrom magick image_data
#' @order 5
#' @export
setMethod("vrImages<-", "vrAssayV2", vrImagesReplacevrAssay)

vrImagesvrImage <- function(object, channel = NULL, as.raster = FALSE, scale.perc = 100){
  
  # check channels
  if(is.null(channel)){
    channel <- object@main_channel
  } else {
    if(!channel %in% vrImageChannelNames(object)){
      warning("'", channel, "' is not among any channel in this vrImage object!")
      return(NULL)
    }
  }
  
  # correct image scale
  if(!is.numeric(scale.perc)){
    stop("scale.perc should be between 0 and 1")
  }
  if(scale.perc <= 0 || scale.perc > 100){
    stop("scale.perc should be between 0 and 100")
  }
  
  # return image
  if(channel!=""){
    
    # get image
    img <- object@image[[channel]]
    if(as.raster){
      
      # return raster image format
      return(img)
      
    } else {
      
      # get image as array if image is stored as a DelayedArray
      if(inherits(img, "Image_Array")){
        # img <- as.array(img@seed)
        img <- as.array(img)
        img <- array(as.raw(img), dim = dim(img))
      }
      
      # read image
      img <- magick::image_read(img)
      
      # scale image if needed
      if(scale.perc < 100){
        img <- image_resize(img, geometry = magick::geometry_size_percent(scale.perc))
      }
      
      # return regular image
      return(img)
    }
  } else{
    warning("No image was found!")
    return(NULL)
  }
}

#' @rdname vrImages
#' @order 4
#' @importFrom magick image_read geometry_size_percent
#'
#' @export
setMethod("vrImages", "vrImage", vrImagesvrImage)

#' @rdname vrImages
#' @order 4
#' @importFrom magick image_read geometry_size_percent
#'
#' @export
setMethod("vrImages", "vrSpatial", vrImagesvrImage)

vrImagesReplacevrImage <- function(object, channel = NULL, value){
  
  if(channel %in% vrImageChannelNames(object)){
    warning("A channel with name '", channel, "' already exists in this vrImage object. \n Overwriting ...")
  }
  
  if(inherits(value, "bitmap")){
    object@image[[channel]] <- value
  } else if(inherits(value, "magick-image")){
    object@image[[channel]] <- magick::image_data(value)
  } else if(inherits(value, "Image_Array")){
    object@image[[channel]] <- value
  } else {
    stop("Please provide either a magick-image or bitmap class image object!")
  }
  
  # return
  object
}

#' @rdname vrImages
#'
#' @importFrom magick image_read
#' @order 6
#' @export
setMethod("vrImages<-", "vrImage", vrImagesReplacevrImage)

#' @rdname vrImages
#'
#' @importFrom magick image_read
#' @order 6
#' @export
setMethod("vrImages<-", "vrSpatial", vrImagesReplacevrImage)

vrMainImageVoltRon <- function(object, assay = NULL){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if assay = all, give a summary
  if(!is.null(assay)){
    if(assay == "all"){
      spatial_names <- unlist(lapply(rownames(SampleMetadata(object)), function(x) paste(vrMainSpatial(object[[x]]), collapse = ",")))
      spatial_names <- data.frame(Assay = assay_names, Spatial = spatial_names)
      return(spatial_names)
    }
  }
  
  # get assay types
  spatial_names <- unlist(lapply(assay_names, function(x) vrMainSpatial(object[[x]])))
  
  # return data
  spatial_data <- data.frame(Assay = assay_names, Spatial = spatial_names)
  
  # return
  return(spatial_data)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' If NULL, the default assay will be used, see \link{vrMainAssay}. If given as "all", then provides a summary of spatial systems across all assays.
#'
#' @rdname vrMainImage
#' @order 2
#' @export
setMethod("vrMainImage", "VoltRon", vrMainImageVoltRon)

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrMainSpatial
#' @order 2
#' @export
setMethod("vrMainSpatial", "VoltRon", vrMainImageVoltRon)

vrMainImageReplaceVoltRon <- function(object, assay = NULL, value){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get spatial metadata
  sample.metadata <- SampleMetadata(object)
  assayclass <- unique(sample.metadata[assay_names,"Assay"])
  
  # check for assay number
  if(length(assayclass) == 1){
    for(assy in assay_names)
      vrMainSpatial(object[[assy]], ignore = TRUE) <- value
  } else {
    stop("You can only set the main spatial system of a single assay")
  }
  
  return(object)
}

#' @param value the name of main image
#'
#' @rdname vrMainImage
#' @order 4
#' @export
setMethod("vrMainImage<-", "VoltRon", vrMainImageReplaceVoltRon)

#' @param value the name of main image
#'
#' @rdname vrMainSpatial
#' @order 4
#' @export
setMethod("vrMainSpatial<-", "VoltRon", vrMainImageReplaceVoltRon)

vrMainImagevrAssay <- function(object) return(object@main_image)

#' @rdname vrMainImage
#' @order 3
#' @export
setMethod("vrMainImage", "vrAssay", vrMainImagevrAssay)

#' @rdname vrMainImage
#' @order 3
#' @export
setMethod("vrMainImage", "vrAssayV2", vrMainImagevrAssay)

#' @rdname vrMainSpatial
#' @order 3
#' @export
setMethod("vrMainSpatial", "vrAssay", vrMainImagevrAssay)

#' @rdname vrMainSpatial
#' @order 3
#' @export
setMethod("vrMainSpatial", "vrAssayV2", vrMainImagevrAssay)

#' @noRd
.replaceMainSpatial <- function(object, ignore = FALSE, value){
  
  if(length(value) %in% c(1,2)){
    
    # get channel name if exists in the value
    if(length(value) == 2){
      channel <- value[2]
      value <- value[1]
    } else {
      channel <- NULL
    }
    
    # set main spatial/image
    if(value %in% vrSpatialNames(object)){
      object@main_image <- value
      
      # set channel
      if(!is.null(channel))
        vrMainChannel(object@image[[value]]) <- channel
      
    } else {
      if(ignore){
        warning("'",value,"' is not a spatial coordinate system in '", vrAssayNames(object),"'. Main system is still set to '", vrMainSpatial(object), "'")
      } else {
        stop("'",value,"' is not a spatial coordinate system in '", vrAssayNames(object),"'. Use ignore = TRUE for ignoring this message")
      }
    }
    
  } else {
    stop("The Main image is set by either: \n    vrMainSpatial(object) <- c('<spatial name>', '<channel name>')\n or vrMainSpatial(object) <- '<spatial name>'")
  }
  
  return(object)
}

#' @param ignore if TRUE, the non-existing spatial coordinate system will be ignored.
#' 
#' @rdname vrMainImage
#' @order 5
#' @export
setMethod("vrMainImage<-", "vrAssay", .replaceMainSpatial)

#' @param ignore if TRUE, the non-existing spatial coordinate system will be ignored.
#' 
#' @rdname vrMainImage
#' @order 5
#' @export
setMethod("vrMainImage<-", "vrAssayV2", .replaceMainSpatial)

#' @param ignore if TRUE, the non-existing spatial coordinate system will be ignored.
#' 
#' @rdname vrMainSpatial
#' @order 5
#' @export
setMethod("vrMainSpatial<-", "vrAssay", .replaceMainSpatial)

#' @param ignore if TRUE, the non-existing spatial coordinate system will be ignored.
#'
#' @rdname vrMainSpatial
#' @order 5
#' @export
setMethod("vrMainSpatial<-", "vrAssayV2", .replaceMainSpatial)

vrImageNamesVoltRon <- function(object, assay = NULL){
  
  # sample metadata
  sample.metadata <- SampleMetadata(object)
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if assay = all, give a summary
  if(!is.null(assay)){
    if(assay == "all"){
      spatial_names <- unlist(lapply(assay_names, function(x) paste(vrSpatialNames(object[[x]]), collapse = ",")))
      main_spatial_names <- unlist(lapply(assay_names, function(x) vrMainSpatial(object[[x]])))
      spatial_names <- data.frame(sample.metadata[assay_names,], Spatial = spatial_names, Main = main_spatial_names)
      return(spatial_names)
    }
  }
  
  # unique names
  spatial_names <- unique(unlist(lapply(assay_names, function(x) vrSpatialNames(object[[x]]))))
  
  return(spatial_names)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' If NULL, the default assay will be used, see \link{vrMainAssay}. If equals to "all", then provides a summary of spatial systems across all assays
#'
#' @rdname vrImageNames
#'
#' @export
setMethod("vrImageNames", "VoltRon", vrImageNamesVoltRon)

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}. If equals to "all", then provides a summary of spatial systems across all assays
#'
#' @rdname vrSpatialNames
#'
#' @export
setMethod("vrSpatialNames", "VoltRon", vrImageNamesVoltRon)

vrImageNamesvrAssay <- function(object) names(object@image)

#' @rdname vrImageNames
#'
#' @export
setMethod("vrImageNames", "vrAssay", vrImageNamesvrAssay)

#' @rdname vrImageNames
#'
#' @export
setMethod("vrImageNames", "vrAssayV2", vrImageNamesvrAssay)

#' @rdname vrSpatialNames
#'
#' @export
setMethod("vrSpatialNames", "vrAssay", vrImageNamesvrAssay)

#' @rdname vrSpatialNames
#'
#' @export
setMethod("vrSpatialNames", "vrAssayV2", vrImageNamesvrAssay)

####
## Channel Methods ####
####

vrMainChannelvrAssay <- function(object, name = NULL){
  if(is.null(name)){
    name <- vrMainSpatial(object)
  }
  return(vrMainChannel(object@image[[name]]))
}

#' @param name the name of the image
#'
#' @rdname vrMainChannel
#' @order 2
#' @export
setMethod("vrMainChannel", "vrAssay", vrMainChannelvrAssay)

#' @param name the name of the image
#'
#' @rdname vrMainChannel
#' @order 2
#' @export
setMethod("vrMainChannel", "vrAssayV2", vrMainChannelvrAssay)

vrMainChannelReplacevrAssay <- function(object, name = NULL, value){
  if(is.null(name)){
    name <- vrMainSpatial(object)
  }
  vrMainChannel(object@image[[name]]) <- value
  return(object)
}

#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @order 4
#' @export
setMethod("vrMainChannel<-", "vrAssay", vrMainChannelReplacevrAssay)

#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @order 4
#' @export
setMethod("vrMainChannel<-", "vrAssayV2", vrMainChannelReplacevrAssay)

#' @rdname vrMainChannel
#' @order 3
#' @export
setMethod("vrMainChannel", "vrImage", function(object){
  return(object@main_channel)
})

#' @rdname vrMainChannel
#' @order 3
#' @export
setMethod("vrMainChannel", "vrSpatial", function(object){
  return(object@main_channel)
})

vrMainChannelReplacevrImage <- function(object, value){
  
  if(value %in% vrImageChannelNames(object)){
    object@main_channel <- value
  } else {
    stop("'",value,"' is not a channel name")
  }
  return(object)
}

#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @method vrMainChannel<- vrImage
#' @order 5
#' @export
setMethod("vrMainChannel<-", "vrImage", vrMainChannelReplacevrImage)

#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @method vrMainChannel<- vrSpatial
#' @order 5
#' @export
setMethod("vrMainChannel<-", "vrSpatial", vrMainChannelReplacevrImage)

vrImageChannelNamesVoltRon <- function(object, assay = NULL){
  
  # get assay names
  if(is.null(assay)){
    assay_names <- vrAssayNames(object, assay = "all")
  } else {
    assay_names <- vrAssayNames(object, assay = assay)
  }
  
  # sample metadata
  sample.metadata <- SampleMetadata(object)
  
  # get image names
  spatial_names <- unlist(lapply(assay_names, function(x) vrMainSpatial(object[[x]])))
  
  # get channel names
  image_channels <- unlist(lapply(assay_names, function(x) paste(vrImageChannelNames(object[[x]]), collapse = ",")))
  
  # return data
  image_data <- data.frame(sample.metadata[assay_names,], Spatial = spatial_names, Channels = image_channels)
  
  # return
  return(image_data)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrImageChannelNames
#'
#' @export
setMethod("vrImageChannelNames", "VoltRon", vrImageChannelNamesVoltRon)

vrImageChannelNamesvrAssay <- function(object, name = NULL){
  
  if(is.null(name)){
    name <- vrMainSpatial(object)
  } else {
    if(!name %in% vrSpatialNames(object))
      stop(name, " is not among any image in this vrAssay object")
  }
  
  return(vrImageChannelNames(object@image[[name]]))
}

#' @param name the key of the image
#'
#' @rdname vrImageChannelNames
#'
#' @export
setMethod("vrImageChannelNames", "vrAssay", vrImageChannelNamesvrAssay)

#' @param name the key of the image
#'
#' @rdname vrImageChannelNames
#'
#' @export
setMethod("vrImageChannelNames", "vrAssayV2", vrImageChannelNamesvrAssay)

vrImageChannelNamesvrImage <- function(object){
  if(is.null(names(object@image))){
    return("No Channels or Images are found!")
  } else{
    return(names(object@image))
  }
}

#' @rdname vrImageChannelNames
#'
#' @export
setMethod("vrImageChannelNames", "vrImage", vrImageChannelNamesvrImage)

#' @rdname vrImageChannelNames
#'
#' @export
setMethod("vrImageChannelNames", "vrSpatial", vrImageChannelNamesvrImage)

####
## Managing Images ####
####

resizeImageVoltRon <- function(object, assay = NULL, name = NULL, reg = FALSE, size = NULL){
  
  assay_names <- vrAssayNames(object, assay = assay)
  
  for(assy in assay_names){
    object[[assy]] <- resizeImage(object[[assy]], name = name, reg = reg, size = size)
  }
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param name the name of the image
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainSpatial}) is requested
#' @param size the width of the resized image
#'
#' @rdname resizeImage
#'
#' @export
setMethod("resizeImage", "VoltRon", resizeImageVoltRon)

resizeImagevrAssay <- function(object, name = NULL, reg = FALSE, size = NULL){
  
  # get main image is main_image is null
  if(is.null(name)) {
    name <- object@main_image
  }
  
  # check registered image
  if(reg){
    if(!paste0(name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", name, "!")
    } else {
      name <- paste0(name, "_reg")
    }
  }
  
  # check main image
  if(!name %in% vrSpatialNames(object)){
    stop(name, " is not among any image in this vrAssay object")
  }
  
  object@image[[name]] <- resizeImage(object@image[[name]], size = size)
  
  # return
  return(object)
}

#' @rdname resizeImage
#'
#' @export
setMethod("resizeImage", "vrAssay", resizeImagevrAssay)

#' @rdname resizeImage
#'
#' @export
setMethod("resizeImage", "vrAssayV2", resizeImagevrAssay)

resizeImagevrImage <- function(object, size = NULL){
  
  # sizefactor
  sizefactor <- image_info(vrImages(object))$width
  
  # check size
  if(is.null(size))
    size = sizefactor
  if(!is.numeric(size))
    stop("width size should be numeric")
  if(!all.equal(size, as.integer(size)) & size > 0)
    stop("width size should be a positive integer")
  if(size < 100)
    stop("width size cannot be less than 100px")
  
  # resize coordinates
  vrCoordinates(object) <- (vrCoordinates(object)*size)/sizefactor
  
  # resize segments
  vrSegments(object) <- lapply(vrSegments(object), function(x) {
    x[,c("x", "y")] <- x[,c("x", "y")]*size/sizefactor
    if(any(colnames(x) %in% c("rx", "ry"))){
      x[,c("rx", "ry")] <- x[,c("rx", "ry")]*size/sizefactor
    }
    return(x)
  })
  
  # resize images
  size <- paste0(size,"x")
  image_names <- vrImageChannelNames(object)
  for(img in image_names){
    img_data <- object@image[[img]]
    if(inherits(img_data, "Image_Array")){
      stop("Currently modulateImage only works on in-memory images!")
    } else {
      img_data <- magick::image_read(img_data)
      img_data <- magick::image_resize(img_data, geometry = size)
      object@image[[img]] <- magick::image_data(img_data) 
    }
  }
  
  # return
  return(object)
}

#' @rdname resizeImage
#'
#' @importFrom magick image_info image_resize image_read image_data
#' @export
setMethod("resizeImage", "vrImage", resizeImagevrImage)

#' @rdname resizeImage
#'
#' @importFrom magick image_info image_resize image_read image_data
#' @export
setMethod("resizeImage", "vrSpatial", resizeImagevrImage)

modulateImageVoltRon <- function(object, assay = NULL, name = NULL, reg = FALSE, channel = NULL, 
                                  brightness = 100, saturation = 100, hue = 100, force = FALSE){
  
  assay_names <- vrAssayNames(object, assay = assay)
  
  for(assy in assay_names){
    object[[assy]] <- modulateImage(object[[assy]], name = name, reg = reg, channel = channel, brightness = brightness, 
                                    saturation = saturation, hue = hue, force = force)
  }
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param name the name of the image
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainSpatial}) is requested
#' @param channel the name of the channel associated with the image
#' @param brightness modulation of brightness as percentage of the current value (100 for no change)
#' @param saturation modulation of saturation as percentage of the current value (100 for no change)
#' @param hue modulation of hue is an absolute rotation of -180 degrees to +180 degrees from the current position corresponding to an argument range of 0 to 200 (100 for no change)
#' @param force if TRUE, all channels will be modulated given no specific channel name
#'
#' @rdname modulateImage
#'
#' @export
setMethod("modulateImage", "VoltRon", modulateImageVoltRon)

modulateImagevrAssay <- function(object,  name = NULL, reg = FALSE, channel = NULL, 
                                  brightness = 100, saturation = 100, hue = 100, force = FALSE){
  
  # check name
  if(is.null(name)) {
    name <- object@main_image
  }
  
  # get registered image
  if(reg){
    if(!paste0(name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", name, "!")
    } else {
      name <- paste0(name, "_reg")
    }
  }
  
  # check main image
  if(!name %in% vrSpatialNames(object)){
    stop(name, " is not among any image in this vrAssay object")
  }
  
  object@image[[name]] <- modulateImage(object@image[[name]], channel = channel, brightness = brightness, 
                                        saturation = saturation, hue = hue, force = force)
  
  # return
  return(object)
}

#' @rdname modulateImage
#'
#' @export
setMethod("modulateImage", "vrAssay", modulateImagevrAssay)

#' @rdname modulateImage
#'
#' @export
setMethod("modulateImage", "vrAssayV2", modulateImagevrAssay)

modulateImagevrImage <- function(object, channel = NULL, brightness = 100, saturation = 100, hue = 100, force = FALSE){
  
  # check main_channels
  if(is.null(channel) && (length(vrImageChannelNames(object)) > 1 && !force)){
    stop("No channel name was specified. \n It is not advised to modulate multiple channels in the same time. \n Please type force = TRUE to allow this behaviour!")
  }
  
  # get channel names
  if(is.null(channel)){
    channel <- vrImageChannelNames(object)
  }
  
  # modulate image
  for(img in channel){
    img_data <- object@image[[img]]
    if(inherits(img_data, "Image_Array")){
      stop("Currently modulateImage only works on in-memory images!")
    } else {
      img_data <- magick::image_read(img_data)
      # img_data <- getImage(object, name = img)
      img_data <- magick::image_modulate(img_data, brightness = brightness, saturation = saturation, hue = hue)
      object@image[[img]] <- magick::image_data(img_data) 
    }
  }
  
  # return
  return(object)
}

#' @rdname modulateImage
#'
#' @importFrom magick image_info image_modulate
#' @export
setMethod("modulateImage", "vrImage", modulateImagevrImage)

#' @rdname modulateImage
#'
#' @importFrom magick image_info image_modulate
#' @export
setMethod("modulateImage", "vrSpatial", modulateImagevrImage)

combineChannelsVoltRon <- function(object, assay = NULL, name = NULL, reg = FALSE, 
                                    channels = NULL, colors = NULL, channel_key = "combined"){
  
  assay_names <- vrAssayNames(object, assay = assay)
  
  for(assy in assay_names){
    object[[assy]] <- combineChannels(object[[assy]], name = name, reg = reg, 
                                      channels = channels, colors = colors, channel_key = channel_key)
  }
  return(object)
}


#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param name the name of the image
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainSpatial}) is requested
#' @param channels the name of the channel associated with the image
#' @param colors the colors associated with each channel
#' @param channel_key the name of the new channel name
#'
#' @rdname combineChannels
#'
#' @export
setMethod("combineChannels", "VoltRon", combineChannelsVoltRon)

combineChannelsvrAssay <- function(object,  name = NULL, reg = FALSE, channels = NULL, colors = NULL, channel_key = "combined"){
  
  # check name
  if(is.null(name)) {
    name <- object@main_image
  }
  
  # get registered image
  if(reg){
    if(!paste0(name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered images with name ", name, "!")
    } else {
      name <- paste0(name, "_reg")
    }
  }
  
  # check main image
  if(!name %in% vrSpatialNames(object)){
    stop(name, " is not among any image in this vrAssay object")
  }
  
  object@image[[name]] <- combineChannels(object@image[[name]], channels = channels, colors = colors, channel_key = channel_key)
  
  # return
  return(object)
}

#' @rdname combineChannels
#'
#' @export
setMethod("combineChannels", "vrAssay", combineChannelsvrAssay)

#' @rdname combineChannels
#'
#' @export
setMethod("combineChannels", "vrAssayV2", combineChannelsvrAssay)

combineChannelsvrImage <- function(object, channels = NULL, colors = NULL, channel_key = "combined"){
  
  # check channel names
  if(is.null(channels)){
    stop("No channel names were given")
  } else {
    if(any(!channels %in% vrImageChannelNames(object))){
      warning("Some channel names do not match with the existing channels.")
    }
  }
  
  # check colors
  if(is.null(colors)){
    stop("No colors were given")
  }
  if(length(colors) != length(channels)){
    stop("The length of colors do not match with the length of channels.")
  }
  
  # configure channel and color names
  colors <- colors[channels %in% vrImageChannelNames(object)]
  channels <- channels[channels %in% vrImageChannelNames(object)]
  names(colors) <- channels
  
  # get images and colorize
  channel_list <- list()
  composite_image <- NULL
  for(img in channels){
    channel_img <- vrImages(object, channel = img)
    color_rgb <-  grDevices::col2rgb(colors[img])[,1]
    imagedata <- as.numeric(magick::image_data(channel_img, channels = "rgb"))
    imagedata[,,1] <- imagedata[,,1] * (color_rgb[1]/255)
    imagedata[,,2] <- imagedata[,,2] * (color_rgb[2]/255)
    imagedata[,,3] <- imagedata[,,3] * (color_rgb[3]/255)
    channel_img <- magick::image_read(imagedata)
    if(is.null(composite_image)){
      composite_image <- channel_img
    } else{
      composite_image <- magick::image_composite(channel_img, composite_image, operator = "Plus")
    }
  }
  
  # combine channels
  vrImages(object, channel = channel_key) <- composite_image
  
  # return
  return(object)
}

#' @rdname combineChannels
#'
#' @importFrom magick image_read image_data image_composite
#' @importFrom grDevices col2rgb
#'
#' @export
setMethod("combineChannels", "vrImage", combineChannelsvrImage)

#' @rdname combineChannels
#'
#' @export
setMethod("combineChannels", "vrSpatial", combineChannelsvrImage)

####
# Other Methods ####
####

#' @rdname vrSpatialPoints
#' @order 4
#'
#' @export
setMethod("vrSpatialPoints", "vrImage", function(object) {
  return(rownames(vrCoordinates(object)))
})

#' @rdname vrSpatialPoints
#' @order 4
#'
#' @export
setMethod("vrSpatialPoints", "vrSpatial", function(object) {
  return(rownames(vrCoordinates(object)))
})

vrSpatialPointsReplacevrImage <- function(object, value) {
  
  # coordinates
  if(length(rownames(object@coords)) != length(value)){
    stop("The number of spatial points is not matching with the input")
  } else {
    rownames(object@coords)  <- value
  }
  
  # segments
  if(length(object@segments) > 0){
    if(length(names(object@segments)) != length(value)){
      stop("The number of spatial points is not matching with the input")
    } else {
      names(object@segments) <- value
    }
  }
  
  # return
  return(object)
}

#' @param value new spatial points
#'
#' @rdname vrSpatialPoints
#' @order 9
#' @export
setMethod("vrSpatialPoints<-", "vrImage", vrSpatialPointsReplacevrImage)

#' @param value new spatial points
#'
#' @rdname vrSpatialPoints
#' @order 9
#' @export
setMethod("vrSpatialPoints<-", "vrSpatial", vrSpatialPointsReplacevrImage)

#' @rdname vrCoordinates
#' @order 3
#' @export
setMethod("vrCoordinates", "vrImage", function(object) {
    return(object@coords)
})

#' @rdname vrCoordinates
#' @order 3
#' @export
setMethod("vrCoordinates", "vrSpatial", function(object) {
  return(object@coords)
})

vrCoordinatesRepkacevrImage <- function(object, value) {
  
  # get coordinates
  coords <- vrCoordinates(object)
  
  # stop if the rownames are not matching
  if(any(sapply(rownames(value),is.null)))
    stop("Provided coordinates data does not have cell/spot/ROI names")
  
  if(!all(rownames(value) %in% rownames(coords)))
    stop("Cant overwrite coordinates, non-existing cells/spots/ROIs!")
  
  # stop if the colnames there are more than two columns
  if(ncol(value) == 2){
    value <- cbind(value, 0)
    colnames(value) <- c("x", "y", "z")
  } else if(ncol(value) == 3){
    colnames(value) <- c("x", "y", "z")
  } else {
    stop("Please make sure that the coordinates matrix have only two or three columns: for x and y coordinates")
  }
  
  methods::slot(object = object, name = 'coords') <- value
  return(object)
}
    
#' @rdname vrCoordinates
#' @order 6
#' @importFrom methods slot
#'
#' @export
setMethod("vrCoordinates<-", "vrImage", vrCoordinatesRepkacevrImage)

#' @rdname vrCoordinates
#' @order 6
#' @importFrom methods slot
#'
#' @export
setMethod("vrCoordinates<-", "vrSpatial", vrCoordinatesRepkacevrImage)

#' @rdname vrSegments
#' @order 4
#' @export
setMethod("vrSegments", "vrImage", function(object) {
  return(object@segments)
})

#' @rdname vrSegments
#' @order 4
#' @export
setMethod("vrSegments", "vrSpatial", function(object) {
  return(object@segments)
})

vrSegmentsReplacevrImage <- function(object, value) {
  
  # get coordinates
  segts <- vrSegments(object)
  
  # stop if the names are not matching
  if(any(sapply(names(value),is.null)))
    stop("Provided coordinates data does not have cell/spot/ROI names")
  
  if(!all(names(value) %in% names(segts)))
    stop("Cant overwrite coordinates, non-existing cells/spots/ROIs!")
  
  methods::slot(object = object, name = 'segments') <- value
  return(object)
}

#' @rdname vrSegments
#' @order 7
#' @importFrom methods slot
#' @export
setMethod("vrSegments<-", "vrImage", vrSegmentsReplacevrImage)

#' @rdname vrSegments
#' @order 7
#' @importFrom methods slot
#' @export
setMethod("vrSegments<-", "vrSpatial", vrSegmentsReplacevrImage)

####
# Demultiplex Images ####
####

#' demuxVoltRon
#'
#' Subsetting/demultiplexing of the VoltRon Object using interactive shiny app
#'
#' @param object a VoltRon object
#' @param max.pixel.size the initial width of the object image
#' @param use.points.only use spatial points instead of the reference image
#' @param shiny.options a list of shiny options (launch.browser, host, port etc.) passed \code{options} arguement of \link{shinyApp}. For more information, see \link{runApp}
#'
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom magick image_scale image_info image_ggplot
#' @importFrom ggplot2 geom_rect
#' @importFrom dplyr filter add_row tibble
#' @importFrom ggrepel geom_label_repel
#'
demuxVoltRon <- function(object, max.pixel.size = 1200, use.points.only = FALSE, shiny.options = list(launch.browser = getOption("shiny.launch.browser", interactive())))
{
  # check if there are only one assay in the object
  sample.metadata <- SampleMetadata(object)
  
  if(length(unique(sample.metadata$Layer)) > 1)
    stop("You can only subset a single VoltRon layer at a time")
  
  # get image
  images <- vrImages(object[[vrAssayNames(object)]], as.raster = TRUE)
  if(!inherits(images, "Image_Array")){
    images <- magick::image_read(images)
  }
  
  # scale 
  imageinfo <- getImageInfo(images)
  scale_factor <- 1
  if(imageinfo$width > max.pixel.size){
    scale_factor <- imageinfo$width/max.pixel.size
  }
  if(use.points.only){
    object_small <- resizeImage(object, size = max.pixel.size)
    image_info_small <- magick::image_info(vrImages(object_small))
    coords <- as.data.frame(vrCoordinates(object_small, reg = FALSE))
    pl <- ggplot() + geom_point(aes_string(x = "x", y = "y"), coords, size = 1.5, color = "black") +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
            axis.line=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
            legend.margin = margin(0,0,0,0), plot.margin = unit( c(0,0,0,0),"in")) +
      coord_fixed()
  } else {
    pl <- plotImage(images, max.pixel.size = max.pixel.size)
  }

  # get the ui and server
  
  # UI ####
  ui <- fluidPage(
    
    # use javascript extensions for Shiny
    shinyjs::useShinyjs(),
    
    # sidebar
    sidebarLayout(position = "left",
                  
                  # Side bar
                  sidebarPanel(
                    tags$style(make_css(list('.well', 'margin', '7%'))),
                    
                    # Interface
                    fluidRow(
                      column(12,h4("Interactive Subsetting"))
                    ),
                    
                    # Buttons
                    fluidRow(
                      column(12,shiny::actionButton("resetpoints", "Remove Box")),
                      br(),
                      column(12,shiny::actionButton("addbox", "Add Box")),
                      br()
                    ),
                    
                    # instructions
                    h4("How to use"),
                    p(style="font-size: 12px;", strong("Single-L-hold-drag:"), "Select area"),
                    p(style="font-size: 12px;", strong("Add Box"), " to set a new subset"),
                    p(style="font-size: 12px;", strong("Remove Box"), " to reset the box"),
                    br(),
                    
                    # Subsets
                    fluidRow(
                      column(12,h4("Selected Subsets")),
                      uiOutput("textbox_ui"),
                      br()
                    ),
                    
                    # Subsets
                    fluidRow(
                      column(12,shiny::actionButton("done", "Done"))
                    ),
                    br(),
                    
                    # panel options
                    width = 3,
                  ),
                  
                  mainPanel(
                    
                    # main image
                    br(),
                    br(),
                    fluidRow(
                      plotOutput("cropped_image",
                                 height = "1000px",
                                 brush = brushOpts(
                                   id = "plot_brush", fill = "green",
                                   resetOnNew = TRUE
                                 )),
                    ),
                    
                    # panel options
                    width = 9
                  )
    )
  )
  
  # Server ####
  server <- function(input, output, session) {
    
    ## Importing images and variables ####
    
    # selected corner list
    selected_corners_list_image <- reactiveVal(dplyr::tibble(box = character()))
    selected_corners_list <- reactiveVal(list())
    
    ## Region Annotators ####
    
    ### Text Box Management ####
    
    # Reactive value to store the number of textboxes
    textboxes <- reactiveVal(numeric(0))
    
    # Initialize textbox values if n > 0, get already existing segments
    textbox_values <- reactiveValues()
    
    # Dynamically generate UI for textboxes and remove buttons
    output$textbox_ui <- renderUI({
      lapply(textboxes(), function(i) {
        column(12,
               textInputwithButton(textinputId = paste0("sample", i), label = paste0("Subset ", i),
                                   buttoninputId = paste0("remove", i), value = isolate(textbox_values[[paste0("sample", i)]]), 
                                   onclick = sprintf('Shiny.setInputValue("remove", %d)', i))
        )
      })
    })
    
    # Observe changes in each textbox to update their values
    observe({
      lapply(textboxes(), function(i) {
        observeEvent(input[[paste0("sample", i)]], {
          textbox_values[[paste0("sample", i)]] <- isolate(input[[paste0("sample", i)]])
        }, ignoreNULL = FALSE)
      })
    })
    
    ### Reset box ####
    observeEvent(input$resetpoints, {
      session$resetBrush("plot_brush")
    })
    
    ### Remove box ####
    
    # Observe event to remove textbox when the button is clicked
    observeEvent(input$remove, {
      
      # remove one point
      selected_corners_list(selected_corners_list()[!(textboxes() == as.numeric(isolate(input$remove)))])
      
      # Update the reactive value to remove the textbox
      textboxes(setdiff(textboxes(), as.numeric(isolate(input$remove))))
      
      # Remove the value from textbox_values
      textbox_values[[paste0("sample", as.numeric(input$remove))]] <- NULL
      
    }, ignoreInit = TRUE)
    
    ### Add box ####
    observeEvent(input$addbox, {
      
      # get corners
      brush <- input$plot_brush
      
      # add a box if brush is active
      if(!is.null(brush)){
        
        # corners 
        corners <- data.frame(x = c(brush$xmin, brush$xmax), 
                              y = c(brush$ymax, brush$ymin))
        
        # record corners
        selected_corners_list(c(selected_corners_list(), list(corners)))
        
        # adjust corners
        corners <- corners*scale_factor
        corners <- FromBoxToCrop(corners, imageinfo)

        # add to box list
        selected_corners_list_image() %>%
          dplyr::add_row(box = corners) %>%
          selected_corners_list_image()
        
        # reset box
        session$resetBrush("plot_brush")
        
        # add buttons
        new_id <- if (length(textboxes()) == 0) 1 else max(textboxes()) + 1
        textboxes(c(textboxes(), new_id))
        textbox_values[[paste0("sample", new_id)]] <- ""
      }
    })
    
    ## Main observable ####
    observe({
      
      # output image
      output[["cropped_image"]] <- renderPlot({
        
        # visualize already selected boxes
        if(length(selected_corners_list()) > 0){
          for (i in 1:length(selected_corners_list())){
            corners <- apply(as.matrix(selected_corners_list()[[i]]),2,as.numeric)
            if(nrow(corners) > 1){
              corners <- as.data.frame(rbind(cbind(corners[1,1], corners[1:2,2]), cbind(corners[2,1], corners[2:1,2])))
              colnames(corners) <- c("x", "y")
              pl <- pl + ggplot2::geom_polygon(aes(x = x, y = y), data = corners, alpha = 0.3, fill = "green", color = "black")
              
            }
          }
        }
        
        # put labels of the already selected polygons
        if(length(selected_corners_list()) > 0){
          for (i in 1:length(selected_corners_list())){
            corners <- selected_corners_list()[[i]]
            corners <- as.data.frame(rbind(cbind(corners[1,1], corners[1:2,2]), cbind(corners[2,1], corners[2:1,2])))
            corners <- data.frame(x = mean(corners[,1]), y = max(corners[,2]), sample = paste("Subset ", isolate(textboxes()[i])))
            pl <- pl +
              ggrepel::geom_label_repel(mapping = aes(x = x, y = y, label = sample), data = corners,
                                        size = 5, direction = "y", nudge_y = 6, box.padding = 0, label.padding = 1, seed = 1, color = "red")
          }
        }
        
        # return graph
        pl
      })
    })
    
    ## Done ####
    
    # show "Done" if a region is selected already
    observe({
      if(nrow(selected_corners_list_image()) > 0){
        shinyjs::show(id = "done")
      } else {
        shinyjs::hide(id = "done")
      }
    })
    
    # observe for done and return the list of objects
    observeEvent(input$done, {
      if(nrow(selected_corners_list_image()) > 0){
        subsets <- list()
        box_list <- selected_corners_list_image()
        
        # collect labels
        sample_names <- sapply(1:length(box_list$box), function(i) input[[paste0("sample",i)]])

        # check if sample names are present
        if(any(sample_names == "")) {
          showNotification("Some subsets have blank (empty!) sample names.")
        } else{
          for(i in 1:length(box_list$box)){
            temp <- subset(object, image = box_list$box[i])
            temp$Sample <- sample_names[i]
            subsets[[sample_names[i]]] <- temp
          }
          stopApp(list(subsets = subsets, subset_info_list = box_list))
        }
        
      } else {
        showNotification("You have not selected a subset yet!")
      }
    })
  }
  
  # configure options
  shiny.options <- configure_shiny_options(shiny.options)
  
  # run app
  shiny::runApp(
    shiny::shinyApp(ui, server, options = list(host = shiny.options[["host"]], port = shiny.options[["port"]], launch.browser = shiny.options[["launch.browser"]]),
                    onStart = function() {
                      onStop(function() {
                      })
                    })
  )
}


#' FromBoxToCrop
#'
#' get magick crop information from a dataframe of box corners
#'
#' @param corners topleft and bottomright coordinates of bounding box
#' @param imageinfo info of the image
#' 
#' @noRd
FromBoxToCrop <- function(corners, imageinfo){
  
  corners <- apply(corners,2,ceiling)
  
  # fix for limits
  corners[1,1] <- ifelse(corners[1,1] < 0, 0, corners[1,1])
  corners[1,1] <- ifelse(corners[1,1] > imageinfo$width, imageinfo$width, corners[1,1])
  corners[2,1] <- ifelse(corners[2,1] < 0, 0, corners[2,1])
  corners[2,1] <- ifelse(corners[2,1] > imageinfo$width, imageinfo$width, corners[2,1])
  corners[1,2] <- ifelse(corners[1,2] < 0, 0, corners[1,2])
  corners[1,2] <- ifelse(corners[1,2] > imageinfo$height, imageinfo$height, corners[1,2])
  corners[2,2] <- ifelse(corners[2,2] < 0, 0, corners[2,2])
  corners[2,2] <- ifelse(corners[2,2] > imageinfo$height, imageinfo$height, corners[2,2])

  # get crop info
  corners <- paste0(abs(corners[2,1]-corners[1,1]), "x",
                    abs(corners[2,2]-corners[1,2]), "+",
                    min(corners[,1]), "+", imageinfo$height - max(corners[,2]))

  # corners 
  return(corners)
}

#' FromSegmentToCrop
#'
#' get magick crop information from coordinates of a segment
#'
#' @param segment coordinates of a segment
#' @param imageinfo info of the image
#' 
#' @export
FromSegmentToCrop <- function(segment, imageinfo){
  
  # make box from segment coordinates
  corners <- matrix(c(0,0,0,0), nrow = 2, ncol = 2)
  corners[1,1] <- min(segment[,1])
  corners[2,1] <- max(segment[,1])
  corners[1,2] <- max(segment[,2])
  corners[2,2] <- min(segment[,2])

  # get crop from box
  corners <- FromBoxToCrop(corners, imageinfo)

  # corners 
  return(corners)
}

