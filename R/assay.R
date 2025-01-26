#' @include zzz.R
#' @importClassesFrom Matrix dgCMatrix dgRMatrix dgeMatrix
#' @importClassesFrom S4Arrays Array
NULL

####
# Objects and Classes ####
####

## vrAssay and vrAssayV2####

## UpdateAssay ####

updateAssayvrAssay <- function(object){
  
  # data matrix and feature data
  data_list <- list(main = object@rawdata, main_norm = object@normdata)
  featuredata_list <- list(main = object@featuredata)
  
  # create assay v2
  methods::new("vrAssayV2",
               data = data_list,
               featuredata = featuredata_list,
               embeddings = object@embeddings,
               image = object@image,
               params = object@params,
               type = object@type,
               name = object@name,
               main_image = object@main_image,
               main_featureset = "main")
}

#' @param object a vrAssay object to be converted to vrAssayV2
#' @rdname updateAssay
#' @method updateAssay vrAssay
#' @importFrom methods new
setMethod("updateAssay", "vrAssay", updateAssayvrAssay)

updateAssayvrAssayV2 <- function(object){
  message("The assay is of version 2, nothing to change!")
  return(object)
}

#' @param object a vrAssayV2 object to be converted to vrAssayV2
#' @rdname updateAssay
#' @method updateAssay vrAssayV2
setMethod("updateAssay", "vrAssayV2", updateAssayvrAssayV2)

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
#' @param main_featureset the name of the main_featureset
#' @param assay_version the assay version
#' @param ... additional arguements passed to \link{formImage}
#'
#' @importFrom methods new
#'
#' @export
#'
formAssay <- function(data = NULL, 
                      coords, 
                      segments = list(), 
                      image = NULL, 
                      params = list(), 
                      type = "ROI", 
                      name = "Assay1", 
                      main_image = "image_1", 
                      main_featureset = NULL, 
                      assay_version = "v2", 
                      ...){

  # get data
  if(is.null(data)){
    data <- matrix(nrow = 0, ncol = nrow(coords))
    colnames(data) <- rownames(coords)
  }

  # get image object
  image <- formImage(coords = coords, segments = segments, image = image, ...)
  image <- list(image)
  names(image) <- main_image

  # check feature
  if(is.null(main_featureset))
    main_featureset <- "main"
  
  # make vrAssay object
  data_list <- list(main = data, main_norm = data)
  names(data_list) <- c(main_featureset, paste0(main_featureset, "_norm"))
  if(assay_version == "v2"){
    object <-   methods::new("vrAssayV2", 
                             data = data_list,
                             image = image, params = params, type = type, name = name, 
                             main_image = main_image, main_featureset = main_featureset)
  } else {
    object <-   methods::new("vrAssay", 
                             rawdata = data, normdata = data,
                             image = image, params = params, type = type, name = name, 
                             main_image = main_image)
  }
  return(object)
}

### Subset vrAssay objects ####

subsetvrAssay <- function(x, subset, spatialpoints = NULL, features = NULL, image = NULL) {
  
  # start 
  object <- x
  
  if (!missing(x = subset)) {
    subset <- rlang::enquo(arg = subset)
  }
  
  # subseting on samples, layers and assays
  if(!is.null(features)){
    
    # select features
    nonmatching_features <- setdiff(features, vrFeatures(object))
    features <- intersect(vrFeatures(object), features)
    
    if(length(features) > 0){
      # object@rawdata <- object@rawdata[rownames(object@rawdata) %in% features,, drop = FALSE]
      # object@normdata <- object@normdata[rownames(object@normdata) %in% features,, drop = FALSE]
      object <- subsetData(object, features = features)
      object <- subsetData(object, features = features)
      
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
      # object@rawdata  <- object@rawdata[,spatialpoints, drop = FALSE]
      # object@normdata  <- object@normdata[,spatialpoints, drop = FALSE]
      object <- subsetData(object, spatialpoints = spatialpoints)
      object <- subsetData(object, spatialpoints = spatialpoints)
      
      # embeddings
      for(embed in vrEmbeddingNames(object)){
        embedding <- vrEmbeddings(object, type = embed)
        vrEmbeddings(object, type = embed) <- embedding[spatialpoints[spatialpoints %in% rownames(embedding)],, drop = FALSE]
      }
      
      # image
      # for(img in vrImageNames(object))
      for(img in vrSpatialNames(object))
        object@image[[img]] <- subsetvrImage(object@image[[img]], spatialpoints = spatialpoints)
        # object@image[[img]] <- subset.vrImage(object@image[[img]], spatialpoints = spatialpoints)
      
    } else if(!is.null(image)) {
      
      # images
      img <- vrMainSpatial(object)
      object@image <- object@image[img]
      object@image[[img]] <- subsetvrImage(object@image[[img]], image = image)
      # object@image[[img]] <- subset.vrImage(object@image[[img]], image = image)
      spatialpoints <- rownames(vrCoordinates(object@image[[img]]))
      
      # data
      # object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% spatialpoints, drop = FALSE]
      # object@normdata  <- object@normdata[,colnames(object@normdata) %in% spatialpoints, drop = FALSE]
      object <- subsetData(object, spatialpoints = spatialpoints)
      object <- subsetData(object, spatialpoints = spatialpoints)
      
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

#' Subsetting vrAssay objects
#'
#' Given a vrAssay object, subset the object given one of the attributes
#'
#' @param x a vrAssay object
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
setMethod("subset", "vrAssay", subsetvrAssay)

#' Subsetting vrAssayV2 objects
#'
#' Given a vrAssayV2 object, subset the object given one of the attributes
#'
#' @param x a vrAssayV2 object
#' @param subset Logical statement for subsetting
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrAssayV2
#' @order 4
#'
#' @export
setMethod("subset", "vrAssayV2", subsetvrAssay)

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
  if(inherits(coords, "IterableMatrix")){
    # BPCells only accepts e1 > e2 ## S4 method for signature 'IterableMatrix,numeric'
    inside <- (!!as.vector(as(coords[,1] > xlim[1], "dgCMatrix")) & 
                 !!!as.vector(as(coords[,1] > xlim[2], "dgCMatrix"))) & 
      (!!as.vector(as(coords[,2] > ylim[1], "dgCMatrix")) & 
         !!!as.vector(as(coords[,2] > ylim[2], "dgCMatrix"))) 
  } else {
    inside <- (coords[,1] > xlim[1] & coords[,1] < xlim[2]) & (coords[,2] > ylim[1] & coords[,2] < ylim[2])
  }
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
#' @importFrom dplyr bind_rows
subsetSegments <- function(segments, image, crop_info){

  # get segments
  segment_names <- names(segments)
  segments <- do.call(dplyr::bind_rows, segments)
  rownames(segments) <- 1:nrow(segments)
  segments <- data.frame(segments, row_id = rownames(segments))
  
  # subset
  cropped_segments <- subsetCoordinates(segments[,c("x","y")], image, crop_info)
  if(any(colnames(segments) %in% c("rx", "ry"))){
    cropped_segments_extra <- segments[rownames(cropped_segments), c("rx", "ry")]
    cropped_segments <- cbind(cropped_segments, cropped_segments_extra)
  }
  cropped_segments <- data.frame(cropped_segments, id = segments[rownames(cropped_segments),1], row_id = rownames(cropped_segments))
  cropped_segments <- cropped_segments %>% right_join(segments[,c(colnames(segments)[1], "row_id")], by = c("row_id" = "row_id"))
  if(any(colnames(segments) %in% c("rx", "ry"))){
    cropped_segments <- cropped_segments[,c(colnames(cropped_segments)[which(grepl(colnames(segments)[1], colnames(cropped_segments)))[1]], "x", "y", "rx", "ry")]
    colnames(cropped_segments) <- c("id", "x", "y", "rx", "ry")
    
  } else {
    cropped_segments <- cropped_segments[,c(colnames(cropped_segments)[which(grepl(colnames(segments)[1], colnames(cropped_segments)))[1]], "x", "y")]
    colnames(cropped_segments) <- c("id", "x", "y")
  }
  # split back to segments
  segments <- split(cropped_segments, cropped_segments[,1])
  segments <- lapply(segments, function(df){
    df[,colSums(is.na(df))<nrow(df), drop = FALSE]
  })
  names(segments) <- segment_names
  
  # return
  return(segments)
}

#' subsetData
#'
#' subsetting data matrices given spatialpoints, features etc.
#'
#' @param object a vrAssay object
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#'
#' @noRd
subsetData <- function(object, spatialpoints = NULL, features = NULL){
  
  # features
  if(!is.null(features)){
    
    if(inherits(object, "vrAssay")){
      if(nrow(object@rawdata) > 0){
        object@rawdata <- object@rawdata[rownames(object@rawdata) %in% features,, drop = FALSE]
        object@normdata <- object@normdata[rownames(object@normdata) %in% features,, drop = FALSE]
      }
    } else {
      main <- vrMainFeatureType(object)
      if(nrow(object@data[[main]]) > 0){
        object@data[[main]] <- object@data[[main]][rownames(object@data[[main]]) %in% features,, drop = FALSE]
        object@data[[paste0(main, "_norm")]] <- object@data[[paste0(main, "_norm")]][rownames(object@data[[paste0(main, "_norm")]]) %in% features,, drop = FALSE]
      }
    }
  }
  
  # spatialpoints
  if(!is.null(spatialpoints)){
    
    if(inherits(object, "vrAssay")){
      # if(nrow(object@rawdata) > 0){
      if(ncol(object@rawdata) > 0){
        object@rawdata  <- object@rawdata[,colnames(object@rawdata) %in% spatialpoints, drop = FALSE]
        object@normdata  <- object@normdata[,colnames(object@normdata) %in% spatialpoints, drop = FALSE]
      }
    } else {
      for(nm in vrFeatureTypeNames(object)){
        # if(nrow(object@data[[nm]]) > 0){
        if(ncol(object@data[[nm]]) > 0){
          object@data[[nm]] <- object@data[[nm]][,colnames(object@data[[nm]]) %in% spatialpoints, drop = FALSE]
          object@data[[paste0(nm, "_norm")]] <- object@data[[paste0(nm, "_norm")]][,colnames(object@data[[paste0(nm, "_norm")]]) %in% spatialpoints, drop = FALSE]
        }
      }
    }
  }
  
  # return
  return(object)
}

#' getData
#'
#' get data matrix
#'
#' @param object a vrAssay object
#'
#' @noRd
getData <- function(object){
  
  if(inherits(object, "vrAssay")){
    data <- object@rawdata
  } else {
    data <- object@data[[vrMainFeatureType(object)]]
  }
  
  return(data)
}

#' updateData
#'
#' update data matrix
#'
#' @param object a vrAssay object
#' @param value the new column names
#'
#' @noRd
updateData <- function(object, value){
  
  if(inherits(object, "vrAssay")){
    if(ncol(object@rawdata) > 0){
      colnames(object@rawdata) <- value
      colnames(object@normdata) <- value 
    }
  } else {
    for(nm in vrFeatureTypeNames(object)){
      if(ncol(object@data[[nm]] > 0)){
        colnames(object@data[[nm]]) <- value
        colnames(object@data[[paste0(nm, "_norm")]]) <- value
      }
    }
  }
  
  return(object)
}

### Feature Methods ####

vrMainFeatureTypevrAssayV2 <- function(object){
  if(inherits(object, "vrAssayV2")){
    return(object@main_featureset)
  } else {
    return(NULL)
  }
}

#' @rdname vrMainFeatureType
#' @order 3
#' @export
setMethod("vrMainFeatureType", "vrAssayV2", vrMainFeatureTypevrAssayV2)

#' @rdname vrMainFeatureType
#' @order 3
#' @export
setMethod("vrMainFeatureType", "vrAssay", vrMainFeatureTypevrAssayV2)

vrMainFeatureTypeReplacevrAssayV2 <- function(object, ignore = FALSE, value){
  if(value %in% names(object@data)){
    object@main_featureset <- value
  } else {
    if(ignore){
      warning("The feature type '", value, "' is not found in '", vrAssayNames(object),"'. Main feature type is still set to '", vrMainFeatureType(object), "'")
    } else {
      stop("The feature type '", value, "' is not found in '", vrAssayNames(object),"'. Use ignore = TRUE for ignoring this message")
    }
  }
  
  return(object)
}

#' @param ignore ignore if some assays dont have the feature set name
#' 
#' @rdname vrMainFeatureType
#' @order 5
#' @export
setMethod("vrMainFeatureType<-", "vrAssayV2", vrMainFeatureTypeReplacevrAssayV2)

#' @param ignore ignore if some assays dont have the feature set name
#' 
#' @rdname vrMainFeatureType
#' @order 5
#' @export
setMethod("vrMainFeatureType<-", "vrAssay", function(object, ignore = FALSE, value){
  stop("vrAssay V1 objects do not have multiple feature types!")
})

vrFeatureTypeNamesvrAssayV2 <- function(object){
  names_data <- names(object@data)
  return(names_data[!grepl("_norm$", names_data)])
}

#' @rdname vrFeatureTypeNames
#'
#' @export
setMethod("vrFeatureTypeNames", "vrAssayV2", vrFeatureTypeNamesvrAssayV2)

#' @rdname vrFeatureTypeNames
#'
#' @export
setMethod("vrFeatureTypeNames", "vrAssay", function(object){
  stop("vrAssay V1 objects do not have multiple feature types!")
})

addFeaturevrAssayV2 <- function(object, data, feature_name){
  
  # get feature name
  featuresets <- vrFeatureTypeNames(object)
  if(feature_name %in% featuresets){
    stop(paste0("Feature type '", feature_name, "' already exists in the assay."))
  }
  
  # check spatial point names in the object
  colnames_data <- colnames(data)
  colnames_data <- stringr::str_remove(colnames_data, pattern = "_Assay[0-9]+$")
  colnames(data) <- paste0(colnames_data, "_", vrAssayNames(object))
  
  # check spatial points
  spatialpoints <- vrSpatialPoints(object)
  if(length(setdiff(colnames(data), vrSpatialPoints(object))) > 0){
    stop("The number of spatial points is not matching with number of points in the input data")
  } 

  # add new features
  feature_list_name <- names(object@data)
  feature_list_name <- c(feature_list_name, feature_name, paste0(feature_name, "_norm"))
  object@data <- c(object@data, list(data,data))
  names(object@data) <- feature_list_name
  
  # return
  return(object)
}

#' @rdname addFeature
#' @method addFeature vrAssayV2
#' 
#' @importFrom stringr str_remove
#' 
#' @export
setMethod("addFeature", "vrAssayV2", addFeaturevrAssayV2)

### Other Methods ####

#' @rdname vrSpatialPoints
#' @order 4
#' 
#' @export
setMethod("vrSpatialPoints", "vrAssay", function(object) {
  return(rownames(vrCoordinates(object)))
})

#' @rdname vrSpatialPoints
#' @order 4
#' 
#' @export
setMethod("vrSpatialPoints", "vrAssayV2", function(object) {
  return(rownames(vrCoordinates(object)))
})

vrSpatialPointsReplacevrAssayV2 <- function(object, value) {
  
  # spatial points 
  spatialpoints <- vrSpatialPoints(object)
  
  # data
  if(length(vrSpatialPoints(object)) != length(value)){
    stop("The number of spatial points is not matching with the input")
  } else {
    if(ncol(getData(object)) > 0){
      object <- updateData(object, value)
    }
  }
  
  # images
  for(img in vrSpatialNames(object)){
    vrSpatialPoints(object@image[[img]]) <- value
  }
  
  # embeddings
  embeddings <- object@embeddings
  embed_names <- names(embeddings)
  if(length(embed_names) > 0){
    for(type in embed_names){
      if(nrow(embeddings[[type]]) > 0){
        rownames(embeddings[[type]]) <- value[match(rownames(embeddings[[type]]), spatialpoints)]
        object@embeddings[[type]] <- embeddings[[type]]
      }
    }
  }
  
  # return
  return(object)
}

#' @rdname vrSpatialPoints
#' @order 8
#' @export
setMethod("vrSpatialPoints<-", "vrAssay", vrSpatialPointsReplacevrAssayV2)

#' @rdname vrSpatialPoints
#' @order 8
#' @export
setMethod("vrSpatialPoints<-", "vrAssayV2", vrSpatialPointsReplacevrAssayV2)

vrFeaturesvrAssay <- function(object) {
  return(rownames(getData(object)))
}

#' @rdname vrFeatures
#' @order 3
#' @export
setMethod("vrFeatures", signature = "vrAssay", definition = vrFeaturesvrAssay)
 
#' @rdname vrFeatures
#' @method vrFeatures vrAssayV2
#' @order 3
#' @export
setMethod("vrFeatures", "vrAssayV2", vrFeaturesvrAssay)

vrFeatureDatavrAssay <- function(object) {
  return(object@featuredata)
}

#' @rdname vrFeatureData
#' @order 3
#' @export
setMethod("vrFeatureData", "vrAssay", vrFeatureDatavrAssay)

vrFeatureDatavrAssayV2 <- function(object, feat_type = NULL){
  if(is.null(feat_type))
    feat_type <- vrMainFeatureType(object)
  return(object@featuredata[[feat_type]])
}

#' @param feat_type the feature set type
#'
#' @rdname vrFeatureData
#' @order 3
#' @export
setMethod("vrFeatureData", "vrAssayV2", vrFeatureDatavrAssayV2)

vrFeatureDataRreplacevrAssay <- function(object, value) {
  object@featuredata <- value
  return(object)
}

#' @rdname vrFeatureData
#' @order 5
#' @export
setMethod("vrFeatureData<-", "vrAssay", vrFeatureDataRreplacevrAssay)

vrFeatureDataReplacevrAssayV2 <- function(object, feat_type = NULL, value) {
  if(is.null(feat_type))
    feat_type <- vrMainFeatureType(object)
  object@featuredata[[feat_type]] <- value
  return(object)
}

#' @rdname vrFeatureData
#' @order 5
#' @export
setMethod("vrFeatureData<-", "vrAssayV2", vrFeatureDataReplacevrAssayV2)

vrAssayNamesvrAssay <- function(object) {
  
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

#' @rdname vrAssayNames
#' @order 4
#' @export
setMethod("vrAssayNames", "vrAssay", vrAssayNamesvrAssay)

#' @rdname vrAssayNames
#' @order 4
#' @export
setMethod("vrAssayNames", "vrAssayV2", vrAssayNamesvrAssay)

vrAssayNamesReplacevrAssay <- function(object, value){
  
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

#' @param value assay name
#' 
#' @rdname vrAssayNames
#' @order 5
#' @importFrom stringr str_replace
setMethod("vrAssayNames<-", "vrAssay", vrAssayNamesReplacevrAssay)

vrAssayNamesReplacevrAssayV2 <- function(object, value){
  
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

#' @param value assay name
#' 
#' @rdname vrAssayNames
#' @order 5
#' @importFrom stringr str_replace
setMethod("vrAssayNames<-", "vrAssayV2", vrAssayNamesReplacevrAssayV2)

vrAssayTypesvrAssay <- function(object) {
  return(object@type)
}

#' @rdname vrAssayTypes
#' @order 3
#' @export
setMethod("vrAssayTypes", "vrAssay", vrAssayTypesvrAssay)

#' @rdname vrAssayTypes
#' @order 3
#' @export
setMethod("vrAssayTypes", "vrAssayV2", vrAssayTypesvrAssay)

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
      message(param, " not found in the param list")
      return(NULL)
    }
  } else {
    return(object@params)
  }
}

vrDatavrAssay <- function(object, features = NULL, feat_type = NULL, norm = FALSE, ...) {
  
  # get assay types
  assay.type <- vrAssayTypes(object)
  
  # for ROIs, cells and spots
  if(assay.type %in% c("ROI", "cell", "spot")){
    
    # check if there are features
    if(!is.null(features)){
      if(!all(features %in% vrFeatures(object))){
        stop("Some features are not available in the assay!")
      }
      
      if(inherits(object, "vrAssay")){
        if(norm){
          return(object@normdata[features,,drop = FALSE])
        } else {
          return(object@rawdata[features,,drop = FALSE])
        }
      } else {
        if(is.null(feat_type))
          feat_type <- vrMainFeatureType(object)
        if(norm){
          return(object@data[[paste0(feat_type, "_norm")]][features,,drop = FALSE])
        } else {
          return(object@data[[feat_type]][features,,drop = FALSE])
        }
      }
      
      # if there are no features requested, return the data
    } else {
      
      if(inherits(object, "vrAssay")){
        if(norm){
          return(object@normdata)
        } else {
          return(object@rawdata)
        }
      } else {
        if(is.null(feat_type))
          feat_type <- vrMainFeatureType(object)
        if(norm){
          return(object@data[[paste0(feat_type, "_norm")]])
        } else {
          return(object@data[[feat_type]])
        }
      }
    }
    
    # for tiles and molecules
  } else {
    
    # check if features are requested
    if(!is.null(features)){
      stop("No features are available for tile and molecule assays!")
    } else{
      
      if(inherits(object, "vrAssay")){
        if(norm){
          return(object@normdata)
        } else {
          return(object@rawdata)
        }
      } else {
        if(is.null(feat_type))
          feat_type <- vrMainFeatureType(object)
        if(norm){
          return(object@data[[paste0(feat_type, "_norm")]])
        } else {
          return(object@data[[feat_type]])
        }
      }
    }
  }
}

#' @rdname vrData
#' @order 3
#'
#' @importFrom magick image_raster
#'
#' @export
setMethod("vrData", "vrAssay", vrDatavrAssay)

#' @rdname vrData
#' @order 3
#'
#' @export
setMethod("vrData", "vrAssayV2", vrDatavrAssay)

generateTileDatavrAssay <- function(object, name = NULL, reg = FALSE, channel = NULL) {
  
  if(vrAssayTypes(object) != "tile"){
    stop("generateTileData can only be used for tile-based assays")
  } else {
    image_data <- as.numeric(vrImages(object, name = name, reg = reg, channel = channel, as.raster = TRUE))
    image_data <- (0.299 * image_data[,,1] + 0.587 * image_data[,,2] + 0.114 * image_data[,,3])
    image_data <- split_into_tiles(image_data, tile_size = vrAssayParams(object, param = "tile.size"))
    image_data <- sapply(image_data, function(x) return(as.vector(x)))
    image_data <- image_data*255
    rownames(image_data) <- paste0("pixel", 1:nrow(image_data))
    colnames(image_data) <- vrSpatialPoints(object)
    feat_type <- vrMainFeatureType(object)
    
    if(inherits(object, "vrAssay")){
      object@rawdata <- object@normdata <- image_data
    } else{
      object@data[[feat_type]] <- image_data
      object@data[[paste0(feat_type, "_norm")]] <- image_data
    }
  }
  return(object)
}

#' @param name the name of the main spatial system
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainSpatial}) is requested
#' @param channel the name of the channel associated with the image
#' 
#' @rdname generateTileData
#' @order 3
#'
#' @export
setMethod("generateTileData", "vrAssay", generateTileDatavrAssay)

#' @rdname generateTileData
#' @order 3
#'
#' @export
setMethod("generateTileData", "vrAssayV2", generateTileDatavrAssay)

vrCoordinatesvrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object)
  }
  
  # check registered coordinates
  if(reg){
    if(!paste0(image_name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered spatial systems with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }
  
  # check coordinates
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any spatial system in this vrAssay object")
  }
  
  # return coordinates
  return(vrCoordinates(object@image[[image_name]]))
}

#' @rdname vrCoordinates
#' @order 3
#' @export
#'
setMethod("vrCoordinates", "vrAssay", vrCoordinatesvrAssay)

#' @rdname vrCoordinates
#' @order 3
#' @export
#'
setMethod("vrCoordinates", "vrAssayV2", vrCoordinatesvrAssay)

vrCoordinatesReplacevrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object)
  }
  
  # check registered coordinates
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }
  
  # check coordinates
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any spatial system in this vrAssay object")
  }
  
  vrCoordinates(object@image[[image_name]]) <- value
  return(object)
}

#' @rdname vrCoordinates
#' @order 5
#' @importFrom methods slot
#'
#' @export
setMethod("vrCoordinates<-", "vrAssay", vrCoordinatesReplacevrAssay)

#' @rdname vrCoordinates
#' @order 5
#' @importFrom methods slot
#'
#' @export
setMethod("vrCoordinates<-", "vrAssayV2", vrCoordinatesReplacevrAssay)

flipCoordinatesvrAssay <- function(object, image_name = NULL, spatial_name = NULL, ...) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # get coordinates
  coords <- vrCoordinates(object, image_name = image_name, ...)
  
  # get image info
  image <- vrImages(object, name = image_name)
  if(!is.null(image)){
    imageinfo <- magick::image_info(vrImages(object, name = image_name))
    height <- imageinfo$height
  } else{
    height <- max(coords[,"y"])
  }
  
  # flip coordinates
  coords[,"y"] <- height - coords[,"y"]
  vrCoordinates(object, image_name = image_name, ...) <- coords
  
  # flip segments
  segments <- vrSegments(object, image_name = image_name, ...)
  if(length(segments) > 0){
    name_segments <- names(segments)
    segments <- do.call("rbind", segments)
    segments[,"y"] <- height - segments[,"y"]
    segments <- split(segments, segments[,1])
    names(segments) <- name_segments
    vrSegments(object, image_name = image_name, ...) <- segments
  }
  
  # return
  return(object)
}

#' @rdname flipCoordinates
#' @order 3
#'
#' @importFrom magick image_info
#'
#' @export
setMethod("flipCoordinates", "vrAssay", flipCoordinatesvrAssay)

#' @rdname flipCoordinates
#' @order 3
#'
#' @export
setMethod("flipCoordinates", "vrAssayV2", flipCoordinatesvrAssay)

vrSegmentsvrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object)
  }
  
  # check registered segments
  if(reg){
    if(!paste0(image_name, "_reg") %in% vrSpatialNames(object)){
      warning("There are no registered spatial systems with name ", image_name, "!")
    } else {
      image_name <- paste0(image_name, "_reg")
    }
  }
  
  # check coordinates
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any spatial system in this vrAssay object")
  }
  
  # return coordinates
  return(vrSegments(object@image[[image_name]]))
}

#' @rdname vrSegments
#' @order 3
#' @export
setMethod("vrSegments", "vrAssay", vrSegmentsvrAssay)

#' @rdname vrSegments
#' @order 3
#' @export
setMethod("vrSegments", "vrAssayV2", vrSegmentsvrAssay)

vrSegmentsReplacevrAssay <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # check main image
  if(is.null(image_name)){
    image_name <- vrMainSpatial(object)
  }
  
  # check registered segments
  if(reg){
    image_name <- paste0(image_name, "_reg")
  }
  
  # check coordinates
  if(!image_name %in% vrSpatialNames(object)){
    stop(image_name, " is not among any spatial system in this vrAssay object")
  }
  
  vrSegments(object@image[[image_name]]) <- value
  return(object)
}

#' @rdname vrSegments
#' @order 6
#' @importFrom methods slot
#' @export
setMethod("vrSegments<-", "vrAssay", vrSegmentsReplacevrAssay)

#' @rdname vrSegments
#' @order 6
#' @importFrom methods slot
#' @export
setMethod("vrSegments<-", "vrAssayV2", vrSegmentsReplacevrAssay)

vrEmbeddingsvrAssay <- function(object, type = "pca", dims = 1:30) {
  
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
#' @order 3
#' @export
setMethod("vrEmbeddings", "vrAssay", vrEmbeddingsvrAssay)

#' @rdname vrEmbeddings
#' @order 3
#' @export
#'
setMethod("vrEmbeddings", "vrAssayV2", vrEmbeddingsvrAssay)

vrEmbeddingsReplacevrAssay <- function(object, type = "pca", value) {
  object@embeddings[[type]] <- value
  return(object)
}

#' @rdname vrEmbeddings
#' @order 4
#' @export
setMethod("vrEmbeddings<-", "vrAssay", vrEmbeddingsReplacevrAssay)

vrEmbeddingsReplacevrAssayV2 <- function(object, type = "pca", value) {
  object@embeddings[[type]] <- value
  return(object)
}

#' @rdname vrEmbeddings
#' @order 4
#' @export
setMethod("vrEmbeddings<-", "vrAssayV2", vrEmbeddingsReplacevrAssayV2)

vrEmbeddingNamesvrAssay <- function(object){
  return(names(object@embeddings))
}

#' @rdname vrEmbeddingNames
#' @order 3
#'
#' @export
setMethod("vrEmbeddingNames", "vrAssay", vrEmbeddingNamesvrAssay)

#' @rdname vrEmbeddingNames
#' @order 3
#'
#' @export
setMethod("vrEmbeddingNames", "vrAssayV2", vrEmbeddingNamesvrAssay)