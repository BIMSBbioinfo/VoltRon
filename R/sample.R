## vrSample ####

### subset ####

#' Methods for vrSample objects
#'
#' Methods for \code{\link{vrSample}} objects for generics defined in other
#' packages
#'
#' @param x A vrSample object
#' @param i the name of layer associated with the sample, see \link{SampleMetadata}
#' @param value a vrLayer object, see \link{vrLayer}
#' 
#' @name vrSample-methods
#' @rdname vrSample-methods
#'
#' @concept vrsample
#'
NULL

#' @describeIn vrSample-methods Accessing vrLayer objects from \code{vrSample} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[',
  signature = c('vrSample', "character"),
  definition = function(x, i){

    # sample names
    layer_names <- names(methods::slot(x, "layer"))

    # check query sample name
    if(!i %in% layer_names){
      stop("There are no layers named ", i, " in this sample")
    }

    # return samples
    return(x@layer[[i]])
  }
)


#' @describeIn vrSample-methods Accessing vrLayer objects from \code{vrSample} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[<-',
  signature = c('vrSample', "character"),
  definition = function(x, i, value){

    # check if value if vrLayer
    if(!inherits(value, "vrLayer")){
      stop("The provided object is not of class vrLayer")
    }

    # sample names
    layer_names <- names(methods::slot(x, "layer"))

    # check query sample name
    if(!i %in% layer_names){
      stop("There are no layers named ", i, " in this sample")
    }

    # change layer
    x@layer[[i]] <- value

    # return
    return(x)
  }
)

## vrBlock ####

### subset ####

#' @describeIn vrSample-methods (deprecated) Accessing vrLayer objects from \code{vrBlock} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[',
  signature = c('vrBlock', "character"),
  definition = function(x, i){
    
    # sample names
    layer_names <- names(methods::slot(x, "layer"))
    
    # check query sample name
    if(!i %in% layer_names){
      stop("There are no layers named ", i, " in this sample")
    }
    
    # return samples
    return(x@layer[[i]])
  }
)

#' @describeIn vrSample-methods (deprecated) Overwriting vrLayer objects from \code{vrBlock} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[<-',
  signature = c('vrBlock', "character"),
  definition = function(x, i, value){
    
    # check if value if vrLayer
    if(!inherits(value, "vrLayer")){
      stop("The provided object is not of class vrLayer")
    }
    
    # sample names
    layer_names <- names(methods::slot(x, "layer"))
    
    # check query sample name
    if(!i %in% layer_names){
      stop("There are no layers named ", i, " in this sample")
    }
    
    # change layer
    x@layer[[i]] <- value
    
    # return
    return(x)
  }
)

## vrLayer ####

### subset ####

#' Methods for vrLayer objects
#'
#' Methods for \code{\link{vrLayer}} objects for generics defined in other
#' packages
#'
#' @param x A vrLayer object
#' @param i the name of assay associated with the layer, see \link{SampleMetadata}
#' @param value a vrAssayV2 object, see \link{vrAssayV2}
#' 
#' @name vrLayer-methods
#' @rdname vrLayer-methods
#'
#' @concept vrlayer
#'
NULL

#' @describeIn vrLayer-methods Accessing vrAssay objects from \code{vrLayer} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[',
  signature = c('vrLayer', "character"),
  definition = function(x, i){

    # if no assay were found, check sample names
    assay_names <- names(methods::slot(x, "assay"))

    # check query sample name
    if(!i %in% assay_names){
      stop("There are no assays named ", i, " in this object")
    } else {
      return(x@assay[[i]])
    }
  }
)

#' @describeIn vrLayer-methods Overwriting vrAssay objects from \code{vrLayer} objects
#' 
#' @importFrom methods slot
setMethod(
  f = '[[<-',
  signature = c('vrLayer', "character"),
  definition = function(x, i, value){

    # if no assay were found, check sample names
    assay_names <- names(methods::slot(x, "assay"))

    # check query sample name
    if(!i %in% assay_names){
      stop("There are no assays named ", i, " in this object")
    }

    x@assay[[i]] <- value
    return(x)
  }
)

####
# Methods ####
####

### vrSample Methods ####

mergevrSample <- function(x, y, samples = NULL){
  
  # start
  object <- x
  object_list <- y
   
  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)
  names(object_list) <- samples
  
  # set VoltRon class
  return(object_list)
}

#' Merging vrSample objects
#'
#' Given a vrSample object, and a list of vrSample objects, merge all.
#'
#' @param x a vrSample object
#' @param y a list of vrSample objects
#' @param samples the sample names
#'
#' @method merge vrSample
setMethod("merge", "vrSample", mergevrSample)

#' Merging vrBlock objects
#'
#' Given a vrBlock object, and a list of vrSample objects, merge all.
#' 
#' @param x a vrSample object
#' @param y a list of vrSample objects
#' @param samples the sample names
#' 
#' @method merge vrBlock
setMethod("merge", "vrBlock", mergevrSample)
# merge.vrBlock <- function(object, object_list, samples = NULL){
#   merge.vrSample(object, object_list = object_list, samples = samples)
# }

subsetvrSample <- function(x, subset, assays = NULL, spatialpoints = NULL, image = NULL) {
  
  # start
  object <- x
  
  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }
  
  # subseting on samples, layers and assays
  layers <- object@layer
  if(!is.null(assays)){
    object@layer <- sapply(layers, function(lay) {
      subsetvrLayer(lay, assays = assays)
    }, USE.NAMES = TRUE, simplify = TRUE)
  } else if(!is.null(spatialpoints)){
    object@layer <- sapply(layers, function(lay) {
      subsetvrLayer(lay, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE, simplify = TRUE)
  } else if(!is.null(image)){
    object@layer <- sapply(layers, function(lay) {
      subsetvrLayer(lay, image = image)
    }, USE.NAMES = TRUE, simplify = TRUE)
  }
  
  # remove NULL assays
  ind <- which(vapply(object@layer, function(x) !is.null(x), logical(1)))
  object@layer <- object@layer[ind]
  
  # check if there are layers
  if(length(object@layer) > 0){
    
    # get updated adjaceny and distance
    catch_connect <- try(slot(object, name = "zlocation"), silent = TRUE)
    if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
      object@zlocation <- object@zlocation[ind]
      object@adjacency <- object@adjacency[ind, ind, drop = FALSE]
    }
    
    # return object
    return(object)
  } else {
    return(NULL)
  }
}

#' Subsetting vrSample objects
#'
#' Given a vrSample object, subset the object given one of the attributes
#'
#' @param x a vrSample object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrSample
#' @order 6
#'
#' @importFrom rlang enquo
setMethod("subset", "vrSample", subsetvrSample)

#' Subsetting vrBlock objects
#'
#' Given a vrBlock object, subset the object given one of the attributes
#' 
#' @param x a vrSample object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrBlock
#' @order 6
setMethod("subset", "vrBlock", subsetvrSample)

# subset.vrBlock <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL){
#   subset.vrSample(object, subset = subset, assays = assays, spatialpoints = spatialpoints, image = image)
# }

#' @rdname vrSpatialPoints
#' @order 5
#' @export
setMethod("vrSpatialPoints", "vrSample", function(object) {
  do.call("c", lapply(object@layer, function(lay) {
    vrSpatialPoints(lay)
  }))
})

#' @rdname vrSpatialPoints
#' @order 5
#' @export
setMethod("vrSpatialPoints", "vrBlock", function(object) {
  do.call("c", lapply(object@layer, function(lay) {
    vrSpatialPoints(lay)
  }))
})

changeAssayNamesvrSample <- function(object, sample.metadata = NULL){
  
  if(is.null(sample.metadata))
    stop("Please provide a sample.metadata")
  
  if(!"NewAssayNames" %in% colnames(sample.metadata))
    stop("Please provide a sample.metadata with NewAssayNames column which includes the new assay names")
  
  # change the assay names of the layers
  layer_names <- names(object@layer)
  for(lyr in layer_names)
    object[[lyr]] <- changeAssayNames(object[[lyr]], sample.metadata = sample.metadata[sample.metadata$Layer == lyr,])
  
  # return
  return(object)
}

#' changeAssayNames.vrSample
#'
#' Change the assay names of assays within a vrSample object
#'
#' @param sample.metadata the sample metadata with NewAssayNames column which includes the new assay names
#' 
#' @rdname changeAssayNames
#'
#' @noRd
setMethod("changeAssayNames", "vrSample", changeAssayNamesvrSample)

changeAssayNamesvrBlock <- function(object, sample.metadata = NULL) {
  object <- changeAssayNamesvrSample(object, sample.metadata = sample.metadata)
  return(object)
}

#' changeAssayNames.vrBlock
#'
#' Change the assay names of assays within a vrBlock object
#' 
#' @param sample.metadata the sample metadata with NewAssayNames column which includes the new assay names
#' 
#' @rdname changeAssayNames
#'
#' @noRd
setMethod("changeAssayNames", "vrBlock", changeAssayNamesvrBlock)

### vrLayer Methods ####

subsetvrLayer <- function(x, subset, assays = NULL, spatialpoints = NULL, image = NULL) {
  
  # start
  object <- x
  
  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }
  
  # subseting on samples, layers and assays
  if(!is.null(assays)){
    
    # get assay names of all assays
    assay_names <- vapply(object@assay, vrAssayNames, character(1))
    if(any(assays %in% assay_names)) {
      assays <- intersect(assays, assay_names)
      object@assay  <- object@assay[which(assay_names %in% assays)]
    } else if(any(assays %in% names(object@assay))) {
      object@assay  <- object@assay[names(object@assay) %in% assays]
    } else {
      return(NULL)
    }
    
  } else if(!is.null(spatialpoints)){
    
    # get points connected to queried spatialpoints
    catch_connect <- try(slot(object, name = "connectivity"), silent = TRUE)
    if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
      if(igraph::vcount(object@connectivity) > 0){
        spatialpoints <- getConnectedSpatialPoints(object, spatialpoints)
        object@connectivity <- subset.Connectivity(object@connectivity, spatialpoints)
      }
    }
    
    # subset assays
    object@assay <- sapply(object@assay, function(assy) {
      if(inherits(assy, "vrAssay")){
        # return(subset.vrAssay(assy, spatialpoints = spatialpoints))
        return(subsetvrAssay(assy, spatialpoints = spatialpoints))
      } else {
        # return(subset.vrAssayV2(assy, spatialpoints = spatialpoints))
        return(subsetvrAssay(assy, spatialpoints = spatialpoints))
      }
    }, USE.NAMES = TRUE, simplify = TRUE)
    
  } else if(!is.null(image)){
    object@assay <- sapply(object@assay, function(assy) {
      if(inherits(assy, "vrAssay")){
        # return(subset.vrAssay(assy, image = image))
        return(subsetvrAssay(assy, image = image))
      } else {
        return(subsetvrAssay(assy, image = image))
      }
    }, USE.NAMES = TRUE, simplify = TRUE)
  }
  
  # remove NULL assays
  object@assay <- object@assay[which(vapply(object@assay, function(x) !is.null(x), logical(1)))]
  
  # set VoltRon class
  if(length(object@assay) > 0){
    return(object)
  } else {
    return(NULL)
  }
}

#' Subsetting vrLayer objects
#'
#' Given a vrLayer object, subset the object given one of the attributes
#'
#' @param x a vrLayer object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrLayer
#' @order 7
#'
#' @importFrom rlang enquo
#' @importFrom methods is
setMethod("subset", "vrLayer", subsetvrLayer)

#' @rdname vrSpatialPoints
#' @order 6
#' @export
setMethod("vrSpatialPoints", "vrLayer", function(object) {
  do.call("c", lapply(object@assay, function(assy) {
      vrSpatialPoints(assy)
  }))
})

#' subset.Connectivity
#'
#' Subsetting the connectivity graph of vrLayer using spatial points
#'
#' @param object the connectivity graph of the vrLayer
#' @param spatialpoints the set of spatial points
#'
#' @importFrom igraph induced_subgraph
#'
#' @noRd
subset.Connectivity <- function(object, spatialpoints = NULL){
  return(igraph::induced_subgraph(object, spatialpoints))
}

#' getConnectedSpatialPoints
#'
#' get spatial points connected to other spatial points in the connectivity graph of vrLayer
#'
#' @param object A vrLayer object
#' @param spatialpoints the set of spatial points
#'
#' @importFrom igraph neighborhood V vcount
#'
#' @noRd
getConnectedSpatialPoints <- function(object, spatialpoints = NULL){
  if(igraph::vcount(object@connectivity) > 0){
    spatialpoints <- intersect(spatialpoints, igraph::V(object@connectivity)$name)
    return(names(unlist(igraph::neighborhood(object@connectivity, nodes = spatialpoints))))
  } else {
    return(spatialpoints)
  }
}

changeAssayNamesvrLayer <- function(object, sample.metadata = NULL){

  if(is.null(sample.metadata))
    stop("Please provide a sample.metadata")

  if(!"NewAssayNames" %in% colnames(sample.metadata))
    stop("Please provide a sample.metadata with NewAssayNames column which includes the new assay names")

  # change the assay names of the connectivity graph if exists
  catch_connect <- try(slot(object, name = "connectivity"), silent = TRUE)
  if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
    if(igraph::vcount(object@connectivity) > 0){
      spatialpoints <- igraph::V(object@connectivity)$name
      old_assay_names <- vapply(object@assay, vrAssayNames, character(1))
      new_assay_names <- sample.metadata$NewAssayNames
      cur_spatialpoints <- spatialpoints
      for(i in seq_len(length(old_assay_names))){
        if(old_assay_names[i]!=new_assay_names[i]){
          ind <- grepl(paste0(old_assay_names[i],"$"), spatialpoints)
          cur_spatialpoints[ind] <- gsub(paste0(old_assay_names[i],"$"), new_assay_names[i], spatialpoints[ind])
        }
      }
      igraph::V(object@connectivity)$name <- cur_spatialpoints
    }
  }

  # change the assay names of vrAssays
  assay_names <- names(object@assay)
  for(assy in assay_names)
    vrAssayNames(object[[assy]]) <- rownames(sample.metadata[sample.metadata$Assay == assy,])

  # return
  return(object)
}

#' changeAssayNamesvrLayer
#'
#' Change the assay names of assays within a vrSample object
#'
#' @rdname changeAssayNames
#'
#' @importFrom igraph V V<- vcount
#' @importFrom methods is
#'
#' @noRd
setMethod("changeAssayNames", "vrLayer", changeAssayNamesvrLayer)
