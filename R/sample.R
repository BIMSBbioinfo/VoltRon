####
# Objects and Classes ####
####

## Auxiliary ####

# Set classes
setOldClass(Classes = c('igraph'))

## vrSample ####

#' The vrSample (VoltRon Sample) Class
#'
#' @slot layer A list of layers (vrLayer)
#' @slot zlocation a vector of z coordinates of layers
#' @slot adjacency an adjacency matrix of connected layers within a block
#'
#' @name vrSample-class
#' @rdname vrSample-class
#' @exportClass vrSample
#'
vrSample <- setClass(
  Class = 'vrSample',
  slots = c(
    layer = 'list',
    zlocation = 'numeric',
    adjacency = "matrix"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrSample',
  definition = function(object) {
    cat(class(x = object), "(VoltRon Block) Object \n")
    layers <- names(unlist(object@layer))
    cat("Layer(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset ####

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[',
  signature = 'vrSample',
  definition = function(x, i, j){

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

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[<-',
  signature = c('vrSample'),
  definition = function(x, i, j, ..., value){

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

#' The vrBlock (VoltRon Block) Class
#'
#' @slot layer A list of layers (vrLayer)
#' @slot zlocation a vector of z coordinates of layers
#' @slot adjacency an adjacency matrix of connected layers within a block
#'
#' @name vrBlock-class
#' @rdname vrBlock-class
#' @exportClass vrBlock
#'
vrBlock <- setClass(
  Class = 'vrBlock',
  slots = c(
    layer = 'list',
    zlocation = 'numeric', 
    adjacency = "matrix"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrBlock',
  definition = function(object) {
    cat(class(x = object), "(VoltRon Block) Object \n")
    layers <- names(unlist(object@layer))
    cat("Layer(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset ####

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[',
  signature = 'vrBlock',
  definition = function(x, i, j){
    
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

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[<-',
  signature = c('vrBlock'),
  definition = function(x, i, j, ..., value){
    
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

#' The vrLayer (VoltRon Layer) Class
#'
#' @slot assay A list of assays (vrAssay)
#' @slot connectivity the connectivity graph
#'
#' @name vrLayer-class
#' @rdname vrLayer-class
#' @exportClass vrLayer
#' 
vrLayer <- setClass(
  Class = 'vrLayer',
  slots = c(
    assay = 'list',
    connectivity = 'igraph'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrLayer',
  definition = function(object) {
    cat(class(x = object), "(VoltRon Layer) Object \n")
    layers <- names(unlist(object@assay))
    cat("Assay(s):", paste(layers, collapse = " "), "\n")
    return(invisible(x = NULL))
  }
)

### subset of assays ####

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ...){

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

#' @importFrom methods slot
#' @noRd
setMethod(
  f = '[[<-',
  signature = c('vrLayer', "character"),
  definition = function(x, i, ..., value){

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

#' Merging vrSample objects
#'
#' Given a vrSample object, and a list of vrSample objects, merge all.
#'
#' @param object a vrSample object
#' @param object_list a list of vrSample objects
#' @param samples the sample names
#'
#' @method merge vrSample
merge.vrSample <- function(object, object_list, samples = NULL){

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)
  names(object_list) <- samples

  # set VoltRon class
  return(object_list)
}

#' Merging vrBlock objects
#'
#' Given a vrBlock object, and a list of vrSample objects, merge all.
#' 
#' @param object a vrSample object
#' @param object_list a list of vrSample objects
#' @param samples the sample names
#' 
#' @method merge vrBlock
merge.vrBlock <- function(object, object_list, samples = NULL){
  merge.vrSample(object, object_list = object_list, samples = samples)
}

#' Subsetting vrSample objects
#'
#' Given a vrSample object, subset the object given one of the attributes
#'
#' @param object a vrSample object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrSample
#' @order 6
#'
#' @importFrom rlang enquo
#'
subset.vrSample <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  layers <- object@layer
  if(!is.null(assays)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, assays = assays)
    }, USE.NAMES = TRUE, simplify = TRUE)
  } else if(!is.null(spatialpoints)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE, simplify = TRUE)
  } else if(!is.null(image)){
    object@layer <- sapply(layers, function(lay) {
      subset.vrLayer(lay, image = image)
    }, USE.NAMES = TRUE, simplify = TRUE)
  }

  # remove NULL assays
  ind <- which(sapply(object@layer, function(x) !is.null(x)))
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

#' Subsetting vrBlock objects
#'
#' Given a vrBlock object, subset the object given one of the attributes
#' 
#' @param object a vrSample object
#' @param subset the subset statement
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param image the subseting string passed to \link{image_crop}
#'
#' @method subset vrBlock
#' @order 6
#'
subset.vrBlock <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL){
  subset.vrSample(object, subset = subset, assays = assays, spatialpoints = spatialpoints, image = image)
}

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

#' Subsetting vrLayer objects
#'
#' Given a vrLayer object, subset the object given one of the attributes
#'
#' @param object a vrLayer object
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
#'
subset.vrLayer <- function(object, subset, assays = NULL, spatialpoints = NULL, image = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subseting on samples, layers and assays
  if(!is.null(assays)){

    # get assay names of all assays
    assay_names <- sapply(object@assay, vrAssayNames)
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
        return(subset.vrAssay(assy, spatialpoints = spatialpoints))
      } else {
        return(subset.vrAssayV2(assy, spatialpoints = spatialpoints))
      }
    }, USE.NAMES = TRUE, simplify = TRUE)
    
  } else if(!is.null(image)){
    object@assay <- sapply(object@assay, function(assy) {
      if(inherits(assy, "vrAssay")){
        return(subset.vrAssay(assy, image = image))
      } else {
        return(subset.vrAssayV2(assy, image = image))
      }
    }, USE.NAMES = TRUE, simplify = TRUE)
  }

  # remove NULL assays
  object@assay <- object@assay[which(sapply(object@assay, function(x) !is.null(x)))]

  # set VoltRon class
  if(length(object@assay) > 0){
    return(object)
  } else {
    return(NULL)
  }
}

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
      old_assay_names <- sapply(object@assay, vrAssayNames)
      new_assay_names <- sample.metadata$NewAssayNames
      cur_spatialpoints <- spatialpoints
      for(i in 1:length(old_assay_names)){
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
