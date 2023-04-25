####
# SpaceRover classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('igraph'))

## SpaceRover ####

#' The SpaceRover Class
#'
#' @slot samples A list of samples for the this project
#' @slot integrated.datasets A list of integrated data objects that indicate integrated spatial layers
#' @slot meta.data Contains meta-information about each sample
#' @slot project Name of the project
#'
#' @name SpaceRover-class
#' @rdname SpaceRover-class
#' @exportClass SpaceRover
#'
SpaceRover <- setClass(
  Class = 'SpaceRover',
  slots = c(
    samples = 'list',
    metadata = "srMetadata",
    sample.metadata = "data.frame",
    zstack = "igraph",
    main.assay = "character",
    project = 'character'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'SpaceRover',
  definition = function(object) {

    # print class
    cat(class(x = object), "Object \n")

    # print samples and layers
    for(samp in names(object@samples)){
      cat(samp, ": \n", sep = "")
      layers <- names(unlist(object@samples[[samp]]@layer))
      cat("  Layers:", paste(layers, collapse = " "), "\n")
    }

    # return invisible
    return(invisible(x = NULL))
  }
)

### subset of samples ####
setMethod(
  f = '[[',
  signature = c('SpaceRover', "character", "missing"),
  definition = function(x, i, j, ...){

    # sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){
      stop("There are no samples named ", i, " in this object")
    }

    return(x@samples[[i]])
  }
)

setMethod(
  f = '[[<-',
  signature = c('SpaceRover', "character", "missing"),
  definition = function(x, i, j, ..., value){

    # sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){
      stop("There are no samples named ", i, " in this object")
    }

    if(!class(value) == "srSample"){
      stop("The provided object is not of class srSample")
    }

    x@samples[[i]] <- value
    return(x)
  }
)

### subset of samples and layers ####
setMethod(
  f = '[[',
  signature = c('SpaceRover', "character", "character"),
  definition = function(x, i, j, ...){
    return(x[[i]]@layer[[j]])
  }
)

setMethod(
  f = '[[<-',
  signature = c('SpaceRover', "character", "character"),
  definition = function(x, i, j, ..., value){

    if(!class(value) == "srLayer"){
      stop("The provided object is not of class srLayer")
    }

    x[[i]]@layer[[j]] <- value
    return(x)
  }
)

####
# Methods ####
####

### Create SpaceRover object ####

#' CreateSpaceRover
#'
#' Create a SpaceRover object
#'
#' @param data the count table
#' @param metadata a metadata object of class \code{srMetadata}
#' @param image the image of the data
#' @param coord the coordinates of the spatial entities
#' @param sample.metadata a data frame of the sample metadata
#' @param zstack the zstack graph to determine the adjacency of spatial entities across layers
#' @param main.assay the name of the main assay of the object
#' @param assay_name the name of the assay
#' @param sample_name the name of the sample
#' @param layer_name the name of the layer
#' @param assay.type the type of the assay (cells, spots, ROIs)
#' @param project project name
#'
#' @export
#' @import igraph
#'
CreateSpaceRover <- function(data, metadata = NULL, image = NULL, coords,
                             sample.metadata = NULL, zstack = NULL,
                             main.assay = NULL, assay.type = "cell",
                             sample_name = NULL, layer_name = NULL,
                             project = NULL, ...){

  # set project name
  if(is.null(project))
    project <- "SpaceRover"

  # labels
  layer_name <- ifelse(is.null(layer_name), "Section1", layer_name)
  sample_name <- ifelse(is.null(sample_name), "Sample1", sample_name)
  if(is.null(colnames(data))){
    entityID <- paste0(assay.type,1:ncol(data))
    colnames(data) <- entityID
  } else {
    entityID <- colnames(data)
  }
  entityID <- paste(entityID, paste(c(main.assay, layer_name, sample_name), collapse = "_"), sep = "_")

  # set meta data
  if(is.null(metadata)){
    metadata <- setSRMetadata(cell = data.frame(), spot = data.frame(), ROI = data.frame())
    slot(metadata, name = assay.type) <- data.frame(Count = colSums(data), row.names = entityID)
  }

  # Coordinates
  if(!is.null(coords)){
    colcoords <- colnames(coords)
    if(all(colcoords %in% c("x","y"))){
      rownames(coords) <- entityID
      coords <- coords[,c("x", "y")]
    } else {
      stop("The colnames of the coordinates matrix should be 'x' and 'y'")
    }
  } else {
    stop("There are no coordinates matrix provided!")
  }

  # set zgraph
  if(is.null(zstack)){
    spatial_entities <- Entities(metadata)
    zstack <- igraph::make_empty_graph(n = length(spatial_entities), directed = FALSE)
    igraph::V(zstack)$name <- spatial_entities
  }

  # create srAssay
  Xenium_assay <- new("srAssay", rawdata = data, normdata = data, coords = coords, image = image, type = assay.type)
  listofAssays <- list(Xenium_assay)
  names(listofAssays) <- main.assay

  # create layers and samples
  listofLayers <- list(new("srLayer", assay = listofAssays))
  names(listofLayers) <- layer_name
  listofSamples <- list(new("srSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setSRSampleMetadata(listofSamples)
  }

  # set SpaceRover class
  new("SpaceRover", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata, zstack = zstack, main.assay = main.assay, project = project)
}

#' CreateSpaceRover
#'
#' Create a SpaceRover object
#'
#' @param samples a list of samples list of class \code{srSamples}
#' @param metadata a metadata object of class \code{srMetadata}
#' @param sample.metadata a data frame of the sample metadata
#' @param zstack the zstack graph to determine the adjacency of spatial entities across layers
#' @param main.assay the name of the main assay of the object
#' @param project project name
#'
#' @export
#' @import igraph
#'
CreateSpaceRover_old <- function(samples, metadata = NULL, sample.metadata = NULL, zstack = NULL, main.assay = NULL, project = NULL){

  # set project name
  if(is.null(project))
    project <- "SpaceRover"

  # check for samples
  if(!is.list(samples))
    stop("Please introduce a list of samples")

  # set meta data
  if(is.null(metadata)){
    metadata <- setSRMetadata(cell = data.frame(), spot = data.frame(), ROI = data.frame())
  }

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setSRSampleMetadata(samples)
  }

  # set zgraph
  if(is.null(zstack)){
    spatial_entities <- Entities(metadata)
    zstack <- igraph::make_empty_graph(n = length(spatial_entities), directed = FALSE)
    igraph::V(zstack)$name <- spatial_entities
  }

  # # labels
  # layer_name <- ifelse(is.null(layer_name), "slide1", layer_name)
  # sample_name <- ifelse(is.null(sample_name), "sample1", sample_name)
  # cellID <- paste(paste0("Cell",1:ncol(rawdata)),
  #                 paste(c(assay_name, layer_name, sample_name), collapse = "_"),
  #                 sep = "_")
  # colnames(rawdata) <- cellID

  # set SpaceRover class
  new("SpaceRover", samples = samples, metadata = metadata, sample.metadata = sample.metadata, zstack = zstack, main.assay = main.assay, project = project)
}

### Merge SpaceRover objects ####

#' @method merge SpaceRover
#'
#' @import igraph
#'
#' @export
#'
merge.SpaceRover <- function(object, object_list, sample_name = NULL, main.assay = NULL) {

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)

  # check if all are spaceRover
  if(!all(lapply(object_list, class) == "SpaceRover"))
     stop("All arguements have to be of SpaceRover class")

  # merge metadata
  metadata_list <- lapply(object_list, function(x) slot(x, name = "metadata"))
  metadata <- merge(metadata_list[[1]], metadata_list[-1])

  # merge sample metadata
  sample.metadata <- NULL
  for(i in 1:length(object_list)){
    sample.metadata <- rbind(sample.metadata, slot(object_list[[i]], "sample.metadata"))
  }
  rownames(sample.metadata) <- paste0("Assay", 1:nrow(sample.metadata))
  if(!is.null(sample_name)){
    sample.metadata$Sample <- sample_name
    sample.metadata$Layer <- paste0("Section", 1:nrow(sample.metadata))
    unique_assay <- unique(sample.metadata$Assay)
    if(nrow(sample.metadata) != length(unique_assay)){
      for(cur_assay in unique_assay){
        cur_assay_ind <- which(sample.metadata$Assay %in% cur_assay)
        sample.metadata$Assay[cur_assay_ind] <- paste0(sample.metadata$Assay[cur_assay_ind], "_", 1:length(cur_assay_ind))
      }
    }
  }

  # combine samples and rename layers
  if(!is.null(sample_name)){
    listofLayers <- NULL
    for(i in 1:length(object_list)){
      cur_object <- object_list[[i]]
      listofLayers <- c(listofLayers, cur_object@samples[[1]]@layer)
    }
    names(listofLayers) <- sample.metadata$Layer
    listofSamples <- list(new("srSample", layer = listofLayers))
    names(listofSamples) <- sample_name
  } else {
    listofSamples <- NULL
    for(i in 1:length(object_list)){
      cur_object <- object_list[[i]]@samples
      listofSamples <- c(listofSamples, cur_object)
    }
  }

  # get main assay
  if(is.null(main.assay))
    main.assay <- "Assay1"

  # merge graphs
  zstack_list <- lapply(object_list, function(x) slot(x, name = "zstack"))
  zstack <- igraph::disjoint_union(zstack_list[1], zstack_list[-1])

  # project
  project <- slot(object_list[[1]], "project")

  # set SpaceRover class
  new("SpaceRover", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata,
      zstack = zstack, main.assay = main.assay, project = project)
}

### Get main assay ####

#' @rdname MainAssay
#' @method MainAssay SpaceRover
#'
#' @export
#'
MainAssay.SpaceRover <- function(object, ...) {

  # get first sample and first layer with the main assay
  main.assay <- object@main.assay
  sample.info <- object@sample.metadata
  sample.info <- sample.info[sample.info$Assay == main.assay,, drop = FALSE]
  sample.info <- sample.info[1,] # GET FIRST FOR NOW!!!!

  # return assay
  return(object[[sample.info$Sample, sample.info$Layer]]@assay[[main.assay]])
}

#' @rdname MainAssay
#' @method MainAssay<- SpaceRover
#'
#' @export
#'
"MainAssay<-.SpaceRover" <- function(object, ..., value) {

  # check class
  if(class(value) != "srAssay") {
    stop("The provided object is not of srAssay class")
  } else {
    # get first sample and first layer with the main assay
    main.assay <- object@main.assay
    sample.info <- object@sample.metadata
    sample.info <- sample.info[sample.info$Assay == main.assay,, drop = FALSE]
    sample.info <- sample.info[1,] # GET FIRST FOR NOW!!!!
    object[[sample.info$Sample, sample.info$Layer]]@assay[[main.assay]] <- value
  }

  return(object)
}


### Get entities ####

#' @rdname Entities
#' @method Entities SpaceRover
#'
#' @export
#'
Entities.SpaceRover <- function(object, ...) {
  return(Entities(object@metadata))
}

### Get coordinates ####

#' @rdname Coordinates
#' @method Coordinates SpaceRover
#'
#' @export
#'
Coordinates.SpaceRover <- function(object, reg = FALSE, ...) {

  # check existing images in the spacerover object
  assay <- MainAssay(object)

  # get image from the assay
  if(reg){
    coords <- assay@coords_reg
  } else {
    coords <- assay@coords
  }

  # return image
  return(coords)
}

#' @rdname Coordinates
#' @method Coordinates<- SpaceRover
#'
#' @export
#'
"Coordinates<-.SpaceRover" <- function(object, ..., value) {

  # get assay
  assay <- MainAssay(object)
  Coordinates(assay) <- value
  MainAssay(object) <- assay

  return(object)
}

