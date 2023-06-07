#' @include zzz.R
#' @include generics.R
#'
#' @useDynLib spaceRover
NULL

####
# VoltRon classes ####
####

## Auxiliary ####

# Set magick-image as an S4 class
setOldClass(Classes = c('igraph'))

## VoltRon ####

#' The VoltRon Class
#'
#' @slot samples A list of samples for the this project
#' @slot integrated.datasets A list of integrated data objects that indicate integrated spatial layers
#' @slot meta.data Contains meta-information about each sample
#' @slot project Name of the project
#'
#' @name VoltRon-class
#' @rdname VoltRon-class
#' @exportClass VoltRon
#'
VoltRon <- setClass(
  Class = 'VoltRon',
  slots = c(
    samples = 'list',
    metadata = "vrMetadata",
    sample.metadata = "data.frame",
    zstack = "igraph",
    main.assay = "character",
    project = 'character'
  )
)

### show ####
setMethod(
  f = 'show',
  signature = 'VoltRon',
  definition = function(object) {

    # print class
    cat(class(x = object), "Object \n")

    # print samples and layers
    all_assays <- NULL
    for(samp in names(object@samples)){
      cat(samp, ": \n", sep = "")
      layers <- names(unlist(object@samples[[samp]]@layer))
      cat("  Layers:", paste(layers, collapse = " "), "\n")
      assays <- sapply(names(object@samples[[samp]]@layer), function(x) names(object[[samp, x]]@assay))
      all_assays <- c(all_assays, assays)
    }

    # print assays
    unique_assays <- unique(all_assays)
    unique_assays <- unique_assays[c(which(unique_assays == object@main.assay),which(unique_assays != object@main.assay))]
    unique_assays[1] <- paste0(unique_assays[1], " (Main)")
    cat("Assays:", paste(unique_assays, collapse = " "), "\n")

    # return invisible
    return(invisible(x = NULL))
  }
)

### $ method ####

#' @export
#' @method $ VoltRon
#'
"$.VoltRon" <- function(x, i, ...) {
  return(SampleMetadata(x)[[i]])
}

#' @export
#' @method $<- VoltRon
#'
"$<-.VoltRon" <- function(x, i, ..., value) {
  if(nrow(SampleMetadata(x)) > 1)
    stop("You can only change the name of a single name")

  # update sample names
  if(i == "Sample"){
    names(x@samples) <- value
  }

  # update sample metadata and metadata
  x@sample.metadata[[i]] <- value
  x@metadata[[i]] <- value

  return(x)
}

### subset of samples and layers ####

setMethod(
  f = '[[',
  signature = c('VoltRon', "character", "missing"),
  definition = function(x, i, j, ...){

    # if no assay were found, check sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){

      # check assays
      sample.metadata <- SampleMetadata(x)
      assay_names <- rownames(sample.metadata)
      if(i %in% assay_names){
        cur_assay <- sample.metadata[i,]
        return(x@samples[[cur_assay$Sample]]@layer[[cur_assay$Layer]]@assay[[cur_assay$Assay]])
      } else {
        stop("There are no samples or assays named ", i, " in this object")
      }

    } else {
      return(x@samples[[i]])
    }
  }
)

setMethod(
  f = '[[<-',
  signature = c('VoltRon', "character", "missing"),
  definition = function(x, i, j, ..., value){

    # sample names
    sample_names <- names(slot(x, "samples"))

    # check query sample name
    if(!i %in% sample_names){

      # check assays
      sample.metadata <- SampleMetadata(x)
      assay_names <- rownames(sample.metadata)
      if(i %in% assay_names){
        cur_assay <- sample.metadata[i,]
        x@samples[[cur_assay$Sample]]@layer[[cur_assay$Layer]]@assay[[cur_assay$Assay]] <- value
      } else {
        stop("There are no samples named ", i, " in this object")
      }
    } else {
      if(!class(value) == "vrSample"){
        stop("The provided object is not of class vrSample")
      }
      x@samples[[i]] <- value
    }
    return(x)
  }
)

setMethod(
  f = '[[',
  signature = c('VoltRon', "character", "character"),
  definition = function(x, i, j, ...){
    return(x[[i]]@layer[[j]])
  }
)

setMethod(
  f = '[[<-',
  signature = c('VoltRon', "character", "character"),
  definition = function(x, i, j, ..., value){

    if(!class(value) == "vrLayer"){
      stop("The provided object is not of class vrLayer")
    }

    x[[i]]@layer[[j]] <- value
    return(x)
  }
)

####
# Methods ####
####

### Create VoltRon object ####

#' CreateVoltRon
#'
#' Create a VoltRon object
#'
#' @param data the count table
#' @param metadata a metadata object of class \code{vrMetadata}
#' @param image the image of the data
#' @param coord the coordinates of the spatial entities
#' @param segments the segments of the spatial entities, optional
#' @param sample.metadata a data frame of the sample metadata
#' @param zstack the zstack graph to determine the adjacency of spatial entities across layers
#' @param main.assay the name of the main assay of the object
#' @param assay_name the name of the assay
#' @param assay.type the type of the assay (cells, spots, ROIs)
#' @param params additional parameters of the object
#' @param sample_name the name of the sample
#' @param layer_name the name of the layer
#' @param project project name
#' @param ... additional parameters passed to VoltRon object
#'
#' @export
#' @import igraph
#'
CreateVoltRon <- function(data, metadata = NULL, image = NULL,
                             coords, segments = list(),
                             sample.metadata = NULL, zstack = NULL,
                             main.assay = "Custom_cell", assay.type = "cell", params = list(),
                             sample_name = NULL, layer_name = NULL,
                             project = NULL, ...){

  # set project name
  if(is.null(project))
    project <- "VoltRon"

  # layer and sample names
  layer_name <- ifelse(is.null(layer_name), "Section1", layer_name)
  if(main.assay == layer_name)
      stop(paste0("'", layer_name, "' cannot be a layer name, since main assay is named '", main.assay, "'."))
  sample_name <- ifelse(is.null(sample_name), "Sample1", sample_name)
  if(main.assay == sample_name)
    stop(paste0("'", sample_name, "' cannot be a sample name, since main assay is named '", main.assay, "'."))

  # entity IDs
  if(is.null(colnames(data))){
    entityID <- paste0(assay.type,1:ncol(data))
    colnames(data) <- entityID
  } else {
    entityID <- colnames(data)
  }
  entityID <- paste(entityID, "Assay1", sep = "_")
  colnames(data) <- entityID

  # set meta data
  if(is.null(metadata)){
    sr_metadata <- setVRMetadata(cell = data.frame(), spot = data.frame(), ROI = data.frame())
    slot(sr_metadata, name = assay.type) <- data.frame(Count = colSums(data), Assay = main.assay, Layer = layer_name, Sample = sample_name, row.names = entityID)
  } else {
    if(any(class(metadata) %in% c("data.frame", "matrix"))){
      sr_metadata <- setVRMetadata(cell = data.frame(), spot = data.frame(), ROI = data.frame())
      if(any(!rownames(metadata) %in% gsub("_Assay1$", "", entityID))){
        stop("Entity IDs are not matching")
      } else {
        metadata <- metadata[gsub("_Assay1$", "", entityID),]
        slot(sr_metadata, name = assay.type) <- data.frame(metadata, Count = colSums(data), Assay = main.assay, Layer = layer_name, Sample = sample_name, row.names = entityID)
      }
    }
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

  # Segments
  if(length(segments) > 0) names(segments) <- entityID

  # set zgraph
  if(is.null(zstack)){
    spatial_entities <- Entities(sr_metadata)
    zstack <- igraph::make_empty_graph(n = length(spatial_entities), directed = FALSE)
    igraph::V(zstack)$name <- spatial_entities
  }

  # create vrAssay
  Xenium_assay <- new("vrAssay", rawdata = data, normdata = data, coords = coords, segments = segments, image = image, params = params, type = assay.type)
  listofAssays <- list(Xenium_assay)
  names(listofAssays) <- main.assay

  # create layers and samples
  listofLayers <- list(new("vrLayer", assay = listofAssays))
  names(listofLayers) <- layer_name
  listofSamples <- list(new("vrSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setVRSampleMetadata(listofSamples)
  }

  # set VoltRon class
  new("VoltRon", samples = listofSamples, metadata = sr_metadata, sample.metadata = sample.metadata, zstack = zstack, main.assay = main.assay, project = project)
}

### Assay Methods ####

#' @rdname MainAssay
#' @method MainAssay VoltRon
#'
#' @export
#'
MainAssay.VoltRon <- function(object, ...) {
  object@main.assay
}

#' @rdname MainAssay
#' @method MainAssay<- VoltRon
#'
#' @export
#'
"MainAssay<-.VoltRon" <- function(object, ..., value) {
  assay_names <- unique(object@sample.metadata$Assay)
  if(!value %in% assay_names){
    stop("There is no assay names '", value, "' in this object")
  } else {
    object@main.assay <- value
  }
  return(object)
}

#' @rdname AddAssay
#' @method AddAssay VoltRon
#'
#' @importfrom igraph union
#' @export
#'
AddAssay.VoltRon <- function(object, assay, assay_name, sample = "Sample1", layer = "Section1"){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay id
  assay_ids <- as.numeric(gsub("Assay", "", rownames(sample.metadata)))
  assay_id <- paste0("Assay", max(assay_ids)+1)
  assay_names <- c(rownames(sample.metadata), assay_id)

  # update sample.metadata and metadata
  object@sample.metadata <- rbind(sample.metadata, c(assay_name, layer, sample))
  rownames(object@sample.metadata) <- assay_names
  object@metadata <- AddAssay(object@metadata,
                              assay = assay, assay_name = assay_name,
                              sample = sample, layer = layer)

  # update sample and layer
  assay_list <- object[[sample, layer]]@assay
  AssayNames(assay) <- assay_id
  new_assay_list <- list(assay)
  names(new_assay_list) <- assay_name
  assay_list <- c(assay_list, new_assay_list)
  object[[sample, layer]]@assay <- assay_list

  # update graph
  newgraph <- igraph::make_empty_graph(n = length(Entities(assay)), directed = FALSE)
  igraph::V(newgraph)$name <- Entities(assay)
  object@zstack <- igraph::union(object@zstack, newgraph)

  # return
  return(object)
}

#' @rdname AssayNames
#' @method AssayNames VoltRon
#'
#' @export
#'
AssayNames.VoltRon <- function(object, assay = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check assays
  if(is.null(assay))
    assay <- MainAssay(object)

  # get assay names
  if(assay %in% sample.metadata$Assay){
    assay_names <- rownames(sample.metadata)[sample.metadata$Assay %in% assay]
  } else {
    if(assay %in% rownames(sample.metadata)) {
      assay_names <- assay
    } else {
      stop("Assay name or type is not found in the object")
    }
  }

  return(assay_names)
}

#' @rdname AssayTypes
#' @method AssayTypes VoltRon
#'
#' @export
#'
AssayTypes.VoltRon <- function(object, assay = NULL){

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get assay types
  assay_types <- sapply(assay_names, function(x) AssayTypes(object[[x]]))

  return(assay_types)
}

### Object Methods ####

#' @method subset VoltRon
#'
#' @importFrom rlang enquo eval_tidy quo_get_expr
#' @import igraph
#' @importFrom stringr str_extract
#' @import shiny
#' @import shinyjs
#'
#' @export
#'
subset.VoltRon <- function(object, subset, samples = NULL, assays = NULL, entities = NULL, features = NULL, image = NULL, interactive = FALSE) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  if(interactive){
    results <- demuxVoltRon(object)
    return(results)
  }

  # subseting on samples
  if(!missing(subset)){

    metadata <- Metadata(GeoMxR1, type = "ROI")
    entities <- rownames(metadata)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
    object <- subset(object, entities = entities)
    return(object)

  } else if(!is.null(samples)){

    sample.metadata <- subset.sampleMetadata(object@sample.metadata, samples = samples)
    metadata <- subset.vrMetadata(object@metadata, samples = samples) # CAN WE CHANGE THIS TO ONLY SUBSET LATER ????
    listofSamples <- object@samples[samples]

  # subsetting on assays name
  } else if(!is.null(assays)) {

    sample.metadata <- subset.sampleMetadata(object@sample.metadata, assays = assays)
    metadata <- subset.vrMetadata(object@metadata, assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, assays = assays)
    }, USE.NAMES = TRUE)

  # subsetting on entity names
  } else if(!is.null(entities)) {

    metadata <- subset.vrMetadata(object@metadata, entities = entities)
    assays <- unique(stringr::str_extract(Entities(metadata), "Assay[0-9]+"))
    sample.metadata <- subset.sampleMetadata(object@sample.metadata, assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, entities = entities)
    }, USE.NAMES = TRUE)

  # subsetting on features
  } else if(!is.null(features)){

    sample.metadata <- SampleMetadata(object)
    assay_names <- AssayNames(object)
    for(assy in assay_names){
      cur_assay <- sample.metadata[assy,]
      vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
      vrassay <- vrlayer[[cur_assay$Assay]]
      vrassay <- subset.vrAssay(vrassay, features = features)
      vrlayer[[cur_assay$Assay]] <- vrassay
      object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer
    }
    metadata <- object@metadata
    listofSamples <- object@samples

  # subsetting on image
  } else if(!is.null(image)) {

    # subsetting based on image magick parameters
    if(class(image) == "character") {

      # check if there are only one image and one assay
      if(nrow(object@sample.metadata) > 1){
        stop("Subseting on images can only be performed on VoltRon objects with a single assay")
      } else {
        sample.metadata <- object@sample.metadata
        samples <- unique(sample.metadata$Sample)
        listofSamples <- sapply(object@samples[samples], function(samp) {
          subset.vrSample(samp, image = image)
        }, USE.NAMES = TRUE)
        entities <-  do.call(c, lapply(listofSamples, Entities.vrSample))
        metadata <- subset.vrMetadata(object@metadata, entities = entities)
      }
    } else if(image){
      results <- demuxVoltRon(object)
      return(results)
    }
  }

  # other attributes
  main.assay <- unique(sample.metadata$Assay)[unique(sample.metadata$Assay) == names(table(sample.metadata$Assay))[which.max(table(sample.metadata$Assay))]]
  zstack <- subgraph(object@zstack, V(object@zstack)[names(V(object@zstack)) %in% Entities(metadata)])
  project <- object@project

  # set VoltRon class
  new("VoltRon", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata, zstack = zstack, main.assay = main.assay, project = project)
}

#' Merging VoltRon objects
#'
#' Given a VoltRon object, and a list of VoltRon object, merge all.
#'
#' @param object a VoltRon Object
#' @param object_list a list of VoltRon objects
#' @param sample_name a single sample name if objects are of the same sample
#' @param main.assay name of the assay
#'
#' @export
#' @method merge VoltRon
#'
#' @import igraph
#'
merge.VoltRon <- function(object, object_list, sample_name = NULL, main.assay = NULL) {

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)

  # check if all are spaceRover
  if(!all(lapply(object_list, class) == "VoltRon"))
     stop("All arguements have to be of VoltRon class")

  # merge metadata and sample metadata
  metadata_list <- lapply(object_list, function(x) slot(x, name = "metadata"))
  metadata <- merge(metadata_list[[1]], metadata_list[-1])
  sample.metadata_list <- lapply(object_list, function(x) slot(x, name = "sample.metadata"))
  sample.metadata <- merge.sampleMetadata(sample.metadata_list, sample_name = sample_name)

  # combine samples and rename layers
  if(!is.null(sample_name)){
    listofLayers <- NULL
    for(i in 1:length(object_list)){
      cur_object <- object_list[[i]]
      listofLayers <- c(listofLayers, cur_object@samples[[1]]@layer)
    }
    names(listofLayers) <- sample.metadata$Layer
    listofSamples <- list(new("vrSample", layer = listofLayers))
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
    main.assay <- names(sort(table(sample.metadata$Assay), decreasing = TRUE))[1]

  # merge graphs
  zstack_list <- lapply(object_list, function(x) slot(x, name = "zstack"))
  zstack <- igraph::disjoint_union(zstack_list[1], zstack_list[-1])

  # project
  project <- slot(object_list[[1]], "project")

  # set VoltRon class
  new("VoltRon", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata,
      zstack = zstack, main.assay = main.assay, project = project)
}

#' @rdname Metadata
#' @method Metadata VoltRon
#'
#' @export
#'
Metadata.VoltRon <- function(object, type = "cell") {
  slot(object@metadata, name = type)
}

#' @rdname Metadata
#' @method Metadata<- VoltRon
#'
#' @export
#'
"Metadata<-.VoltRon" <- function(object, type = "cell", ..., value) {
  slot(object@metadata, name = type) <- value
  return(object)
}

#' @rdname SampleMetadata
#' @method SampleMetadata VoltRon
#'
#' @export
#'
SampleMetadata.VoltRon <- function(object, ...) {
  object@sample.metadata
}

#' @rdname Entities
#' @method Entities VoltRon
#'
#' @export
#'
Entities.VoltRon <- function(object, ...) {
  return(Entities(object@metadata))
}

#' @rdname Features
#' @method Features VoltRon
#'
#' @export
#'
Features.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get all features
  features <- NULL
  for(assy in assay_names)
    features <- c(features, Features(object[[assy]]))

  return(unique(features))
}

#' @rdname Data
#' @method Data VoltRon
#'
#' @export
#'
Data.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get all coordinates
  returndata_list <- list()
  for(i in 1:length(assay_names))
    returndata_list[[i]] <- Data(object[[assay_names[i]]], ...)

  return(do.call(cbind, returndata_list))
}

#' @rdname Graph
#' @method Graph VoltRon
#'
#' @export
#'
Graph.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)
  assay_pattern <- paste0(assay_names, collapse = "|")
  node_names <- Entities(object)[grepl(assay_pattern, Entities(object))]

  returngraph <- induced_subgraph(object@zstack, node_names)
  return(returngraph)
}

#' @rdname Coordinates
#' @method Coordinates VoltRon
#'
#' @export
#'
Coordinates.VoltRon <- function(object, reg = FALSE, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get all coordinates
  coords <- NULL
  for(assy in assay_names)
    coords <- rbind(coords, Coordinates(object[[assy]], reg = reg))

  # return image
  return(coords)
}

#' @rdname Coordinates
#' @method Coordinates<- VoltRon
#'
#' @export
#'
"Coordinates<-.VoltRon" <- function(object, reg = FALSE, ..., value) {

  # check the number of assays in the object
  if(nrow(object@sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays are not permitted!")

  # get assay
  cur_assay <- object@sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # change coordinates
  Coordinates(vrassay, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @rdname Segments
#' @method Segments VoltRon
#'
#' @export
#'
Segments.VoltRon <- function(object, reg = FALSE, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get all coordinates
  segts <- NULL
  for(assy in assay_names)
    segts <- rbind(segts, Segments(object[[assy]], reg = reg))

  # return image
  return(segts)
}

#' @rdname Segments
#' @method Segments<- VoltRon
#'
#' @export
#'
"Segments<-.VoltRon" <- function(object, reg = FALSE, ..., value) {

  # check the number of assays in the object
  if(nrow(object@sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays are not permitted!")

  # get assay
  cur_assay <- object@sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # change coordinates
  Segments(vrassay, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @rdname Embeddings
#' @method Embeddings VoltRon
#'
#' @export
#'
Embeddings.VoltRon <- function(object, assay = NULL, type = "pca", ..., value) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get all coordinates
  returndata_list <- list()
  for(i in 1:length(assay_names))
    returndata_list[[i]] <- Embeddings(object[[assay_names[i]]], type = type, ...)

  return(do.call(rbind, returndata_list))
}

#' @rdname Embeddings
#' @method Embeddings<- VoltRon
#'
#' @export
#'
"Embeddings<-.VoltRon" <- function(object, assay = NULL, type = "pca", ..., value) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names){
    assayobject <- object[[assy]]
    Embeddings(assayobject, type = type) <- value[grepl(paste0(assy, "$"), rownames(value)),]
    object[[assy]] <- assayobject
  }

  return(object)
}

