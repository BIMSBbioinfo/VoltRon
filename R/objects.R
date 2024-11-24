#' @include zzz.R
#' @include generics.R
#' @useDynLib VoltRon
NULL

####
# VoltRon classes ####
####

## VoltRon ####

#' The VoltRon Class
#'
#' @slot samples A list of samples (vrSample)
#' @slot metadata A vrMetadata object that includes metadata of ROIs, spots, and cells
#' @slot sample.metadata Contains meta-information about each sample, layer and assay
#' @slot graph A igraph object
#' @slot main.assay The type of the main assay (i.e. Visium, Xenium, GeoMx etc.)
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
    graph = "list",
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

    # sample metadata
    sample.metadata <- SampleMetadata(object)

    # get sample and layer names
    sample_names <- unique(sample.metadata$Sample)
    show_length <- min(5,length(sample_names))
    for(samp in sample_names[1:show_length]){
      cat(samp, ": \n", sep = "")
      layers <- unique(sample.metadata$Layer[sample.metadata$Sample == samp])
      # cat("  Layers:", paste(layers, collapse = " "), "\n")
      layers <- split(layers, ceiling(seq_along(layers)/5))
      cat("  Layers:", paste(layers[[1]], collapse = " "), "\n")
      if(length(layers) > 1){
        for(i in 2:length(layers)){
          cat("         ", paste(layers[[i]], collapse = " "), "\n")
        } 
      }
    }

    # get assay names
    unique_assays <- unique(sample.metadata$Assay)

    # print
    if(length(sample_names) > 5){
      cat("...", "\n")
      cat("There are", length(sample_names), "samples in total", "\n")
    }

    # print assays
    main.assay <- vrMainAssay(object)
    unique_assays <- unique_assays[c(which(unique_assays == main.assay),which(unique_assays != main.assay))]
    unique_assays[1] <- paste0(unique_assays[1], "(Main)")
    cat("Assays:", paste(unique_assays, collapse = " "), "\n")
    
    # print features
    main.feat <- vrMainFeatureType(object)
    if(!is.null(main.feat)){
      main.feat <- unique(vrMainFeatureType(object)$Feature)
      unique_features <- vrFeatureTypeNames(object)
      unique_features <- unique_features[c(which(unique_features == main.feat),which(unique_features != main.feat))]
      unique_features[1] <- paste0(unique_features[1], "(Main)")
      cat("Features:", paste(unique_features, collapse = " "), "\n") 
    }

    # return invisible
    return(invisible(x = NULL))
  }
)

####
# Methods ####
####

#' Methods for VoltRon
#'
#' Methods for \code{\link{VoltRon}} objects for generics defined in other
#' packages
#'
#' @param x A VoltRon object
#' @param i,value Depends on the usage
#' \describe{
#'  \item{\code{$}, \code{$<-}}{Name (\code{i}) of a single metadata column from the main assay, see \link{vrMainAssay}}
#'  \item{\code{[[}, \code{[[<-}}{
#'    If only \code{i} is given, either a vrSample object or a vrAssay for \code{i} (and \code{value}) being name of the sample or assay.
#'    If both \code{i} and \code{j} are given, vrLayer with layer name \code{j} (and \code{value}) of vrSample with same name \code{i}.
#'  }
#' }
#' @param j Depends on the usage, see \code{i}.
#' @param ... Arguments passed to other methods
#'
#' @name VoltRon-methods
#' @rdname VoltRon-methods
#'
#' @concept voltron
#'
NULL

## $ method ####

#' @describeIn VoltRon-methods Metadata access for \code{VoltRon} objects
#' 
#' @export
#' @method $ VoltRon
#'
"$.VoltRon" <- function(x, i, ...) {

  # get assay names
  assay_names <- vrAssayNames(x)

  # metadata
  metadata <- Metadata(x, assay = assay_names)

  # get metadata column
  # return(metadata[[i]])
  # return(metadata[,i, drop = TRUE])
  return(as.vector(metadata[[i]]))
}

#' @describeIn VoltRon-methods Metadata overwrite for \code{VoltRon} objects
#'
#' @export
#' @method $<- VoltRon
#'
"$<-.VoltRon" <- function(x, i, value) {

  # sample metadata
  sample.metadata <- SampleMetadata(x)

  # get assay names
  assay_names <- vrAssayNames(x)

  # metadata
  metadata <- Metadata(x, assay = assay_names)

  # dont change Assays or Layers
  if(i %in% c("Assay", "Layer")){
    stop("Changing names of assay types or layers aren't allowed!")
  }

  # change/insert either sample names of metadata columns of main assays
  if(i == "Sample"){
    if(!any(length(value) %in% c(1,nrow(sample.metadata)))){
      stop("New sample names should of length 1 or the same number of assays!")
    } else {
      sample.metadata[[i]] <- value
      x <- changeSampleNames(x, samples = value)
    }
  } else {
    if(length(value) == 1 | nrow(metadata) == length(value)){
      metadata[[i]] <- value
      Metadata(x, assay = assay_names) <- metadata
    } else {
      stop("The new or the existing column should of length 1 or the same as the number of rows")
    }
  }

  return(x)
}

#' @describeIn VoltRon-methods Autocompletion for \code{$} access for \code{VoltRon} objects
#'
#' @inheritParams utils::.DollarNames
#'
#' @importFrom utils .DollarNames
#' @method .DollarNames VoltRon
#'
".DollarNames.VoltRon" <- function(x, pattern = '') {
  meta.data <- as.list(x = Metadata(x))
  return(.DollarNames(x = meta.data, pattern = pattern))
}

### subset of samples and layers ####

#' @describeIn VoltRon-methods Accessing vrAssay or vrSample objects from \code{VoltRon} objects
#'
#' @aliases [[,VoltRon-methods
#' @docType methods
#'
#' @export
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
        assay_list <- x@samples[[cur_assay$Sample]]@layer[[cur_assay$Layer]]@assay
        assay_names <- sapply(assay_list, vrAssayNames)
        return(assay_list[[which(assay_names == rownames(cur_assay))]])
      } else {
        stop("There are no samples or assays named ", i, " in this object")
      }

    } else {
      return(x@samples[[i]])
    }
  }
)

#' @describeIn VoltRon-methods Overwriting vrAssay or vrSample objects from \code{VoltRon} objects
#'
#' @aliases [[<-,VoltRon-methods
#' @docType methods
#'
#' @return \code{[[<-}: \code{x} with the metadata or associated objects added
#' as \code{i}; if \code{value} is \code{NULL}, removes metadata or associated
#' object \code{i} from object \code{x}
#'
#' @export
#'
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
      if(!inherits(value, "vrSample") & !inherits(value, "vrBlock")  ) {
        stop("The provided object is not of class vrSample")
      }
      x@samples[[i]] <- value
    }
    return(x)
  }
)

#' @describeIn VoltRon-methods Accessing vrLayer objects from \code{VoltRon} objects
#'
#' @aliases [[,VoltRon-methods
#' @docType methods
#'
#' @export
#'
setMethod(
  f = '[[',
  signature = c('VoltRon', "character", "character"),
  definition = function(x, i, j, ...){
    return(x[[i]]@layer[[j]])
  }
)

#' @describeIn VoltRon-methods Overwriting vrLayer objects from \code{VoltRon} objects
#'
#' @aliases [[<-,VoltRon-methods
#' @docType methods
#'
#' @return \code{[[<-}: \code{x} with the metadata or associated objects added
#' as \code{i}; if \code{value} is \code{NULL}, removes metadata or associated
#' object \code{i} from object \code{x}
#'
#' @export
#'
setMethod(
  f = '[[<-',
  signature = c('VoltRon', "character", "character"),
  definition = function(x, i, j, ..., value){

    if(!inherits(value, "vrLayer")){
      stop("The provided object is not of class vrLayer")
    }

    x[[i]]@layer[[j]] <- value
    return(x)
  }
)

### Create VoltRon object ####

#' formVoltRon
#'
#' Create a VoltRon object
#'
#' @param data the feature matrix of spatialpoints
#' @param metadata a metadata object of class \link{vrMetadata}
#' @param image a singelton or list of images as magick-image objects
#' @param coords the coordinates of the spatial points
#' @param segments the list of segments each associated with a spatial point
#' @param sample.metadata a data frame of the sample metadata, see \link{SampleMetadata}
#' @param main.assay the name of the main assay
#' @param assay.type the type of the assay (tile, molecule, cell, spot or ROI)
#' @param params additional parameters
#' @param sample_name the name of the sample
#' @param layer_name the name of the layer
#' @param image_name the name/key of the image
#' @param feature_name the name/key of the feature set
#' @param assay_version the VoltRon object version
#' @param project project name
#' @param ... additional parameters passed to \link{formAssay}
#'
#' @importFrom igraph make_empty_graph V vertices
#' @importFrom methods new
#' @importFrom data.table data.table
#' @importFrom rlang %||%
#' @importFrom ids random_id
#' @importFrom Matrix colSums
#'
#' @export
#'
formVoltRon <- function(data = NULL, 
                        metadata = NULL, 
                        image = NULL,
                        coords,
                        segments = list(),
                        sample.metadata = NULL,
                        main.assay = NULL, 
                        assay.type = "cell", 
                        params = list(),
                        sample_name = NULL, 
                        layer_name = NULL, 
                        image_name = NULL, 
                        feature_name = NULL,
                        project = NULL, 
                        version = "v2", ...){

  # set project name
  if(is.null(project))
    project <- "VoltRon"

  # check VoltRon object version
  if(!version %in% c("v1", "v2")){
    stop("'version' has to be set to either 'v1' or 'v2'")
  }
  
  # layer and sample names
  if(is.null(main.assay))
    main.assay <- paste0("Custom_", assay.type)
  layer_name <- ifelse(is.null(layer_name), "Section1", layer_name)
  if(main.assay == layer_name)
      stop(paste0("'", layer_name, "' cannot be a layer name, since main assay is named '", main.assay, "'."))
  sample_name <- ifelse(is.null(sample_name), "Sample1", sample_name)
  if(main.assay == sample_name)
    stop(paste0("'", sample_name, "' cannot be a sample name, since main assay is named '", main.assay, "'."))
  image_name <- ifelse(is.null(image_name), "image_1", image_name)

  # entity IDs from either the data or metadata
  if(!is.null(data)){

    # check for colnames of the raw data
    if(is.null(colnames(data))){
      entityID_nopostfix <- paste0(assay.type,1:ncol(data))
    } else {
      entityID_nopostfix <- colnames(data)
    }

  } else{

    # make empty data if data is missing
    data <- matrix(nrow = 0, ncol = nrow(metadata))

    # check for metadata
    if(!is.null(metadata)) {

      # check row names if exists
      if(is.null(rownames(metadata)) && is.null(metadata$id)){
        entityID_nopostfix <- paste0(assay.type,1:nrow(metadata))
        rownames(metadata) <- entityID
      } else {
        entityID_nopostfix <- metadata$id %||% rownames(metadata)
      }
    } else {
      stop("Either data or metadata has to be provided to build a VoltRon object")
    }
  }

  # Metadata
  vr_metadata_list <- setVRMetadata(metadata, 
                              data, 
                              entityID_nopostfix, 
                              main.assay, 
                              assay.type,
                              sample_name, 
                              layer_name, 
                              version)
  vr_metadata <- vr_metadata_list$vr_metadata
  entityID <- vr_metadata_list$entityID
  colnames(data) <- entityID
  
  # Coordinates
  if(!is.null(coords)){
    if(inherits(coords, "data.frame")){
      coords <- as.matrix(coords)
    }
    if(!inherits(coords, "matrix")){
      stop("Coordinates table should either of a matrix or data.frame class!")
    }
    if(ncol(coords) == 2){
      coords <- cbind(coords,0)
    } else if(ncol(coords) == 3){
      rownames(coords) <- entityID
    } else {
      stop("The length of colnames of the coordinates matrix should be either two or three!")
    }
    rownames(coords) <- entityID
    colnames(coords) <- c("x", "y", "z")
  } else {
    stop("There are no coordinate matrix provided!")
  }

  # create vrAssay
  Assay <- formAssay(data = data, 
                     coords = coords, 
                     segments = segments, 
                     image = image, 
                     params = params, 
                     type = assay.type, 
                     name = "Assay1", 
                     main_image = image_name, 
                     main_featureset = feature_name, 
                     ...)
  listofAssays <- list(Assay)
  names(listofAssays) <- main.assay

  # create layers
  listofLayers <- list(methods::new("vrLayer",
                                    assay = listofAssays,
                                    connectivity = igraph::make_empty_graph(directed = FALSE) + igraph::vertices(entityID)))
  names(listofLayers) <- layer_name
  
  # create samples
  listofSamples <- list(methods::new("vrBlock", 
                                     layer = listofLayers, 
                                     zlocation = c("Section1" = 0),
                                     adjacency = matrix(0, nrow = 1, ncol = 1,
                                                        dimnames = list("Section1", "Section1"))))
                        

  names(listofSamples) <- sample_name

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setVRSampleMetadata(listofSamples)
  }

  # set VoltRon class
  methods::new("VoltRon", samples = listofSamples, metadata = vr_metadata, sample.metadata = sample.metadata, main.assay = main.assay, project = project)
}

### Assay Methods ####

#' Main Assay
#'
#' Get and set the main assay of a VoltRon object
#'
#' @param object a VoltRon object
#' @rdname vrMainAssay
#'
#' @export
#'
vrMainAssay <- function(object) {
  object@main.assay
}

#' @param value new assay name
#'
#' @rdname vrMainAssay
#'
#' @export
"vrMainAssay<-" <- function(object, value) {
  sample.metadata <- SampleMetadata(object)
  assay_names <- unique(sample.metadata$Assay)
  if(!value %in% assay_names){
    stop("There is no assay names '", value, "' in this object")
  } else {
    object@main.assay <- value
  }
  return(object)
}


#' @rdname addAssay
#' @method addAssay VoltRon
#'
#' @importFrom igraph make_empty_graph add_edges vertices
#'
#' @export
#' @noRd
addAssay.VoltRon <- function(object, assay, metadata = NULL, assay_name, sample = "Sample1", layer = "Section1"){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay id
  assay_ids <- as.numeric(gsub("Assay", "", rownames(sample.metadata)))
  assay_id <- paste0("Assay", max(assay_ids)+1)
  assay_names <- c(rownames(sample.metadata), assay_id)

  # update sample.metadata and metadata
  object@sample.metadata <- rbind(sample.metadata, c(assay_name, layer, sample))
  rownames(object@sample.metadata) <- assay_names
  object@metadata <- addAssay(object@metadata, metadata = metadata,
                              assay = assay, assay_name = assay_name,
                              sample = sample, layer = layer)

  # get sample and layer
  curlayer <- object[[sample, layer]]
  assay_list <- curlayer@assay

  # change assay name and add to the layer
  vrAssayNames(assay) <- assay_id
  new_assay_list <- list(assay)
  names(new_assay_list) <- assay_name
  assay_list <- c(assay_list, new_assay_list)
  object[[sample, layer]]@assay <- assay_list

  # add connectivities of assay to the layer
  catch_connect <- try(slot(curlayer, name = "connectivity"), silent = TRUE)
  if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
    g_assay <- igraph::make_empty_graph(directed = FALSE) + igraph::vertices(vrSpatialPoints(object, assay = assay_id))
    g_layer <- curlayer@connectivity + g_assay
    object[[sample, layer]]@connectivity <- g_layer 
  }

  # return
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrAssayNames
#' @order 2
#' @export
#' 
vrAssayNames.VoltRon <- function(object, assay = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check assays
  if(is.null(assay))
    assay <- vrMainAssay(object)

  # get assay names
  if(any(assay == "all")){
    assay_names <- rownames(sample.metadata)
  } else {
    if(all(assay %in% sample.metadata$Assay)){
      assay_names <- rownames(sample.metadata)[sample.metadata$Assay %in% assay]
    } else {
      if(all(assay %in% rownames(sample.metadata))) {
        assay_names <- assay
      } else {
        stop("Assay name or type is not found in the object")
      }
    }
  }

  return(assay_names)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' 
#' @rdname vrAssayTypes
#' @order 2
#'
#' @export
vrAssayTypes.VoltRon <- function(object, assay = NULL){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assay types
  assay_types <- sapply(assay_names, function(x) vrAssayTypes(object[[x]]))

  return(assay_types)
}


#' changeSampleNames.VoltRon
#'
#' Change the sample names of the VoltRon object and reorient layers if needed
#'
#' @param samples a single or a set of sample names
#' 
#' @rdname changeSampleNames
#'
#' @importFrom dplyr n_distinct %>% distinct select mutate group_by
#' @importFrom methods new
#'
#' @noRd
changeSampleNames.VoltRon <- function(object, samples = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # old to new samples table
  samples_table <- data.frame(sample.metadata, AssayID = rownames(sample.metadata), NewSample = samples)

  # check if multiple new sample names are associated with the same section of one sample
  check_samples_table <- samples_table %>%
    dplyr::group_by(Assay, Sample) %>% dplyr::mutate(n = dplyr::n_distinct(NewSample)) %>%
    select(c("Assay", "Sample", "n")) %>% distinct()
  if(any(check_samples_table$n > 1)){
    message("Overwriting the sample names of assays that were original from a single layer of a sample arent allowed")
    message("Check Sample Metadata for the correct Sample reassignment")
    print(sample.metadata)
    stop()
  }

  # assign new sample names to samples and sample metadata
  new_sample.metadata <- NULL
  new_listofSamples <- list()
  for(cur_sample in unique(samples)){

    # current sample and sample table
    cur_sample.metadata <- samples_table[samples_table$NewSample == cur_sample,]

    # for each unique sample names, combine layers and multiple samples into one
    listofLayers <- NULL
    uniq_old_samples <- unique(cur_sample.metadata$Sample)
    for(i in 1:length(uniq_old_samples)){
      listofLayers <- c(listofLayers, object[[uniq_old_samples[i]]]@layer)
    }
    cur_sample.metadata$comb <- paste(cur_sample.metadata$Sample, cur_sample.metadata$Layer, sep = "_")
    cur_sample.metadata$NewLayer <- paste0("Section", as.numeric(factor(cur_sample.metadata$comb, levels = unique(cur_sample.metadata$comb))))
    # names(listofLayers) <- cur_sample.metadata$NewLayer
    names(listofLayers) <- unique(cur_sample.metadata$NewLayer) ## CHANGE THIS LATER IF NEEDED ####
    
    # make layer adjacency and get distance
    adjacency <- matrix(0, nrow = length(listofLayers), ncol = length(listofLayers),
                  dimnames = list(names(listofLayers), names(listofLayers)))
    diag(adjacency) <- 1
    # distance <- matrix(NA, nrow = length(listofLayers), ncol = length(listofLayers),
    #                    dimnames = list(names(listofLayers), names(listofLayers)))
    # diag(distance) <- 0
    zlocation <- rep(0,length(listofLayers))
    names(zlocation) <- names(listofLayers)
    
    # make new block
    # listofSamples <- list(methods::new("vrBlock", 
    #                                    layer = listofLayers, adjacency = adjacency, distance = distance))
    listofSamples <- list(methods::new("vrBlock",
                                       layer = listofLayers, zlocation = zlocation, adjacency = adjacency))
    names(listofSamples) <- cur_sample
    new_listofSamples <- c(new_listofSamples, listofSamples)
    new_sample.metadata <- rbind(new_sample.metadata, cur_sample.metadata)
  }

  # assign new samples and layers to metadata
  metadata <- changeSampleNames(Metadata(object, type = "all"), sample_metadata_table = new_sample.metadata)

  # sample metadata
  new_sample.metadata <- new_sample.metadata[,c("Assay", "NewLayer", "NewSample")]
  colnames(new_sample.metadata) <- c("Assay", "Layer", "Sample")

  # reinsert object elements
  object@sample.metadata <- new_sample.metadata
  object@samples <- new_listofSamples
  object@metadata <- metadata

  # return
  return(object)
}

#' changeAssayNames.VoltRon
#'
#' Change the sample names of the VoltRon object and reorient layers if needed
#'
#' @rdname changeAssayNames
#' @method changeAssayNames VoltRon
#'
#' @param object a VoltRon object
#' @param assays a set of assay names
#'
#' @noRd
changeAssayNames.VoltRon <- function(object, assays = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the length of the new assay names
  if(nrow(sample.metadata) != length(assays))
    stop("The set of new assay names should be of the number of assays in the VoltRon object.")

  # check the uniqueness of the assay names
  if(length(unique(assays)) != length(assays))
    stop("Each new assay name should be unique")

  # attach new names of sample.metadata
  sample.metadata$NewAssayNames <- assays

  # change assay names in layers
  samples <- unique(sample.metadata$Sample)
  for(samp in samples){
    object[[samp]] <- changeAssayNames(object[[samp]], sample.metadata = sample.metadata[sample.metadata$Sample == samp,])
  }

  # return
  return(object)
}

#' addLayerConnectivity
#'
#' add connectivity information to the assays (vrAssay) of the same layer (vrLayer)
#'
#' @param object a VoltRon object
#' @param connectivity a metadata of edges representing connected spatial points across assays
#' @param sample sample name
#' @param layer layer name
#'
#' @importFrom igraph add_edges
#'
#' @noRd
addLayerConnectivity <- function(object, connectivity, sample, layer){

  # get sample and layer
  curlayer <- object[[sample, layer]]

  # make edges from connectivity matrix
  connectivity <- as.vector(t(as.matrix(connectivity)))

  # add edges
  object[[sample, layer]]@connectivity <- igraph::add_edges(curlayer@connectivity, edges = connectivity)

  # return
  return(object)
}

### Layer Methods ####

#' addBlockConnectivity
#'
#' add connectivity information to the layers (vrLayer) of the same block (Block)
#'
#' @param object a VoltRon object
#' @param connectivity a metadata of edges representing connected layers within a block
#' @param zlocation 
#' @param sample sample name
#'
#' @noRd
addBlockConnectivity <- function(object, connectivity, zlocation = NULL, sample){
  
  # get sample and layer
  cursample <- object[[sample]]
  
  # update z location/coordinates
  if(!is.null(zlocation)){
    cursample@zlocation[names(cursample@zlocation)] <- zlocation
  }
  
  # update adjacency
  adjacency <- cursample@adjacency
  for(i in 1:nrow(connectivity)){
    adjacency[connectivity[i,1], connectivity[i,2]] <- 
      adjacency[connectivity[i,2], connectivity[i,1]] <- 1
  }
  cursample@adjacency <- adjacency
  
  # return sample
  object[[sample]] <- cursample
  
  # return
  return(object)
}

#' getBlockConnectivity
#'
#' get connected assays
#'
#' @param object a VoltRon object
#' @param connectivity a metadata of edges representing connected layers within a block
#' @param zlocation 
#' @param sample sample name
#'
#' @importFrom igraph components graph_from_adjacency_matrix
#' 
#' @noRd
getBlockConnectivity <- function(object, assay){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get samples
  sample_metadata <- SampleMetadata(object)
  samples <- unique(sample_metadata[assay_names, "Sample"])
  
  # get list of connected assays
  assay_list <- list()
  for(samp in samples){
    cur_sample_metadata <- sample_metadata[sample_metadata$Sample == samp,]
    cur_assaynames <- assay_names[assay_names %in% rownames(cur_sample_metadata)]
    cur_sections <- cur_sample_metadata[cur_assaynames, "Layer"]
    
    catch_connect <- try(slot(object[[samp]], name = "adjacency"), silent = TRUE)
    if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
      adjacency <- object[[samp]]@adjacency
      adjacency <- adjacency[match(cur_sections,rownames(adjacency)), match(cur_sections,rownames(adjacency))]
      colnames(adjacency) <- rownames(adjacency) <- cur_assaynames
      components <- igraph::components(igraph::graph_from_adjacency_matrix(adjacency))
      assay_list <- c(assay_list, split(names(components$membership), components$membership))
    } else {
      assay_list <- c(assay_list, cur_assaynames)
    }
  }
  
  # return list
  assay_list
}

### Object Methods ####

#' Subsetting VoltRon objects
#'
#' Given a VoltRon object, subset the object given one of the attributes
#'
#' @param object a vrAssay object
#' @param subset Logical statement for subsetting
#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#' @param image the subseting string passed to \link{image_crop}
#' @param interactive TRUE if interactive subsetting on the image is demanded
#' @param use.points.only if \code{interactive} is \code{TRUE}, use spatial points instead of the reference image
#' @param shiny.options a list of shiny options (launch.browser, host, port etc.) passed \code{options} arguement of \link{shinyApp}. For more information, see \link{runApp}
#' 
#' @rdname subset
#' @aliases subset
#' @method subset VoltRon
#'
#' @importFrom rlang enquo eval_tidy quo_get_expr quo_text
#' @importFrom stringr str_extract
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#' # example data
#' data("visium_data")
#' 
#' # subset based on assay
#' subset(visium_data, assays = "Assay1")
#' subset(visium_data, assays = "Visium")
#' 
#' # subset based on samples
#' subset(visium_data, samples = "Anterior1")
#' 
#' # subset based on assay
#' subset(visium_data, spatialpoints = c("GTTATATTATCTCCCT-1_Assay1", "GTTTGGGTTTCGCCCG-1_Assay1"))
#' 
#' # subset based on features
#' subset(visium_data, features = c("Map3k19", "Rab3gap1"))
#' 
#' # interactive subsetting
#' \dontrun{
#' visium_subset_data <- subset(visium_data, interactive = TRUE)
#' visium_subset <- visium_subset_data$subsets[[1]]
#' }
subset.VoltRon <- function(object, subset, samples = NULL, assays = NULL, spatialpoints = NULL, features = NULL, image = NULL, interactive = FALSE, use.points.only = FALSE, 
                           shiny.options = list(launch.browser = getOption("shiny.launch.browser", interactive()))) {

  # subseting based on subset argument
  if (!missing(x = subset)) {
    # subset_data <- subset
    subset <- rlang::enquo(arg = subset)
  }
  if(!missing(subset)){
    metadata <- Metadata(object)
    name <- strsplit(rlang::quo_text(subset), split = " ")[[1]][1]
    if(name %in% colnames(metadata)){
      if(inherits(metadata, "data.table")){
        spatialpoints <- metadata$id[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
      } else if(inherits(metadata, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame"))){
        stop("Direct subsetting for Ondisk VoltRon objects are currently not possible!")
        # spatialpoints <- as.vector(metadata$id)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
      } else {
        if(!is.null(rownames(metadata))){
          cur_data <- rownames(metadata)
        } else {
          cur_data <- metadata$id
        }
        spatialpoints <- rownames(metadata)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
      }
    } else {
      stop("Column '", name, "' is not found in the metadata")
    }
    object <- subset(object, spatialpoints = spatialpoints)
    return(object)
  }

  # subseting on other attributes
  attrinfo <- c(sapply(list(samples, assays, spatialpoints, features), function(x) length(x) > 0), interactive)
  if(sum(attrinfo) > 1){
    stop("Please choose only one of the subsetting attributes: 'samples', 'assays', 'spatialpoints', 'features' or 'interactive'")
  }

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # subsetting
  if(!is.null(samples)){

    # check assays associated with samples and subset for assays
    if(all(samples %in% sample.metadata$Sample)){
      assays <- rownames(sample.metadata)[sample.metadata$Sample %in% samples]
      return(subset.VoltRon(object, assays = assays))
    } else {
      stop("Some requested samples are not found in this VoltRon object!")
    }

  } else if(!is.null(assays)){

    # subset for assays
    sample.metadata <- subset_sampleMetadata(sample.metadata, assays = assays)
    metadata <- subset.vrMetadata(Metadata(object, type = "all"), assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, assays = assays)
    }, USE.NAMES = TRUE)

  } else if(!is.null(spatialpoints)) {

    # subsetting on entity names
    metadata <- subset.vrMetadata(Metadata(object, type = "all"), spatialpoints = spatialpoints)
    samples <- vrSampleNames(metadata)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE)
    spatialpoints <-  do.call("c", lapply(listofSamples, vrSpatialPoints.vrSample))
    metadata <- subset.vrMetadata(Metadata(object, type = "all"), spatialpoints = spatialpoints)
    sample.metadata <- subset_sampleMetadata(sample.metadata, assays = vrAssayNames.vrMetadata(metadata))

  } else if(!is.null(features)){
    
    # subsetting on features
    assay_names <- vrAssayNames(object)
    for(assy in assay_names){
      if(inherits(object[[assy]], "vrAssay")){
        object[[assy]] <- subset.vrAssay(object[[assy]], features = features) 
      } else {
        object[[assy]] <- subset.vrAssayV2(object[[assy]], features = features)
      }
    } 
    metadata <- Metadata(object, type = "all")
    listofSamples <- object@samples

  } else if(!is.null(image)) {

    # subsetting on image
    if(inherits(image, "character")){

      # check if there are only one image and one assay
      numlayers <- paste0(sample.metadata$Layer, sample.metadata$Sample)
      if(length(unique(numlayers)) > 1){
        stop("Subseting on images can only be performed on VoltRon objects with a single layer")
      } else {
        samples <- unique(sample.metadata$Sample)
        listofSamples <- sapply(object@samples[samples], function(samp) {
          subset.vrSample(samp, image = image)
        }, USE.NAMES = TRUE)
        spatialpoints <-  do.call(c, lapply(listofSamples, vrSpatialPoints.vrSample))
        metadata <- subset.vrMetadata(Metadata(object, type = "all"), spatialpoints = spatialpoints)
      }
    } else {
      stop("Please provide a character based subsetting notation, see magick::image_crop documentation")
    }
  } else if(interactive){
    
    # interactive subsetting
    results <- demuxVoltRon(object, use.points.only = use.points.only, shiny.options = shiny.options)
    return(results)
  }

  # main.assay
  main.assay <- unique(sample.metadata$Assay)[unique(sample.metadata$Assay) == names(table(sample.metadata$Assay))[which.max(table(sample.metadata$Assay))]]

  # project
  project <- object@project

  # subset graphs
  graph_list <- subset_graphs(object, 
                              spatialpoints = vrSpatialPoints(metadata, assay = vrAssayNames(object)))

  # set VoltRon class
  methods::new("VoltRon",
               samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata,
               graph = graph_list, main.assay = main.assay, project = project)
}

#' Merging VoltRon objects
#'
#' Given a VoltRon object, and a list of VoltRon objects, merge all.
#'
#' @param object a VoltRon Object
#' @param object_list a list of VoltRon objects
#' @param samples a single sample name or multiple sample names of the same size as the given VoltRon objects
#' @param main.assay the name of the main assay
#'
#' @rdname merge
#' @aliases merge
#' @method merge VoltRon
#' @importFrom methods new
#'
#' @export
merge.VoltRon <- function(object, object_list, samples = NULL, main.assay = NULL) {

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)

  # check if all are VoltRon
  if(!all(lapply(object_list, class) == "VoltRon"))
    stop("All arguements have to be of VoltRon class")

  # sample metadata list
  sample.metadata_list <- lapply(object_list, function(x) slot(x, name = "sample.metadata"))

  # old assay names
  old_assay_names <- do.call(c, lapply(sample.metadata_list, rownames))

  # merge sample metadata
  sample.metadata <- merge_sampleMetadata(sample.metadata_list)

  # merge metadata and sample metadata
  message("Merging metadata ...")
  metadata_list <- lapply(object_list, function(x) slot(x, name = "metadata"))
  metadata <- merge(metadata_list[[1]], metadata_list[-1])

  # combine samples and rename layers
  message("Merging blocks and layers ...")
  listofSamples <- NULL
  for(i in 1:length(object_list)){
    cur_object <- object_list[[i]]@samples
    listofSamples <- c(listofSamples, cur_object)
  }

  # get main assay
  if(is.null(main.assay))
      main.assay <- names(sort(table(sample.metadata$Assay), decreasing = TRUE))[1]

  # project
  project <- slot(object_list[[1]], "project")

  # set VoltRon class
  object <- methods::new("VoltRon", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata, main.assay = main.assay, project = project)

  # change assay names and sample names
  object <- changeAssayNames(object, assays = rownames(sample.metadata))

  # change sample names
  if(!is.null(samples))
    object$Sample <- samples

  # return
  object
}

#' @rdname vrSpatialPoints
#' @order 2
#' 
#' @export
vrSpatialPoints.VoltRon <- function(object, assay = NULL) {

  # get assays
  assay <- vrAssayNames(object, assay = assay)

  # return
  return(vrSpatialPoints(object@metadata, assay = assay))
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrFeatures
#' @order 2
#' @export
vrFeatures.VoltRon <- function(object, assay = NULL) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all features
  features <- NULL
  for(assy in assay_names)
    features <- c(features, vrFeatures(object[[assy]]))

  return(unique(features))
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrFeatureData
#' @order 2
#' @export
vrFeatureData.VoltRon <- function(object, assay = NULL, feat_type = NULL) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all features
  features <- vrFeatureData(object[[assay_names[1]]], feat_type = feat_type)

  # return
  return(features)
}

#' @param value new feature metadata
#' 
#' @rdname vrFeatureData
#' @order 4
#' @export
"vrFeatureData<-.VoltRon" <- function(object, assay = NULL, value) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names)
    vrFeatureData(object[[assy]]) <- value

  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param features the set of features
#' @param feat_type the feature set type 
#' @param norm TRUE if normalized data should be returned
#' @param ... additional parameters passed to other methods and \link{vrImages}
#'
#' @rdname vrData
#' @order 2
#' 
#' @importFrom dplyr full_join mutate_all coalesce
#'
#' @export
vrData.VoltRon <- function(object, assay = NULL, features = NULL, feat_type = NULL, norm = FALSE, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  data <- NULL
  for(i in 1:length(assay_names)){
    cur_data <- vrData(object[[assay_names[i]]], features = features, feat_type = feat_type, norm = norm, ...)
    if(inherits(cur_data, c("data.frame", "Matrix", "matrix"))){
      cur_data <- data.frame(cur_data, feature.ID = rownames(cur_data), check.names = FALSE) 
    } 
    if(i == 1){
      data <- cur_data
    } else {
      data <- merge_data(data, cur_data, by = "feature.ID")
    }
  }
  if("feature.ID" %in% colnames(data)){
    rownames(data) <- data$feature.ID
    data <- data[,!colnames(data) %in% "feature.ID"] 
    data <- as.matrix(data)
    data <- replaceNaMatrix(data, 0)
    colnames(data) <- gsub("\\.","-", colnames(data))
  }

  return(data)
}

#' @importFrom Matrix Matrix
merge_data <- function(data1, data2, by = "feature.ID"){
  if(inherits(data1, c("data.frame", "Matrix"))){
    
    # merge
    data1 <- dplyr::full_join(data1, data2, by = "feature.ID")
    
  } else if(inherits(data1, c("IterableMatrix"))) {
    rownames_all <- unique(c(rownames(data1), rownames(data2)))
    
    # first data
    m <- Matrix::Matrix(nrow = length(rownames_all) - length(rownames(data1)), ncol = ncol(data1), data = 0, sparse = TRUE)
    data1_new <- rbind(data1, m)
    rownames(data1_new) <- c(rownames(data1), setdiff(rownames_all, rownames(data1)))
    data1_new <- data1_new[rownames_all,]
    
    # second data
    m <- Matrix::Matrix(nrow = length(rownames_all) - length(rownames(data2)), ncol = ncol(data2), data = 0, sparse = TRUE)
    data2_new <- rbind(data2, m)
    rownames(data2_new) <- c(rownames(data2), setdiff(rownames_all, rownames(data2)))
    data2_new <- data2_new[rownames_all,]
   
    # merge 
    data1 <- cbind(data1_new, data2_new)
  }
  return(data1)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param ... additional parameters passed to \link{generateTileData.vrAssay}
#' 
#' @rdname generateTileData
#' @order 2
#'
#' @export
generateTileData.VoltRon <- function(object, assay = NULL, ...) {
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # check if assay types are all tiles
  assay_types <- vrAssayTypes(object, assay = assay)
  if(!all(assay_types == "tile"))
    stop("generateTileData can only be used for tile-based assays")
  
  # get tile data for all assays
  for(assy in assay_names)
    object[[assy]] <- generateTileData(object[[assy]], ...)
}

#' @rdname vrEmbeddings
#' @order 2
#'
#' @export
vrEmbeddings.VoltRon <- function(object, assay = NULL, type = "pca", dims = 1:30) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  returndata_list <- list()
  for(i in 1:length(assay_names))
    returndata_list[[i]] <- vrEmbeddings(object[[assay_names[i]]], type = type, dims = dims)

  return(do.call(rbind, returndata_list))
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param type the key name for the embedding
#' @param overwrite Whether the existing embedding with name 'type' should be overwritten
#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @order 4
#'
#' @export
#'
"vrEmbeddings<-.VoltRon" <- function(object, assay = NULL, type = "pca", overwrite = FALSE, value) {

  # check if the embedding exists
  if(type %in% vrEmbeddingNames(object) && !overwrite)
    stop("An embedding named '", type, "' already exists in this object. Do overwrite = TRUE for replacing with the existing one.")

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names){
    assayobject <- object[[assy]]
    if(vrAssayTypes(assayobject) %in% c("ROI", "cell", "spot")){
      vrEmbeddings(assayobject, type = type) <- value[grepl(paste0(assy, "$"), rownames(value)),, drop = FALSE]
    } else {
      vrEmbeddings(assayobject, type = type) <- value[vrSpatialPoints(assayobject),, drop = FALSE]
    }
    object[[assy]] <- assayobject
  }

  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrEmbeddingNames
#' @order 2
#'
#' @export
vrEmbeddingNames.VoltRon <- function(object, assay = NULL){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assay types
  embed_names <- unique(unlist(lapply(assay_names, function(x) vrEmbeddingNames(object[[x]]))))

  return(embed_names)
}

#### Feature ####

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' If NULL, the default assay will be used, see \link{vrMainAssay}. If given as "all", then provides a summary of spatial systems across all assays.
#' @param data new data matrix for new feature set
#' @param feature_name the name of the new feature set
#' 
#' @rdname addFeature
#' @method addFeature VoltRon
#' 
#' @importFrom stringr str_replace
#' 
#' @export
addFeature.VoltRon <- function(object, assay = NULL, data, feature_name){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  if(length(assay_names) > 1){
    stop("You cannot add new features to multiple assays at once!")
  }
  
  # add assay
  object[[assay_names]] <- addFeature(object[[assay_names]], data = data, feature_name = feature_name)
  
  # return
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' If NULL, the default assay will be used, see \link{vrMainAssay}. If given as "all", then provides a summary of spatial systems across all assays.
#'
#' @rdname vrMainFeatureType
#' @order 2
#' @export
vrMainFeatureType.VoltRon <- function(object, assay = NULL){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if assay = all, give a summary
  if(!is.null(assay)){
    if(assay == "all"){
      featuretype_names <- unlist(lapply(rownames(SampleMetadata(object)), function(x) paste(vrMainFeatureType(object[[x]]), collapse = ",")))
      featuretype_names <- data.frame(Assay = assay_names, Feature = featuretype_names)
      return(featuretype_names)
    }
  }
  
  # get assay types
  featuretype_names <- unlist(lapply(assay_names, function(x) vrMainFeatureType(object[[x]])))
  
  # return data
  if(!is.null(featuretype_names)){
    featuretype_data <- data.frame(Assay = assay_names, Feature = featuretype_names)
    return(featuretype_data)
  } else {
    return(NULL)
  }
}

#' @param value the name of main feature set
#'
#' @rdname vrMainFeatureType
#' @order 4
#' @export
"vrMainFeatureType<-.VoltRon" <- function(object, assay = NULL, value){
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # assays 
  assay_names <- vrAssayNames(object, assay = assay)
  unique_assays <- unique(sample_metadata[assay_names, "Assay"])
  if(length(unique_assays) > 1){
    stop("You can only set the main feature type of a single assay type")
  } else {
    for(assy in assay_names){
      vrMainFeatureType(object[[assy]]) <- value
    }
  }
  
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' If NULL, the default assay will be used, see \link{vrMainAssay}. If given as "all", then provides a summary of spatial systems across all assays
#'
#' @rdname vrFeatureTypeNames
#'
#' @export
vrFeatureTypeNames.VoltRon <- function(object, assay = NULL){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # if assay = all, give a summary
  if(!is.null(assay)){
    if(assay == "all"){
      feature_names <- unlist(lapply(assay_names, function(x) paste(vrFeatureTypeNames(object[[x]]), collapse = ",")))
      feature_names <- data.frame(Assay = assay_names, Feature = feature_names)
      return(feature_names)
    }
  }
  
  feature_names <- unique(unlist(lapply(assay_names, function(x) vrFeatureTypeNames(object[[x]]))))
  
  return(feature_names)
}

#### Metadata ####

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param type the assay type: ROI, spot or cell, or all for the entire metadata object
#' 
#' @rdname Metadata
#'
#' @importFrom methods slotNames
#' @export
Metadata.VoltRon <- function(object, assay = NULL, type = NULL) {

  # check type
  if(!is.null(type)){
    
    if(type == "all"){
      return(object@metadata)
    } else {
      if(!is.null(assay)){
        stop("Please specify either assay or type, not both!")
      }
      if(type %in% methods::slotNames(object@metadata)){
        return(slot(object@metadata, name = type))
      }
    }
  } else{
    type <- unique(vrAssayTypes(object, assay = assay))
    if(length(type) > 1)
      stop("Select only metadata with a single assay type!")
  }

  # get assay metadata from matching type
  if(type %in% methods::slotNames(object@metadata)){

    # sample metadata
    sample.metadata <- SampleMetadata(object)

    # get assay names
    assay_names <- vrAssayNames(object, assay = assay)

    # get metadata
    metadata <- slot(object@metadata, name = type)
    if(inherits(metadata, "data.table")){
      metadata <- subset(metadata, assay_id %in% assay_names)
    } else if(inherits(metadata, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame"))){
      if("assay_id" %in% colnames(metadata)){
        metadata_list <- list()
        for(assy in assay_names){
          metadata_list[[assy]] <- metadata[metadata$assay_id == assy,]
        }
        metadata <- do.call("rbind", metadata_list)
      } else {
        ind <- stringr::str_extract(as.vector(metadata$id), "Assay[0-9]+") %in% assay_names
        metadata <- metadata[ind,]
      }
    } else {
      metadata <- metadata[stringr::str_extract(rownames(metadata), "Assay[0-9]+") %in% assay_names, ]
    }
    return(metadata)
  } else {
    stop("Please provide one of five assay types: 'ROI', 'cell', 'spot', 'molecule' or 'tile'.")
  }
}

#' @param value new metadata
#'
#' @rdname Metadata
#' @method Metadata<- VoltRon
#'
#' @export
"Metadata<-.VoltRon" <- function(object, assay = NULL, type = NULL, value) {

  if(!is.data.frame(value) && !inherits(value, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame")))
    stop("The new or updated metadata has to be a data frame")

  `%notin%` <- Negate(`%in%`)
  if(is.null(rownames(value)) && "id" %notin% colnames(value))
    stop("The new metadata should have row names or a column called 'id' to match its rows with the existing one")

  if(is.null(type)){
    type <- unique(vrAssayTypes(object, assay = assay))
  }

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get metadata
  metadata <- slot(object@metadata, name = type)

  # if(type %in% c("ROI", "cell", "spot")){
  if(!is.null(rownames(metadata))){

    # replace the metadata (or some part of it) with the new value
    if(length(setdiff(rownames(value), rownames(metadata))) == 0){

      # check columns of the new table
      new_columns <- setdiff(colnames(value), colnames(metadata))

      # current metadata shouldnt have columns that value doesnt have
      if(length(setdiff(colnames(metadata), colnames(value))) > 0)
        stop("Some columns of new data frame are not available in the metadata")

      # if new columns appear, update the column names of the metadata'
      if(length(new_columns) > 0){
        value <- value[,c(colnames(metadata), new_columns)]
        for(cur_col in new_columns){
          if(is.numeric(value[[cur_col]])){
            metadata[[cur_col]] <- NA
          } else {
            metadata[[cur_col]] <- ""
          }
        }
      }

      # replace data
      metadata[rownames(value), ] <- value
      slot(object@metadata, name = type) <- metadata

    } else {
      stop("Some rows of new data frame are not available in the metadata")
    }

  # } else if(type %in% c("tile", "molecule")){
  } else if("id" %in% colnames(metadata)){

    # replace the metadata (or some part of it) with the new value
    if(length(setdiff(value$id, metadata$id)) == 0){

      # check columns of the new table
      new_columns <- setdiff(colnames(value), colnames(metadata))

      # current metadata shouldnt have columns that value doesnt have
      if(length(setdiff(colnames(metadata), colnames(value))) > 0)
        stop("Some columns of new data frame are not available in the metadata")

      # if new columns appear, update the column names of the metadata'
      if(length(new_columns) > 0){
        value <- value[,colnames(value)[colnames(value) %in% c(colnames(metadata), new_columns)], with = FALSE]
        # value <- value[,c(colnames(metadata), new_columns), with = FALSE]
        for(cur_col in new_columns){
          if(is.numeric(value[[cur_col]])){
            metadata[[cur_col]] <- NA
          } else {
            metadata[[cur_col]] <- ""
          }
        }
      }

      # replace data
      metadata <- value
      slot(object@metadata, name = type) <- metadata

    } else {
      stop("Some rows of new data frame are not available in the metadata")
    }

  } else {
    # stop("Please provide one of three assay types: 'ROI', 'cell', 'spot'.")
    stop("The metadata should either have rownames or a column called 'id'!")
  }
  
  return(object)
}

#' SampleMetadata
#'
#' Get the sample metadata of a VoltRon object
#'
#' @param object a VoltRon object
#'
#' @export
SampleMetadata <- function(object) {
  object@sample.metadata
}

#### Spatial ####

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param image_name (deprecated, use \code{spatial_name}) the name/key of the image associated with the coordinates
#' @param spatial_name the name/key of the spatial system associated with the coordinates
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainImage}) is requested
#'
#' @rdname vrCoordinates
#' @order 2
#' @export
#'
vrCoordinates.VoltRon <- function(object, assay = NULL, image_name = NULL, spatial_name = NULL, reg = FALSE) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get sample metadata
  sample_metadata <- SampleMetadata(object)

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # get all coordinates
  coords <- NULL
  for(assy in assay_names){
    
    # get coordinates
    cur_coords <- vrCoordinates(object[[assy]], image_name = image_name, reg = reg)
    
    # update zlocation
    sample_name <- sample_metadata[assy, "Sample"]
    
    catch_connect <- try(slot(object, name = "zlocation"), silent = TRUE)
    if(!is(catch_connect, 'try-error') && !methods::is(catch_connect,'error')){
      zlocation <- object[[sample_name]]@zlocation 
      cur_coords[,"z"] <- rep(zlocation[sample_metadata[assy, "Layer"]], nrow(cur_coords)) 
    }
    
    # merge coordinates
    if(!is.null(coords)){
      coords <- rbind(coords, cur_coords)
    } else {
      coords <- cur_coords
    }
  }

  # return image
  return(coords)
}

#' @param value new coordinates of spatial points
#' 
#' @rdname vrCoordinates
#' @order 4
#' @export
"vrCoordinates<-.VoltRon" <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the number of assays in the object
  if(nrow(sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays in the same time are not permitted!")

  # get assay
  cur_assay <- sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # change coordinates
  # vrCoordinates(vrassay, image_name = image_name, reg = reg) <- value
  vrCoordinates(vrassay, spatial_name = image_name, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param image_name (deprecated, use \code{spatial_name}) the name/key of the image associated with the coordinates
#' @param spatial_name the name/key of the spatial system associated with the coordinates
#' @param reg TRUE if registered coordinates of the main image (\link{vrMainImage}) is requested
#'
#' @rdname vrSegments
#' @order 2
#' @export
vrSegments.VoltRon <- function(object, assay = NULL, image_name = NULL, spatial_name = NULL, reg = FALSE) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # get all coordinates
  segts <- NULL
  for(assy in assay_names)
    segts <- c(segts, vrSegments(object[[assy]], spatial_name = image_name, reg = reg))
    # segts <- c(segts, vrSegments(object[[assy]], image_name = image_name, reg = reg))

  # return image
  return(segts)
}

#' @param value new segment coordinates of spatial points
#' 
#' @rdname vrSegments
#' @order 5
#' @export
"vrSegments<-.VoltRon" <- function(object, image_name = NULL, spatial_name = NULL, reg = FALSE, value) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the number of assays in the object
  if(nrow(sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays are not permitted!")

  # get assay
  cur_assay <- sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # change coordinates
  # vrSegments(vrassay, image_name = image_name, reg = reg) <- value
  vrSegments(vrassay, spatial_name = image_name, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param image_name (deprecated, use \code{spatial_name}) the name/key of the image
#' @param spatial_name the name/key of the spatial system associated with the coordinates
#' @param ... additional parameters passed to \link{vrCoordinates} and \link{vrSegments}
#' 
#' @rdname flipCoordinates
#' @order 2
#'
#' @export
flipCoordinates.VoltRon <- function(object, assay = NULL, image_name = NULL, spatial_name = NULL, ...){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get spatial name
  if(!is.null(spatial_name)) 
    image_name <- spatial_name
  
  # flip coordinates
  for(assy in assay_names){
    object[[assy]] <- flipCoordinates(object[[assy]], spatial_name = image_name, ...)
    # object[[assy]] <- flipCoordinates(object[[assy]], image_name = image_name, ...)
  }
  return(object)
}

#### Graphs ####

#' vrGraph
#'
#' Get graph of a VoltRon object
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param graph.type the type of the graph, either custom or given by \link{getProfileNeighbors} or \link{getSpatialNeighbors} functions
#'
#' @rdname vrGraph
#'
#' @importFrom igraph induced_subgraph
#' @export
vrGraph <- function(object, assay = NULL, graph.type = "kNN") {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  node_names <- vrSpatialPoints(object, assay = assay_names)

  # check if there exists graphs
  if(length(names(object@graph)) == 0)
    stop("There are no graphs in this VoltRon object!")

  # check graph type
  if(!graph.type %in% names(object@graph))
    stop("The graph name '", graph.type, "' can't be found in this VoltRon object!")

  # return graph
  if(length(vrGraphNames(object)) > 0){
    returngraph <- igraph::induced_subgraph(object@graph[[graph.type]], node_names)
    return(returngraph)
  } else {
    warning("This VoltRon object does not have any graphs yet!")
    return(NULL)
  }
}

#' @param value new graph
#' 
#' @rdname vrGraph
#'
#' @importFrom igraph disjoint_union induced_subgraph
#' @export
"vrGraph<-" <- function(object, graph.type = "kNN", value) {

  # check value
  if(!inherits(value, "igraph"))
    stop("The 'value' should be of an igraph class!")

  # all vertices
  spobject <- vrSpatialPoints(object)

  # check if there exists graphs
  graph <- object@graph
  if(length(names(object@graph)) == 0 || !graph.type %in% names(object@graph)){
    graph[[graph.type]] <- make_empty_graph(directed = FALSE) + vertices(spobject)
  }

  # vertices
  new_vert <- V(value)$name

  # edges
  subg_inv <- igraph::induced_subgraph(graph[[graph.type]], spobject[!spobject%in%new_vert])
  graph[[graph.type]] <- igraph::disjoint_union(value, subg_inv)

  # update object
  object@graph <- graph

  # return
  return(object)
}

#' vrGraphNames
#'
#' Get names of all graphs
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrGraphNames
#'
#' @export
vrGraphNames <- function(object, assay = NULL){
  return(names(object@graph))
}

#' subset_graphs
#'
#' Given a VoltRon object and a vrMetadata, subset the graph
#'
#' @param object a VoltRon Object
#' @param spatialpoints a set of spatial points
#'
#' @importFrom igraph subgraph V
#'
#' @noRd
subset_graphs <- function(object, spatialpoints){

  # graph names
  graphnames <- vrGraphNames(object)

  # # get spatialpoints
  # spatialpoints <- vrSpatialPoints(metadata, assay = vrAssayNames(object))

  # for all graphs
  if(!is.null(graphnames)){
    graph_list <- object@graph
    for(g in vrGraphNames(object)){
      cur_graph <- graph_list[[g]]
      cur_graph<- igraph::subgraph(cur_graph, igraph::V(cur_graph)[names(igraph::V(cur_graph)) %in% spatialpoints])
      graph_list[[g]] <- cur_graph
    }
  } else {
    graph_list <- list()
  }

  return(graph_list)
}

#' merge_graphs
#'
#' Given a VoltRon object, and a list of VoltRon objects, merge their graphs.
#'
#' @param object a VoltRon Object
#' @param object_list a list of VoltRon objects
#'
#' @importFrom igraph disjoint_union
#'
#' @noRd
merge_graphs <- function(object, object_list){

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  if(inherits(object, "VoltRon")){
    object_list <- c(object, object_list)
  } else {
    object_list <- c(list(object), object_list)
  }

  # choose objects
  obj1 <- object_list[[1]]
  obj2 <- object_list[[2]]

  # initial combination
  if(length(object_list) > 2){
    combined_graph <- merge_graphs(obj1, obj2)
    for(i in 3:(length(object_list))){
      combined_graph <- merge_graphs(combined_graph, object_list[[i]])
    }
  } else {
    updateobjects <- updateGraphAssay(obj1, obj2)
    obj1 <- updateobjects$object1
    obj2 <- updateobjects$object2
    combined_graph <- igraph::disjoint_union(obj1, obj2)
  }

  return(combined_graph)
}

#' updateGraphAssay
#'
#' @param object1 VoltRon object
#' @param object2 VoltRon object
#'
#' @importFrom igraph V
#' @importFrom stringr str_extract
#'
#' @noRd
updateGraphAssay <- function(object1, object2){

  if(inherits(object1, "VoltRon"))
    object1 <- vrGraph(object1, assay = "all")
  if(inherits(object2, "VoltRon"))
    object2 <- vrGraph(object2, assay = "all")

  # get assay types
  assaytype <- unique(stringr::str_extract(igraph::V(object1)$name, "Assay[0-9]+$"))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", 1:length(assaytype))
  vertex_names <- igraph::V(object1)$name
  temp <- vertex_names
  for(i in 1:length(assaytype))
    temp[grepl(paste0(assaytype[i],"$"), vertex_names)] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                vertex_names[grepl(paste0(assaytype[i],"$"), vertex_names)])
  igraph::V(object1)$name <- temp

  # get assay types
  assaytype <- unique(stringr::str_extract(igraph::V(object2)$name, "Assay[0-9]+$"))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  vertex_names <- igraph::V(object2)$name
  temp <- vertex_names
  for(i in 1:length(assaytype))
    temp[grepl(paste0(assaytype[i],"$"), vertex_names)] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                vertex_names[grepl(paste0(assaytype[i],"$"), vertex_names)])
  igraph::V(object2)$name <- temp

  # return
  return(list(object1 = object1, object2 = object2))
}

#' combineGraphs
#'
#' Combining the edges of multiple graphs
#'
#' @param object a VoltRon Object
#' @param graph.names a vector of graph names
#' @param graph.weights the weights for edges of each graph.
#' @param graph.key the name of the combined graph
#'
#' @importFrom igraph union edge_attr_names as_adjacency_matrix graph_from_adjacency_matrix
#'
#' @export
combineGraphs <- function(object, graph.names = NULL, graph.weights = NULL, graph.key = "combined"){

  if(!inherits(object, "VoltRon"))
    stop("Object must be of VoltRon class!")

  if(length(graph.names) == 0)
    stop("Please provide graph names")

  if(any(!graph.names %in% vrGraphNames(object))){
    graph.names <- setdiff(graph.names, vrGraphNames(object))
    stop("The following graphs are not included in the VoltRon object: ",
         paste(graph.names, sep = ",", collapse = TRUE))
  }

  # check weights
  if(is.null(graph.weights)){
    graph.weights <- rep(0.5, length(graph.names))
  }
  if(length(graph.weights) != length(graph.names)){
    stop("The weights should be of the length of graph names")
  }
  if(any(!is.numeric(graph.weights))){
    stop("Weights should be numeric")
  }
  if(sum(graph.weights) != 1){
    stop("Weights should sum up to 1!")
  }
  names(graph.weights) <- graph.names

  # collect graphs
  allmat <- NULL
  # gr_list <- list()
  for(gr in graph.names){
    # gr_list[[gr]] <- vrGraph(object, graph.type = gr)
    # weights <- E(gr_list[[gr]])$weight
    # if(is.null(weights)){
    #   E(gr_list[[gr]])$weight <- graph.weights[gr]
    # } else {
    #   E(gr_list[[gr]])$weight <- weights*graph.weights[gr]
    # }
    cur_graph <- vrGraph(object, graph.type = gr)
    if("weight" %in% igraph::edge_attr_names(cur_graph)){
      adjmat <- igraph::as_adjacency_matrix(cur_graph, attr = "weight")
    } else {
      adjmat <- igraph::as_adjacency_matrix(cur_graph)
    }
    adjmat <- adjmat*graph.weights[gr]
    if(is.null(allmat)){
      allmat <- adjmat
    } else {
      allmat <- allmat + adjmat
    }
  }

  # union of graphs
  # combined_gr <- igraph::union(gr_list[[1]], gr_list[-1])
  # combined_gr <- simplify(combined_gr, edge.attr.comb=list(weight="sum"))
  # vrGraph(object, graph.type = graph.key) <- combined_gr
  vrGraph(object, graph.type = graph.key) <- igraph::graph_from_adjacency_matrix(allmat, mode = "undirected", weighted = TRUE, diag = FALSE)

  # return
  return(object)
}

