#' @include zzz.R
#' @include generics.R
#'
#' @useDynLib VoltRon
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
    graph = "igraph",
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
    sample_names <- names(object@samples)
    show_length <- min(5,length(sample_names))
    for(samp in sample_names[1:show_length]){
      cat(samp, ": \n", sep = "")
      layers <- names(unlist(object@samples[[samp]]@layer))
      cat("  Layers:", paste(layers, collapse = " "), "\n")
      assays <- sapply(names(object@samples[[samp]]@layer), function(x) names(object[[samp, x]]@assay))
      all_assays <- c(all_assays, assays)
    }

    if(length(sample_names) > 5){
      cat("...", "\n")
      cat("There are", length(sample_names), "samples in total", "\n")
    }

    # print assays
    main.assay <- vrMainAssay(object)
    unique_assays <- unique(all_assays)
    unique_assays <- unique_assays[c(which(unique_assays == main.assay),which(unique_assays != main.assay))]
    unique_assays[1] <- paste0(unique_assays[1], "(Main)")
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

  # get assay names
  assay_names <- vrAssayNames(x)

  # metadata
  metadata <- Metadata(x, assay = assay_names)

  # get metadata column
  return(metadata[[i]])
}

#' @export
#' @method $<- VoltRon
#'
"$<-.VoltRon" <- function(x, i, ..., value) {

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
    if(!nrow(metadata) %in% c(1,length(value))){
      stop("The new or the existing column should of length 1 or the same as the number of rows")
    } else {
      metadata[[i]] <- value
      Metadata(x, assay = assay_names) <- metadata
    }
  }

  return(x)
}

#' #' @export
#' #' @method $<- VoltRon
#' #'
#' "$<-.VoltRon" <- function(x, i, ..., value) {
#'
#'   # update sample metadata and metadata
#'   x@sample.metadata[[i]] <- value
#'   Metadata(x)[[i]] <- value
#'
#'   # update if sample names are provided
#'   if(i == "Sample"){
#'     if(!length(x@sample.metadata) %in% c(1,length(value))){
#'       stop("Sample names should be of length 1 or the same length as the number of assays!")
#'     }
#'     names(x@samples) <- value
#'   }
#'
#'   return(x)
#' }

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

#' formVoltRon
#'
#' Create a VoltRon object
#'
#' @param data the count table
#' @param metadata a metadata object of class \code{vrMetadata}
#' @param image the image of the data
#' @param coord the coordinates of the spatial points
#' @param segments the segments of the spatial points, optional
#' @param sample.metadata a data frame of the sample metadata
#' @param graph the graph to determine the adjacency of spatial points across layers
#' @param main.assay the name of the main assay of the object
#' @param assay_name the name of the assay
#' @param assay.type the type of the assay (cells, spots, ROIs)
#' @param params additional parameters of the object
#' @param sample_name the name of the sample
#' @param layer_name the name of the layer
#' @param project project name
#'
#' @import igraph
#'
#' @export
#'
formVoltRon <- function(data, metadata = NULL, image = NULL,
                             coords, segments = list(),
                             sample.metadata = NULL, graph = NULL,
                             main.assay = "Custom_cell", assay.type = "cell", params = list(),
                             sample_name = NULL, layer_name = NULL,
                             project = NULL){

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
        slot(sr_metadata, name = assay.type) <- data.frame(Count = colSums(data), Assay = main.assay, Layer = layer_name, Sample = sample_name, metadata, row.names = entityID)
      }
    }
  }

  # Coordinates
  if(!is.null(coords)){
    colcoords <- colnames(coords)
    if(length(colnames(coords)) == 2){
      rownames(coords) <- entityID
      colnames(coords) <- c("x", "y")
    } else {
      stop("The length of colnames of the coordinates matrix should two!")
    }
  } else {
    stop("There are no coordinates matrix provided!")
  }

  # Segments
  if(length(segments) > 0) names(segments) <- entityID

  # set zgraph
  if(is.null(graph)){
    spatial_points <- vrSpatialPoints(sr_metadata)
    graph <- igraph::make_empty_graph(n = length(spatial_points), directed = FALSE)
    igraph::V(graph)$name <- spatial_points
  }

  # create vrAssay
  Xenium_assay <- new("vrAssay", rawdata = data, normdata = data, coords = coords, segments = segments, image = as.raster(image), params = params, type = assay.type)
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
  new("VoltRon", samples = listofSamples, metadata = sr_metadata, sample.metadata = sample.metadata, graph = graph, main.assay = main.assay, project = project)
}

### Assay Methods ####

#' @rdname vrMainAssay
#' @method vrMainAssay VoltRon
#'
#' @export
#'
vrMainAssay.VoltRon <- function(object, ...) {
  object@main.assay
}

#' @param value new assay name
#'
#' @rdname vrMainAssay
#' @method vrMainAssay<- VoltRon
#'
#' @export
#'
"vrMainAssay<-.VoltRon" <- function(object, ..., value) {
  sample.metadata <- SampleMetadata(object)
  assay_names <- unique(sample.metadata$Assay)
  if(!value %in% assay_names){
    stop("There is no assay names '", value, "' in this object")
  } else {
    object@main.assay <- value
  }
  return(object)
}

#' @param assay assay
#' @param assay_name assay name
#' @param sample sample name
#' @param layer layer name
#'
#' @rdname addAssay
#' @method addAssay VoltRon
#'
#' @importFrom igraph union V
#' @export
#'
addAssay.VoltRon <- function(object, assay, assay_name, sample = "Sample1", layer = "Section1"){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay id
  assay_ids <- as.numeric(gsub("Assay", "", rownames(sample.metadata)))
  assay_id <- paste0("Assay", max(assay_ids)+1)
  assay_names <- c(rownames(sample.metadata), assay_id)

  # update sample.metadata and metadata
  object@sample.metadata <- rbind(sample.metadata, c(assay_name, layer, sample))
  rownames(object@sample.metadata) <- assay_names
  object@metadata <- addAssay(object@metadata,
                              assay = assay, assay_name = assay_name,
                              sample = sample, layer = layer)

  # update sample and layer
  assay_list <- object[[sample, layer]]@assay
  vrAssayNames(assay) <- assay_id
  new_assay_list <- list(assay)
  names(new_assay_list) <- assay_name
  assay_list <- c(assay_list, new_assay_list)
  object[[sample, layer]]@assay <- assay_list

  # update graph
  newgraph <- igraph::make_empty_graph(n = length(vrSpatialPoints(assay)), directed = FALSE)
  igraph::V(newgraph)$name <- vrSpatialPoints(assay)
  object@graph <- igraph::union(object@graph, newgraph)

  # return
  return(object)
}

#' @param assay assay
#'
#' @rdname vrAssayNames
#' @method vrAssayNames VoltRon
#'
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

#' @param assay assay
#'
#' @rdname vrAssayTypes
#' @method vrAssayTypes VoltRon
#'
#' @export
#'
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
#' @rdname changeSampleNames
#' @method changeSampleNames VoltRon
#'
#' @param object a VoltRon object
#' @param samples a single or a set of sample names
#'
#' @import dplyr
#'
changeSampleNames.VoltRon <- function(object, samples = NULL){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # old to new samples table
  samples_table <- data.frame(sample.metadata, AssayID = rownames(sample.metadata), NewSample = samples)

  # check if multiple new sample names are associated with the same section of one sample
  check_samples_table <- samples_table %>%
    group_by(Assay, Sample) %>% mutate(n = dplyr::n_distinct(NewSample)) %>%
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
    names(listofLayers) <- cur_sample.metadata$NewLayer
    listofSamples <- list(new("vrSample", layer = listofLayers))
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

  return(object)
}


### Object Methods ####

#' Subsetting VoltRon objects
#'
#' Given a VoltRon object, subset the object given one of the attributes
#'
#' @param object A vrAssay object
#' @param subset Logical statement for subsetting
#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#' @param features the set of features to subset the object
#' @param image the subseting string passed to \code{magick::image_crop}
#' @param interactive TRUE if interactive subsetting on the image is demanded
#'
#' @rdname subset
#' @method subset VoltRon
#'
#' @importFrom rlang enquo eval_tidy quo_get_expr
#' @importFrom stringr str_extract
#' @import igraph
#'
#' @export
#'
subset.VoltRon <- function(object, subset, samples = NULL, assays = NULL, spatialpoints = NULL, features = NULL, image = NULL, interactive = FALSE) {

  # subseting based on subset argument
  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }
  if(!missing(subset)){
    metadata <- Metadata(object)
    spatialpoints <- rownames(metadata)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
    object <- subset(object, spatialpoints = spatialpoints)
    return(object)
  }

  # subseting on other attributes
  attrinfo <- c(sapply(list(samples, assays, spatialpoints, features), function(x) length(x) > 0), interactive)
  if(sum(attrinfo) > 1){
    stop("Please choose only one of the subsetting attributes: 'samples', 'assays', 'spatialpoints', 'features' or 'interactive'")
  }

  if(!is.null(samples)){

    sample.metadata <- subset.sampleMetadata(SampleMetadata(object), samples = samples)
    metadata <- subset.vrMetadata(object@metadata, samples = samples) # CAN WE CHANGE THIS TO ONLY SUBSET LATER ????
    listofSamples <- object@samples[samples]

  # subsetting on assays name
  } else if(!is.null(assays)) {

    sample.metadata <- subset.sampleMetadata(SampleMetadata(object), assays = assays)
    metadata <- subset.vrMetadata(object@metadata, assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, assays = assays)
    }, USE.NAMES = TRUE)

  # subsetting on entity names
  } else if(!is.null(spatialpoints)) {

    metadata <- subset.vrMetadata(object@metadata, spatialpoints = spatialpoints)
    assays <- unique(stringr::str_extract(vrSpatialPoints(metadata), "Assay[0-9]+"))
    sample.metadata <- subset.sampleMetadata(SampleMetadata(object), assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE)

  # subsetting on features
  } else if(!is.null(features)){

    sample.metadata <- SampleMetadata(object)
    assay_names <- vrAssayNames(object)
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
      sample.metadata <- SampleMetadata(object)
      if(nrow(sample.metadata) > 1){
        stop("Subseting on images can only be performed on VoltRon objects with a single assay")
      } else {
        samples <- unique(sample.metadata$Sample)
        listofSamples <- sapply(object@samples[samples], function(samp) {
          subset.vrSample(samp, image = image)
        }, USE.NAMES = TRUE)
        spatialpoints <-  do.call(c, lapply(listofSamples, vrSpatialPoints.vrSample))
        metadata <- subset.vrMetadata(object@metadata, spatialpoints = spatialpoints)
      }
    } else {
      stop("Please provide a character based subsetting notation, see magick documentation")
    }
  } else if(interactive){
    results <- demuxVoltRon(object)
    return(results)
  }

  # other attributes
  main.assay <- unique(sample.metadata$Assay)[unique(sample.metadata$Assay) == names(table(sample.metadata$Assay))[which.max(table(sample.metadata$Assay))]]
  graph <- igraph::subgraph(object@graph, V(object@graph)[names(V(object@graph)) %in% vrSpatialPoints(metadata)])
  project <- object@project

  # set VoltRon class
  new("VoltRon", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata, graph = graph, main.assay = main.assay, project = project)
}

#' Merging VoltRon objects
#'
#' Given a VoltRon object, and a list of VoltRon object, merge all.
#'
#' @param object a VoltRon Object
#' @param object_list a list of VoltRon objects
#' @param samples a single sample name or multiple sample names of the same size as the given VoltRon objects
#' @param main.assay name of the assay
#'
#' @method merge VoltRon
#'
#' @import igraph
#'
#' @export
#'
#'
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
  metadata_list <- lapply(object_list, function(x) slot(x, name = "metadata"))
  metadata <- merge(metadata_list[[1]], metadata_list[-1])

  # combine samples and rename layers
  listofSamples <- NULL
  for(i in 1:length(object_list)){
    cur_object <- object_list[[i]]@samples
    listofSamples <- c(listofSamples, cur_object)
  }

  # merge graphs
  graph <- merge_graphs(object_list[[1]], object_list[-1])

  # get main assay
  if(is.null(main.assay))
      main.assay <- names(sort(table(sample.metadata$Assay), decreasing = TRUE))[1]

  # project
  project <- slot(object_list[[1]], "project")

  # set VoltRon class
  object <- new("VoltRon", samples = listofSamples, metadata = metadata, sample.metadata = sample.metadata,
                graph = graph, main.assay = main.assay, project = project)

  # change assay names and sample names
  for(assy in rownames(sample.metadata))
    vrAssayNames(object[[assy]]) <- assy

  # change sample names
  if(!is.null(samples))
    object$Sample <- samples

  # return
  object
}

merge_graphs <- function(object, object_list){

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  if(class(object) == "VoltRon"){
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

updateGraphAssay <- function(object1, object2){

  if(class(object1) == "VoltRon")
    object1 <- vrGraph(object1, assay = "all")
  if(class(object2) == "VoltRon")
    object2 <- vrGraph(object2, assay = "all")

  # get assay types
  assaytype <- unique(stringr::str_extract(V(object1)$name, "Assay[0-9]+$"))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", 1:length(assaytype))
  vertex_names <- V(object1)$name
  temp <- vertex_names
  for(i in 1:length(assaytype))
    temp[grepl(paste0(assaytype[i],"$"), vertex_names)] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                vertex_names[grepl(paste0(assaytype[i],"$"), vertex_names)])
  V(object1)$name <- temp

  # get assay types
  assaytype <- unique(stringr::str_extract(V(object2)$name, "Assay[0-9]+$"))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  vertex_names <- V(object2)$name
  temp <- vertex_names
  for(i in 1:length(assaytype))
    temp[grepl(paste0(assaytype[i],"$"), vertex_names)] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                vertex_names[grepl(paste0(assaytype[i],"$"), vertex_names)])
  V(object2)$name <- temp

  # return
  return(list(object1 = object1, object2 = object2))
}





#' @param assay assay
#' @param type the assay type: ROI, spot or cell
#'
#' @rdname Metadata
#' @method Metadata VoltRon
#'
#' @export
#'
Metadata.VoltRon <- function(object, assay = NULL, type = NULL) {
  if(is.null(type)){
    type <- unique(vrAssayTypes(object, assay = assay))
  } else{
    if(type == "all")
      return(object@metadata)
  }
  return(slot(object@metadata, name = type))
}

#' @param type the assay type: ROI, spot or cell
#' @param assay assay
#' @param value new metadata
#'
#' @rdname Metadata
#' @method Metadata<- VoltRon
#'
#' @export
#'
"Metadata<-.VoltRon" <- function(object, type = NULL, assay = NULL, ..., value) {

  if(is.null(type)){
    type <- unique(vrAssayTypes(object, assay = assay))
  }

  slot(object@metadata, name = type) <- value

  if(type == "all"){
    all_types <- names(slotToList(object@metadata))
    for(cur_type in all_types){
      cur_metadata <- slot(object@metadata, name = cur_type)
      if(nrow(cur_metadata) > 0){
        slot(object@metadata, name = cur_type) <- value
      }
      slot(object@metadata, name = cur_type)
    }
  } else {
    slot(object@metadata, name = type) <- value
  }

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

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints VoltRon
#'
#' @export
#'
vrSpatialPoints.VoltRon <- function(object, ...) {
  return(vrSpatialPoints(object@metadata))
}

#' @param assay assay
#'
#' @rdname vrFeatures
#' @method vrFeatures VoltRon
#'
#' @export
#'
vrFeatures.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all features
  features <- NULL
  for(assy in assay_names)
    features <- c(features, vrFeatures(object[[assy]]))

  return(unique(features))
}

#' @param assay assay
#'
#' @rdname vrFeatureData
#' @method vrFeatureData VoltRon
#'
#' @export
#'
vrFeatureData.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all features
  features <- vrFeatureData(object[[assay_names[1]]], ...)

  # return
  return(features)
}

#' @param assay assay
#' @param value new Feature Data
#'
#' @rdname vrFeatureData
#' @method vrFeatureData<- VoltRon
#'
#' @export
#'
"vrFeatureData<-.VoltRon" <- function(object, assay = NULL, ..., value) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names)
    vrFeatureData(object[[assy]]) <- value

  return(object)
}

#' @param assay assay
#'
#' @rdname vrData
#' @method vrData VoltRon
#'
#' @export
#'
vrData.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  data <- NULL
  for(i in 1:length(assay_names)){
    cur_data <- vrData(object[[assay_names[i]]], ...)
    cur_data <- data.frame(cur_data, feature.ID = rownames(cur_data), check.names = FALSE)
    if(i == 1){
      data <- cur_data
    } else {
      data <- merge(data, cur_data, by = "feature.ID", all = TRUE)
    }
  }
  rownames(data) <- data$feature.ID
  data <- data[,!colnames(data) %in% "feature.ID"]
  data[is.na(data)] <- 0
  data <- as.matrix(data)
  colnames(data) <- gsub("\\.","-", colnames(data))

  return(data)
}

#' @param assay assay
#'
#' @rdname vrGraph
#' @method vrGraph VoltRon
#'
#' @export
#'
vrGraph.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  assay_pattern <- paste0(assay_names, collapse = "|")
  node_names <- vrSpatialPoints(object)[grepl(assay_pattern, vrSpatialPoints(object))]

  returngraph <- induced_subgraph(object@graph, node_names)
  return(returngraph)
}

#' @param assay assay
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrCoordinates
#' @method vrCoordinates VoltRon
#'
#' @export
#'
vrCoordinates.VoltRon <- function(object, reg = FALSE, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  coords <- NULL
  for(assy in assay_names)
    coords <- rbind(coords, vrCoordinates(object[[assy]], reg = reg))

  # return image
  return(coords)
}

#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D coordinates
#'
#' @rdname vrCoordinates
#' @method vrCoordinates<- VoltRon
#'
#' @export
#'
"vrCoordinates<-.VoltRon" <- function(object, reg = FALSE, ..., value) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the number of assays in the object
  if(nrow(sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays are not permitted!")

  # get assay
  cur_assay <- sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # change coordinates
  vrCoordinates(vrassay, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @rdname flipCoordinates
#' @method flipCoordinates VoltRon
#'
#' @export
#'
flipCoordinates.VoltRon <- function(object, ...){
  sample.metadata <- SampleMetadata(object)
  assay_names <- rownames(sample.metadata)
  for(assy in assay_names){
    object[[assy]] <- flipCoordinates(object[[assy]], ...)
  }
  return(object)
}

#' @param assay assay
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrSegments
#' @method vrSegments VoltRon
#'
#' @export
#'
vrSegments.VoltRon <- function(object, reg = FALSE, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  segts <- NULL
  for(assy in assay_names)
    segts <- c(segts, vrSegments(object[[assy]], reg = reg))

  # return image
  return(segts)
}

#' @param assay assay
#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D segments for each spatial point
#'
#' @rdname vrSegments
#' @method vrSegments<- VoltRon
#'
#' @export
#'
"vrSegments<-.VoltRon" <- function(object, reg = FALSE, ..., value) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the number of assays in the object
  if(nrow(sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays are not permitted!")

  # get assay
  cur_assay <- sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # change coordinates
  vrSegments(vrassay, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#' @param assay assay
#' @param dims the set of dimensions of the embedding data
#' @param type the key name for the embedding
#'
#' @rdname vrEmbeddings
#' @method vrEmbeddings VoltRon
#'
#' @export
#'
vrEmbeddings.VoltRon <- function(object, assay = NULL, type = "pca", dims = 1:30, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  returndata_list <- list()
  for(i in 1:length(assay_names))
    returndata_list[[i]] <- vrEmbeddings(object[[assay_names[i]]], type = type, dims = dims, ...)

  return(do.call(rbind, returndata_list))
}

#' @param type the key name for the embedding
#' @param assay assay
#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @method vrEmbeddings<- VoltRon
#'
#' @export
#'
"vrEmbeddings<-.VoltRon" <- function(object, assay = NULL, type = "pca", ..., value) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names){
    assayobject <- object[[assy]]
    vrEmbeddings(assayobject, type = type) <- value[grepl(paste0(assy, "$"), rownames(value)),]
    object[[assy]] <- assayobject
  }

  return(object)
}

