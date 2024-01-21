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
      cat("  Layers:", paste(layers, collapse = " "), "\n")
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
#'  \item{\code{$}, \code{$<-}}{Name (\code{i}) of a single metadata column from the main assay, see \code{vrMainAssay}}
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
  return(metadata[[i]])
}

#' @describeIn VoltRon-methods Metadata overwrite for \code{VoltRon} objects
#'
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
#'
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
      if(!inherits(value, "vrSample")) {
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
#' @param data the count table
#' @param metadata a metadata object of class \code{vrMetadata}
#' @param image the image of the data
#' @param coords the coordinates of the spatial points
#' @param segments the segments of the spatial points, optional
#' @param sample.metadata a data frame of the sample metadata
#' @param main.assay the name of the main assay of the object
#' @param assay.type the type of the assay (cells, spots, ROIs)
#' @param params additional parameters of the object
#' @param sample_name the name of the sample
#' @param layer_name the name of the layer
#' @param image_name the name/key of the image
#' @param project project name
#'
#' @importFrom igraph make_empty_graph V vertices
#' @importFrom methods new
#' @importFrom data.table data.table
#' @importFrom rlang %||%
#' @importFrom ids random_id
#'
#' @export
#'
formVoltRon <- function(data = NULL, metadata = NULL, image = NULL,
                             coords,
                             segments = list(),
                             sample.metadata = NULL,
                             main.assay = "Custom_cell", assay.type = "cell", params = list(),
                             sample_name = NULL, layer_name = NULL, image_name = NULL,
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

  # set meta data if its empty
  if(is.null(metadata)){

    # set metadata
    vr_metadata <- setVRMetadata(molecule = data.table::data.table(), cell = data.frame(), spot = data.frame(), ROI = data.frame(), tile = data.table::data.table())

    # create entity IDs using Assay index, make it colnames
    entityID <- stringr::str_replace(entityID_nopostfix, pattern = "$", paste0("_Assay1"))
    colnames(data) <- entityID

    # create metadata
    slot(vr_metadata, name = assay.type) <- data.frame(Count = colSums(data), Assay = main.assay, Layer = layer_name, Sample = sample_name, row.names = entityID)

  } else {
    if(any(class(metadata) %in% c("data.table", "data.frame", "matrix"))){
      vr_metadata <- setVRMetadata(molecule = data.table::data.table(), cell = data.frame(), spot = data.frame(), ROI = data.frame(), tile = data.table::data.table())

      # if metadata is a data.table
      if(inherits(metadata, "data.table")){

        # check ID names
        if(length(setdiff(metadata$id, entityID_nopostfix)) > 0){
          stop("Entity IDs are not matching")
        } else {

          # entity IDs
          metadata <- subset(metadata, subset = entityID_nopostfix %in% id)

          # create entity IDs using Assay index, make it colnames
          set.seed(nrow(metadata$id))
          entityID <- paste0(metadata$id, "_", ids::random_id(bytes = 3, use_openssl = FALSE))
          colnames(data) <- entityID

          if(nrow(data) > 0){
            slot(vr_metadata, name = assay.type) <- data.table::data.table(id = entityID, assay_id = "Assay1", Count = colSums(data), Assay = main.assay,
                                                                           Layer = layer_name, Sample = sample_name, metadata)
          } else{
            slot(vr_metadata, name = assay.type) <- data.table::data.table(id = entityID, assay_id = "Assay1", Assay = main.assay,
                                                                           Layer = layer_name, Sample = sample_name, metadata)
          }
        }

      # if metadata is a regular data.frame
      } else if(inherits(metadata, "data.frame")){

        # check row names
        if(length(setdiff(rownames(metadata), entityID_nopostfix)) > 0){
          stop("Entity IDs are not matching")
        } else {

          # entity IDs
          metadata <- metadata[entityID_nopostfix,]

          # create entity IDs using Assay index, make it colnames
          entityID <- stringr::str_replace(entityID_nopostfix, pattern = "$", paste0("_Assay1"))
          colnames(data) <- entityID

          # create metadata
          if(nrow(data) > 0){
            slot(vr_metadata, name = assay.type) <- data.frame(Count = colSums(data), Assay = main.assay, Layer = layer_name, Sample = sample_name, metadata, row.names = entityID)
          } else{
            slot(vr_metadata, name = assay.type) <- data.frame(Assay = main.assay, Layer = layer_name, Sample = sample_name, metadata, row.names = entityID)
          }
        }
      }
    }
  }

  # Coordinates
  if(!is.null(coords)){
    if(length(colnames(coords)) == 2){
      rownames(coords) <- entityID
      colnames(coords) <- c("x", "y")
    } else {
      stop("The length of colnames of the coordinates matrix should two!")
    }
  } else {
    stop("There are no coordinates matrix provided!")
  }

  # create vrAssay
  Assay <- formAssay(data = data, coords = coords, segments = segments, image = image, params = params, type = assay.type, name = "Assay1", main_image = image_name)
  listofAssays <- list(Assay)
  names(listofAssays) <- main.assay

  # create layers and samples
  listofLayers <- list(methods::new("vrLayer",
                                    assay = listofAssays,
                                    connectivity = igraph::make_empty_graph(directed = FALSE) + igraph::vertices(entityID)))
  names(listofLayers) <- layer_name
  listofSamples <- list(methods::new("vrSample", layer = listofLayers))
  names(listofSamples) <- sample_name

  # set sample meta data
  if(is.null(sample.metadata)){
    sample.metadata <- setVRSampleMetadata(listofSamples)
  }

  # set VoltRon class
  methods::new("VoltRon", samples = listofSamples, metadata = vr_metadata, sample.metadata = sample.metadata, main.assay = main.assay, project = project)
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
#' @param metadata a predefined metadata
#' @param assay_name assay name
#' @param connectivity a metadata of edges representing connected spatial points across assays
#' @param sample sample name
#' @param layer layer name
#'
#' @rdname addAssay
#' @method addAssay VoltRon
#'
#' @importFrom igraph make_empty_graph add_edges vertices
#'
#' @export
#'
addAssay.VoltRon <- function(object, assay, metadata = NULL, assay_name, connectitivity = NULL, sample = "Sample1", layer = "Section1"){

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
  g_assay <- igraph::make_empty_graph(directed = FALSE) + igraph::vertices(vrSpatialPoints(assay))
  g_layer <- curlayer@connectivity + g_assay
  object[[sample, layer]]@connectivity <- g_layer

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
    names(listofLayers) <- cur_sample.metadata$NewLayer
    listofSamples <- list(methods::new("vrSample", layer = listofLayers))
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

#' configureConnectivity
#'
#' add connectivity information to the assays of the same layer
#'
#' @param assay assay
#' @param metadata a predefined metadata
#' @param assay_name assay name
#' @param connectivity a metadata of edges representing connected spatial points across assays
#' @param sample sample name
#' @param layer layer name
#'
#' @importFrom igraph add_edges
#'
addConnectivity <- function(object, connectivity, sample, layer){

  # get sample and layer
  curlayer <- object[[sample, layer]]

  # make edges from connectivity matrix
  connectivity <- as.vector(t(as.matrix(connectivity)))

  # add edges
  curlayer@connectivity <- igraph::add_edges(curlayer@connectivity, edges = connectivity)

  # return
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
#' @param use_points if \code{interactive} is \code{TRUE}, use spatial points instead of the reference image
#'
#' @rdname subset
#' @method subset VoltRon
#'
#' @importFrom rlang enquo eval_tidy quo_get_expr quo_text
#' @importFrom stringr str_extract
#' @importFrom methods new
#'
#' @export
#'
subset.VoltRon <- function(object, subset, samples = NULL, assays = NULL, spatialpoints = NULL, features = NULL, image = NULL, interactive = FALSE, use_points = FALSE) {

  # subseting based on subset argument
  if (!missing(x = subset)) {
    subset <- rlang::enquo(arg = subset)
  }
  if(!missing(subset)){
    metadata <- Metadata(object)
    name <- strsplit(rlang::quo_text(subset), split = " ")[[1]][1]
    if(name %in% colnames(metadata)){
      spatialpoints <- rownames(metadata)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
    } else {
      stop("Column '", name, "' is not found in the metadata")
    }
    spatialpoints <- rownames(metadata)[eval_tidy(rlang::quo_get_expr(subset), data = metadata)]
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
    metadata <- subset.vrMetadata(object@metadata, assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, assays = assays)
    }, USE.NAMES = TRUE)

  # subsetting on entity names
  } else if(!is.null(spatialpoints)) {

    metadata <- subset.vrMetadata(object@metadata, spatialpoints = spatialpoints)
    assays <- unique(stringr::str_extract(vrSpatialPoints(metadata), "Assay[0-9]+"))
    sample.metadata <- subset_sampleMetadata(sample.metadata, assays = assays)
    samples <- unique(sample.metadata$Sample)
    listofSamples <- sapply(object@samples[samples], function(samp) {
      subset.vrSample(samp, spatialpoints = spatialpoints)
    }, USE.NAMES = TRUE)

  # subsetting on features
  } else if(!is.null(features)){
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
    if(inherits(image, "character")){

      # check if there are only one image and one assay
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
      stop("Please provide a character based subsetting notation, see magick::image_crop documentation")
    }
  } else if(interactive){
    results <- demuxVoltRon(object, use_points = use_points)
    return(results)
  }

  # main.assay
  main.assay <- unique(sample.metadata$Assay)[unique(sample.metadata$Assay) == names(table(sample.metadata$Assay))[which.max(table(sample.metadata$Assay))]]

  # project
  project <- object@project

  # subset graphs
  graph_list <- subset_graphs(object, metadata)

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
#' @param main.assay name of the assay
#'
#' @method merge VoltRon
#' @importFrom methods new
#'
#' @export
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
  message("Merging metadata ...")
  metadata_list <- lapply(object_list, function(x) slot(x, name = "metadata"))
  metadata <- merge(metadata_list[[1]], metadata_list[-1])

  # combine samples and rename layers
  message("Merging samples and layers ...")
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
  for(assy in rownames(sample.metadata))
    vrAssayNames(object[[assy]]) <- assy

  # change sample names
  if(!is.null(samples))
    object$Sample <- samples

  # return
  object
}

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints VoltRon
#'
#' @export
#'
vrSpatialPoints.VoltRon <- function(object, ...) {
  return(vrSpatialPoints(object@metadata, ...))
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
#' @param norm TRUE if normalized data should be returned
#' @param ... additional parameters passed to \code{vrData.vrAssay}
#'
#' @rdname vrData
#' @method vrData VoltRon
#'
#' @export
#'
vrData.VoltRon <- function(object, assay = NULL, norm = FALSE, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  data <- NULL
  for(i in 1:length(assay_names)){
    cur_data <- vrData(object[[assay_names[i]]], norm = norm, ...)
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
#' @param dims the set of dimensions of the embedding data
#' @param type the key name for the embedding, i.e. "pca" or "umap"
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

#' @param assay assay
#' @param type the key name for the embedding
#' @param overwrite Whether the existing embedding with name 'type' should be overwritten
#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @method vrEmbeddings<- VoltRon
#'
#' @export
#'
"vrEmbeddings<-.VoltRon" <- function(object, assay = NULL, type = "pca", overwrite = FALSE, ..., value) {

  # check if the embedding exists
  if(type %in% vrEmbeddingNames(object) && !overwrite)
    stop("An embedding named '", type, "' already exists in this object. Do overwrite = TRUE for replacing with the existing one.")

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # set embeddings
  for(assy in assay_names){
    assayobject <- object[[assy]]
    if(vrAssayTypes(assayobject) %in% c("ROI", "cell", "spot")){
      vrEmbeddings(assayobject, type = type) <- value[grepl(paste0(assy, "$"), rownames(value)),]
    } else {
      vrEmbeddings(assayobject, type = type) <- value[vrSpatialPoints(assayobject),]
    }
    object[[assy]] <- assayobject
  }

  return(object)
}

#' @param assay assay
#'
#' @rdname vrEmbeddingNames
#' @method vrEmbeddingNames VoltRon
#'
#' @export
#'
vrEmbeddingNames.VoltRon <- function(object, assay = NULL){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get assay types
  embed_names <- unique(unlist(lapply(assay_names, function(x) vrEmbeddingNames(object[[x]]))))

  return(embed_names)
}

#### Metadata ####

#' @param assay assay
#' @param type the assay type: ROI, spot or cell, or all for the entire metadata object
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
  if(type %in% slotNames(object@metadata)){

    # sample metadata
    sample.metadata <- SampleMetadata(object)

    # get assay names
    assay_names <- vrAssayNames(object, assay = assay)

    # get metadata
    metadata <- slot(object@metadata, name = type)
    if(inherits(metadata, "data.table")){
      metadata <- subset(metadata, assay_id %in% assay_names)
    } else {
      metadata <- metadata[stringr::str_extract(rownames(metadata), "Assay[0-9]+") %in% assay_names, ]
    }
    return(metadata)
  } else {
    stop("Please provide one of three assay types: 'ROI', 'cell', 'spot' and 'molecules'.")
  }
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

  if(!is.data.frame(value))
    stop("The new or updated metadata has to be a data frame")

  if(is.null(rownames(value)))
    stop("The new metadata should have row names to match its rows with the existing one")

  if(is.null(type)){
    type <- unique(vrAssayTypes(object, assay = assay))
  }

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get metadata
  metadata <- slot(object@metadata, name = type)
  # cur_metadata <- metadata[stringr::str_extract(rownames(metadata), "Assay[0-9]+") %in% assay_names, ]

  # if(type %in% slotNames(object@metadata)){
  if(type %in% c("ROI", "cell", "spot")){

    # replace the metadata (or some part of it) with the new value
    if(length(setdiff(rownames(values), rownames(metadata))) == 0){

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

  } else if(type %in% c("tile", "molecule")){

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
      # metadata[id %in% value$id, names(metadata):=value]
      metadata <- value
      slot(object@metadata, name = type) <- metadata

    } else {
      stop("Some rows of new data frame are not available in the metadata")
    }

  } else {
    stop("Please provide one of three assay types: 'ROI', 'cell', 'spot'.")
  }


  # if(type == "all"){
  #   all_types <- names(slotToList(object@metadata))
  #   for(cur_type in all_types){
  #     cur_metadata <- slot(object@metadata, name = cur_type)
  #     if(nrow(cur_metadata) > 0){
  #       slot(object@metadata, name = cur_type) <- value
  #     }
  #     slot(object@metadata, name = cur_type)
  #   }
  # } else {
  #   slot(object@metadata, name = type) <- value
  # }

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

#### Spatial ####

#' @param assay assay
#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrCoordinates
#' @method vrCoordinates VoltRon
#'
#' @export
#'
vrCoordinates.VoltRon <- function(object, assay = NULL, image_name = NULL, reg = FALSE, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  coords <- NULL
  for(assy in assay_names)
    coords <- rbind(coords, vrCoordinates(object[[assy]], image_name = image_name, reg = reg))

  # return image
  return(coords)
}

#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D coordinates
#'
#' @rdname vrCoordinates
#' @method vrCoordinates<- VoltRon
#'
#' @export
#'
"vrCoordinates<-.VoltRon" <- function(object, image_name = NULL, reg = FALSE, ..., value) {

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check the number of assays in the object
  if(nrow(sample.metadata) > 1)
    stop("Changing the coordinates of multiple assays in the same time are not permitted!")

  # get assay
  cur_assay <- sample.metadata[1,]
  vrlayer <- object[[cur_assay$Sample, cur_assay$Layer]]
  vrassay <- vrlayer[[cur_assay$Assay]]

  # change coordinates
  vrCoordinates(vrassay, image_name = image_name, reg = reg) <- value
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
#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#'
#' @rdname vrSegments
#' @method vrSegments VoltRon
#'
#' @export
#'
vrSegments.VoltRon <- function(object, assay = NULL, image_name = NULL, reg = FALSE, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get all coordinates
  segts <- NULL
  for(assy in assay_names)
    segts <- c(segts, vrSegments(object[[assy]], image_name = image_name, reg = reg))

  # return image
  return(segts)
}

#' @param image_name the key of the image associated with the coordinates
#' @param reg TRUE if registered segments are being updated
#' @param value the new set of 2D segments for each spatial point
#'
#' @rdname vrSegments
#' @method vrSegments<- VoltRon
#'
#' @export
#'
"vrSegments<-.VoltRon" <- function(object, image_name = NULL, reg = FALSE, ..., value) {

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
  vrSegments(vrassay, image_name = image_name, reg = reg) <- value
  vrlayer[[cur_assay$Assay]] <- vrassay
  object[[cur_assay$Sample, cur_assay$Layer]] <- vrlayer

  return(object)
}

#### Graphs ####

#' @param assay assay
#' @param graph.type the type of the graph, either custom or given by \code{getProfileNeighbors} or \code{getSpatialNeighbors} functions
#'
#' @rdname vrGraph
#' @method vrGraph VoltRon
#'
#' @importFrom igraph induced_subgraph
#' @export
#'
vrGraph.VoltRon <- function(object, assay = NULL, graph.type = "kNN", ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  # assay_pattern <- paste0(assay_names, "$", collapse = "|")
  # node_names <- vrSpatialPoints(object)[grepl(assay_pattern, vrSpatialPoints(object))]
  node_names <- vrSpatialPoints(object, assay = assay_names)

  # check if there exists graphs
  if(length(names(object@graph)) == 0)
    stop("There are no graphs in this VoltRon object!")

  # check graph type
  if(!graph.type %in% names(object@graph))
    stop("The graph name '", graph.type, "' can't be found in this VoltRon object!")

  # return graph
  if(length(object@graph[[graph.type]]) > 0){
    returngraph <- igraph::induced_subgraph(object@graph[[graph.type]], node_names)
    return(returngraph)
  } else {
    warning("This VoltRon object does not have any graphs yet!")
    return(NULL)
  }
}

#' @param assay assay
#' @param value new Feature Data
#'
#' @rdname vrGraph
#' @method vrGraph<- VoltRon
#'
#' @export
#'
"vrGraph<-.VoltRon" <- function(object, assay = NULL, graph.type = "kNN", ..., value) {

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

#' @param assay assay
#'
#' @rdname vrGraphNames
#' @method vrGraphNames VoltRon
#'
#' @export
#'
vrGraphNames.VoltRon <- function(object, assay = NULL){
  return(names(object@graph))
}

#' subset_graphs
#'
#' Given a VoltRon object and a vrMetadata, subset the graph
#'
#' @param object a VoltRon Object
#' @param metadata a vrMetadata Object
#'
#' @importFrom igraph subgraph V
#'
#' @noRd
subset_graphs <- function(object, metadata){

  # graph names
  graphnames <- vrGraphNames(object)

  # for all graphs
  if(!is.null(graphnames)){
    graph_list <- object@graph
    for(g in vrGraphNames(object)){
      cur_graph <- graph_list[[g]]
      cur_graph<- igraph::subgraph(cur_graph, igraph::V(cur_graph)[names(igraph::V(cur_graph)) %in% vrSpatialPoints(metadata)])
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
#' @param graph.types a vector of graph names
#' @param weights the weights for edges of each graph.
#' @param graph.key the name of the combined graph
#'
#' @importFrom igraph union
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
    if("weight" %in% edge_attr_names(cur_graph)){
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

