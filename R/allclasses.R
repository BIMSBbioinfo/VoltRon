#' @include zzz.R
#' @include allgenerics.R
#' @useDynLib VoltRon
NULL

## Auxiliary ####

# pseudo IterableMatrix for BPCells
if(!requireNamespace("BPCells", quietly = TRUE)){
  setClass("IterableMatrix")
} 

## vrImage ####

# Set class union
suppressWarnings({
  setClassUnion("image_matrix", 
                members = c("matrix", 
                            "data.frame",
                            "dgRMatrix", 
                            "dgeMatrix", 
                            "Array", 
                            "IterableMatrix"))
})

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

## vrAssay ####

# Set class union
suppressWarnings({
  setClassUnion("data_matrix", 
                members = c("matrix", 
                            "dgCMatrix", "dgRMatrix", "dgeMatrix", 
                            "Array", 
                            "IterableMatrix"))
})

#' The vrAssay (VoltRon Assay) Class
#'
#' @slot rawdata raw data
#' @slot normdata normalized data
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

## vrAssayV2 ####

#' The vrAssayV2 (VoltRon Assay) Class
#'
#' @slot data list of count/normalized datasets
#' @slot featuredata list of feature metadata
#' @slot embeddings list of embeddings
#' @slot image a list of vrImage objects
#' @slot params additional parameters used by different assay types
#' @slot type the type of the assay (tile, molecule, cell, spot, ROI)
#' @slot name the assay name
#' @slot main_image the key of the main image
#' @slot main_featureset the key of the main feature set
#'
#' @name vrAssayV2-class
#' @rdname vrAssayV2-class
#' @exportClass vrAssayV2
vrAssayV2 <- setClass(
  Class = 'vrAssayV2',
  slots = c(
    data = "list",
    featuredata = 'list',
    embeddings = "list",
    image = "list",
    params = "list",
    type = "character",
    name = "character",
    main_image = "character",
    main_featureset = "character"
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrAssayV2',
  definition = function(object) {
    
    # check if there is a data or rawdata slot in assay object
    cat(
      paste0("vrAssayV2 (VoltRon Assay V2) of ", 
             nrow(vrCoordinates(object)), " spatial points and ", 
             nrow(object@data[[vrMainFeatureType(object)]]), " features (", vrMainFeatureType(object), "). \n")
    )
    
    return(invisible(x = NULL))
  }
)

## vrLayer ####

# Set classes
setOldClass(Classes = c('igraph'))

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

## vrMetadata ####

suppressWarnings({
  setClassUnion("metadata_data",
                members = c("data.table", 
                            "data.frame",
                            if (requireNamespace("S4Vectors", quietly = TRUE)) "DataFrame" else NULL,
                            if (requireNamespace("HDF5DataFrame", quietly = TRUE)) "HDF5DataFrame" else NULL, 
                            if (requireNamespace("ZarrDataFrame", quietly = TRUE)) "ZarrDataFrame" else NULL))
})

#' The vrMetadata (VoltRon Metadata) Class
#'
#' @slot tile the metadata of tiles
#' @slot molecule the metadata of molecules
#' @slot cell the metadata of cells
#' @slot spot the metadata of spot
#' @slot ROI the metadata of ROI
#'
#' @name vrMetadata-class
#' @rdname vrMetadata-class
#' @exportClass vrMetadata
#'
vrMetadata <- setClass(
  Class = 'vrMetadata',
  slots = c(
    molecule = 'metadata_data',
    cell = 'metadata_data',
    spot = 'metadata_data',
    ROI = 'metadata_data',
    tile = 'metadata_data'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrMetadata',
  definition = function(object) {
    cat("VoltRon Metadata Object \n")
    cat("This object includes: \n")
    lapply(methods::slotNames(object), function(x){
      if(nrow(slot(object, name = x))){
        cat("  ", nrow(slot(object, name = x)), paste0(x, "s"), "\n")
      }
    })
    return(invisible(x = NULL))
  }
)

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
    for(samp in sample_names[seq_len(show_length)]){
      cat(samp, ": \n", sep = "")
      layers <- unique(sample.metadata$Layer[sample.metadata$Sample == samp])
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
      if(length(main.feat) == 1){
        unique_features <- unique_features[c(which(unique_features == main.feat),which(unique_features != main.feat))]
        unique_features[1] <- paste0(unique_features[1], "(Main)") 
      }
      cat("Features:", paste(unique_features, collapse = " "), "\n") 
    }
    
    # return invisible
    return(invisible(x = NULL))
  }
)