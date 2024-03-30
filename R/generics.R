#' @include zzz.R
#'
NULL

####
# General ####
####

#' vrFeatures
#'
#' Get features from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrFeatures
#' @export vrFeatures
#'
vrFeatures <- function(object, ...) {
  UseMethod(generic = 'vrFeatures', object = object)
}

#' Feature Data
#'
#' Get and set feature data from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrFeatureData
#' @export vrFeatureData
#'
vrFeatureData <- function(object, ...) {
  UseMethod(generic = 'vrFeatureData', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value feature data
#'
#' @rdname vrFeatureData
#' @export vrFeatureData<-
#'
"vrFeatureData<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrFeatureData<-', object = object)
}

#' Get spatially variable feature
#'
#' get spatially variable features in an assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname getFeatures
#' @export getFeatures
#'
getFeatures <- function(object, ...) {
  UseMethod(generic = 'getFeatures', object = object)
}

#' vrData
#'
#' Get data from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrData
#' @export vrData
#'
vrData <- function(object, ...) {
  UseMethod(generic = 'vrData', object = object)
}

#' changeSampleNames
#'
#' change sample names of VoltRon or other objects
#'
#' @param object a VoltRon or vrMetadata object
#'
#' @rdname changeSampleNames
#' @export changeSampleNames
#'
#' @noRd
changeSampleNames <- function(object, samples = NULL, sample_metadata_table) {
  UseMethod(generic = 'changeSampleNames', object = object)
}

#' changeAssayNames
#'
#' change assay names of VoltRon or other objects
#'
#' @param object a VoltRon, vrSample or vrLayer object
#'
#' @rdname changeAssayNames
#' @export changeAssayNames
#'
#' @noRd
changeAssayNames <- function(object, ...) {
  UseMethod(generic = 'changeAssayNames', object = object)
}

####
# Assay ####
####

#' Add Assay
#'
#' add assay to the object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname addAssay
#' @export addAssay
#'
addAssay <- function(object, ...) {
  UseMethod(generic = 'addAssay', object = object)
}

#' Get Assay names
#'
#' Given a VoltRon, vrMetadata or vrAssay object, get/set names of assays.
#'
#' @param object An object 
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#'
#' @rdname vrAssayNames
#' @export vrAssayNames
#'
vrAssayNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrAssayNames', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value assay name
#'
#' @rdname vrAssayNames
#' @export vrAssayNames<-
#'
"vrAssayNames<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrAssayNames<-', object = object)
}

#' Get Assay types
#' 
#' Given a VoltRon or vrAssay object, get types of assays. 
#' Here, an assay type is of either tile, molecule, cell, spot or ROI.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrAssayTypes
#' @export vrAssayTypes
#'
vrAssayTypes <- function(object, ...) {
  UseMethod(generic = 'vrAssayTypes', object = object)
}

#' Get Assay types
#'
#' Given a vrAssay object, if there are any, get a list of parameters of the assay(s)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrAssayParams
#' @export vrAssayParams
#'
vrAssayParams <- function(object, ...) {
  UseMethod(generic = 'vrAssayParams', object = object)
}

#' Main Assay
#'
#' Get and set the main assay of a VoltRon object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrMainAssay
#' @export vrMainAssay
#'
vrMainAssay <- function(object, ...) {
  UseMethod(generic = 'vrMainAssay', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value assay name
#'
#' @rdname vrMainAssay
#' @export vrMainAssay<-
#'
"vrMainAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrMainAssay<-', object = object)
}

####
# Sample ####
####

#' Get Sample names
#'
#' Given a vrMetadata object, give names of samples
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrSampleNames
#' @export vrSampleNames
#'
vrSampleNames <- function(object, ...) {
  UseMethod(generic = 'vrSampleNames', object = object)
}

####
# Metadata ####
####

#' Metadata
#'
#' Get the metadata of a VoltRon object.
#'
#' @param object a VoltRon object
#'
#' @rdname Metadata
#' @export Metadata
#'
Metadata <- function(object, assay = NULL, type = NULL) {
  UseMethod(generic = 'Metadata', object = object)
}

#' @param value new metadata
#' 
#' @rdname Metadata
#' @export Metadata<-
#'
"Metadata<-" <- function(object, assay = NULL, type = NULL, value) {
  UseMethod(generic = 'Metadata<-', object = object)
}

####
# Processing ####
####

#' Normalize Data
#' 
#' Given a VoltRon or vrAssay object, normalize the raw count data.
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname normalizeData
#' @export normalizeData
#'
normalizeData <- function(object, assay = NULL, method = "LogNorm", desiredQuantile = 0.9, 
                          scale = 0.2, sizefactor = 10000) {
  UseMethod(generic = 'normalizeData', object = object)
}

####
# Embedding ####
####

#' Get Embedding names
#'
#' Given a VoltRon or vrAssay object, give names of embeddings
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrEmbeddingNames
#' @export vrEmbeddingNames
#'
vrEmbeddingNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrEmbeddingNames', object = object)
}

#' vrEmbeddings
#'
#' Given a VoltRon or vrAssay object, get embeddings of spatial points
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings
#'
vrEmbeddings <- function(object, assay = NULL, type = "pca", dims = 1:30) {
  UseMethod(generic = 'vrEmbeddings', object = object)
}

#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings<-
#'
"vrEmbeddings<-" <- function(object, assay = NULL, type = "pca", overwrite = FALSE, value) {
  UseMethod(generic = 'vrEmbeddings<-', object = object)
}

####
# Spatial ####
####

#' vrSpatialPoints
#'
#' Get and set spatial entities (tile, molecules, cells, spots and ROI).
#'
#' @param object An object
#' @param ... arguments passed to other methods
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints
#'
vrSpatialPoints <- function(object, assay = NULL) {
  UseMethod(generic = 'vrSpatialPoints', object = object)
}

#' @param object An object
#' @param value names for spatial points
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints<-
#'
"vrSpatialPoints<-" <- function(object, value) {
  UseMethod(generic = 'vrSpatialPoints<-', object = object)
}

#' vrCoordinates
#'
#' Given a VoltRon, vrAssay or vrImage object, get and set the coordinates of assays and coordinate systems
#'
#' @param object a VoltRon, vrAssay or vrImage object
#'
#' @rdname vrCoordinates
#' @export vrCoordinates
#'
vrCoordinates <- function(object, assay = NULL, image_name = NULL, reg = FALSE) {
  UseMethod(generic = 'vrCoordinates', object = object)
}

#' @param value new coordinates of spatial points
#'
#' @rdname vrCoordinates
#' @export vrCoordinates<-
#'
"vrCoordinates<-" <- function(object, image_name = NULL, reg = FALSE, value) {
  UseMethod(generic = 'vrCoordinates<-', object = object)
}

#' vrSegments
#'
#' Given a VoltRon, vrAssay or vrImage object, get and set the segment coordinates of assays and coordinate systems
#'
#' @param object a VoltRon, vrAssay or vrImage object
#'
#' @rdname vrSegments
#' @export vrSegments
#'
vrSegments <- function(object, assay = NULL, image_name = NULL, reg = FALSE) {
  UseMethod(generic = 'vrSegments', object = object)
}

#' @param value new segment coordinates of spatial points
#'
#' @rdname vrSegments
#' @export vrSegments<-
#'
"vrSegments<-" <- function(object, image_name = NULL, reg = FALSE, value) {
  UseMethod(generic = 'vrSegments<-', object = object)
}

#' flipCoordinates
#'
#' Flip the coordinates of the spatial points in the y axis direction
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname flipCoordinates
#' @export flipCoordinates
#'
flipCoordinates <- function(object, ...) {
  UseMethod(generic = 'flipCoordinates', object = object)
}

####
# Graph ####
####

#' vrGraph
#'
#' Get graph of a VoltRon object
#'
#' @param object A VoltRon object
#'
#' @rdname vrGraph
#' @export vrGraph
#'
vrGraph <- function(object, assay = NULL, graph.type = "kNN") {
  UseMethod(generic = 'vrGraph', object = object)
}

#' @param value new graph
#'
#' @rdname vrGraph
#' @export vrGraph<-
#'
"vrGraph<-" <- function(object, graph.type = "kNN", value) {
  UseMethod(generic = 'vrGraph<-', object = object)
}

#' vrGraphNames
#'
#' Get names of all graphs
#'
#' @param object a VoltRon object
#'
#' @rdname vrGraphNames
#' @export vrGraphNames
#'
vrGraphNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrGraphNames', object = object)
}

####
# Image ####
####

#' vrImages
#'
#' Get images of VoltRon objects
#'
#' @param object A VoltRon, vrAssay or vrImage object
#'
#' @rdname vrImages
#' @export vrImages
#'
vrImages <- function(object, assay = NULL, name = NULL, reg = FALSE, channel = NULL, as.raster = FALSE, scale.perc = 100) {
  UseMethod(generic = 'vrImages', object = object)
}

#' @param value new image
#'
#' @rdname vrImages
#' @export vrImages<-
#'
"vrImages<-" <- function(object, name = NULL, reg = FALSE, channel = NULL, value) {
  UseMethod(generic = 'vrImages<-', object = object)
}

#' vrImageNames
#'
#' Get names of all images
#'
#' @param object A VoltRon or vrAssay object
#'
#' @rdname vrImageNames
#' @export vrImageNames
#'
vrImageNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrImageNames', object = object)
}

#' vrImageChannelNames
#'
#' Get names of all image channels
#'
#' @param object A VoltRon, vrAssay or vrImage object
#'
#' @rdname vrImageChannelNames
#' @export vrImageChannelNames
#'
vrImageChannelNames <- function(object, assay = NULL, name = NULL) {
  UseMethod(generic = 'vrImageChannelNames', object = object)
}

#' vrMainImage
#'
#' Get the main image
#'
#' @param object A VoltRon or vrAssay object
#'
#' @rdname vrMainImage
#' @export vrMainImage
#'
vrMainImage <- function(object, assay = NULL) {
  UseMethod(generic = 'vrMainImage', object = object)
}

#' vrMainImage
#'
#' Set the main image
#'
#' @param value the name of main image
#'
#' @rdname vrMainImage
#' @export vrMainImage<-
#'
"vrMainImage<-" <- function(object, value) {
  UseMethod(generic = 'vrMainImage<-', object = object)
}

#' vrMainChannel
#'
#' Get the main channel name of the image
#'
#' @param object a vrAssay or vrImage object
#'
#' @rdname vrMainChannel
#' @export vrMainChannel
#'
vrMainChannel <- function(object, name = NULL) {
  UseMethod(generic = 'vrMainChannel', object = object)
}

#' vrMainChannel
#'
#' Set the main channel name of the image
#'
#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @export vrMainChannel<-
#'
"vrMainChannel<-" <- function(object, name = NULL, value) {
  UseMethod(generic = 'vrMainChannel<-', object = object)
}


#' resizeImage
#'
#' Resizing Magick images
#'
#' @param object A VoltRon, vrAssay or vrImage object
#'
#' @rdname resizeImage
#' @export resizeImage
#'
resizeImage <- function(object, assay = NULL, name = NULL, reg = FALSE, size) {
  UseMethod(generic = 'resizeImage', object = object)
}

#' modulateImage
#'
#' Modulating Magick images
#'
#' @param object A VoltRon, vrAssay or vrImage object
#'
#' @rdname modulateImage
#' @export modulateImage
#'
modulateImage <- function(object, assay = NULL, name = NULL, reg = FALSE, channel = NULL, 
                          brightness = 100, saturation = 100, hue = 100, force = FALSE) {
  UseMethod(generic = 'modulateImage', object = object)
}

#' combineChannels
#'
#' Combining channels into novel channels of the same image
#'
#' @param object A VoltRon, vrAssay or vrImage object
#'
#' @rdname combineChannels
#' @export combineChannels
#'
combineChannels <- function(object, assay = NULL, name = NULL, reg = FALSE, channels = NULL, 
                            colors = NULL, channel_key = "combined") {
  UseMethod(generic = 'combineChannels', object = object)
}

####
# Conversion ####
####

#' as.Seurat
#'
#' Generic methods for conversion into a Seurat object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Seurat
#' @export as.Seurat
#'
as.Seurat <- function(object, ...) {
  UseMethod(generic = 'as.Seurat', object = object)
}

#' as.AnnData
#'
#' Generic methods for conversion into a AnnData object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.AnnData
#' @export as.AnnData
#'
as.AnnData <- function(object, ...) {
  UseMethod(generic = 'as.AnnData', object = object)
}

#' as.Zarr
#'
#' Generic methods for conversion into a Seurat object
#'
#' @param object An object
#' @param out_path output path to ome.zarr
#' @param image_id image name
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Zarr
#' @export as.Zarr
#'
as.Zarr <- function(object, out_path, image_id, ...) {
  UseMethod(generic = 'as.Zarr', object = object)
}


#' as.Giotto
#'
#' Generic methods for conversion into a Giotto object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.Giotto
#' @export as.Giotto
#'
as.Giotto <- function(object, ...) {
  UseMethod(generic = 'as.Giotto', object = object)
}

#' as.VoltRon
#'
#' Generic methods for conversion into a VoltRon object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.VoltRon
#' @export as.VoltRon
#'
as.VoltRon <- function(object, ...) {
  UseMethod(generic = 'as.VoltRon', object = object)
}
