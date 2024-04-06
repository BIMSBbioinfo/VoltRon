#' @include zzz.R
#'
NULL

####
# General ####
####

#' vrFeatures
#'
#' Get names of the features.
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrFeatures
#' @export vrFeatures
#' @order 1
vrFeatures <- function(object, assay = NULL) {
  UseMethod(generic = 'vrFeatures', object = object)
}

#' Feature Data
#'
#' Get and set feature data from the main.assay
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrFeatureData
#' @export vrFeatureData
#' @order 1
vrFeatureData <- function(object, assay = NULL) {
  UseMethod(generic = 'vrFeatureData', object = object)
}

#' @rdname vrFeatureData
#' @export vrFeatureData<-
#' @noRd
"vrFeatureData<-" <- function(object, assay = NULL, value) {
  UseMethod(generic = 'vrFeatureData<-', object = object)
}

#' getFeatures
#'
#' Get variable features of assays
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname getFeatures
#' @export getFeatures
#'
getFeatures <- function(object, assay = NULL, max.count = 1, n = 3000) {
  UseMethod(generic = 'getFeatures', object = object)
}

#' vrData
#'
#' Get data of assays
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrData
#' @export vrData
#' @order 1
vrData <- function(object, assay = NULL, features = NULL, norm = FALSE) {
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
#' add a new assay (vrAssay object) to a VoltRon object
#'
#' @param object a VoltRon object
#' @param assay a vrAssay object
#' @param metadata a predefined metadata
#' @param assay_name assay name of the new added assay
#' @param sample sample name
#' @param layer layer name
#'
#' @rdname addAssay
#' @export addAssay
#'
addAssay <- function(object, assay, metadata = NULL, assay_name, sample = "Sample1", layer = "Section1") {
  UseMethod(generic = 'addAssay', object = object)
}

#' Get Assay names
#'
#' Given a VoltRon, vrMetadata or vrAssay object, get/set names of assays.
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrAssayNames
#' @export vrAssayNames
#' @order 1
vrAssayNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrAssayNames', object = object)
}

#' @rdname vrAssayNames
#' @export vrAssayNames<-
#' @noRd
"vrAssayNames<-" <- function(object, value) {
  UseMethod(generic = 'vrAssayNames<-', object = object)
}

#' Get assay types
#' 
#' Given a VoltRon or vrAssay object, get types of assays. 
#' Here, an assay type is of either tile, molecule, cell, spot or ROI.
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrAssayTypes
#' @export vrAssayTypes
#' @order 1
vrAssayTypes <- function(object, assay = NULL) {
  UseMethod(generic = 'vrAssayTypes', object = object)
}

####
# Sample ####
####

#' Get Sample names
#'
#' Given a vrMetadata object, give names of samples
#'
#' @param object a vrMetadata object
#'
#' @rdname vrSampleNames
#' @export vrSampleNames
#'
vrSampleNames <- function(object) {
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
#' @noRd
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
#' @order 1
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
#' @order 1
vrEmbeddings <- function(object, assay = NULL, type = "pca", dims = 1:30) {
  UseMethod(generic = 'vrEmbeddings', object = object)
}

#' @param value new embedding data
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings<-
#' @noRd
"vrEmbeddings<-" <- function(object, assay = NULL, type = "pca", overwrite = FALSE, value) {
  UseMethod(generic = 'vrEmbeddings<-', object = object)
}

####
# Spatial ####
####

#' vrSpatialPoints
#'
#' Get and set spatial entities.
#'
#' @param object a VoltRon, vrSample, vrLayer, vrAssay or vrImage object
#' @param ... arguments passed to other methods
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints
#' @order 1
#'
#' @examples
#' vrSpatialPoints(visium_data)
#' vrSpatialPoints(visium_data, assay = "Visium")
#' vrSpatialPoints(visium_data, assay = "Assay1")
vrSpatialPoints <- function(object, assay = NULL) {
  UseMethod(generic = 'vrSpatialPoints', object = object)
}

#' @param value names for spatial points
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints<-
#' @noRd
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
#' @order 1
#'
vrCoordinates <- function(object, assay = NULL, image_name = NULL, reg = FALSE) {
  UseMethod(generic = 'vrCoordinates', object = object)
}

#' @rdname vrCoordinates
#' @export vrCoordinates<-
#' @noRd
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
#' @order 1
vrSegments <- function(object, assay = NULL, image_name = NULL, reg = FALSE) {
  UseMethod(generic = 'vrSegments', object = object)
}

#' @rdname vrSegments
#' @export vrSegments<-
#' @noRd
"vrSegments<-" <- function(object, image_name = NULL, reg = FALSE, value) {
  UseMethod(generic = 'vrSegments<-', object = object)
}

#' flipCoordinates
#'
#' Flip the coordinates of spatial points in the y axis direction. 
#'
#' @param object a VoltRon, vrAssay or vrImage object
#'
#' @rdname flipCoordinates
#' @export flipCoordinates
#' @order 1
flipCoordinates <- function(object, assay = NULL, image_name = NULL, ...) {
  UseMethod(generic = 'flipCoordinates', object = object)
}

####
# Image ####
####

#' vrImages
#'
#' Get images of VoltRon objects
#'
#' @param object a VoltRon, vrAssay or vrImage object
#'
#' @rdname vrImages
#' @export vrImages
#' @order 1
vrImages <- function(object, assay = NULL, name = NULL, reg = FALSE, channel = NULL, as.raster = FALSE, scale.perc = 100) {
  UseMethod(generic = 'vrImages', object = object)
}

#' @rdname vrImages
#' @export vrImages<-
#' @noRd
"vrImages<-" <- function(object, name = NULL, reg = FALSE, channel = NULL, value) {
  UseMethod(generic = 'vrImages<-', object = object)
}

#' vrImageNames
#'
#' Get names of all images
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrImageNames
#' @export vrImageNames
#' 
vrImageNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrImageNames', object = object)
}

#' vrSpatialNames
#'
#' Get names of all spatial systems
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrSpatialNames
#' @export vrSpatialNames
#' 
vrSpatialNames <- function(object, assay = NULL) {
  UseMethod(generic = 'vrSpatialNames', object = object)
}

#' vrImageChannelNames
#'
#' Get names of all image channels
#'
#' @param object a VoltRon, vrAssay or vrImage object
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
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrMainImage
#' @export vrMainImage
#' @order 1
vrMainImage <- function(object, assay = NULL) {
  UseMethod(generic = 'vrMainImage', object = object)
}

#' vrMainSpatial
#'
#' Get the main spatial system name
#'
#' @param object a VoltRon or vrAssay object
#'
#' @rdname vrMainSpatial
#' @export vrMainSpatial
#' @order 1
vrMainSpatial <- function(object, assay = NULL) {
  UseMethod(generic = 'vrMainSpatial', object = object)
}

#' vrMainImage
#'
#' Set the main image
#'
#' @param value the name of main image
#'
#' @rdname vrMainImage
#' @export vrMainImage<-
#' @noRd
"vrMainImage<-" <- function(object, value) {
  UseMethod(generic = 'vrMainImage<-', object = object)
}

#' vrMainSpatial
#'
#' Set the main image
#'
#' @param value the name of main image
#'
#' @rdname vrMainSpatial
#' @export vrMainSpatial<-
#' @noRd
"vrMainSpatial<-" <- function(object, value) {
  UseMethod(generic = 'vrMainSpatial<-', object = object)
}

#' vrMainChannel
#'
#' Get and set the main channel name of the image
#'
#' @param object a vrAssay or vrImage object
#'
#' @rdname vrMainChannel
#' @export vrMainChannel
#' @order 1
vrMainChannel <- function(object, name = NULL) {
  UseMethod(generic = 'vrMainChannel', object = object)
}

#' @param value the name of main channel
#'
#' @rdname vrMainChannel
#' @export vrMainChannel<-
#' @noRd
"vrMainChannel<-" <- function(object, name = NULL, value) {
  UseMethod(generic = 'vrMainChannel<-', object = object)
}


#' resizeImage
#'
#' Resizing Magick images
#'
#' @param object a VoltRon, vrAssay or vrImage object
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
#' @param object a VoltRon, vrAssay or vrImage object
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
#' @param object a VoltRon, vrAssay or vrImage object
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

#' as.VoltRon
#'
#' Generic methods for conversion into a VoltRon object
#'
#' @param object a VoltRon object
#' @param ... Arguments passed to other methods
#'
#' @rdname as.VoltRon
#' @export as.VoltRon
as.VoltRon <- function(object, ...) {
  UseMethod(generic = 'as.VoltRon', object = object)
}

#' as.Zarr
#'
#' Generic methods to save VoltRon or magick-image objects as zarr files
#'
#' @param object a VoltRon or magick-image object
#' @param out_path output path to zarr file
#' @param image_id image name
#'
#' @rdname as.Zarr
#' @export as.Zarr
#'
as.Zarr <- function(object, out_path, image_id) {
  UseMethod(generic = 'as.Zarr', object = object)
}
