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


#' Get profile specific neighbors
#'
#' get neighbors in an assay given omic profiles
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname getProfileNeighbors
#' @export getProfileNeighbors
#'
getProfileNeighbors <- function(object, ...) {
  UseMethod(generic = 'getProfileNeighbors', object = object)
}

#' changeSampleNames
#'
#' change sample names of VoltRon or other objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname changeSampleNames
#' @export changeSampleNames
#'
#' @noRd
#'
changeSampleNames <- function(object, ...) {
  UseMethod(generic = 'changeSampleNames', object = object)
}

#' changeAssayNames
#'
#' change assay names of VoltRon or other objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname changeAssayNames
#' @export changeAssayNames
#'
#' @noRd
#'
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
#' Given a VoltRon object, give names of assays of some type, or main assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrAssayNames
#' @export vrAssayNames
#'
vrAssayNames <- function(object, ...) {
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
#' Given a VoltRon object, give types of assays of some type, or of the main assay.
#' Here, an assay type is of either cell, spot or ROI.
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
#' Given a VoltRon object, if there are any, get a list of parameters of the assay(s)
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
#' Get and set the main assay
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
# Metadata ####
####

#' Metadata
#'
#' Get the metadata
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname Metadata
#' @export Metadata
#'
Metadata <- function(object, ...) {
  UseMethod(generic = 'Metadata', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value metadata
#'
#' @rdname Metadata
#' @export Metadata<-
#'
"Metadata<-" <- function(object, ..., value) {
  UseMethod(generic = 'Metadata<-', object = object)
}

#' SampleMetadata
#'
#' Get the sample metadata
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname SampleMetadata
#' @export SampleMetadata
#'
SampleMetadata <- function(object, ...) {
  UseMethod(generic = 'SampleMetadata', object = object)
}

####
# Processing ####
####

#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname normalizeData
#' @export
#'
normalizeData <- function(object, ...) {
  UseMethod(generic = 'normalizeData', object = object)
}

####
# Embedding ####
####

#' Get Embedding names
#'
#' Given a VoltRon object, give names of embeddings
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrEmbeddingNames
#' @export vrEmbeddingNames
#'
vrEmbeddingNames <- function(object, ...) {
  UseMethod(generic = 'vrEmbeddingNames', object = object)
}

#' vrEmbeddings
#'
#' Get embeddings of spatial points
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings
#'
vrEmbeddings <- function(object, ...) {
  UseMethod(generic = 'vrEmbeddings', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value embeddings
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings<-
#'
"vrEmbeddings<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrEmbeddings<-', object = object)
}

#' getPCA
#'
#' calculate getPCA of the VoltRon objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname getPCA
#' @export getPCA
#'
getPCA <- function(object, ...) {
  UseMethod(generic = 'getPCA', object = object)
}

#' getUMAP
#'
#' calculate getUMAP of the VoltRon objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname getUMAP
#' @export getUMAP
#'
getUMAP <- function(object, ...) {
  UseMethod(generic = 'getUMAP', object = object)
}

####
# Spatial ####
####

#' vrSpatialPoints
#'
#' Get and set spatial entities (cells, spots, ROI)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints
#'
vrSpatialPoints <- function(object, ...) {
  UseMethod(generic = 'vrSpatialPoints', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value names for spatial points
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints<-
#'
"vrSpatialPoints<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrSpatialPoints<-', object = object)
}

#' Get spatial neighbors
#'
#' get neighbors in an assay given spatial coordinates
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname getSpatialNeighbors
#' @export getSpatialNeighbors
#'
getSpatialNeighbors <- function(object, ...) {
  UseMethod(generic = 'getSpatialNeighbors', object = object)
}

#' vrCoordinates
#'
#' Get and set the coordinates of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrCoordinates
#' @export vrCoordinates
#'
vrCoordinates <- function(object, ...) {
  UseMethod(generic = 'vrCoordinates', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value coordinates
#'
#' @rdname vrCoordinates
#' @export vrCoordinates<-
#'
"vrCoordinates<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrCoordinates<-', object = object)
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

#' vrSegments
#'
#' Get and set the segments of spatian points of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrSegments
#' @export vrSegments
#'
vrSegments <- function(object, ...) {
  UseMethod(generic = 'vrSegments', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value segments
#'
#' @rdname vrSegments
#' @export vrSegments<-
#'
"vrSegments<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrSegments<-', object = object)
}

####
# Graph ####
####

#' vrGraph
#'
#' Get graph from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrGraph
#' @export vrGraph
#'
vrGraph <- function(object, ...) {
  UseMethod(generic = 'vrGraph', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value an igraph object
#'
#' @rdname vrGraph
#' @export vrGraph<-
#'
"vrGraph<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrGraph<-', object = object)
}

#' vrGraphNames
#'
#' Get names of all graphs
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrGraphNames
#' @export vrGraphNames
#'
vrGraphNames <- function(object, ...) {
  UseMethod(generic = 'vrGraphNames', object = object)
}

####
# Image ####
####

#' vrImages
#'
#' Get images of a VoltRon assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrImages
#' @export vrImages
#'
vrImages <- function(object, ...) {
  UseMethod(generic = 'vrImages', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value image
#'
#' @rdname vrImages
#' @export vrImages<-
#'
"vrImages<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrImages<-', object = object)
}

#' vrImageNames
#'
#' Get names of all images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrImageNames
#' @export vrImageNames
#'
vrImageNames <- function(object, ...) {
  UseMethod(generic = 'vrImageNames', object = object)
}

#' vrImageChannelNames
#'
#' Get names of all image channels
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrImageChannelNames
#' @export vrImageChannelNames
#'
vrImageChannelNames <- function(object, ...) {
  UseMethod(generic = 'vrImageChannelNames', object = object)
}

#' vrMainImage
#'
#' Get the main image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrMainImage
#' @export vrMainImage
#'
vrMainImage <- function(object, ...) {
  UseMethod(generic = 'vrMainImage', object = object)
}

#' vrMainImage
#'
#' Set the main image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrMainImage
#' @export vrMainImage<-
#'
"vrMainImage<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrMainImage<-', object = object)
}

#' vrMainChannel
#'
#' Get the main channel name of the image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrMainChannel
#' @export vrMainChannel
#'
vrMainChannel <- function(object, ...) {
  UseMethod(generic = 'vrMainChannel', object = object)
}

#' vrMainChannel
#'
#' Set the main channel name of the image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname vrMainChannel
#' @export vrMainChannel<-
#'
"vrMainChannel<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrMainChannel<-', object = object)
}


#' resizeImage
#'
#' Resizing Magick images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @rdname resizeImage
#' @export resizeImage
#'
resizeImage <- function(object, ...) {
  UseMethod(generic = 'resizeImage', object = object)
}

#' modulateImage
#'
#' Modulating Magick images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#'
#' @rdname modulateImage
#' @export modulateImage
#'
modulateImage <- function(object, ...) {
  UseMethod(generic = 'modulateImage', object = object)
}

#' combineChannels
#'
#' Combining channels into novel channels of the same image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#'
#' @rdname combineChannels
#' @export combineChannels
#'
combineChannels <- function(object, ...) {
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
