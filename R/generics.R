#' @include zzz.R
#'
NULL

#' Normalize Data
#'
#' Normalize the count data present in a given assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname normalizeData
#' @export normalizeData
#'
normalizeData <- function(object, ...) {
  UseMethod(generic = 'normalizeData', object = object)
}

#' Metadata
#'
#' Get the metadata
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Metadata}: The name of the default assay
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
#' @return \code{Metadata<-}: The coordinates updated
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
#' @return \code{SampleMetadata}: The name of the default assay
#'
#' @rdname SampleMetadata
#' @export SampleMetadata
#'
SampleMetadata <- function(object, ...) {
  UseMethod(generic = 'SampleMetadata', object = object)
}

#' Add Assay
#'
#' add assay to the object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{addAssay}: The name of the default assay
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
#' @return \code{vrAssayNames}: The name of the default assay
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
#' @return \code{vrAssayNames<-}: The coordinates updated
#'
#' @rdname vrAssayNames
#' @export vrAssayNames<-
#'
"vrAssayNames<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrAssayNames<-', object = object)
}

#' changeSampleNames
#'
#' change sample names of VoltRon or other objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{changeSampleNames}: The name of the default assay
#'
#' @rdname changeSampleNames
#' @export changeSampleNames
#'
changeSampleNames <- function(object, ...) {
  UseMethod(generic = 'changeSampleNames', object = object)
}

#' Get Assay types
#'
#' Given a VoltRon object, give types of assays of some type, or of the main assay.
#' Here, an assay type is of either cell, spot or ROI.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrAssayTypes}: The name of the default assay
#'
#' @rdname vrAssayTypes
#' @export vrAssayTypes
#'
vrAssayTypes <- function(object, ...) {
  UseMethod(generic = 'vrAssayTypes', object = object)
}

#' Main Assay
#'
#' Get and set the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrMainAssay}: The name of the default assay
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
#' @return \code{vrMainAssay<-}: The name of the default assay updated
#'
#' @rdname vrMainAssay
#' @export vrMainAssay<-
#'
"vrMainAssay<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrMainAssay<-', object = object)
}

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
#' @return \code{vrFeatureData}: The name of the default assay
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

#' vrData
#'
#' Get data from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrData}: The name of the default assay
#'
#' @rdname vrData
#' @export vrData
#'
vrData <- function(object, ...) {
  UseMethod(generic = 'vrData', object = object)
}

#' vrGraph
#'
#' Get graph from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrGraph}: The name of the default assay
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
#' @return \code{vrGraph<-}: graph updated
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
#' @return \code{vrGraphNames}: The names of all graphs
#'
#' @rdname vrGraphNames
#' @export vrGraphNames
#'
vrGraphNames <- function(object, ...) {
  UseMethod(generic = 'vrGraphNames', object = object)
}


#' Get profile specific neighbors
#'
#' get neighbors in an assay given omic profiles
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname getProfileNeighbors
#' @export getProfileNeighbors
#'
getProfileNeighbors <- function(object, ...) {
  UseMethod(generic = 'getProfileNeighbors', object = object)
}

#' Get spatial neighbors
#'
#' get neighbors in an assay given spatial coordinates
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname getSpatialNeighbors
#' @export getSpatialNeighbors
#'
getSpatialNeighbors <- function(object, ...) {
  UseMethod(generic = 'getSpatialNeighbors', object = object)
}

#' Get spatially variable feature
#'
#' get spatially variable features in an assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname getFeatures
#' @export getFeatures
#'
getFeatures <- function(object, ...) {
  UseMethod(generic = 'getFeatures', object = object)
}


#' vrCoordinates
#'
#' Get and set the coordinates of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrCoordinates}: The name of the default assay
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
#' @return \code{vrCoordinates<-}: The coordinates updated
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
#' @return \code{flipCoordinates}: The name of the default assay
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
#' @return \code{vrSegments}: The name of the default assay
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
#' @return \code{vrSegments<-}: The coordinates updated
#'
#' @rdname vrSegments
#' @export vrSegments<-
#'
"vrSegments<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrSegments<-', object = object)
}

#' vrSubcellular
#'
#' Get and set the subcellular data
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrSubcellular}: The name of the default assay
#'
#' @rdname vrSubcellular
#' @export vrSubcellular
#'
vrSubcellular <- function(object, ...) {
  UseMethod(generic = 'vrSubcellular', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#' @param value subcellular data
#'
#' @return \code{vrSubcellular<-}: The coordinates updated
#'
#' @rdname vrSubcellular
#' @export vrSubcellular<-
#'
"vrSubcellular<-" <- function(object, ..., value) {
  UseMethod(generic = 'vrSubcellular<-', object = object)
}

#' vrDistances
#'
#' Get distances between spatial points using coordinates
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrDistances}: The name of the default assay
#'
#' @rdname vrDistances
#' @export vrDistances
#'
vrDistances <- function(object, ...) {
  UseMethod(generic = 'vrDistances', object = object)
}

#' vrEmbeddings
#'
#' Get embeddings of spatial points
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrEmbeddings}: The name of the default assay
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

#' vrImages
#'
#' Get images of a VoltRon assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrImages}: The name of the default assay
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
#' @return \code{vrImageNames}: The name of the default assay
#'
#' @rdname vrImageNames
#' @export vrImageNames
#'
vrImageNames <- function(object, ...) {
  UseMethod(generic = 'vrImageNames', object = object)
}

#' vrMainImage
#'
#' Get the main image
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrMainImage}: The name of the default assay
#'
#' @rdname vrMainImage
#' @export vrMainImage
#'
vrMainImage <- function(object, ...) {
  UseMethod(generic = 'vrMainImage', object = object)
}

#' resizeImage
#'
#' Resizing Magick images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{resizeImage}: The name of the default assay
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
#' @return \code{modulateImage}: The name of the default assay
#'
#' @rdname modulateImage
#' @export modulateImage
#'
modulateImage <- function(object, ...) {
  UseMethod(generic = 'modulateImage', object = object)
}

#' as.Seurat
#'
#' Generic methods for conversion into a Seurat object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{as.Seurat}: The name of the default assay
#'
#' @rdname as.Seurat
#' @export as.Seurat
#'
as.Seurat <- function(object, ...) {
  UseMethod(generic = 'as.Seurat', object = object)
}

#' as.Giotto
#'
#' Generic methods for conversion into a Giotto object
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{as.Giotto}: The name of the default assay
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
#' @return \code{as.VoltRon}: The name of the default assay
#'
#' @rdname as.VoltRon
#' @export as.VoltRon
#'
as.VoltRon <- function(object, ...) {
  UseMethod(generic = 'as.VoltRon', object = object)
}
