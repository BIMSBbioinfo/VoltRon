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
#'
#' @return \code{Metadata<-}: The coordinates updated
#'
#' @rdname Metadata
#' @export Metadata<-
#'
"Metadata<-" <- function(object, ...) {
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
#'
#' @return \code{vrAssayNames<-}: The coordinates updated
#'
#' @rdname vrAssayNames
#' @export vrAssayNames<-
#'
"vrAssayNames<-" <- function(object, ...) {
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
#'
#' @return \code{vrMainAssay<-}: The name of the default assay updated
#'
#' @rdname vrMainAssay
#' @export vrMainAssay<-
#'
"vrMainAssay<-" <- function(object, ...) {
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
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints<-
#'
"vrSpatialPoints<-" <- function(object, ...) {
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
#'
#' @return \code{vrFeatureData<-}: The coordinates updated
#'
#' @rdname vrFeatureData
#' @export vrFeatureData<-
#'
"vrFeatureData<-" <- function(object, ...) {
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

#' Get neighbors
#'
#' get neighbors in an assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname getNeighbors
#' @export getNeighbors
#'
getNeighbors <- function(object, ...) {
  UseMethod(generic = 'getNeighbors', object = object)
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
#'
#' @return \code{vrCoordinates<-}: The coordinates updated
#'
#' @rdname vrCoordinates
#' @export vrCoordinates<-
#'
"vrCoordinates<-" <- function(object, ...) {
  UseMethod(generic = 'vrCoordinates<-', object = object)
}

#' flipCoords
#'
#' Get and set the coordinates of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{flipCoords}: The name of the default assay
#'
#' @rdname flipCoords
#' @export flipCoords
#'
flipCoords <- function(object, ...) {
  UseMethod(generic = 'flipCoords', object = object)
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
#'
#' @return \code{vrSegments<-}: The coordinates updated
#'
#' @rdname vrSegments
#' @export vrSegments<-
#'
"vrSegments<-" <- function(object, ...) {
  UseMethod(generic = 'vrSegments<-', object = object)
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
#'
#' @return \code{vrEmbeddings}: The name of the default assay
#'
#' @rdname vrEmbeddings
#' @export vrEmbeddings<-
#'
"vrEmbeddings<-" <- function(object, ...) {
  UseMethod(generic = 'vrEmbeddings<-', object = object)
}

#' getPCA
#'
#' calculate getPCA of the VoltRon objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{getPCA}: The name of the default assay
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
#' @return \code{getUMAP}: The name of the default assay
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
#'
#' @return \code{vrImages<-}: The coordinates updated
#'
#' @rdname vrImages
#' @export vrImages<-
#'
"vrImages<-" <- function(object, ...) {
  UseMethod(generic = 'vrImages<-', object = object)
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
