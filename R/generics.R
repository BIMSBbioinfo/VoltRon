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
#' @concept data-access
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
#' @concept data-access
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
#' @concept data-access
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
#' @return \code{AddAssay}: The name of the default assay
#'
#' @rdname AddAssay
#' @export AddAssay
#'
#' @concept data-access
#'
AddAssay <- function(object, ...) {
  UseMethod(generic = 'AddAssay', object = object)
}

#' Get Assay names
#'
#' Given a VoltRon object, give names of assays of some type, or main assay.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{AssayNames}: The name of the default assay
#'
#' @rdname AssayNames
#' @export AssayNames
#'
#' @concept data-access
#'
AssayNames <- function(object, ...) {
  UseMethod(generic = 'AssayNames', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{AssayNames<-}: The coordinates updated
#'
#' @rdname AssayNames
#' @export AssayNames<-
#'
#' @concept data-access
#'
"AssayNames<-" <- function(object, ...) {
  UseMethod(generic = 'AssayNames<-', object = object)
}

#' Get Assay types
#'
#' Given a VoltRon object, give types of assays of some type, or of the main assay.
#' Here, an assay type is of either cell, spot or ROI.
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{AssayTypes}: The name of the default assay
#'
#' @rdname AssayTypes
#' @export AssayTypes
#'
#' @concept data-access
#'
AssayTypes <- function(object, ...) {
  UseMethod(generic = 'AssayTypes', object = object)
}

#' Main Assay
#'
#' Get and set the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{MainAssay}: The name of the default assay
#'
#' @rdname MainAssay
#' @export MainAssay
#'
#' @concept data-access
#'
MainAssay <- function(object, ...) {
  UseMethod(generic = 'MainAssay', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{MainAssay<-}: The name of the default assay updated
#'
#' @rdname MainAssay
#' @export MainAssay<-
#'
#' @concept data-access
#'
"MainAssay<-" <- function(object, ...) {
  UseMethod(generic = 'MainAssay<-', object = object)
}

#' vrSpatialPoints
#'
#' Get and set spatial entities (cells, spots, ROI)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{vrSpatialPoints}: The name of the default assay
#'
#' @rdname vrSpatialPoints
#' @export vrSpatialPoints
#'
#' @concept data-access
#'
vrSpatialPoints <- function(object, ...) {
  UseMethod(generic = 'vrSpatialPoints', object = object)
}

#' Features
#'
#' Get features from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Features}: The name of the default assay
#'
#' @rdname Features
#' @export Features
#'
#' @concept data-access
#'
Features <- function(object, ...) {
  UseMethod(generic = 'Features', object = object)
}

#' Feature Data
#'
#' Get and set feature data from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{FeatureData}: The name of the default assay
#'
#' @rdname FeatureData
#' @export FeatureData
#'
#' @concept data-access
#'
FeatureData <- function(object, ...) {
  UseMethod(generic = 'FeatureData', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{FeatureData<-}: The coordinates updated
#'
#' @rdname FeatureData
#' @export FeatureData<-
#'
#' @concept data-access
#'
"FeatureData<-" <- function(object, ...) {
  UseMethod(generic = 'FeatureData<-', object = object)
}

#' Data
#'
#' Get data from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Data}: The name of the default assay
#'
#' @rdname Data
#' @export Data
#'
#' @concept data-access
#'
Data <- function(object, ...) {
  UseMethod(generic = 'Data', object = object)
}

#' Graph
#'
#' Get graph from the main.assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Graph}: The name of the default assay
#'
#' @rdname Graph
#' @export Graph
#'
#' @concept data-access
#'
Graph <- function(object, ...) {
  UseMethod(generic = 'Graph', object = object)
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
#' @rdname getSpatialFeatures
#' @export getSpatialFeatures
#'
getSpatialFeatures <- function(object, ...) {
  UseMethod(generic = 'getSpatialFeatures', object = object)
}


#' Coordinates
#'
#' Get and set the coordinates of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Coordinates}: The name of the default assay
#'
#' @rdname Coordinates
#' @export Coordinates
#'
#' @concept data-access
#'
Coordinates <- function(object, ...) {
  UseMethod(generic = 'Coordinates', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Coordinates<-}: The coordinates updated
#'
#' @rdname Coordinates
#' @export Coordinates<-
#'
#' @concept data-access
#'
"Coordinates<-" <- function(object, ...) {
  UseMethod(generic = 'Coordinates<-', object = object)
}

#' Segments
#'
#' Get and set the segments of spatian points of the main assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Segments}: The name of the default assay
#'
#' @rdname Segments
#' @export Segments
#'
#' @concept data-access
#'
Segments <- function(object, ...) {
  UseMethod(generic = 'Segments', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Segments<-}: The coordinates updated
#'
#' @rdname Segments
#' @export Segments<-
#'
#' @concept data-access
#'
"Segments<-" <- function(object, ...) {
  UseMethod(generic = 'Segments<-', object = object)
}

#' Distances
#'
#' Get distances between spatial points using coordinates
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Distances}: The name of the default assay
#'
#' @rdname Distances
#' @export Distances
#'
#' @concept data-access
#'
Distances <- function(object, ...) {
  UseMethod(generic = 'Distances', object = object)
}

#' Embeddings
#'
#' Get embeddings of spatial points
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Embeddings}: The name of the default assay
#'
#' @rdname Embeddings
#' @export Embeddings
#'
#' @concept data-access
#'
Embeddings <- function(object, ...) {
  UseMethod(generic = 'Embeddings', object = object)
}

#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Embeddings}: The name of the default assay
#'
#' @rdname Embeddings
#' @export Embeddings<-
#'
#' @concept data-access
#'
"Embeddings<-" <- function(object, ...) {
  UseMethod(generic = 'Embeddings<-', object = object)
}

#' PCA
#'
#' calculate PCA of the VoltRon objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{PCA}: The name of the default assay
#'
#' @rdname PCA
#' @export PCA
#'
#' @concept data-access
#'
PCA <- function(object, ...) {
  UseMethod(generic = 'PCA', object = object)
}

#' UMAP
#'
#' calculate UMAP of the VoltRon objects
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{UMAP}: The name of the default assay
#'
#' @rdname UMAP
#' @export UMAP
#'
#' @concept data-access
#'
UMAP <- function(object, ...) {
  UseMethod(generic = 'UMAP', object = object)
}

#' Image
#'
#' Get images of a VoltRon assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Image}: The name of the default assay
#'
#' @rdname Image
#' @export Image
#'
#' @concept data-access
#'
Image <- function(object, ...) {
  UseMethod(generic = 'Image', object = object)
}

#' ResizeImage
#'
#' Resizing Magick images
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{ResizeImage}: The name of the default assay
#'
#' @rdname ResizeImage
#' @export ResizeImage
#'
#' @concept data-access
#'
ResizeImage <- function(object, ...) {
  UseMethod(generic = 'ResizeImage', object = object)
}
