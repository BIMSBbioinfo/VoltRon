#' @include zzz.R
#'
NULL

#' Find neighbors
#'
#' find neighbors in an assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return Returns object after normalization
#'
#' @rdname findNeighbors
#' @export findNeighbors
#'
findNeighbors <- function(object, ...) {
  UseMethod(generic = 'findNeighbors', object = object)
}

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
#' Given a SpaceRover object, give names of assays of some type, or main assay.
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
#' Given a SpaceRover object, give types of assays of some type, or of the main assay.
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

#' Entities
#'
#' Get and set spatial entities (cells, spots, ROI)
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{Entities}: The name of the default assay
#'
#' @rdname Entities
#' @export Entities
#'
#' @concept data-access
#'
Entities <- function(object, ...) {
  UseMethod(generic = 'Entities', object = object)
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

#' Image
#'
#' Get images of a spaceRover assay
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
