#' #' merge
#' #'
#' #' Merge multiple spaceRover objects
#' #'
#' #' @param object An object
#' #' @param ... Arguments passed to other methods
#' #'
#' #' @return \code{merge}: The name of the default assay
#' #'
#' #' @rdname merge
#' #' @export merge
#' #'
#' merge <- function(object, ...) {
#'   UseMethod(generic = 'merge', object = object)
#' }

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

#' Image
#'
#' Get and set images of a spaceRover assay
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
