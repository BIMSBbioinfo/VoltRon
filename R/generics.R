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

#' getImage
#'
#' Get and set images of a spaceRover assay
#'
#' @param object An object
#' @param ... Arguments passed to other methods
#'
#' @return \code{getImage}: The name of the default assay
#'
#' @rdname getImage
#' @export getImage
#'
#' @concept data-access
#'
getImage <- function(object, ...) {
  UseMethod(generic = 'getImage', object = object)
}
