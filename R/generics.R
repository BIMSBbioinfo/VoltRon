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

#' #' Get Image
#' #'
#' #' Get and set the image of a spatial assay
#' #'
#' #' @param object An object
#' #' @param ... Arguments passed to other methods
#' #'
#' #' @return \code{getImage}: The name of the default assay
#' #'
#' #' @rdname getImage
#' #' @export getImage
#' #'
#' #' @concept data-access
#' #'
#' getImage <- function(object, ...) {
#'   UseMethod(generic = 'getImage', object = object)
#' }
