fill.na <- function(x, i = 5) {
  if (is.na(x)[i]) {
    return(round(mean(x, na.rm = TRUE), 0))
  }
  else {
    return(round(x[i], 0))
  }
}

#' slotApply
#'
#' apply to slots
#'
#' @param x object
#' @param FUN function
#' @param ... arguments passed to \code{FUN}
#'
#' @importFrom methods slot
#' @export
#'
slotApply <- function(x,FUN,...){
  cl <- class(x)
  result <- list()
  for(i in methods::slotNames(cl)){
    result[[i]] <- FUN(methods::slot(x,i),...)
  }
  result
}

#' slotToList
#'
#' slot to list
#'
#' @param x object
#'
#' @importFrom methods slot slotNames
#' @export
#'
slotToList <- function(x){
  returnlist <- list()
  namesslot <- methods::slotNames(x)
  for(cur_slot in namesslot)
    returnlist[[cur_slot]] <- methods::slot(x, name = cur_slot)
  returnlist
}

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
}

jaccard_similarity <- function(mat) {
  matinv <- 1 - mat
  return((mat %*% mat)/(nrow(mat) - (matinv %*% matinv)))
}

#' make_css
#'
#' make_css from \code{tableHTML} package
#'
#' @param ... css style definitions. Each object you provide must be a list of three elements. The first element will be a vector of the selectors to be styled (e.g. table, th, an id or html class). If the first element is a vector of length greater than one then the selectors will be comma separated in the css. The second element will be a vector of the css definitions and the third element will a vector of the values of those definitions.
#' @param file Character sting. If a file name is provided then the css code will be printed into that file. If the argument is NULL (default) then a string will be returned.
#'
#' @importFrom htmltools HTML
#' @export
#'
make_css <- function (..., file = NULL)
{
  css_defs <- list(...)
  for (x in css_defs) {
    if ((!is.list(x)) | (length(x) != 3L)) {
      stop("Each element in ... needs to be a list of three elements")
    }
    if (length(x[[2]]) != length(x[[3]])) {
      stop("The second and third elements of each list need to have the same length")
    }
  }
  all_css <- vapply(css_defs, function(x) {
    css_comp <- paste0(x[[2]], ": ", x[[3]], ";")
    style <- paste(css_comp, collapse = "\n  ")
    to_be_styled <- paste(x[[1]], collapse = ",\n")
    paste0(to_be_styled, " {\n  ", style, "\n}\n")
  }, FUN.VALUE = character(1))
  css_string <- htmltools::HTML(paste(all_css, collapse = "\n"))
  if (is.null(file)) {
    css_string
  }
  else {
    cat(css_string, file = file)
    invisible(NULL)
  }
}
