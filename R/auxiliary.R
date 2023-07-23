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
