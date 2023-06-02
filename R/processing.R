#' @include generics.R
#'
NULL

#' @rdname normalizeData
#' @concept preprocessing
#' @method normalizeData SpaceRover
#'
#' @export
normalizeData.SpaceRover <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- normalizeData(cur_assay, ...)
  }
  return(object)
}

#' @rdname normalizeData
#' @concept preprocessing
#' @method normalizeData srAssay
#'
#' @export
normalizeData.srAssay <- function(object, method = "LogNorm", desiredQuantile = 0.9) {

  # size factor
  rawdata <- object@rawdata
  sizefactor <- colSums(rawdata)

  if(method == "LogNorm"){
    # log transformation
    sizefactor <- matrix(rep(sizefactor, nrow(rawdata)), byrow = T, nrow = nrow(rawdata))
    normdata <- (rawdata/sizefactor)*1000000
    normdata <- log(normdata + 1)
  } else if(method == "QuanNorm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
  } else if(method == "LogQuanNorm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
    normdata <- log(normdata + 1)
  } else if(method == "CLR") {
    normdata <- apply(rawdata, 2, function(x) {
      log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    })
  }

  # get normalized data
  object@normdata <- normdata

  # return
  return(object)
}
