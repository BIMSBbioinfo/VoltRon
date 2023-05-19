#' @include generics.R
#'
NULL

#' @rdname NormalizeData
#' @concept preprocessing
#' @method NormalizeData SpaceRover
#'
#' @export
NormalizeData.SpaceRover <- function(object, assay = NULL, ...) {

  # # check assays
  # if(is.null(assay))
  #   assay <- object@main.assay
  #
  # # get all assays that are of main assay
  # assay_names <- unique(object@sample.metadata$Assay)
  # if(!assay %in% assay_names)
  #   stop("There are no assays named '", assay, "' in this object!")
  # sample.metadata <- object@sample.metadata[object@sample.metadata == assay,]
  # assays <- rownames(sample.metadata)

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- NormalizeData(cur_assay, ...)
  }
  return(object)
}

#' @rdname NormalizeData
#' @concept preprocessing
#' @method NormalizeData srAssay
#'
#' @export
NormalizeData.srAssay <- function(object, method = "LogNorm", desiredQuantile = 0.9) {

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
  }

  # get normalized data
  object@normdata <- normdata

  # return
  return(object)
}

#' @rdname NormalizeData
#' @concept preprocessing
#'
#' @export
NormalizeData.default <- function(object, ...) {
  Seurat::NormalizeData(object, ...)
}
