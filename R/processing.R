#' @rdname NormalizeData
#' @method NormalizeData SpaceRover
#'
#' @export
#'
NormalizeData.SpaceRover <- function(object, ...) {

  # check assays
  sample.metadata <- SampleMetadata(object)
  assay_names <- rownames(sample.metadata)
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- NormalizeData(cur_assay)
  }
  return(object)
}

#' @rdname NormalizeData
#' @method NormalizeData srAssay
#'
#' @export
#'
NormalizeData.srAssay <- function(object, ...) {

  # size factor
  rawdata <- object@rawdata
  sizefactor <- colSums(rawdata)

  # log transformation
  sizefactor <- matrix(rep(sizefactor, nrow(rawdata)), byrow = T, nrow = nrow(rawdata))
  normdata <- (rawdata/sizefactor)*1000000
  normdata <- log(normdata + 1)

  # get normalized data
  object@normdata <- normdata

  # return
  return(object)
}
