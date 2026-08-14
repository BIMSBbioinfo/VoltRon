#' @docType package
#' @name VoltRon-package
#' @rdname VoltRon-package
#' @aliases VoltRon-package
#'
"_PACKAGE"

utils::globalVariables(
  c(
    "Assay",
    "KeyPoint",
    "NewSample",
    "Sample",
    "rx",
    "ry",
    "score",
    "segment",
    "str_pad",
    "value",
    "variable",
    "x",
    "y",
    "RTS_ID",
    "Module",
    "assoc_test",
    "from_value",
    "mean_value",
    "n",
    "p.adjust",
    "p_assoc",
    "p_segreg",
    "segreg_test",
    "to_value",
    "type",
    "region",
    "sample.metadata",
    "cell_id",
    "from.x",
    "from.y",
    "to.x",
    "to.y",
    ".SD",
    ".hasSlot",
    ":=",
    "assay_id",
    "id",
    "slideID",
    "J",
    "majortype",
    "addImg",
    "X1",
    "V1",
    "in_tissue",
    "fill",
    "colour",
    "color",
    "color_group",
    "group",
    ".N",
    "assay_id",
    "postfix",
    "cell_assay_id",
    "obs",
    "feat",
    "osbm"
  ),
  package = "VoltRon",
  add = FALSE
)

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  
  optional <- c(
    BPCells        = "IterableMatrix",
    DelayedArray   = "DelayedMatrix"
  )
  
  for (pkg in names(optional)) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      methods::setIs(optional[[pkg]], "data_matrix", where = ns)
      # conditional methods go here too, same reason:
      # methods::setMethod("vrData", optional[[pkg]], function(object, ...) ...)
    }
  }
  invisible()
}