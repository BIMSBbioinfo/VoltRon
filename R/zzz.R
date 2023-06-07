#' @docType package
#' @name VoltRon-package
#' @rdname VoltRon-package
#'
#' @importFrom conflicted conflict_prefer
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad = function(libname, pkgname) {

  # Normalize Data
  # conflicted::conflict_prefer("NormalizeData", "VoltRon")
  # registerS3method("NormalizeData", "VoltRon", "NormalizeData", envir = getNamespace("VoltRon"))
  # registerS3method('NormalizeData', 'VoltRon', "NormalizeData.VoltRon", envir = getNamespace("VoltRon"))
  # registerS3method('NormalizeData', 'srAssay', "NormalizeData.srAssay", envir = getNamespace("VoltRon"))
}
