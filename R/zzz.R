#' @docType package
#' @name spaceRover-package
#' @rdname spaceRover-package
#'
#' @importFrom conflicted conflict_prefer
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad = function(libname, pkgname) {

  # Normalize Data
  # conflicted::conflict_prefer("NormalizeData", "spaceRover")
  # registerS3method("NormalizeData", "SpaceRover", "NormalizeData", envir = getNamespace("spaceRover"))
  # registerS3method('NormalizeData', 'SpaceRover', "NormalizeData.SpaceRover", envir = getNamespace("spaceRover"))
  # registerS3method('NormalizeData', 'srAssay', "NormalizeData.srAssay", envir = getNamespace("spaceRover"))
}
