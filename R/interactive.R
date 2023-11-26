####
# Basilisk Environment ####
####

#' The Python environment
#'
#' Defines a conda environment via Basilisk, which is used
#' to convert R objects to Zarr stores.
#' This environment has been adapted from zellkonverter::.AnnDataDependencies.
#' Reference: https://bioconductor.org/packages/release/bioc/vignettes/basilisk/inst/doc/motivation.html
#'
#' @importFrom basilisk BasiliskEnvironment
#'
#' @keywords internal
py_env <- basilisk::BasiliskEnvironment(
  envname="VoltRon_basilisk_env",
  pkgname="VoltRon",
  packages=c(
    "numpy==1.*",
    "pandas==1.*",
    "anndata==0.7.*",
    "h5py==3.*",
    "hdf5==1.*",
    "natsort==7.*",
    "packaging==20.*",
    "scipy==1.*",
    "sqlite==3.*",
    "zarr==2.*",
    "numcodecs==0.*"
  ),
  pip=c(
    "ome-zarr==0.2.1"
  )
)

####
# Conversion into Zarr for Vitessce ####
####

#' Title
#'
#' @param vrimage VoltRon image
#' @param out_path output path to ome.zarr
#' @param image_id image name
#'
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
#' @importFrom reticulate import
#' @importFrom magick image_raster
#' @importFrom grDevices col2rgb
#'
#' @export
#'
vrImage_to_zarr <- function (vrimage, out_path, image_id = "main_image")
{
  img_arr <- apply(as.matrix(magick::image_raster(vrimage, tidy = FALSE)), c(1, 2), col2rgb)
  proc <- basilisk::basiliskStart(py_env)
  on.exit(basilisk::basiliskStop(proc))
  success <- basilisk::basiliskRun(proc, function(img_arr, image_id, out_path) {
    zarr <- reticulate::import("zarr")
    ome_zarr <- reticulate::import("ome_zarr")
    z_root <- zarr$open_group(out_path, mode = "w")
    obj_list <- function(...) {
      retval <- stats::setNames(list(), character(0))
      param_list <- list(...)
      for (key in names(param_list)) {
        retval[[key]] = param_list[[key]]
      }
      retval
    }
    default_window <- obj_list(start = 0, min = 0, max = 255, end = 255)
    ome_zarr$writer$write_image(image = img_arr,
                                group = z_root,
                                axes = "cyx",
                                omero = obj_list(name = image_id, version = "0.3",
                                                 rdefs = obj_list(),
                                                 channels = list(obj_list(label = "r", color = "FF0000", window = default_window),
                                                                 obj_list(label = "g", color = "00FF00", window = default_window),
                                                                 obj_list(label = "b", color = "0000FF", window = default_window))))
    return(TRUE)
  }, img_arr = img_arr, image_id = image_id, out_path = out_path)
  return(success)
}
