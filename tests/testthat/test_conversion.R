test_that("as.AnnData", {
  
  # library
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # file
  h5ad_file <- tempfile(fileext = ".h5ad")
  zarr_file <- tempfile(fileext = ".zarr")
  
  # xenium to anndata
  as.AnnData(xenium_data, file = h5ad_file)
  as.AnnData(xenium_data, file = h5ad_file, assay = "Assay1")
  as.AnnData(xenium_data, file = h5ad_file, flip_coordinates = TRUE)
  
  # visium to anndata
  as.AnnData(visium_data, file = h5ad_file)
  as.AnnData(visium_data, file = h5ad_file, assay = "Assay1")
  as.AnnData(visium_data, file = h5ad_file, flip_coordinates = TRUE)
  
  # xenium to anndata
  as.AnnData(xenium_data, file = zarr_file)
  as.AnnData(xenium_data, file = zarr_file, assay = "Assay1")
  as.AnnData(xenium_data, file = zarr_file, flip_coordinates = TRUE)
  
  # visium to anndata
  as.AnnData(visium_data, file = zarr_file)
  as.AnnData(visium_data, file = zarr_file, assay = "Assay1")
  as.AnnData(visium_data, file = zarr_file, flip_coordinates = TRUE)
  
  # clean file
  file.remove(h5ad_file)
  unlink(zarr_file, recursive = TRUE)
  expect_equal(1,1L)
})