# library
skip_on_ci() # dont test on GitHub actions

# file
h5ad_file <- tempfile(fileext = ".h5ad")
zarr_file <- tempfile(fileext = ".zarr")
ometiff_file <- tempfile(fileext = ".ome.tiff")

test_that("as.AnnData", {

  # get data
  data("visium_data")
  data("xenium_data")

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

test_that("as.AnnData, python path", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # python.path
  expect_error(as.AnnData(visium_data, file = h5ad_file, python.path = ""))
  
  # python.path
  python.path <- system("which python", intern = TRUE)
  expect_error(as.AnnData(visium_data, file = zarr_file, python.path = python.path))
  expect_error(as.AnnData(visium_data, file = zarr_file, python.path = ""))
  
  # options path
  options(voltron.python.path = python.path)
  expect_error(as.AnnData(visium_data, file = zarr_file))
  options(voltron.python.path = NULL)
  expect_true(as.AnnData(visium_data, file = zarr_file))
  
  # clean file
  expect_equal(1,1L)
})

test_that("as.ometiff", {
  
  # get image
  data.file <- system.file(file.path('extdata', 'DAPI.tif'), package='VoltRon')
  
  # save image
  as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file)
  
  # clean file
  file.remove(ometiff_file)
  expect_equal(1,1L)
})

test_that("as.ometiff, python path", {
  
  # get image
  data.file <- system.file(file.path('extdata', 'DAPI.tif'), package='VoltRon')
  
  # save image
  as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file)
  
  # python.path
  expect_error(as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file, python.path = ""))
  
  # python.path
  python.path <- system("which python", intern = TRUE)
  expect_error(as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file, python.path = python.path))
  
  # options path
  options(voltron.python.path = python.path)
  expect_error(as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file))
  options(voltron.python.path = NULL)
  expect_true(as.OmeTiff(magick::image_read(data.file), out_path = ometiff_file))
  
  # clean file
  file.remove(ometiff_file)
  expect_equal(1,1L)
})

test_that("as.Seurat", {
  
  # check Seurat
  skip_if_not_installed("Seurat")
  
  # one sample
  test1 <- VoltRon::as.Seurat(xenium_data, cell.assay = "Xenium", type = "image")
  test1_Voltron <- as.VoltRon(test1, assay_name = "Xenium")
  
  # multiple samples
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "sample1"
  xenium_data2 <- merge(xenium_data2, xenium_data)
  test1 <- VoltRon::as.Seurat(xenium_data2, cell.assay = "Xenium", type = "image")
  test1_Voltron <- as.VoltRon(test1, assay_name = "Xenium")
  
  expect_equal(1,1L)
  
})