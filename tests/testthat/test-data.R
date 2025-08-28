# Testing functions of assay data ####
test_that("data", {
  
  # get data
  data("visium_data")
  
  # get and set feature types
  vrMainFeatureType(visium_data)
  vrMainFeatureType(visium_data, assay = "Assay1") <- "RNA"
  expect_warning(vrMainFeatureType(visium_data) <- "main")
  
  # get data
  data <- vrData(visium_data, norm = FALSE)
  data <- vrData(visium_data, norm = TRUE)
  data <- vrData(visium_data, assay = "Assay1", norm = TRUE)
  
  # get data with features
  data <- vrData(visium_data, features = vrFeatures(visium_data)[1:2], norm = FALSE)
  expect_error(data <- vrData(visium_data, features = "feature", norm = FALSE))
  data <- vrData(visium_data, features = vrFeatures(visium_data)[1:2], norm = TRUE)
  expect_error(data <- vrData(visium_data, features = "feature", norm = TRUE))
  
  # return
  expect_equal(1,1L)
})

# Testing functions of assay data ####
test_that("molecule/tile data", {
  
  # get data
  data("merged_object")

  # get data of molecule assay
  data <- vrData(merged_object, assay = "MolAssay")
  expect_equal(nrow(data), 0)
  expect_equal(ncol(data), 
               length(vrSpatialPoints(merged_object, assay = "MolAssay")))
  
  # get data of an image data
  img <- system.file(package = "VoltRon", "extdata/DAPI.tif")
  vrimgdata <- importImageData(img, tile.size = 10)
  data <- vrData(vrimgdata)
  expect_equal(nrow(data), 0)
  expect_equal(ncol(data), 
               length(vrSpatialPoints(vrimgdata)))
  
  # return
  expect_equal(1,1L)
})