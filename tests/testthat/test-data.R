# Testing functions of assay data ####
test_that("data", {
  
  # get data
  data("visium_data")
  
  # get and set feature types
  vrMainFeatureType(visium_data)
  vrMainFeatureType(visium_data, assay = "Assay1") <- "RNA"
  expect_error(vrMainFeatureType(visium_data) <- "main")
  
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