# Testing functions of assay data ####
test_that("featuredata", {
  
  # get data
  data("visium_data")
  
  # get and set feature types
  vrMainFeatureType(visium_data)
  vrMainFeatureType(visium_data, assay = "Assay1") <- "RNA"
  expect_error(vrMainFeatureType(visium_data) <- "main")
  
  # get feature data
  expect_equal(vrFeatureTypeNames(visium_data), "RNA")
  vrFeatureData(visium_data, feat_type = "RNA")
  expect_equal(vrFeatureData(visium_data, feat_type = "Prot"), NULL)
  
  # return
  expect_equal(1,1L)
})