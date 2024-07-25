# Testing functions of subsetting ####
test_that("subset", {
  
  # get data
  data("visium_data")
  
  # subset based on assay
  subset(visium_data, assays = "Assay1")
  subset(visium_data, assays = "Visium")
  expect_error(subset(visium_data, assays = "Visium2"))
  
  # subset based on samples
  subset(visium_data, samples = "Anterior1")
  expect_error(subset(visium_data, samples = "Anterior2"))
  
  # subset based on assay
  subset(visium_data, spatialpoints = c("GTTATATTATCTCCCT-1_Assay1", "GTTTGGGTTTCGCCCG-1_Assay1"))
  expect_error(subset(visium_data, spatialpoints = c("point")))
  
  # subset based on features
  subset(visium_data, features = c("Map3k19", "Rab3gap1"))
  expect_error(subset(visium_data, features = c("feature")))
  
  # return
  expect_equal(1,1L)
})