# Testing functions of manipulating spatialpoints ####
test_that("spatialpoints", {
  
  # get data
  data("visium_data")
  
  # get spatial points
  print(head(vrSpatialPoints(visium_data)))
  print(head(vrSpatialPoints(visium_data, assay = "Assay1")))
  
  # subset on spatial points
  spatialpoints <- vrSpatialPoints(visium_data)
  visium_data_sub <- subset(visium_data, spatialpoints = spatialpoints[1:5])
  
  expect_equal(1,1L)
})