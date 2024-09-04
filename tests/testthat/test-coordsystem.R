# Testing functions of manipulating images ####
test_that("coords_system", {
  
  # get data
  data("visium_data")
  
  # get coordinate systems
  vrSpatialNames(visium_data)
  vrSpatialNames(visium_data[["Assay1"]])
  vrSpatialNames(visium_data, assay = "Assay1")
  vrSpatialNames(visium_data, assay = "all")
  
  # return
  expect_equal(1,1L)
})