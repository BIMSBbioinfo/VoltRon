# Testing functions of manipulating images ####
test_that("coords_system", {
  
  # get data
  data("visium_data")
  
  # get coordinate systems
  vrSpatialNames(visium_data)
  vrSpatialNames(visium_data[["Assay1"]])
  vrSpatialNames(visium_data, assay = "Assay1")
  vrSpatialNames(visium_data, assay = "all")
  
  # set coordinate systems
  vrMainSpatial(visium_data[["Assay1"]]) <- "main"
  vrMainSpatial(visium_data[["Assay1"]]) <- c("main", "H&E")
  expect_error(vrMainSpatial(visium_data[["Assay1"]]) <- c("main", "DAPI"))
  expect_error(rMainSpatial(visium_data[["Assay1"]]) <- c("main2", "H&E"))
  vrMainSpatial(visium_data) <- c("main", "H&E")
  vrMainSpatial(visium_data, assay = "Assay1") <- c("main", "H&E")
  
  # return
  expect_equal(1,1L)
})