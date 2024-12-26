# Testing functions of manipulating images ####
test_that("coords_system", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # get coordinate systems
  vrSpatialNames(visium_data)
  vrSpatialNames(visium_data[["Assay1"]])
  vrSpatialNames(visium_data, assay = "Assay1")
  vrSpatialNames(visium_data, assay = "all")
  
  # set coordinate systems
  vrMainSpatial(visium_data[["Assay1"]]) <- "main"
  vrMainSpatial(visium_data[["Assay1"]]) <- c("main", "H&E")
  expect_error(vrMainSpatial(visium_data[["Assay1"]]) <- c("main", "DAPI"))
  expect_error(vrMainSpatial(visium_data[["Assay1"]]) <- c("main2", "H&E"))
  vrMainSpatial(visium_data) <- c("main", "H&E")
  vrMainSpatial(visium_data, assay = "Assay1") <- c("main", "H&E")
  
  # set coordinate system when ignore is specified
  expect_warning(vrMainSpatial(visium_data[["Assay1"]], ignore = TRUE) <- "main2")
  expect_error(vrMainSpatial(visium_data[["Assay1"]], ignore = FALSE) <- "main2")
  
  # mergeddata
  mergeddata <- merge(xenium_data, visium_data)
  expect_true(nrow(vrSpatialNames(mergeddata, assay = "all")) == 2)
  vrMainSpatial(mergeddata) <- c("main", "H&E")
  vrMainAssay(mergeddata) <- "Xenium"
  expect_error(vrMainSpatial(mergeddata) <- c("main", "H&E"))

  # return
  expect_equal(1,1L)
})