# Testing functions of manipulating coordinates ####
test_that("coordinates", {
  
  # get data
  data("visium_data")
  
  # coordinates
  coords <- vrCoordinates(visium_data)
  coords <- vrCoordinates(visium_data, image_name = "main")
  coords <- vrCoordinates(visium_data, spatial_name = "main")
  expect_warning(coords <- vrCoordinates(visium_data, reg = TRUE))
  expect_warning(coords <- vrCoordinates(visium_data, assay = "Assay1", reg = TRUE))
  
  # update coordinates
  vrCoordinates(visium_data) <- coords*2
  expect_error(vrCoordinates(visium_data, reg = TRUE) <- coords*3)
  
  # flip coordinates
  visium_data <- flipCoordinates(visium_data)
  
  # segments
  segments <- vrSegments(visium_data)
  expect_warning(segments <- vrSegments(visium_data, reg = TRUE))
  expect_warning(segments <- vrSegments(visium_data, assay = "Assay1", reg = TRUE))
  
  expect_equal(1,1L)
})

test_that("2d vs 3d", {
  
  # get data
  data("xenium_data")
  
  # change to 2 dimensions, might occur in previous versions of VoltRon
  xenium_data[["Assay1"]]@image$main@coords <- xenium_data[["Assay1"]]@image$main@coords[,1:2]
  
  # coordinates
  coords <- vrCoordinates(xenium_data)
  
  # merge with second object and try again
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "Sample2"
  xenium_data <- merge(xenium_data, xenium_data2)
  
  expect_equal(1,1L)
})