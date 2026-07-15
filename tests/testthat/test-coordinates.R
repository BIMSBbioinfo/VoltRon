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

test_that("segments",{
  
  # get data
  data("visium_data")
  
  # segments
  segments <- vrSegments(visium_data)
  expect_warning(segments <- vrSegments(visium_data, reg = TRUE))
  expect_warning(segments <- vrSegments(visium_data, assay = "Assay1", reg = TRUE))
  
  # get data
  data("merged_object")
  
  # check segments
  segments <- vrSegments(merged_object)
  expect_equal(segments, checkSegments(segments))
  expect_error(checkSegments(list(2,3)), "segments have to be named")
  names(segments) <- NULL
  expect_error(checkSegments(segments), "segments have to be named")

  expect_equal(1,1L)
  
})

test_that("replace coordinates", {
  
  # get data
  data("merged_object")

  # replace coords of selected assays
  coords <- vrCoordinates(merged_object, assay = "Assay1")
  vrCoordinates(merged_object, assay = "Assay1") <- coords[,c("x", "y")] * 2

  # replace coords of multiple assays
  coords <- vrCoordinates(merged_object, assay = "MolAssay")
  expect_error(
    vrCoordinates(merged_object, assay = "MolAssay") <- coords[,c("x", "y")] * 2,
    "Changing the coordinates of multiple assays in the same time are not permitted"
  )
})

test_that("replace segments", {
  
  # get data
  data("merged_object")
  
  # replace segments of selected assays
  segt <- vrSegments(merged_object, assay = "Assay3")
  segt_new <- lapply(segt, function(sg){
    sg[,c("x", "y")] <- sg[,c("x", "y")]*2
    sg
  })
  vrSegments(merged_object, assay = "Assay3") <- segt_new
  
  # replace coords of multiple assays
  segt <- vrSegments(merged_object, assay = "ROIAssay")
  expect_error(
    vrCoordinates(merged_object, assay = "ROIAssay") <- segt,
    "Changing the coordinates of multiple assays in the same time are not permitted"
  )
  
  # replace segments without id works
  segt <- vrSegments(merged_object, assay = "Assay3")
  segt_new <- lapply(segt, function(sg) sg[,c("x", "y")])
  vrSegments(merged_object, assay = "Assay3") <- segt_new
  
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

