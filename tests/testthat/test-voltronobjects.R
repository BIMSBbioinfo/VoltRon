# Testing functions of manipulating coordinates ####
test_that("coordinates", {
  data("visium_data")

  # coordinates
  coords <- vrCoordinates(visium_data)
  coords <- vrCoordinates(visium_data, reg = TRUE)

  # update coordinates
  vrCoordinates(visium_data) <- coords*2
  vrCoordinates(visium_data, reg = TRUE) <- coords*3

  # segments
  segments <- vrSegments(visium_data)
  segments <- vrSegments(visium_data, reg = TRUE)
  expect_equal(1,1L)
})

# Testing functions of manipulating coordinates ####
test_that("image", {
  data("visium_data")

  # get image
  images <- vrImages(visium_data)
  images <- vrImages(visium_data, main_image = "main_image")

  # manipulate image
  visium_data_resize <- resizeImage(visium_data, size = 400)
  visium_data_modulate <- modulateImage(visium_data, brightness = 400)

  # add new image
  visium_data[["Assay1"]]@image[["new_image"]] <- vrImages(visium_data_resize)

  # get image names
  expect_equal(vrImageNames(visium_data), c("main_image", "new_image"))

  # get main image
  expect_equal(vrMainImage(visium_data[["Assay1"]]), "main_image")

  # change main image
  vrMainImage(visium_data[["Assay1"]]) <- "new_image"
  expect_equal(vrMainImage(visium_data[["Assay1"]]), "new_image")

  # return
  expect_equal(1,1L)
})
