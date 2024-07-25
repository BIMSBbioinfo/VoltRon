# Testing functions of manipulating images ####
test_that("image", {
  
  # get data
  data("visium_data")
  
  # get image
  images <- vrImages(visium_data)
  images <- vrImages(visium_data, name = "main")
  expect_error(images <- vrImages(visium_data, name = "main2"))
  images <- vrImages(visium_data, name = "main", channel = "H&E")
  expect_warning(images <- vrImages(visium_data, name = "main", channel = "H&E2"))
  
  # manipulate image
  visium_data_resize <- resizeImage(visium_data, size = 400)
  visium_data_modulate <- modulateImage(visium_data, brightness = 400)
  
  # add new image
  visium_data[["Assay1"]]@image[["new_image"]] <- vrImages(visium_data_resize)
  
  # get image names
  expect_equal(vrImageNames(visium_data), c("main", "new_image"))
  
  # get main image
  expect_equal(vrMainImage(visium_data[["Assay1"]]), "main")
  expect_equal(vrMainSpatial(visium_data[["Assay1"]]), "main")
  
  # change main image
  vrMainImage(visium_data[["Assay1"]]) <- "new_image"
  vrMainSpatial(visium_data[["Assay1"]]) <- "new_image"
  expect_equal(vrMainImage(visium_data[["Assay1"]]), "new_image")
  expect_equal(vrMainSpatial(visium_data[["Assay1"]]), "new_image")
  
  # return
  expect_equal(1,1L)
})