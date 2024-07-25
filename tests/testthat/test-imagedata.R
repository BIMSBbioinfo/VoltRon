# Testing functions of image datasets ####
test_that("imagedata", {
  
  # import image data
  img <- system.file(package = "VoltRon", "extdata/DAPI.tif")
  vrimgdata <- importImageData(img, tile.size = 10, stack.id = 1)
  vrimgdata <- importImageData(img, tile.size = 10, stack.id = 2)
  
  # import image data when image is a stack
  img_magick <- magick::image_read(img)
  img_stack <- c(img_magick, img_magick)
  vrimgdata <- importImageData(img_stack, tile.size = 10, stack.id = 1)
  vrimgdata <- importImageData(img_stack, tile.size = 10, stack.id = 2)
  expect_error(vrimgdata <- importImageData(img_stack, tile.size = 10, stack.id = 3))
  
  # return
  expect_equal(1,1L)
})