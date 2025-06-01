img <- system.file(package = "VoltRon", "extdata/DAPI.tif")

test_that("import", {
  
  # import image data
  vrimgdata <- importImageData(img, tile.size = 10)
  vrimgdata <- importImageData(img, tile.size = 2)
  vrimgdata <- importImageData(img, tile.size = 200)
  vrimgdata <- importImageData(img, tile.size = 300)

  # import image data when image is a stack
  img_magick <- magick::image_read(img)
  img_stack <- c(img_magick, img_magick)
  vrimgdata <- importImageData(img_stack, tile.size = 10)

  # return
  expect_equal(1,1L)
})

test_that("visualize", {
  
  # import image data
  vrimgdata <- importImageData(img, tile.size = 100)
  
  # plot
  vrimgdata$class <- sample(paste(1:6),
                            length(vrSpatialPoints(vrimgdata)), 
                            replace = TRUE) 
  vrSpatialPlot(vrimgdata, alpha = 0.6, group.by = "class")
  vrimgdata$count <- 1:length(vrSpatialPoints(vrimgdata)) 
  vrSpatialFeaturePlot(vrimgdata, alpha = 0.6, features = "count")
  
  # return
  expect_equal(1,1L)
})