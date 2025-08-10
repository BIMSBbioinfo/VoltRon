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

test_that("import image voltron data", {
  
  # get image
  imgfile <- system.file("extdata", "DAPI.tif", package = "VoltRon")

  # tile size
  imgdata <- importImageData(imgfile, tile.size = 4, image_name = "main")
  imgdata <- importImageData(imgfile, tile.size = 1, image_name = "main")
  imgdata <- importImageData(imgfile, tile.size = 200, image_name = "main")
  expect_equal(vrImageChannelNames(imgdata)$Spatial, "main")
  expect_error(imgdata <- importImageData("", tile.size = 200, image_name = "main"))
  
  # channel names
  imgdata <- importImageData(imgfile, tile.size = 200, image_name = "main", channel_names = "DAPI")
  expect_equal(vrImageChannelNames(imgdata)$Channels, "DAPI")
  expect_error(importImageData(imgfile, tile.size = 10, 
                               channel_names = c("ch1", "ch3")))
  # multiple images
  imgfile <- c(imgfile, imgfile)
  imgdata <- importImageData(imgfile, tile.size = 4)
  expect_equal(vrImageChannelNames(imgdata)$Channels, "channel_1,channel_2")
  imgdata <- importImageData(imgfile, tile.size = 4, channel_names = c("DAPI", "DAPI2"))
  expect_equal(vrImageChannelNames(imgdata)$Channels, "DAPI,DAPI2")
  
  # return
  expect_equal(1,1L)
})

test_that("import ome.tiff", {
  
  # skip
  skip_if_not_installed("RBioFormats")
  
  # get image
  img.ometiff <- system.file(package = "VoltRon","extdata/xy_12bit__plant.ome.tiff")
  
  # ome.tiff
  vrimagedata <- importImageData(img.ometiff, series = 1, resolution = 1)
  expect_error(importImageData(img.ometiff, series = 1))
})