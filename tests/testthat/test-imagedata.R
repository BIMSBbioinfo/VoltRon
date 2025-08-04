img <- system.file(package = "VoltRon", "extdata/DAPI.tif")

# Testing functions of image datasets ####
test_that("imagedata", {
  
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

# Testing functions of image datasets ####
test_that("image data visualization", {

  # import image data
  vrimgdata <- importImageData(img, tile.size = 10)
  
  # embeddings
  vrimgdata <- generateTileData(vrimgdata)
  vrimgdata <- getPCA(vrimgdata, dims = 15, overwrite = TRUE)
  vrimgdata <- getUMAP(vrimgdata, dims = 1:15)
  
  # embeddings
  vrEmbeddingPlot(vrimgdata, embedding = "umap")
  vrimgdata$pca <- vrEmbeddings(vrimgdata, type = "pca")[,1]
  vrEmbeddingFeaturePlot(vrimgdata, embedding = "umap", features = "pca")
  
  # spatialplot
  vrimgdata$pca <- vrEmbeddings(vrimgdata, type = "pca")[,1]
  vrSpatialFeaturePlot(vrimgdata, features = "pca")
  
  # return
  expect_equal(1,1L)
})