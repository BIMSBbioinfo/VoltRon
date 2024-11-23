# create dir
dir.create(td <- tempfile())
output_zarr <- paste0(td, "/xenium_data_zarr_test")
output_h5ad <- paste0(td, "/xenium_data_h5_test")
output_h5ad_2 <- paste0(td, "/xenium_data_h5_test")
output_h5ad_merged <- paste0(td, "/xenium_data_merged_h5_test")

test_that("save and load on disk", {
  
  # get data
  data("xenium_data")
  
  # HDF5
  xenium_data2 <- saveVoltRon(xenium_data, 
                              output = output_h5ad, 
                              format = "HDF5VoltRon", 
                              replace = TRUE)
  xenium_data2 <- loadVoltRon(dir = output_h5ad)

  # zarr
  xenium_data2 <- saveVoltRon(xenium_data, 
                              output = output_zarr,
                              format = "ZarrVoltRon", 
                              replace = TRUE)
  xenium_data2 <- loadVoltRon(dir = output_zarr)

  # remove files
  unlink(output_h5ad, recursive = TRUE)
  unlink(output_zarr, recursive = TRUE)

  expect_equal(1,1L)
})

test_that("double write and merging", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # write merged data
  xenium_data2 <- xenium_data 
  xenium_data2$Sample <- "XeniumR2"
  xenium_data2 <- merge(xenium_data, xenium_data2)
  xenium_data3 <- saveVoltRon(xenium_data2, 
                              output = output_zarr, 
                              format = "ZarrVoltRon", 
                              replace = TRUE)
  
  # merged ondisk data
  xenium_data_disk <- saveVoltRon(xenium_data, 
                              output = output_h5ad, 
                              format = "HDF5VoltRon", 
                              replace = TRUE)
  xenium_data2 <- xenium_data 
  xenium_data2$Sample <- "XeniumR2"
  xenium_data2_disk <- saveVoltRon(xenium_data2, 
                                   output = output_h5ad_2, 
                                   format = "HDF5VoltRon", 
                                   replace = TRUE)
  xenium_data_merged <- merge(xenium_data_disk, xenium_data2_disk)
  xenium_data_merged_disk <- saveVoltRon(xenium_data_merged, 
                                         output = output_h5ad_merged, 
                                         format = "HDF5VoltRon", 
                                         replace = TRUE)
  
  # remove files
  unlink(output_h5ad, recursive = TRUE)
  unlink(output_h5ad_2, recursive = TRUE)
  unlink(output_h5ad_merged, recursive = TRUE)
  unlink(output_zarr, recursive = TRUE)
  
  expect_equal(1,1L)
})

test_that("subsetting", {
  
  # get data
  data("xenium_data")
  
  # HDF5
  xenium_data2 <- saveVoltRon(xenium_data, 
                              output = output_h5ad, 
                              format = "HDF5VoltRon", 
                              replace = TRUE)

  # by spatialpoints
  spatpoints <- vrSpatialPoints(xenium_data2)[1:30]
  xenium_data2_subset <- subset(xenium_data2, spatialpoints = spatpoints)
  expect_equal(vrSpatialPoints(xenium_data2_subset), spatpoints)
  
  # visualize
  vrSpatialPlot(xenium_data2_subset)
  vrSpatialFeaturePlot(xenium_data2_subset, features = "KRT15")
  
  # by image
  xenium_data2_subset <- subset(xenium_data2, image = "290x202+98+17")
  expect_equal(length(vrSpatialPoints(xenium_data2_subset)), 392)
  expect_equal(is(vrImages(xenium_data2_subset)), "magick-image")
  expect_equal(is(vrImages(xenium_data2_subset, as.raster = TRUE)), "Image_Array")
  
  # visualize
  vrSpatialPlot(xenium_data2_subset)
  vrSpatialFeaturePlot(xenium_data2_subset, features = "KRT15")
  
  # remove files
  unlink(output_h5ad, recursive = TRUE)
  
  expect_equal(1,1L)
})

test_that("embeddings", {
  
  # get data
  data("xenium_data")
  
  # HDF5
  xenium_data2 <- saveVoltRon(xenium_data, 
                              output = output_h5ad, 
                              format = "HDF5VoltRon", 
                              replace = TRUE)
  xenium_data2 <- getPCA(xenium_data2, features = vrFeatures(xenium_data2), type = "pca_key")
  
  # check embeddings
  emb <- vrEmbeddings(xenium_data2, type = "pca_key")
  expect_true(all(dim(emb) > 1))
  expect_true(all(vrEmbeddingNames(xenium_data2) == c("pca", "umap", "pca_key")))
  
  # remove files
  unlink(output_h5ad, recursive = TRUE)
  
  expect_equal(1,1L)
})

test_that("visualization", {
  
  # get data
  data("visium_data")
  
  # HDF5
  visium_data2 <- saveVoltRon(visium_data, 
                              output = output_h5ad, 
                              format = "HDF5VoltRon", 
                              replace = TRUE)

  # check embeddings
  vrSpatialPlot(visium_data2, group.by = "Sample")
  vrSpatialFeaturePlot(visium_data2, features = "Count")
  
  # remove files
  unlink(output_h5ad, recursive = TRUE)
  
  expect_equal(1,1L)
})