# Testing functions of manipulating assays ####
test_that("assay", {
  data("visium_data")

  # get assay names
  expect_equal(vrAssayNames(visium_data), "Assay1")

  # get assay object
  print(visium_data[["Assay1"]])

  # subset on assay name
  visium_data_sub <- subset(visium_data, assays = "Assay1")

  # subset on assay type
  visium_data_sub <- subset(visium_data, assays = "Visium")

  expect_equal(1,1L)
})

# Testing functions of manipulating samples ####
test_that("sample", {
  data("visium_data")

  # get sample metadata
  print(SampleMetadata(visium_data))

  # change sample name
  visium_data$Sample <- "Test_Sample_Name"

  # check metadata
  expect_equal(unique(visium_data$Sample), "Test_Sample_Name")

  # check list name
  expect_equal(unique(names(visium_data@samples)), "Test_Sample_Name")

  # check metadata name
  sample.metadata <- SampleMetadata(visium_data)
  expect_equal(sample.metadata$Sample == "Test_Sample_Name", TRUE)

  expect_equal(1,1L)
})

# Testing functions of manipulating spatialpoints ####
test_that("spatialpoints", {
  data("visium_data")

  # get spatial points
  print(head(vrSpatialPoints(visium_data)))
  print(head(vrSpatialPoints(visium_data, assay = "Assay1")))

  # subset on spatial points
  spatialpoints <- vrSpatialPoints(visium_data)
  visium_data_sub <- subset(visium_data, spatialpoints = spatialpoints[1:5])

  expect_equal(1,1L)
})

# Testing functions of manipulating coordinates ####
test_that("coordinates", {
  data("visium_data")

  # coordinates
  coords <- vrCoordinates(visium_data)
  coords <- vrCoordinates(visium_data, reg = TRUE)
  coords <- vrCoordinates(visium_data, assays = "Assay1", reg = TRUE)

  # update coordinates
  vrCoordinates(visium_data) <- coords*2
  vrCoordinates(visium_data, reg = TRUE) <- coords*3

  # flip coordinates
  visium_data <- flipCoordinates(visium_data)
  visium_data <- flipCoordinates(visium_data, reg = TRUE)

  # segments
  segments <- vrSegments(visium_data)
  segments <- vrSegments(visium_data, reg = TRUE)
  segments <- vrSegments(visium_data, assays = "Assay1", reg = TRUE)
  expect_equal(1,1L)
})

# Testing functions of manipulating images ####
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

# Testing functions of manipulating images ####
test_that("embeddings", {
  data("visium_data")

  # write embedding
  vrEmbeddings(visium_data, type = "pca") <- vrCoordinates(visium_data)

  # check overwrite embeddings
  expect_error(vrEmbeddings(visium_data, type = "pca") <- vrCoordinates(visium_data))
  vrEmbeddings(visium_data, type = "pca", overwrite = TRUE) <- vrCoordinates(visium_data)

  # return
  expect_equal(1,1L)
})
