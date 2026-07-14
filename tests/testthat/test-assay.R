# Testing functions of manipulating embeddings ####
test_that("add cell assay", {
  
  data("xenium_data")
  assay <- vrAssayNames(xenium_data)
  sample_metadata <- SampleMetadata(xenium_data)

  # form assay
  new_assay <- formAssay(data = vrData(xenium_data),
                         coords = vrCoordinates(xenium_data), 
                         segments = vrSegments(xenium_data),
                         type = "cell",
                         image = vrImages(xenium_data),
                         main_image = vrMainImage(xenium_data[[assay]]),
                         name = assay)
  expect_equal(
    vrSpatialPoints(new_assay),
    vrSpatialPoints(xenium_data)
  )
  
  # add assay
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  assay = new_assay,
                                  assay_name = "Xenium",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_equal(unique(Metadata(xenium_data2)$assay_id), c("Assay1", "Assay2"))
  expect_equal(SampleMetadata(xenium_data2)$Assay, c("Xenium", "Xenium"))
  expect_equal(rownames(SampleMetadata(xenium_data2)), c("Xenium", "Xenium"))
  expect_equal(nrow(SampleMetadata(xenium_data2)), 2)
  
  
  # form assay
  new_assay <- formAssay(data = vrData(xenium_data),
                         coords = vrCoordinates(xenium_data), 
                         segments = vrSegments(xenium_data),
                         type = "cell",
                         image = vrImages(xenium_data),
                         main_image = vrMainImage(xenium_data[[assay]]),
                         name = assay)
  
  # add with metadata
  metadata <- data.frame(points = vrSpatialPoints(xenium_data))
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  metadata = metadata,
                                  assay = new_assay,
                                  assay_name = "Xenium",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_true("points" %in% colnames(Metadata(xenium_data2)))
  expect_identical(vrSpatialPoints(xenium_data), 
                   Metadata(xenium_data2, assay = "Assay2")$points)
  
  # we check if adding an assay with a metadata with rownames will
  # remove the rownames an populate id instead
  metadata <- data.frame(row.names = vrSpatialPoints(xenium_data), 
                         points = vrSpatialPoints(xenium_data))
  expect_identical(rownames(metadata), vrSpatialPoints(xenium_data))
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  metadata = metadata,
                                  assay = new_assay,
                                  assay_name = "Xenium",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_true("points" %in% colnames(Metadata(xenium_data2)))
  expect_identical(vrSpatialPoints(xenium_data), 
                   Metadata(xenium_data2, assay = "Assay2")$points)
  expect_false(any(rownames(metadata) %in%  
                   rownames(Metadata(xenium_data2, assay = "Assay2"))))
})

test_that("add ROI assay", {
  
  data("xenium_data")
  assay <- vrAssayNames(xenium_data)
  sample_metadata <- SampleMetadata(xenium_data)
  
  coords <- vrCoordinates(xenium_data)[1:2,,drop = FALSE]
  new_assay <- formAssay(coords = coords, 
                         segments = vrSegments(xenium_data)[1:2],
                         type = "ROI",
                         image = vrImages(xenium_data, assay = assay),
                         main_image = vrMainImage(xenium_data[[assay]]),
                         name = assay)
  metadata <- data.frame(check.rows = FALSE, 
                         row.names = rownames(coords), 
                         rep("art", nrow(coords)))
  colnames(metadata) <- "label"
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  assay = new_assay,
                                  metadata = metadata,
                                  assay_name = "random_ROI_assay",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_true("label" %in% colnames(Metadata(xenium_data2, type = "ROI")))
  expect_true(nrow(Metadata(xenium_data2, type = "ROI")) == 2)

  # return
  expect_equal(1,1L)
})

test_that("add ROI assay (segments only)", {
  
  data("xenium_data")
  assay <- vrAssayNames(xenium_data)
  sample_metadata <- SampleMetadata(xenium_data)
  
  segments <- vrSegments(xenium_data)[1:2]
  new_assay <- formAssay(segments = segments,
                         type = "ROI",
                         image = vrImages(xenium_data, assay = assay),
                         main_image = vrMainImage(xenium_data[[assay]]),
                         name = assay)
  metadata <- data.frame(check.rows = FALSE, 
                         row.names = names(segments), 
                         rep("art", length(segments)))
  colnames(metadata) <- "label"
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  assay = new_assay,
                                  metadata = metadata,
                                  assay_name = "random_ROI_assay",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_true("label" %in% colnames(Metadata(xenium_data2, type = "ROI")))
  expect_true(nrow(Metadata(xenium_data2, type = "ROI")) == 2)
  
  # return
  expect_equal(1,1L)
})

test_that("segments vs coordinates", {
  
  data("xenium_data")
  assay <- vrAssayNames(xenium_data)
  sample_metadata <- SampleMetadata(xenium_data)
  
  # non-matching segments and coords should throw an error
  segments <- vrSegments(xenium_data)
  segments <- segments[-c(1:2)]
  coords <- vrCoordinates(xenium_data)
  expect_error(
    formAssay(coords = coords,
              segments = segments,
              image = vrImages(xenium_data, assay = assay),
              main_image = vrMainImage(xenium_data[[assay]]),
              name = assay),
    "Number of segments do not match with the number of points!")
  
})