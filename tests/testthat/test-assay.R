# Testing functions of manipulating embeddings ####
test_that("add assay", {
  
  # get data
  data("xenium_data")

  # assay name and metadata
  assay <- vrAssayNames(xenium_data)
  sample_metadata <- SampleMetadata(xenium_data)

  # add cell assay
  new_assay <- formAssay(data = vrData(xenium_data),
                         coords = vrCoordinates(xenium_data), 
                         segments = vrSegments(xenium_data),
                         type = "cell",
                         image = vrImages(xenium_data),
                         main_image = vrMainImage(xenium_data[[assay]]),
                         name = assay)
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                  assay = new_assay,
                                  assay_name = "Xenium",
                                  sample = sample_metadata[assay, "Sample"],
                                  layer = sample_metadata[assay, "Layer"])
  expect_equal(unique(Metadata(xenium_data2)$assay_id), c("Assay1", "Assay2"))
  
  # add cell assay with metadata
  metadata <- data.frame(points = vrSpatialPoints(xenium_data))
  xenium_data2 <- addAssayVoltRon(xenium_data,
                                   metadata = metadata,
                                   assay = new_assay,
                                   assay_name = "Xenium",
                                   sample = sample_metadata[assay, "Sample"],
                                   layer = sample_metadata[assay, "Layer"])
  expect_equal(unique(Metadata(xenium_data2)$assay_id), c("Assay1", "Assay2"))
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
  expect_equal(unique(Metadata(xenium_data2)$assay_id), c("Assay1", "Assay2"))
  expect_true("points" %in% colnames(Metadata(xenium_data2)))
  expect_identical(vrSpatialPoints(xenium_data), 
                   Metadata(xenium_data2, assay = "Assay2")$points)
  expect_false(any(rownames(metadata) %in%  
                   rownames(Metadata(xenium_data2, assay = "Assay2"))))
  
  # add ROI assay
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