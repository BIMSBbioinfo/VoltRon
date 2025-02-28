test_that("assay", {

  # get data
  data("visium_data")

  # get assay names
  expect_equal(vrAssayNames(visium_data), "Assay1")

  # get assay object
  visium_data[["Assay1"]]

  # subset on assay name
  visium_data_sub <- subset(visium_data, assays = "Assay1")

  # subset on assay type
  visium_data_sub <- subset(visium_data, assays = "Visium")

  expect_equal(1,1L)
})

test_that("sample", {

  # get data
  data("visium_data")

  # get sample metadata
  SampleMetadata(visium_data)

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

test_that("merge objects", {
  
  # get data
  data("xenium_data")
  data("visium_data")
  data("melc_data")
  
  # merge two of same types
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "XeniumR2"
  merged_data <- merge(xenium_data, xenium_data2, verbose = FALSE)
  
  # merge two of different types
  merged_data <- merge(xenium_data, visium_data, verbose = FALSE)
  
  # merge multiple
  merge_list <- list(xenium_data, visium_data, melc_data)
  merged_data <- merge(merge_list[[1]], merge_list[-1], verbose = FALSE)
  
  expect_equal(1,1L)
})

