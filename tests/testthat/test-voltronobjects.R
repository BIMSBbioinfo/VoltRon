# Testing functions of manipulating assays ####
test_that("assay", {

  # get data
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

  # get data
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
