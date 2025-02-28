test_that("add column (numeric)", {

  # get data
  data("visium_data")
  data("xenium_data")
  
  # add a numeric column
  xenium_data2 <- addMetadata(xenium_data, value = 1, label = "column")
  old_column <- xenium_data2$column
  
  # replace value of label, and check difference
  xenium_data2 <- addMetadata(xenium_data2, value = 1:nrow(Metadata(xenium_data2)), label = "column")
  expect_false(all(xenium_data2$column == old_column))
  
  # update metadata of only one assay
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "XeniumR2"
  xenium_data2 <- merge(xenium_data, xenium_data2, verbose = FALSE)
  xenium_data2 <- addMetadata(xenium_data2, assay = "Assay1", value = 1:nrow(Metadata(xenium_data2, assay = "Assay1")), label = "column")
  
  # check if Assay2 label is empty
  tmp <- Metadata(xenium_data2, assay = "Assay2")$column
  expect_true(all(is.na(tmp)))
  
  expect_equal(1,1L)
})

test_that("add column (text)", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # add a numeric column
  xenium_data2 <- addMetadata(xenium_data, value = "text", label = "column")
  old_column <- xenium_data2$column
  
  # replace value of label, and check difference
  xenium_data2 <- addMetadata(xenium_data2, value = paste0("text", 1:nrow(Metadata(xenium_data2))), label = "column")
  expect_false(all(xenium_data2$column == old_column))
  
  # update metadata of only one assay
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "XeniumR2"
  xenium_data2 <- merge(xenium_data, xenium_data2, verbose = FALSE)
  xenium_data2 <- addMetadata(xenium_data2, assay = "Assay1", value = paste0("text", 1:nrow(Metadata(xenium_data2, assay = "Assay1"))), label = "column")
  
  # check if Assay2 label is empty text
  tmp <- Metadata(xenium_data2, assay = "Assay2")$column
  expect_true(all(tmp == ""))
  
  expect_equal(1,1L)
})

test_that("add column ($ method)", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # add a text column
  xenium_data2 <- xenium_data
  xenium_data2$column <- "text"
  old_column <- xenium_data2$column
  
  # replace value of label, and check difference
  xenium_data2$column <- paste0("text", 1:nrow(Metadata(xenium_data2)))
  expect_false(all(xenium_data2$column == old_column))
  
  # add a numeric column
  xenium_data2 <- xenium_data
  xenium_data2$column <- 1
  old_column <- xenium_data2$column
  
  # replace value of label, and check difference
  xenium_data2$column <- 1:nrow(Metadata(xenium_data2))
  expect_false(all(xenium_data2$column == old_column))
  
  expect_equal(1,1L)
})