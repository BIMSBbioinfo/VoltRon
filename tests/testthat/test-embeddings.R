# Testing functions of manipulating embeddings ####
test_that("embeddings", {
  
  # get data
  data("visium_data")
  
  # write embedding
  vrEmbeddings(visium_data, type = "pca") <- vrCoordinates(visium_data)
  
  # check overwrite embeddings
  expect_error(vrEmbeddings(visium_data, type = "pca") <- vrCoordinates(visium_data))
  vrEmbeddings(visium_data, type = "pca", overwrite = TRUE) <- vrCoordinates(visium_data)
  
  # return
  expect_equal(1,1L)
})