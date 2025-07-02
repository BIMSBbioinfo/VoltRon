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

# Testing functions of manipulating embeddings ####
test_that("embeddings from embeddings", {
  
  # get data
  data("visium_data")
  
  # test PCA on embeddings associated with obj
  # Get row names
  row_ids <- rownames(vrCoordinates(visium_data))
  
  # Create dummy embeddings with the same rows and 100 columns
  dummy_emb <- matrix(
    rnorm(length(row_ids) * 100),
    nrow = length(row_ids),
    ncol = 100
  )
  rownames(dummy_emb) <- row_ids
  vrEmbeddings(visium_data, type = "embed", overwrite=TRUE) <- dummy_emb
  
  #do PCA
  visium_data <- getPCA(visium_data, data.type = "embed", pca.key = "embedding_PCA", overwrite = TRUE)
  expect_true("embedding_PCA" %in% vrEmbeddingNames(visium_data))
  
  PCAs <- vrEmbeddings(visium_data, type = "embedding_PCA")
  expect_false(is.null(PCAs))
  expect_gt(nrow(PCAs), 0)
  expect_gt(ncol(PCAs), 0)
  
  # return
  expect_equal(1,1L)
})
