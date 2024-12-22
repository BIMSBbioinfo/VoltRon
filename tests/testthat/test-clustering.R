test_that("spatial", {
  
  # get data
  data("xenium_data")
  
  # profile neighbors, kNN
  xenium_data <- getProfileNeighbors(xenium_data, method = "kNN", k = 10, data.type = "pca")
  graphs <- vrGraph(xenium_data, graph.type = "kNN")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  xenium_data <- getClusters(xenium_data, graph = "kNN", label = "cluster_knn")
  expect_true(is.numeric(unique(xenium_data$cluster_knn)))
  
  # profile neighbors, SNN
  xenium_data <- getProfileNeighbors(xenium_data, method = "SNN", k = 10, data.type = "pca")
  graphs <- vrGraph(xenium_data, graph.type = "SNN")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN")
  expect_true(is.numeric(unique(xenium_data$cluster_SNN)))
  
  # check parameters
  expect_error(xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN", resolution = 0))
  expect_error(xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN", resolution = c(0,1)))
  expect_error(xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN", resolution = -0.5))
  expect_error(xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN", resolution = 1, abundance_limit = 0))
  expect_error(xenium_data <- getClusters(xenium_data, graph = "SNN", label = "cluster_SNN", resolution = 1, abundance_limit = 1.1))
  
  # return
  expect_equal(1,1L)
})