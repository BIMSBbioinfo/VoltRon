test_that("spatial graphs", {
  
  # get data
  data("xenium_data")
  
  # spatial graphs
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
})

test_that("overwrite graphs", {
  
  # spatial graphs
  data("xenium_data")
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  
  # overwrite using subset of metadata groups
  data("xenium_data")
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", group.by = "clusters", group.ids = c(1,2,3,4), verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", group.by = "clusters", group.ids = c(3,4,5,6), verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
})