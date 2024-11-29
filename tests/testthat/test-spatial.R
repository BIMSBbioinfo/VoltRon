test_that("spatial", {
  
  # get data
  data("xenium_data")
  
  # spatial neighbors, delaunay
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay")
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  
  # spatial neighbors, spatialkNN
  xenium_data <- getSpatialNeighbors(xenium_data, method = "spatialkNN", k = 5)
  graphs <- vrGraph(xenium_data, graph.type = "spatialkNN")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  
  # spatial neighbors, radius
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10)
  graphs <- vrGraph(xenium_data, graph.type = "radius")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  
  # return
  expect_equal(1,1L)
})