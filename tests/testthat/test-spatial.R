test_that("spatial neighbors", {
  
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

test_that("spatial tests", {
  
  # get data
  data("xenium_data")
  
  # spatial neighbors, radius based
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10)
  
  # getis ord test
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY"))
  expect_true(all(c("GNLY_hotspot_stat", "GNLY_hotspot_pvalue", "GNLY_hotspot_flag") %in% colnames(Metadata(xenium_data))))
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count"))
  expect_error(getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count1")))
  expect_error(getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY2", "Count1")))
  
  # multiple assays
  data("xenium_data")
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "Sample2"
  xenium_data <- merge(xenium_data, xenium_data2)
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10)
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count"))
  expect_true(all(c("GNLY_hotspot_stat", "GNLY_hotspot_pvalue", "GNLY_hotspot_flag") %in% colnames(Metadata(xenium_data))))
  expect_true(all(c("Count_hotspot_stat", "Count_hotspot_pvalue", "Count_hotspot_flag") %in% colnames(Metadata(xenium_data))))
})

test_that("niche clutsering", {
  
  # single assay
  data("xenium_data")
  xenium_data2 <- getSpatialNeighbors(xenium_data, radius = 15, method = "radius")
  xenium_data2 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius")
  expect_equal(vrFeatureTypeNames(xenium_data2), c("RNA", "Niche"))
  xenium_data2 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius", new_feature_name = "Niche2")
  expect_equal(vrFeatureTypeNames(xenium_data2), c("RNA", "Niche", "Niche2"))
  expect_error(getNicheAssay(xenium_data2, label = "clusters1", graph.type = "radius"))
  rm(xenium_data2)

  # multiple assays
  data("xenium_data")
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "Sample2"
  xenium_data2 <- merge(xenium_data2, xenium_data)
  xenium_data2 <- getSpatialNeighbors(xenium_data2, radius = 15, method = "radius")
  xenium_data3 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius")
  metadata <- vrFeatureTypeNames(xenium_data3, assay = "all")
  expect_equal(metadata$Feature, c("RNA,Niche", "RNA,Niche"))
  xenium_data3 <- getNicheAssay(xenium_data2, assay = "Assay1", label = "clusters", graph.type = "radius")
  metadata <- vrFeatureTypeNames(xenium_data3, assay = "all")
  expect_equal(metadata$Feature, c("RNA,Niche", "RNA"))
  
})