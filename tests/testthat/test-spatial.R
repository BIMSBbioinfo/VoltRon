test_that("spatial neighbors", {
  
  # get data
  data("xenium_data")
  
  # spatial neighbors, delaunay
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  vrSpatialPlot(xenium_data, graph.name = "delaunay", group.by = "clusters")
  
  # spatial neighbors, spatialkNN
  xenium_data <- getSpatialNeighbors(xenium_data, method = "spatialkNN", k = 5, verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "spatialkNN")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  vrSpatialPlot(xenium_data, graph.name = "spatialkNN", group.by = "clusters")
  
  # spatial neighbors, radius
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10, verbose = FALSE)
  graphs <- vrGraph(xenium_data, graph.type = "radius")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  vrSpatialPlot(xenium_data, graph.name = "radius", group.by = "clusters")
  
  # return
  expect_equal(1,1L)
})

test_that("spatial neighbors for subsets", {
  
  # get data
  data("xenium_data")
  
  # merge two of same types
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "XeniumR2"
  merged_data <- merge(xenium_data, xenium_data2, verbose = FALSE)
  
  # spatial neighbors, delaunay
  merged_data <- getSpatialNeighbors(merged_data, assay = "Assay1", method = "delaunay", verbose = FALSE)
  graphs <- vrGraph(merged_data, graph.type = "delaunay")
  expect_true(inherits(graphs,"igraph"))
  expect_true(length(igraph::E(graphs)) > 0)
  vrSpatialPlot(merged_data, graph.name = "delaunay", group.by = "clusters")
  
  # return
  expect_equal(1,1L)
})

test_that("spatial neighbors with distances", {
  
  # get data
  data("xenium_data")
  
  # spatial neighbors, spatialkNN
  xenium_data <- getSpatialNeighbors(xenium_data, method = "spatialkNN", k = 5, verbose = FALSE, 
                                     calculate.distances = TRUE)
  graphs <- vrGraph(xenium_data, graph.type = "spatialkNN")
  graphdist <- igraph::get.edge.attribute(graphs, name = "dist")
  expect_true(length(graphdist) > 1)
  
  # spatial neighbors, radius
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10, verbose = FALSE, 
                                     calculate.distances = TRUE)
  graphs <- vrGraph(xenium_data, graph.type = "radius")
  graphdist <- igraph::get.edge.attribute(graphs, name = "dist")
  expect_true(length(graphdist) > 1)
  
  # return
  expect_equal(1,1L)
})

test_that("spatial tests", {
  
  # get data
  data("xenium_data")
  
  # spatial neighbors, radius based
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10, verbose = FALSE)
  
  # getis ord test, replace
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY"), verbose = FALSE)
  tmp <- xenium_data$GNLY_hotspot_flag
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY"), verbose = FALSE, alpha.value = 0.0001)
  expect_false(all(xenium_data$GNLY_hotspot_flag == tmp))
  
  # getis ord test
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY"), verbose = FALSE)
  expect_true(all(c("GNLY_hotspot_stat", "GNLY_hotspot_padj", "GNLY_hotspot_flag") %in% colnames(Metadata(xenium_data))))
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count"), verbose = FALSE)
  expect_error(getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count1"), verbose = FALSE))
  expect_error(getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY2", "Count1"), verbose = FALSE))
  
  # multiple assays
  data("xenium_data")
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "Sample2"
  xenium_data <- merge(xenium_data, xenium_data2, verbose = FALSE)
  xenium_data <- getSpatialNeighbors(xenium_data, method = "radius", radius = 10, verbose = FALSE)
  xenium_data <- getHotSpotAnalysis(xenium_data, graph.type = "radius", features = c("GNLY", "Count"), verbose = FALSE)
  expect_true(all(c("GNLY_hotspot_stat", "GNLY_hotspot_padj", "GNLY_hotspot_flag") %in% colnames(Metadata(xenium_data))))
  expect_true(all(c("Count_hotspot_stat", "Count_hotspot_padj", "Count_hotspot_flag") %in% colnames(Metadata(xenium_data))))
})

test_that("niche clustering", {

  data("xenium_data")
  
  ####
  # single assay
  ####
  
  # build niche assay
  xenium_data2 <- getSpatialNeighbors(xenium_data, radius = 15, method = "radius", verbose = FALSE)
  xenium_data2 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius")
  expect_equal(vrFeatureTypeNames(xenium_data2), c("RNA", "Niche"))
  xenium_data2 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius", new_feature_name = "Niche2")
  expect_equal(vrFeatureTypeNames(xenium_data2), c("RNA", "Niche", "Niche2"))
  expect_error(getNicheAssay(xenium_data2, label = "clusters1", graph.type = "radius"))
  
  # cluster niches
  vrMainFeatureType(xenium_data2) <- "Niche"
  xenium_data2 <- getClusters(xenium_data2, nclus = 3, method = "kmeans", label = "cluster_niches")
  expect_true("cluster_niches" %in% colnames(Metadata(xenium_data2)))
  expect_error(xenium_data2 <- getClusters(xenium_data2, nclus = 0, method = "kmeans", label = "cluster_niches"))
  expect_error(xenium_data2 <- getClusters(xenium_data2, nclus = 1.1, method = "kmeans", label = "cluster_niches"))
  expect_error(xenium_data2 <- getClusters(xenium_data2, nclus = -1, method = "kmeans", label = "cluster_niches"))
  expect_error(xenium_data2 <- getClusters(xenium_data2, nclus = c(1,2), method = "kmeans", label = "cluster_niches"))
  
  # clean
  rm(xenium_data2)

  ####
  # multiple assays
  ####
  
  xenium_data2 <- xenium_data
  xenium_data2$Sample <- "Sample2"
  xenium_data2 <- merge(xenium_data2, xenium_data, verbose = FALSE)
  xenium_data2 <- getSpatialNeighbors(xenium_data2, radius = 15, method = "radius", verbose = FALSE)
  xenium_data3 <- getNicheAssay(xenium_data2, label = "clusters", graph.type = "radius")
  metadata <- vrFeatureTypeNames(xenium_data3, assay = "all")
  expect_equal(metadata$Feature, c("RNA,Niche", "RNA,Niche"))
  xenium_data3 <- getNicheAssay(xenium_data2, assay = "Assay1", label = "clusters", graph.type = "radius")
  metadata <- vrFeatureTypeNames(xenium_data3, assay = "all")
  expect_equal(metadata$Feature, c("RNA,Niche", "RNA"))
  
  # clean
  rm(xenium_data2)
  rm(xenium_data3)
  
})

test_that("neighborhood analysis", {
  
  # data
  data("xenium_data")
  
  # neighborhood test
  xenium_data <- getSpatialNeighbors(xenium_data, method = "delaunay", verbose = FALSE)
  results <- vrNeighbourhoodEnrichment(xenium_data, group.by = "clusters", graph.type = "delaunay")
  expect_error(results <- vrNeighbourhoodEnrichment(xenium_data, group.by = "clusters2", graph.type = "delaunay"))
})

test_that("multi assay neigh. analysis", {
  
  # data
  data("merged_object")
  
  # get spatial neighbors
  merged_object <- getSpatialNeighbors(merged_object, 
                                       assay = c("Assay1", "Assay4"), 
                                       group.by = list(
                                         Assay1 = "clusters",
                                         Assay4 = "gene"
                                         ),
                                       group.ids = list(
                                         Assay4 = c("KRT14", "KRT15")
                                         ),
                                       graph.key = "delaunay", 
                                       verbose = FALSE)
  
  # check if graph has multiple assays
  expect_contains(
    unique(stringr::str_extract(names(V(vrGraph(merged_object, 
                                                assay = c("Assay1", "Assay4"),
                                                graph.type = "delaunay"))), 
                                "Assay[0-9]+$")),
    c("Assay1", "Assay4")
  )
    
  # neighborhood test
  results <- vrNeighbourhoodEnrichment(merged_object, 
                                       assay = c("Assay1", "Assay4"), 
                                       group.by = list(Assay1 = "clusters", 
                                                       Assay4 = "gene"), 
                                       graph.type = "delaunay")

})

test_that("hot spot analysis with tiles", {
  
  # get data
  data("merged_object")
  vrMainAssay(merged_object) <- "MolAssay"
  
  # getis ord test of molecules with tiles
  merged_object <- getHotSpotAnalysis(merged_object, features = "gene", verbose = FALSE, n.tile = 2)
  vrSpatialPlot(merged_object, group.by = "gene_hotspot_flag")
  
  # getis ord test of molecules with no tiles
  expect_error({merged_object <- getHotSpotAnalysis(merged_object, features = "gene", verbose = FALSE)})
  merged_object <- getSpatialNeighbors(merged_object, method = "radius", radius = 30, verbose = FALSE)
  merged_object <- getHotSpotAnalysis(merged_object, features = "gene", graph.type = "radius", verbose = FALSE)
  vrSpatialPlot(merged_object, group.by = "gene_hotspot_flag")
  
  # get data
  data("merged_object")
  vrMainAssay(merged_object) <- "CellAssay"
  expect_error({merged_object <- getHotSpotAnalysis(merged_object, features = "CellType", 
                                      group.ids = c("MyelomaCells", "CD8_TCells"), 
                                      verbose = FALSE)})
  
  # getis ord test of cells with tiles, gets a warning
  merged_object <- getHotSpotAnalysis(merged_object, features = "CellType", 
                                      group.ids = c("MyelomaCells", "CD8_TCells"), 
                                      verbose = FALSE, n.tile = 10)
  vrSpatialPlot(merged_object, group.by = "CellType_hotspot_flag")
  
  # getis ord test of cells with tiles and a graph, gets a warning, should ignore graph
  merged_object <- getSpatialNeighbors(merged_object, method = "radius", radius = 30, verbose = FALSE)
  expect_warning({merged_object <- getHotSpotAnalysis(merged_object, features = "CellType", 
                                      group.ids = c("MyelomaCells", "CD8_TCells"), 
                                      verbose = FALSE, graph.type = "radius", n.tile = 10)})
  vrSpatialPlot(merged_object, group.by = "CellType_hotspot_flag")
  
  # getis ord test of cells with feature and group.ids, gets a warning, should ignore group.ids
  expect_warning({merged_object <- getHotSpotAnalysis(merged_object, features = "Count", 
                                      group.ids = c("MyelomaCells", "CD8_TCells"), 
                                      verbose = FALSE, graph.type = "radius", n.tile = 10)})
  vrSpatialPlot(merged_object, group.by = "Count_hotspot_flag")
  
})