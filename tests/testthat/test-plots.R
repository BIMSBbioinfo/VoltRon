# Testing plotting functions
test_that("plots", {

  # get data
  data("visium_data")
  data("xenium_data")

  # get custom colors
  colors <- scales::hue_pal()(length(unique(xenium_data$clusters)))
  names(colors) <- unique(xenium_data$clusters)

  # embedding plot
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", group.ids = c(1,3,4), label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", colors = colors, label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", group.ids = c(1,3,4), colors = colors[c(1,3,4)], label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "clusters")
  vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "Sample")
  expect_error(vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "art"))

  # spatial plot
  vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", group.ids = c(1,3,4), plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", colors = colors, plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", group.ids = c(1,3,4), colors = colors[c(1,3,4)], plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "black")
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "white")
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "main")
  expect_warning(vrSpatialPlot(xenium_data, group.by = "clusters", background = c("main", "DAPI2")))
  
  # spatial plot without segmentation
  vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = FALSE)

  # spatial plot of visium
  vrSpatialPlot(visium_data)

  # spatial plot of melc data
  vrSpatialPlot(melc_data, group.by = "Clusters")
  expect_error(vrSpatialPlot(melc_data, group.by = "Clusters_new"))

  # feature plots
  vrSpatialFeaturePlot(visium_data, features = "Count")
  vrSpatialFeaturePlot(visium_data, features = "Stat1", norm = TRUE, log = TRUE)
  expect_error(vrSpatialFeaturePlot(visium_data, features = "Count_new"))
  
  # return
  expect_equal(1,1L)
})

# Testing plotting functions
test_that("missing_embedding_values", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  
  # change embeddings 
  vrEmbeddings(xenium_data, type = "new_umap") <- vrEmbeddings(xenium_data, type = "umap")[sample(1:length(vrSpatialPoints(xenium_data)), 500),]
  vrEmbeddingPlot(xenium_data, embedding = "new_umap")
  expect_error(vrEmbeddingPlot(xenium_data, embedding = "new_umap", group.by = "cluster"))
  vrEmbeddingFeaturePlot(xenium_data, embedding = "new_umap", features = "Count")
  expect_error(vrEmbeddingFeaturePlot(xenium_data, embedding = "new_umap", features = "Counts"))
  vrEmbeddingFeaturePlot(xenium_data, embedding = "new_umap", features = "REXO4")
  expect_error(vrEmbeddingFeaturePlot(xenium_data, embedding = "new_umap", features = "REXO4s"))
  
  # return
  expect_equal(1,1L)
})

# Testing plotting functions
test_that("rasterization", {
  
  # get data
  data("xenium_data")
  
  # spatial plot
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "black", n.tile = 100)
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "black", n.tile = 1)
  vrSpatialPlot(xenium_data, group.by = "clusters", background = "black", n.tile = 10)
  expect_warning(vrSpatialPlot(xenium_data, group.by = "clusters", background = c("main", "DAPI2")))
  
  # feature plots
  vrSpatialFeaturePlot(xenium_data, features = "Count", n.tile = 20)
  vrSpatialFeaturePlot(xenium_data, features = "KRT14", norm = TRUE, log = TRUE, n.tile = 10)
  expect_error(vrSpatialFeaturePlot(xenium_data, features = "Count_new"))
  
  # return
  expect_equal(1,1L)
})