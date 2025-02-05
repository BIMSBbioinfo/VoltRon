# Testing plotting functions
test_that("non-spatial plots", {

  # get data
  data("xenium_data")
  
  # get custom colors
  # colors <- scales::hue_pal()(length(unique(xenium_data$clusters)))
  colors <- hue_pal(length(unique(xenium_data$clusters)))
  names(colors) <- unique(xenium_data$clusters)
  
  # embedding plot
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", group.ids = c(1,3,4), label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", colors = colors, label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", embedding = "umap", group.ids = c(1,3,4), colors = colors[c(1,3,4)], label = T)
  vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "clusters")
  vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "Sample")
  expect_error(vrEmbeddingPlot(xenium_data, group.by = "clusters", ncol = 3, split.by = "art"))
  
  # embedding feature plot
  vrEmbeddingFeaturePlot(xenium_data, features = c("ACTA2", "TACSTD2"), embedding = "umap", combine.features = TRUE)
  
  # scatterplot
  vrScatterPlot(xenium_data, feature.1 = "NKG7", feature.2 = "TRAC")
  xenium_data <- normalizeData(xenium_data)
  vrScatterPlot(xenium_data, feature.1 = "NKG7", feature.2 = "TRAC", norm = TRUE)
})


# Testing plotting functions
test_that("spatial plots", {

  # get data
  data("xenium_data")

  # get custom colors
  # colors <- scales::hue_pal()(length(unique(xenium_data$clusters)))
  colors <- hue_pal(length(unique(xenium_data$clusters)))
  names(colors) <- unique(xenium_data$clusters)

  # spatial plot, groups and colors
  vrSpatialPlot(xenium_data, group.by = "clusters", plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", group.ids = c(1,3,4), plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", colors = colors, plot.segments = TRUE)
  vrSpatialPlot(xenium_data, group.by = "clusters", group.ids = c(1,3,4), colors = colors[c(1,3,4)], plot.segments = TRUE)
  
  # spatial plot, background color
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "black")
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "yellow")
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "white")
  
  # spatial plot with spatial and channel arguments
  vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main")
  vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main", background.color = "yellow")
  vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main", channel = "DAPI")
  vrSpatialPlot(xenium_data, group.by = "clusters", channel = "DAPI")
  vrSpatialPlot(xenium_data, group.by = "clusters", channel = "DAPI2")
  vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main", channel = "DAPI2")
  expect_error(vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main2"))
  expect_error(vrSpatialPlot(xenium_data, group.by = "clusters", spatial = "main2", channel = "DAPI2"))

  # spatial plot, old background argument
  expect_warning(
    expect_warning(vrSpatialPlot(xenium_data, group.by = "clusters", background = "main")))
  expect_warning(
      expect_warning(vrSpatialPlot(xenium_data, group.by = "clusters", background = c("main", "DAPI2"))))
  expect_warning(
    expect_warning(
      expect_error(vrSpatialPlot(xenium_data, group.by = "clusters", background = "main2"))))
  
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
  vrSpatialFeaturePlot(xenium_data, features = c("ACTA2", "TACSTD2"), norm = TRUE, combine.features = TRUE)
  
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
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "black", n.tile = 100)
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "black", n.tile = 1)
  vrSpatialPlot(xenium_data, group.by = "clusters", background.color = "black", n.tile = 10)

  # feature plots
  vrSpatialFeaturePlot(xenium_data, features = "Count", n.tile = 20)
  vrSpatialFeaturePlot(xenium_data, features = "KRT14", norm = TRUE, log = TRUE, n.tile = 10)
  expect_error(vrSpatialFeaturePlot(xenium_data, features = "Count_new"))
  vrSpatialFeaturePlot(xenium_data, features = c("ACTA2", "TACSTD2"), norm = TRUE, n.tile = 100, combine.features = TRUE)
  vrSpatialFeaturePlot(xenium_data, features = c("ACTA2", "TACSTD2"), norm = TRUE, n.tile = 2, combine.features = TRUE)
  
  # embedding plots
  vrEmbeddingPlot(xenium_data, n.tile = 1200, group.by = "clusters")
  vrEmbeddingFeaturePlot(xenium_data, n.tile = 1200, features = c("ACTA2", "TACSTD2"))
  vrEmbeddingFeaturePlot(xenium_data, n.tile = 100, features = c("ACTA2", "TACSTD2"), embedding = "umap", combine.features = TRUE)
  vrEmbeddingFeaturePlot(xenium_data, n.tile = 2, features = c("ACTA2", "TACSTD2"), embedding = "umap", combine.features = TRUE)
  vrEmbeddingPlot(xenium_data, n.tile = 2, group.by = "clusters")
  vrEmbeddingFeaturePlot(xenium_data, n.tile = 10, features = c("ACTA2"))
  
  # return
  expect_equal(1,1L)
})

# testing multilayer plots
test_that("multilayer", {
  
  skip_if_not_installed("ggnewscale")
  data("merged_object")
  
  # single
  vrSpatialPlot(merged_object)
  
  # cell vs ROI (without segments)
  vrSpatialPlot(merged_object, plot.segments = FALSE) |>
    addSpatialLayer(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.3)
  
  # cell vs ROI (with segments)
  vrSpatialPlot(merged_object, plot.segments = TRUE) |>
    addSpatialLayer(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue"))

  # ROI vs cell
  vrSpatialPlot(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue")) |>
    addSpatialLayer(merged_object, assay = "Assay1")
  vrSpatialPlot(merged_object, assay = "Assay3", group.by = "Sample", alpha = 1, colors = list(Block = "blue")) |>
    addSpatialLayer(merged_object, assay = "Assay1", plot.segments = TRUE, alpha = 0.4)
  vrSpatialPlot(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue")) |>
    addSpatialLayer(merged_object, assay = "Assay1", n.tile = 100)
  
  # cell vs molecule (without segments)
  vrSpatialPlot(merged_object, plot.segments = FALSE) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"))
  
  # cell vs molecule 
  vrSpatialPlot(merged_object, plot.segments = TRUE) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"))
  
  # molecule vs cell (with segments)
  vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green")) |>
    addSpatialLayer(merged_object, assay = "Assay1")
  vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green")) |>
    addSpatialLayer(merged_object, assay = "Assay1", plot.segments = TRUE, alpha = 0.4)
  
  # cells, ROIs and molecules together
  vrSpatialPlot(merged_object, plot.segments = TRUE) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green")) |>
    addSpatialLayer(merged_object, assay = "Assay3", group.by = "Layer", alpha = 0.4, colors = list(Section3 = "blue"))
  
  expect_equal(1,1)
})

# testing multilayer plots
# TODO: tiling multilayer visualization behavior is not ideal right now
test_that("multilayer (with tiling)", {
  
  skip_if_not_installed("ggnewscale")
  data("merged_object")
  
  # single
  vrSpatialPlot(merged_object)
  
  # cell vs ROI (without segments)
  vrSpatialPlot(merged_object, plot.segments = FALSE, n.tile = 100) |>
    addSpatialLayer(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue"))
  
  # ROI vs cell
  vrSpatialPlot(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue")) |>
    addSpatialLayer(merged_object, assay = "Assay1", n.tile = 100)
  
  # cell vs molecule (without segments)
  vrSpatialPlot(merged_object, plot.segments = FALSE) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"), n.tile = 100)
  vrSpatialPlot(merged_object, plot.segments = FALSE, n.tile = 100) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"), n.tile = 100)
  
  # cell vs molecule 
  vrSpatialPlot(merged_object, plot.segments = TRUE) |>
    addSpatialLayer(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"), n.tile = 100)
  
  # molecule vs cell (with segments)
  vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green")) |>
    addSpatialLayer(merged_object, assay = "Assay1", n.tile = 100)
  vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green")) |>
    addSpatialLayer(merged_object, assay = "Assay1", plot.segments = TRUE)
  
  # molecule vs ROI
  vrSpatialPlot(merged_object, assay = "Assay2", group.by = "gene", alpha = 1, colors = list(KRT15 = "blue", KRT14 = "green"), n.tile = 100) |>
    addSpatialLayer(merged_object, assay = "Assay3", group.by = "Sample", alpha = 0.4, colors = list(Block = "blue"))
  
  expect_equal(1,1)
})