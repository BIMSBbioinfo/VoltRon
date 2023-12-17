test_that("voltronobjects", {
  data("visium_data")
  vrSpatialPlot(visium_data, group.by = "Sample", alpha = 0.3)
})
