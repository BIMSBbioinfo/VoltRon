test_that("subset", {
  
  # get data
  data("visium_data")
  data("xenium_data")
  data("merged_object")
  
  # subset based on assay
  subset(visium_data, assays = "Assay1")
  subset(visium_data, assays = "Visium")
  expect_error(subset(visium_data, assays = "Visium2"))
  
  # subset based on samples
  subset(visium_data, samples = "Anterior1")
  expect_error(subset(visium_data, samples = "Anterior2"))
  
  # subset based on assay
  subset(visium_data, spatialpoints = c("GTTATATTATCTCCCT-1_Assay1", "GTTTGGGTTTCGCCCG-1_Assay1"))
  expect_error(subset(visium_data, spatialpoints = c("point")))
  
  # subset based on features
  subset(visium_data, features = c("Map3k19", "Rab3gap1"))
  expect_error(subset(visium_data, features = c("feature")))
  
  # subset based on image
  tmp <- subset(merged_object, assay = "Assay1")
  tmp <- expect_warning(subset(tmp, image = "199x74+125+179.973327636719"))
  tmp <- subset(xenium_data, image = "240x150+135+188")
  
  # return
  expect_equal(1,1L)
})

test_that("subset image", {
  
  # visium
  data("visium_data")
  
  # return empty object 
  expect_warning(
    expect_null(subset(visium_data, image = "25x20+243+190"))
  )
})