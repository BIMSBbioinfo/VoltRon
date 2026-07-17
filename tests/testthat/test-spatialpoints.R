# Testing functions of manipulating spatialpoints ####
test_that("spatialpoints", {
  
  # get data
  data("visium_data")
  
  # get spatial points
  expect_equal(head(vrSpatialPoints(visium_data)),
               c("CCTTGACCACTTTATT-1_Assay1", "ATTTGTCTTGGGAGCT-1_Assay1", "TCACGCATTGTAGATC-1_Assay1", "CCGAGCTGTGCTTGTC-1_Assay1",
                 "GCATGGGTACTGACGC-1_Assay1", "AGTCGGCCCAAACGAC-1_Assay1"))
  expect_equal(head(vrSpatialPoints(visium_data, assay = "Assay1")),
               c("CCTTGACCACTTTATT-1_Assay1", "ATTTGTCTTGGGAGCT-1_Assay1", "TCACGCATTGTAGATC-1_Assay1", "CCGAGCTGTGCTTGTC-1_Assay1",
                 "GCATGGGTACTGACGC-1_Assay1", "AGTCGGCCCAAACGAC-1_Assay1"))
  
  # subset on spatial points
  spatialpoints <- vrSpatialPoints(visium_data)
  visium_data_sub <- subset(visium_data, spatialpoints = spatialpoints[1:5])
  
  expect_equal(1,1L)
})