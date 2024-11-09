# Testing functions of manipulating spatialpoints ####
test_that("spatialpoints", {
  
  # get data
  data("visium_data")
  
  # get spatial points
  expect_equal(head(vrSpatialPoints(visium_data)),
               c("AAAGGCTCTCGCGCCG-1_Assay1","AAATGGCCCGTGCCCT-1_Assay1","AAATTACACGACTCTG-1_Assay1","AAGACATACGTGGTTT-1_Assay1",
                 "ACCTACTATAAATCTA-1_Assay1", "ACGCGGGCCAAGGACA-1_Assay1"))
  expect_equal(head(vrSpatialPoints(visium_data, assay = "Assay1")),
               c("AAAGGCTCTCGCGCCG-1_Assay1","AAATGGCCCGTGCCCT-1_Assay1","AAATTACACGACTCTG-1_Assay1","AAGACATACGTGGTTT-1_Assay1", 
                 "ACCTACTATAAATCTA-1_Assay1", "ACGCGGGCCAAGGACA-1_Assay1"))
  
  # subset on spatial points
  spatialpoints <- vrSpatialPoints(visium_data)
  visium_data_sub <- subset(visium_data, spatialpoints = spatialpoints[1:5])
  
  expect_equal(1,1L)
})