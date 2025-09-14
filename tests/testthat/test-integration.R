test_that("integration", {
  
  # cell to ROI
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay1", to = "Assay3", 
                                features = "CellType")
  expect_identical(c("main", "main_pseudo"),
                   vrFeatureTypeNames(merged_object, assay = "Assay3"))
  vrMainFeatureType(merged_object, assay = "Assay3") <- "main_pseudo"
  vrMainAssay(merged_object) <- "ROIAssay"
  expect_contains(vrFeatures(merged_object),
                  c("CD4_TCells", "CD8_TCells"))
  
  # cell to ROI, all features
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay1", to = "Assay3", new_feature_name = "main_pseudo2")
  expect_identical(c("main", "main_pseudo", "main_pseudo2"),
                   vrFeatureTypeNames(merged_object, assay = "Assay3"))
  vrMainFeatureType(merged_object, assay = "Assay3") <- "main_pseudo2"
  vrMainAssay(merged_object) <- "ROIAssay"
  expect_identical(vrFeatures(merged_object, assay = "ROIAssay"),
                   vrFeatures(merged_object, assay = "CellAssay"))
  
  # ROI to cell
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay3", to = "Assay1", 
                                features = "annotation")
  expect_contains(merged_object$annotation,
                  c("DCIS_Subtype1", "DCIS_Subtype2", "Immune"))
  
  # ROI to molecules
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay3", to = "Assay2", 
                                features = "annotation")
  vrMainAssay(merged_object) <- "MolAssay"
  expect_contains(merged_object$annotation,
                  c("DCIS_Subtype1", "DCIS_Subtype2", "Immune"))
  
  # return
  expect_equal(1,1L)
})