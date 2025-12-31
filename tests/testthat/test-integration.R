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

test_that("ROI to ROI", {
  
  # update object
  merged_object2 <- merged_object
  merged_object2[["Assay3"]]@image$main@coords <- rbind(merged_object2[["Assay3"]]@image$main@coords,
                                                        merged_object2[["Assay3"]]@image$main@coords[3,])
  rownames(merged_object2[["Assay3"]]@image$main@coords)[4] <- "DCIS_Subtype3_Assay3"
  merged_object2[["Assay3"]]@image$main@coords <- merged_object2[["Assay3"]]@image$main@coords[-1,]
  merged_object2[["Assay3"]]@image$main@segments <- c(merged_object2[["Assay3"]]@image$main@segments,
                                                      list(merged_object2[["Assay3"]]@image$main@segments[[3]]))
  names(merged_object2[["Assay3"]]@image$main@segments)[4] <- "DCIS_Subtype3_Assay3"
  merged_object2[["Assay3"]]@image$main@segments <- merged_object2[["Assay3"]]@image$main@segments[-1]
  merged_object2@metadata@ROI <- 
    rbind(merged_object2@metadata@ROI[merged_object2@metadata@ROI$assay_id == "Assay3",],
          merged_object2@metadata@ROI[3,],
          merged_object2@metadata@ROI[merged_object2@metadata@ROI$assay_id == "Assay5",])
  rownames(merged_object2@metadata@ROI)[4] <- "DCIS_Subtype3_Assay3"
  merged_object2@metadata@ROI$id[4] <- "DCIS_Subtype3_Assay3"
  merged_object2@metadata@ROI$annotation[4] <- "temp"
  merged_object2@metadata@ROI <- merged_object2@metadata@ROI[-1,]
  
  # transfer from ROI to ROI
  merged_object2 <- transferData(merged_object2, from = "Assay3", to = "Assay5", features = "annotation", new_feature_name = "annotation2")
  expect_equal(merged_object2@metadata@ROI$annotation[3], "temp")
  expect_equal(merged_object2@metadata@ROI$annotation[4], "undefined")
  expect_equal(merged_object2@metadata@ROI$annotation[6], "DCIS_Subtype2,temp")
})