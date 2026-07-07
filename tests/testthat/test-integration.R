test_that("cell to ROI", {
  
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
})
  
test_that("cell to ROI, all features", {
  # cell to ROI, all features
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay1", to = "Assay3", new_feature_name = "main_pseudo2")
  expect_identical(c("main", "main_pseudo2"),
                   vrFeatureTypeNames(merged_object, assay = "Assay3"))
  vrMainFeatureType(merged_object, assay = "Assay3") <- "main_pseudo2"
  vrMainAssay(merged_object) <- "ROIAssay"
  expect_identical(vrFeatures(merged_object, assay = "ROIAssay"),
                   vrFeatures(merged_object, assay = "CellAssay"))
})

test_that("ROI to cell", {
  data("merged_object")
  merged_object <- transferData(merged_object, from = "Assay3", to = "Assay1", 
                                features = "annotation")
  expect_contains(merged_object$annotation,
                  c("DCIS_Subtype1", "DCIS_Subtype2", "Immune"))
})
  
test_that("ROI to molecules", {
  
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

test_that("Single-cell to VoltRon", {
  
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("SingleCellExperiment")
  
  data(pbmc_small, package = "SeuratObject")
  expect_error(.getTransferReference(pbmc_small, sc.assay = "RNA2"))
  
  pbmc_small <- Seurat::as.SingleCellExperiment(pbmc_small)
  expect_error(.getTransferReference(pbmc_small, sc.assay = "RNA2"))
  
  # make pseudo voltron object
  datax <- SummarizedExperiment::assay(pbmc_small)
  coords <- matrix(runif(2*ncol(datax)), ncol = 2)
  rownames(coords) <- colnames(datax)
  vr_pseudo <- formVoltRon(data = datax, coords = coords)
  
  # transfer counts 
  vr_pseudo <- transferData(vr_pseudo, from = pbmc_small)
  datax <- vrData(vr_pseudo, feat_type = "counts_import")
  expect_equal(dim(datax), dim(pbmc_small))
  
  # transfer some features
  features <- rownames(pbmc_small)[1:20]
  vr_pseudo <- transferData(vr_pseudo, from = pbmc_small, 
                            features = features,
                            new_feature_name = "temp")
  datax <- vrData(vr_pseudo, feat_type = "temp")
  expect_equal(dim(datax), c(length(features), ncol(pbmc_small)))
  
  # dont mix features and metadata features
  features <- c(rownames(pbmc_small)[1:20], "groups")
  expect_error(
    vr_pseudo <- transferData(vr_pseudo, from = pbmc_small, 
                              features = features,
                              new_feature_name = "temp"), 
    "Data and Metadata features cannot be transfered in the same time"
  )
  
  # dont transfer more than one metadata feature
  features <- c("letter.idents", "groups")
  expect_error(
    vr_pseudo <- transferData(vr_pseudo, from = pbmc_small, 
                              features = features,
                              new_feature_name = "temp"), 
    "Only one metadata feature can be transfered at a time"
  )
})