# Testing functions of assay data ####
test_that("featuredata", {
  
  # get data
  data("visium_data")
  
  # get and set feature types
  vrMainFeatureType(visium_data)
  vrMainFeatureType(visium_data, assay = "Assay1") <- "RNA"
  expect_warning(vrMainFeatureType(visium_data) <- "main")
  
  # get feature data
  expect_equal(vrFeatureTypeNames(visium_data), "RNA")
  vrFeatureData(visium_data, feat_type = "RNA")
  expect_equal(vrFeatureData(visium_data, feat_type = "Prot"), NULL)
  
  # return
  expect_equal(1,1L)
})

# Testing functions of manipulating embeddings ####
test_that("add feature", {
  
  # get data
  data("xenium_data")
  
  # add feature
  xenium_data2 <- addFeature(xenium_data,
                             data = vrData(xenium_data), 
                             feature_name = "Protein")
  expect_equal(vrFeatureTypeNames(xenium_data2), 
               c("RNA", "Protein"))
  
  # get data
  data <- vrData(xenium_data2, feat_type = c("RNA"))
  data <- vrData(xenium_data2, feat_type = c("Protein"))
  data <- vrData(xenium_data2, feat_type = c("RNA", "Protein"))
  expect_equal(rownames(data), 
               c(paste(vrFeatures(xenium_data), "RNA", sep = "_"),
                 paste(vrFeatures(xenium_data), "Protein", sep = "_")))
})