mapping_parameters_file <- system.file(package = "VoltRon", "extdata/xenium_mappingparameters.rds")
mapping_parameters <- readRDS(mapping_parameters_file)
mapping_parameters$nonrigid <- "OpenST"
mapping_parameters_nonrigid_file <- system.file(package = "VoltRon", "extdata/xenium_mappingparameters_nonrigid.rds")
mapping_parameters_nonrigid <- readRDS(mapping_parameters_nonrigid_file)
mapping_parameters_nonrigid$nonrigid <- "OpenST"

test_that("registeration rigid + affine", {
  
  # get data
  data("xenium_data")
  
  # non-interactive
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters, 
                      interactive = FALSE)
  
  # check mapping parameters
  expect_error({
    registerSpatialData(reference_spatdata = xenium_data, 
                        query_spatdata = xenium_data, 
                        interactive = FALSE)
  })
  
  # change parameters
  mapping_parameters$Matcher <- "BRUTE-FORCE"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters, 
                      interactive = FALSE)
  mapping_parameters$Method <- "Homography"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters, 
                      interactive = FALSE)
  mapping_parameters$Method <- "Homography"
  mapping_parameters$Matcher <- "FLANN"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters, 
                      interactive = FALSE)
  
  
  # return
  expect_equal(1,1L)
})

test_that("registeration non-rigid", {
  
  # get data
  data("xenium_data")
  
  # non-interactive
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  
  # change parameters
  mapping_parameters_nonrigid$Matcher <- "BRUTE-FORCE"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  mapping_parameters_nonrigid$Method <- "Homography"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  mapping_parameters_nonrigid$Method <- "Homography"
  mapping_parameters_nonrigid$Matcher <- "FLANN"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  
  # return
  expect_equal(1,1L)
})

test_that("registeration non-rigid simpleitk", {
  
  # get data
  data("xenium_data")
  
  # simpleitk 
  mapping_parameters_nonrigid$nonrigid <- "SimpleITK"
  
  # non-interactive
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  
  # change parameters
  mapping_parameters_nonrigid$Matcher <- "BRUTE-FORCE"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  mapping_parameters_nonrigid$Method <- "Homography"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  mapping_parameters_nonrigid$Method <- "Homography"
  mapping_parameters_nonrigid$Matcher <- "FLANN"
  registerSpatialData(reference_spatdata = xenium_data, 
                      query_spatdata = xenium_data, 
                      mapping_parameters = mapping_parameters_nonrigid, 
                      interactive = FALSE)
  
  # return
  expect_equal(1,1L)
})