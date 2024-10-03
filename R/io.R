####
# save ####
####

saveVoltRon <- function (object, 
                         assay = NULL,
                         format = c("InMemoryVoltRon", "HDF5VoltRon", "ZarrVoltRon"), 
                         output = "my_se", 
                         replace = FALSE, 
                         chunkdim = NULL, 
                         level = NULL, 
                         as.sparse = NA, 
                         verbose = FALSE) 
{
  # check object
  if (!is(object, "VoltRon")) 
    stop("'object' must be a VoltRon object")
  
  # check output
  if (!isSingleString(output)) 
    stop("'output' must be a single string specifying the path ", 
         "to the directory where to save the ", class(object), 
         " object (the directory will be created if needed)")
  
  # check if the object was previously saved on disk
  paths <- .get_unique_links(object)
  paths <- paths[paths != "try-error"]
  if(length(paths) > 1){
    stop("There are multiple paths that this VoltRon object is saved to, cannot write!")
  } else if(length(paths) == 1){
    message("There are existing paths in the object, using those instead of the provided 'ondisk_path'")
    format <- ifelse(grepl(".zarr$", paths), "ZarrVoltRon", "HDF5VoltRon")
    output <- base::dirname(paths)
    replace <- FALSE
  } else {
    if(length(format) > 1){
      message("No paths are found in the object, and no format is chosen, saving as rds only!")
      format <- "InMemoryVoltRon"
    }
  }
  
  # save VoltRon on disk
  if(format != "InMemoryVoltRon"){
    
    # create or replace output folder
    if (!isTRUEorFALSE(replace)) 
      stop("'replace' must be TRUE or FALSE")
    if (!dir.exists(output)) {
      create_dir(output)
    } else if(replace){
      replace_dir(output) 
    }
    
    # determine format
    switch(format,
           HDF5VoltRon = {
             ondisk_path <- file.path(output, "assays.h5")
           }, 
           ZarrVoltRon = {
             ondisk_path <- file.path(output, "assays.zarr")
           })
    rds_path <- file.path(output, "se.rds")
    
    
    # write on disk
    object <- .write_VoltRon(object, assay = assay, format = format, rds_path = rds_path, ondisk_path = ondisk_path, 
                             chunkdim = chunkdim, level = level, as.sparse = as.sparse, 
                             verbose = verbose)
    
    # serialize rds file
    .serialize_VoltRonObject(object, rds_path, verbose = verbose)
    
  # save VoltRon on memory
  } else {
    rds_path <- paste0(output, "_", paste0(prefix, "se.rds"))
    saveRDS(object, file = rds_path)
  }
  
  # return
  object
}

####
# load ####
####

loadVoltRon <- function(dir="my_se")
{
  # check dir
  if (!isSingleString(dir))
    stop(wmsg("'dir' must be a single string specifying the path ",
              "to the directory containing ", .THE_EXPECTED_STUFF))
  if (!dir.exists(dir)) {
    if (file.exists(dir))
      stop(wmsg("\"", dir, "\" is a file, not a directory"))
    stop(wmsg("directory \"", dir, "\" not found"))
  }
  
  # get rds path
  rds_path <- file.path(dir, paste0("se.rds"))
  
  # load h5/zarr store
  ans <- try(.read_VoltRon(rds_path), silent=TRUE)
  if (inherits(ans, "try-error"))
    stop_if_bad_dir(dir, prefix = "")
  ans
}

####
# ondisk Methods ####
####

.serialize_VoltRonObject <- function(object, rds_path, verbose)
{
  # assay_names 
  assay_names <- vrAssayNames(object, assay = "all")
  
  # update all assays
  for(assy in assay_names)
    object[[assy]] <- shorten_assay_links(object[[assy]])
  
  # verbose and save rds
  if (verbose)
    message("Serialize ", class(object), " object to ",
            ifelse(file.exists(rds_path), "existing ", ""),
            "RDS file:\n  ", rds_path)
  saveRDS(object, file=rds_path)
}

shorten_assay_links <- function(object)
{
  # data
  
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      
      object@data[[feat]] <- modify_seeds(object@data[[feat]],
                                     function(x) {
                                       x@filepath <- basename(x@filepath)
                                       x
                                     })
      object@data[[paste0(feat, "_norm")]] <- modify_seeds(object@data[[paste0(feat, "_norm")]],
                                      function(x) {
                                        x@filepath <- basename(x@filepath)
                                        x
                                      })  
      
    }
    
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    object@rawdata <- modify_seeds(object@rawdata,
                                   function(x) {
                                     x@filepath <- basename(x@filepath)
                                     x
                                   })
    object@normdata <- modify_seeds(object@normdata,
                                    function(x) {
                                      x@filepath <- basename(x@filepath)
                                      x
                                    })  
  }
  
  # images
  object <- shorten_assay_links_images(object)
  
  # return
  object
}

shorten_assay_links_images <- function(object){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    for(ch in channels){
      
      img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
      img <- modify_seeds(img,
                          function(x) {
                            ImageArray::filepath(x) <- basename(ImageArray::filepath(x))
                            x
                          })
      vrImages(object, name = spat, channel = ch) <- img 
    } 
  }
  
  # return
  object
}

modify_seeds <- function (x, FUN, ...) 
{
  if (is(x, "DelayedUnaryOp")) {
    x@seed <- modify_seeds(x@seed, FUN, ...)
  }
  else if (is(x, "DelayedNaryOp")) {
    x@seeds <- lapply(x@seeds, modify_seeds, FUN, ...)
  }
  else {
    x <- FUN(x, ...)
  }
  x
}

.write_VoltRon <- function(object, assay = NULL, format, rds_path, ondisk_path, chunkdim=NULL, level=NULL, as.sparse=NA, verbose=FALSE)
{
  # check object
  if (!is(object, "VoltRon"))
    stop("'object' must be a VoltRon object")
  
  # check output and other related arguments
  if (!isSingleString(rds_path) || rds_path == "")
    stop("'rds_path' must be a a non-empty string ",
         "specifying the path to the RDS file ",
         "where to write the ", class(object), " object")
  if (!isSingleString(ondisk_path) || ondisk_path == "")
    stop("'ondisk_path' must be a a non-empty string ",
         "specifying the path to the HDF5 file ",
         "where to write the assays of the ", class(object), " object")
  if (!isTRUEorFALSE(verbose))
    stop("'verbose' must be TRUE or FALSE")
  
  if(format == "HDF5VoltRon"){
    object <- write_h5_samples(object, assay = assay, h5_path = ondisk_path, chunkdim, level, as.sparse, verbose)
  } else if(format == "ZarrVoltRon"){
    object <- write_zarr_samples(object, assay = assay, zarr_path = ondisk_path, chunkdim, level, as.sparse, verbose)
  } else {
    stop("'format' should be either 'HDF5VoltRon' or 'ZarrVoltRon'")
  }
  invisible(object)
}

.read_VoltRon <- function(rds_path)
{
  # check rds file
  if (!file.exists(rds_path))
    stop(wmsg("file not found: ", rds_path))
  if (dir.exists(rds_path))
    stop(wmsg("'", rds_path, "' is a directory, not a file"))
  
  # check VoltRon object
  ans <- readRDS(rds_path)
  if (!is(ans, "VoltRon"))
    stop(wmsg("the object serialized in \"", rds_path, "\" is not ",
              "a VoltRon object"))
  
  # get dir name
  dir <- dirname(rds_path)
  
  # restore assay links
  ans@assays <- restore_absolute_assay2h5_links(ans@assays, dir)
  
  # 
  ans
}

####
## HDF5 Support ####
####

write_h5_samples <- function(object, assay = NULL, h5_path, chunkdim, level,
                             as.sparse, verbose)
{
  # check packages
  if(!requireNamespace('HDF5Array'))
    stop("Please install HDF5Array package!")
  if(!requireNamespace('rhdf5'))
    stop("Please install rhdf5 package!")
  if(!requireNamespace('ImageArray'))
    stop("Please install ImageArray package!")
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # open h5
  # file.h5 <- hdf5r::H5File$new(h5_path, mode="w")
  cat(paste0("HDF5 file: ", h5_path, "\n"))
  rhdf5::h5createFile(h5_path)
  
  # iterate over assays
  for (assy in assay_names) {
    
    # get assay object
    assay_object <- object[[assy]]
    
    # create assay group in h5
    # file.h5$create_group(assy)
    rhdf5::h5createGroup(h5_path, group = assy)
    
    # get data and write
    cat(paste0("  Writing '", assy, "' data \n"))
    assay_object <- writeHDF5ArrayInVrData(object = assay_object, 
                                           h5_path,
                                           name = assy,
                                           chunkdim=chunkdim, 
                                           level=level,
                                           as.sparse=as.sparse,
                                           with.dimnames=TRUE,
                                           verbose=verbose)
    
    # get image data and write
    assay_object <- writeHDF5ArrayInImage(object = assay_object, 
                                          h5_path,
                                          name = assy,
                                          chunkdim=chunkdim, 
                                          level=level,
                                          as.sparse=as.sparse,
                                          with.dimnames=FALSE,
                                          verbose=verbose)
    
    # write assay back
    object[[assy]] <- assay_object
    
  }
  object
}

writeHDF5ArrayInVrData <- function(object, 
                                   h5_path,
                                   name,
                                   chunkdim, 
                                   level,
                                   as.sparse,
                                   with.dimnames=FALSE,
                                   verbose){
 
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      
      # raw data
      a <- vrData(object, feat_type = feat, norm = FALSE)
      if(!inherits(a, "DelayedArray")){
        a <- HDF5Array::writeHDF5Array(a, 
                                       h5_path, 
                                       name = paste0(name, "/", feat),
                                       chunkdim=chunkdim, 
                                       level=level,
                                       as.sparse=as.sparse,
                                       with.dimnames=with.dimnames,
                                       verbose=verbose)
        object@data[[feat]] <- a   
      }
      
      # normalized data
      a <- vrData(object, feat_type = feat, norm = TRUE)
      if(!inherits(a, "DelayedArray")){
        a <- HDF5Array::writeHDF5Array(a, 
                                       h5_path, 
                                       name = paste0(name, "/", feat, "_norm"),
                                       chunkdim=chunkdim, 
                                       level=level,
                                       as.sparse=as.sparse,
                                       with.dimnames=with.dimnames,
                                       verbose=verbose)
        object@data[[paste0(feat, "_norm")]] <- a  
      }
      
    }
    
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    
    # raw data
    a <- vrData(object, norm = FALSE)
    if(!inherits(a, "DelayedArray")){
      a <- HDF5Array::writeHDF5Array(a, 
                                     h5_path, 
                                     name = paste0(name, "/rawdata"),
                                     chunkdim=chunkdim, 
                                     level=level,
                                     as.sparse=as.sparse,
                                     with.dimnames=TRUE,
                                     verbose=verbose)
      object@rawdata <- a 
    }
    
    # normalized data
    a <- vrData(object, norm = TRUE)
    if(!inherits(a, "DelayedArray")){
      a <- HDF5Array::writeHDF5Array(a, 
                                     h5_path, 
                                     name = paste0(name, "/normdata"),
                                     chunkdim=chunkdim, 
                                     level=level,
                                     as.sparse=as.sparse,
                                     with.dimnames=TRUE,
                                     verbose=verbose)
      object@normdata <- a
    }
    
  }

  return(object)
}

writeHDF5ArrayInImage <- function(object, 
                                  h5_path,
                                  name,
                                  chunkdim, 
                                  level=level,
                                  as.sparse,
                                  with.dimnames,
                                  verbose){

  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){

    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    for(ch in channels){
      
      # open group for spatial system
      # file.h5$create_group(paste0(name, "/", spat))
      cat(paste0("  Writing '", name, "' image channel '", ch, "' for spatial system '", spat,"' \n"))
      rhdf5::h5createGroup(h5_path, group = paste0(name, "/", spat))
      
      # get image and write to h5
      img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
      
      # write image
      if(!inherits(img, "DelayedArray")){
        img <- ImageArray::writeImageArray(img,
                                           output = gsub(".h5$", "", h5_path),
                                           name = paste0(name, "/", spat, "/", ch), 
                                           format = "HDF5ImageArray", 
                                           replace = FALSE, 
                                           chunkdim=chunkdim,
                                           level=level,
                                           as.sparse=as.sparse,
                                           verbose=verbose)
        vrImages(object, name = spat, channel = ch) <- img 
      }
    } 
  }
  
  return(object)
}

####
## ZARR Support ####
####

write_zarr_samples <- function(object, assay = NULL, zarr_path, chunkdim, level,
                             as.sparse, verbose)
{
  # check packages
  if(!requireNamespace('ZarrArray'))
    stop("Please install ZarrArray package!")
  if(!requireNamespace('pizzarr'))
    stop("Please install pizzarr package!")
  if(!requireNamespace('ImageArray'))
    stop("Please install ImageArray package!")
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # open h5
  cat(paste0("ZARR array: ", zarr_path, "\n"))
  zarr.array <- pizzarr::zarr_open(store = zarr_path)
  
  # iterate over assays
  for(assy in assay_names){
    
    # get assay object
    assay_object <- object[[assy]]
    
    # create assay group in h5
    zarr.array$create_group(assy)
    
    # get data and write
    cat(paste0("  Writing '", assy, "' data \n"))
    assay_object <- writeZarrArrayInVrData(object = assay_object, 
                                           zarr_path,
                                           name = assy,
                                           chunkdim=chunkdim, 
                                           level=level,
                                           as.sparse=as.sparse,
                                           with.dimnames=TRUE,
                                           verbose=verbose)
    
    # TODO: image zarr conversion does not work for bitmap arrays now
    # get image data and write
    assay_object <- writeZarrArrayInImage(object = assay_object,
                                          zarr_path,
                                          name = assy, 
                                          chunkdim=chunkdim, 
                                          level=level,
                                          as.sparse=as.sparse,
                                          with.dimnames=FALSE,
                                          verbose=verbose)
    
    # write assay back
    object[[assy]] <- assay_object
    
  }
  object
}

writeZarrArrayInVrData <- function(object, 
                                   zarr_path,
                                   name,
                                   chunkdim, 
                                   level,
                                   as.sparse,
                                   with.dimnames=FALSE,
                                   verbose){
  
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      
      # raw data
      a <- vrData(object, feat_type = feat, norm = FALSE)
      if(!inherits(a, "DelayedArray")){
        a <- ZarrArray::writeZarrArray(a, 
                                       zarr_path, 
                                       name = paste0(name, "/", feat),
                                       chunkdim=chunkdim, 
                                       level=level,
                                       as.sparse=as.sparse,
                                       with.dimnames=with.dimnames,
                                       verbose=verbose)
        object@data[[feat]] <- a   
      }
      
      # normalized data
      a <- vrData(object, feat_type = feat, norm = TRUE)
      if(!inherits(a, "DelayedArray")){
        a <- ZarrArray::writeZarrArray(a, 
                                       zarr_path, 
                                       name = paste0(name, "/", feat, "_norm"),
                                       chunkdim=chunkdim, 
                                       level=level,
                                       as.sparse=as.sparse,
                                       with.dimnames=with.dimnames,
                                       verbose=verbose)
        object@data[[paste0(feat, "_norm")]] <- a  
      }
    }
    
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    
    # raw data
    a <- vrData(object, norm = FALSE)
    if(!inherits(a, "DelayedArray")){
      a <- ZarrArray::writeZarrArray(a, 
                                     zarr_path, 
                                     name = paste0(name, "/rawdata"),
                                     chunkdim=chunkdim, 
                                     level=level,
                                     as.sparse=as.sparse,
                                     with.dimnames=TRUE,
                                     verbose=verbose)
      object@rawdata <- a   
    }
    
    # normalized data
    a <- vrData(object, norm = TRUE)
    if(!inherits(a, "DelayedArray")){
      a <- ZarrArray::writeZarrArray(a, 
                                     zarr_path, 
                                     name = paste0(name, "/normdata"),
                                     chunkdim=chunkdim, 
                                     level=level,
                                     as.sparse=as.sparse,
                                     with.dimnames=TRUE,
                                     verbose=verbose)
      object@normdata <- a 
    }
    
  }
  
  return(object)
}

writeZarrArrayInImage <- function(object, 
                                  zarr_path,
                                  name = assy,
                                  chunkdim=chunkdim, 
                                  level=level,
                                  as.sparse=as.sparse,
                                  with.dimnames=FALSE,
                                  verbose=verbose){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    for(ch in channels){
      
      # open group for spatial system
      cat(paste0("  Writing '", name, "' image channel '", ch, "' for spatial system '", spat,"' \n"))
      zarr.array <- pizzarr::zarr_open(store = zarr_path)
      zarr.array$create_group(paste0(name, "/", spat))
      
      # get image and write to h5
      img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
      
      # write image
      if(!inherits(img, "DelayedArray")){
        img <- ImageArray::writeImageArray(img,
                                           output = gsub(".zarr$", "", zarr_path),
                                           name = paste0(name, "/", spat, "/", ch), 
                                           format = "ZarrImageArray", 
                                           replace = FALSE, 
                                           chunkdim=chunkdim,
                                           level=level,
                                           as.sparse=as.sparse,
                                           verbose=verbose)
        vrImages(object, name = spat, channel = ch) <- img 
      }
    } 
  }
  
  return(object)
}

####
# Auxiliary ####
####

isTRUEorFALSE <- function (x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

isSingleString <- function (x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

create_dir <- function (dir){
  if (file.exists(dir)) 
    stop("\"", dir, "\" already exists and is a file, ", 
         "not a directory")
  if (!suppressWarnings(dir.create(dir))) 
    stop("cannot create directory \"", dir, "\"")
}

replace_dir <- function(dir){
  if (unlink(dir, recursive = TRUE) != 0L) 
    stop("failed to delete directory \"", dir, "\"")
  if (!suppressWarnings(dir.create(dir))) 
    stop("cannot create directory \"", dir, "\"")
}

check_and_delete_files <- function (rds_path, h5_path, replace) 
{
  if (dir.exists(rds_path) || dir.exists(h5_path)) 
    stop("\"", rds_path, "\" and/or \"", h5_path, "\" ", 
         "already exist and are directories, not files")
  if (!(file.exists(rds_path) || file.exists(h5_path))) 
    return()
  if (!replace) 
    stop("Files \"", rds_path, "\" and/or \"", h5_path, 
         "\" ", "already exist. Use a different 'prefix' or use ", 
         "'replace=TRUE' if you really want to replace them.")
  if (unlink(rds_path, recursive = TRUE) != 0L) 
    stop("failed to delete file \"", rds_path, "\"")
  if (unlink(h5_path, recursive = TRUE) != 0L) 
    stop("failed to delete file \"", h5_path, "\"")
}

.get_unique_links <- function(object, assay = NULL)
{
  # assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # iterate over assays
  all_links <- NULL
  for(assy in assay_names){
    
    # get data and image links
    all_links <- c(all_links, .get_unique_data_links(object[[assy]]))
    all_links <- c(all_links, .get_unique_data_links(object[[assy]]))
                           
  }
  
  # return
  unique(all_links)
}

.get_unique_data_links <- function(object){
  
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  all_links <- NULL
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      cur_path <- try(path(vrData(object, feat_type = feat, norm = FALSE)), silent = TRUE)
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
      cur_path <- try(path(vrData(object, feat_type = feat, norm = TRUE)), silent = TRUE)
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    }
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    cur_path <- try(path(vrData(object, norm = FALSE), silent = TRUE))
    all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    cur_path <- try(path(vrData(object, norm = TRUE), silent = TRUE))
    all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    
  }
  
  # return
  return(all_links)
}

.get_unique_image_links <- function(object){
 
  # links
  all_links <- NULL
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    for(ch in channels){
      cur_path <- try(path(vrImages(object, name = spat, channel = ch, as.raster = TRUE), silent = TRUE))
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    }
  }
  
  # return
  return(all_links)
}

stop_if_bad_dir <- function(dir, prefix = "")
{
  .THE_EXPECTED_STUFF <- c(
    "an HDF5-based SummarizedExperiment object ",
    "previously saved with saveHDF5SummarizedExperiment",
    "()"
  )
  if (prefix == "") {
    msg <- c("directory \"", dir, "\" does not seem to contain ",
             .THE_EXPECTED_STUFF)
  } else {
    msg <- c("Directory \"", dir, "\" does not seem to contain ",
             head(.THE_EXPECTED_STUFF, n=-1L),
             "(..., prefix=\"", prefix, "\"). ",
             "Make sure you're using the same 'prefix' ",
             "that was used when the object was saved.")
  }
  stop(wmsg(msg))
}