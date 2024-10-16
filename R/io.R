####
# save ####
####

#' saveVoltRon
#'
#' save VoltRon object in memory or on disk
#' 
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' @param format the format the object should be written: InMemoryVoltRon (rds only), HDF5VoltRon (h5), or ZarrVoltRon (zarr).
#' @param output When saving, the directory will be created if it doesn't already exist. If the directory already exists and no prefix is specified and replace is set to TRUE, then it's replaced with an empty directory.
#' @param replace When no prefix is specified, should a pre-existing directory be replaced with a new empty one? The content of the pre-existing directory will be lost!
#' @param chunkdim The dimensions of the chunks to use for writing the assay data to disk.
#' @param level The compression level to use for writing the assay data to disk.
#' @param as.sparse as.sparse
#' @param verbose verbose
#'
#' @export
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
                             chunkdim = chunkdim, level = level, as.sparse = as.sparse, verbose = verbose)
    
    # serialize rds file
    .serialize_VoltRonObject(object, rds_path, verbose = verbose)
    
  # save VoltRon on memory
  } else {
    rds_path <- paste0(output, "_", paste0("se.rds"))
    saveRDS(object, file = rds_path)
  }
  
  # return
  object
}

####
# load ####
####


#' loadVoltRon
#'
#' load VoltRon object from memory or disk
#' 
#' @param dir the directory that VoltRon object is found.
#'
#' @export
loadVoltRon <- function(dir="my_se")
{
  if(!requireNamespace('DelayedArray'))
    stop("Please install DelayedArray package!")
  
  # check dir
  if (!isSingleString(dir))
    stop(paste0("'dir' must be a single string specifying the path ",
              "to the directory containing an rds and/or .h5/.zarr file!"))
  if (!dir.exists(dir)) {
    if (file.exists(dir))
      stop(paste0("\"", dir, "\" is a file, not a directory"))
    stop(paste0("directory \"", dir, "\" not found"))
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

#' .serialize_VoltRonObject
#'
#' @param object a VoltRon object
#' @param rds_path the path to rds
#' @param verbose verbose
#'
#' @noRd
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

#' modify_seeds
#'
#' @noRd
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

#' .write_VoltRon
#'
#' @noRd
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

####
## HDF5 Support ####
####

#' write_h5_samples
#'
#' @noRd
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
  cat(paste0("HDF5 file: ", h5_path, "\n"))
  rhdf5::h5createFile(h5_path)
  
  # iterate over assays
  for (assy in assay_names) {
    
    # get assay object
    assay_object <- object[[assy]]
    
    # create assay group in h5
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

#' writeHDF5ArrayInVrData
#'
#' @noRd
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
      if(!inherits(a, c("DelayedArray", "IterableMatrix"))){
        # a <- HDF5Array::writeHDF5Array(a, 
        #                                h5_path, 
        #                                name = paste0(name, "/", feat),
        #                                chunkdim=chunkdim, 
        #                                level=level,
        #                                as.sparse=as.sparse,
        #                                with.dimnames=with.dimnames,
        #                                verbose=verbose)
        if(!inherits(a, "dgCMatrix"))
          a <- as(a, "dgCMatrix")
        a <- BPCells::write_matrix_hdf5(a, 
                                        path = h5_path, 
                                        group = paste0(name, "/", feat))
        # chunk_size = chunkdim)
        object@data[[feat]] <- a   
        
      }
      
      # normalized data
      a <- vrData(object, feat_type = feat, norm = TRUE)
      if(!inherits(a, c("DelayedArray", "IterableMatrix"))){
        if(!inherits(a, "dgCMatrix"))
          a <- as(a, "dgCMatrix")
        a <- BPCells::write_matrix_hdf5(a, 
                                        path = h5_path, 
                                        group = paste0(name, "/", feat, "_norm"))
        # chunk_size = chunkdim)
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

#' writeHDF5ArrayInImage
#'
#' @noRd
writeHDF5ArrayInImage <- function(object, 
                                  h5_path,
                                  name,
                                  chunkdim, 
                                  level,
                                  as.sparse,
                                  with.dimnames,
                                  verbose){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    if(!grepl("No Channels", channels)){
      for(ch in channels){
        
        # open group for spatial system
        # file.h5$create_group(paste0(name, "/", spat))
        cat(paste0("  Writing '", name, "' image channel '", ch, "' for spatial system '", spat,"' \n"))
        rhdf5::h5createGroup(h5_path, group = paste0(name, "/", spat))
        
        # get image and write to h5
        img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
        
        # write image
        if(!inherits(img, "Image_Array")){
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
  }
  
  return(object)
}

####
## ZARR Support ####
####

#' write_zarr_samples
#'
#' @noRd
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

#' writeZarrArrayInVrData
#'
#' @noRd
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

#' writeZarrArrayInImage
#'
#' @noRd
writeZarrArrayInImage <- function(object, 
                                  zarr_path,
                                  name ,
                                  chunkdim, 
                                  level,
                                  as.sparse,
                                  with.dimnames=FALSE,
                                  verbose=verbose){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    if(!grepl("No Channels", channels)){
      for(ch in channels){
        
        # open group for spatial system
        cat(paste0("  Writing '", name, "' image channel '", ch, "' for spatial system '", spat,"' \n"))
        zarr.array <- pizzarr::zarr_open(store = zarr_path)
        zarr.array$create_group(paste0(name, "/", spat))
        
        # get image and write to h5
        img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
        
        # write image
        if(!inherits(img, "Image_Array")){
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
  }
  
  return(object)
}

#' .read_VoltRon
#'
#' @noRd
.read_VoltRon <- function(rds_path)
{
  # check rds file
  if (!file.exists(rds_path))
    stop(paste0("file not found: ", rds_path))
  if (dir.exists(rds_path))
    stop(paste0("'", rds_path, "' is a directory, not a file"))
  
  # check VoltRon object
  object <- readRDS(rds_path)
  if (!is(object, "VoltRon"))
    stop(paste0("the object serialized in \"", rds_path, "\" is not ",
                "a VoltRon object"))
  
  # get dir name
  dir <- dirname(rds_path)
  
  # assay_names 
  assay_names <- vrAssayNames(object, assay = "all")
  
  # restore assay links
  for(assy in assay_names)
    object[[assy]] <- restore_absolute_assay_links(object[[assy]], dir)
  
  # return object
  object
}

####
## get links ####
####

#' .get_unique_links
#'
#' @noRd
.get_unique_links <- function(object, assay = NULL)
{
  # assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # iterate over assays
  all_links <- NULL
  for(assy in assay_names){
    
    # get data and image links
    all_links <- c(all_links, .get_unique_data_links(object[[assy]]))
    all_links <- c(all_links, .get_unique_image_links(object[[assy]]))
    
  }
  
  # return
  unique(all_links)
}

#' .get_unique_data_links
#'
#' @noRd
.get_unique_data_links <- function(object){
  
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  all_links <- NULL
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      cur_path <- try(getPath(vrData(object, feat_type = feat, norm = FALSE)), silent = TRUE)
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
      cur_path <- try(getPath(vrData(object, feat_type = feat, norm = TRUE)), silent = TRUE)
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    }
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    cur_path <- try(getPath(vrData(object, norm = FALSE)), silent = TRUE)
    all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    cur_path <- try(getPath(vrData(object, norm = TRUE)), silent = TRUE)
    all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    
  }
  
  # return
  return(all_links)
}

#' .get_unique_image_links
#'
#' @noRd
.get_unique_image_links <- function(object){
  
  # links
  all_links <- NULL
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    for(ch in channels){
      cur_path <- try(DelayedArray::path(vrImages(object, name = spat, channel = ch, as.raster = TRUE)), silent = TRUE)
      all_links <- c(all_links, ifelse(is(cur_path, "try-error"), "try-error", cur_path))
    }
  }
  
  # return
  return(all_links)
}

getPath <- function(object){
  if(inherits(object, "DelayedArray")){
    return(DelayedArray::path(object))
  } else if(inherits(object, "IterableMatrix")){
    return(getIterableMatrixPath(object))
  } else {
    stop()
  }
}

getIterableMatrixPath <- function(object){
  if(!inherits(object, "IterableMatrix")){
    stop("object should be an object of IterableMatrix")
  }
  slot_names <- slotNames(object)
  if("path" %in% slot_names){
    return(object@path)
  } else if("matrix" %in% slot_names){
    return(getIterableMatrixPath(object@matrix))
  } else if("matrix_list" %in% slot_names){
    return(unlist(lapply(object@matrix_list, getIterableMatrixPath)))
  }
}

####
## shorten links ####
####

#' shorten_assay_links
#'
#' @noRd
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
                                       shorten_assay_links_data(x)
                                     })
      object@data[[paste0(feat, "_norm")]] <- modify_seeds(object@data[[paste0(feat, "_norm")]],
                                      function(x) {
                                        shorten_assay_links_data(x)
                                      })  
      
    }
    
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    object@rawdata <- modify_seeds(object@rawdata,
                                   function(x) {
                                     shorten_assay_links_data(x)
                                   })
    object@normdata <- modify_seeds(object@normdata,
                                    function(x) {
                                      shorten_assay_links_data(x)
                                    })  
  }
  
  # images
  object <- shorten_assay_links_images(object)
  
  # return
  object
}

#' shorten_assay_links_images
#'
#' @noRd
shorten_assay_links_images <- function(object){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    if(!grepl("No Channels", channels)){
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
  }
  
  # return
  return(object)
}

shorten_assay_links_data <- function(object){
  if(inherits(object, "DelayedArray")){
    object@filepath <- basename(object@filepath)
  } else if(inherits(object, "IterableMatrix")){
    object <- shorten_assay_links_bpcells(object)
  }
  return(object)
}

shorten_assay_links_bpcells <- function(object){
  if(!inherits(object, "IterableMatrix")){
    stop("object should be an object of IterableMatrix")
  }
  slot_names <- slotNames(object)
  if("path" %in% slot_names){
    object@path <- basename(object@path)
  } else if("matrix" %in% slot_names){
    object@matrix <- shorten_assay_links_bpcells(object@matrix)
  } else if("matrix_list" %in% slot_names){
    object_list <- object@matrix_list
    for(i in 1:length(object_list))
      object_list[[i]] <- shorten_assay_links_bpcells(object_list[[i]])
  }
  return(object)
}

####
## restore links ####
####

#' restore_absolute_assay_links
#'
#' @noRd
restore_absolute_assay_links <- function(object, dir){
  
  # check if there is a data or rawdata slot in assay object
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  
  # get data with a specific feature
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    
    feature_types <- vrFeatureTypeNames(object)
    for(feat in feature_types){
      
      object@data[[feat]] <- modify_seeds(object@data[[feat]],
                                          function(x) {
                                            restore_absolute_links(x, dir)
                                          })
      object@data[[paste0(feat, "_norm")]] <- modify_seeds(object@data[[paste0(feat, "_norm")]],
                                                           function(x) {
                                                             restore_absolute_links(x, dir)
                                                           })  
      
    }
    
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    object@rawdata <- modify_seeds(object@rawdata,
                                   function(x) {
                                     restore_absolute_links(x, dir)
                                     x
                                   })
    object@normdata <- modify_seeds(object@normdata,
                                    function(x) {
                                      restore_absolute_links(x, dir)
                                    })  
  }
  
  # images
  object <- restore_absolute_assay_links_images(object, dir)
  
  # return
  object
}

#' restore_absolute_assay_links_images
#'
#' @noRd
restore_absolute_assay_links_images <- function(object, dir){
  
  # for each spatial system
  spatial_names <- vrSpatialNames(object)
  for(spat in spatial_names){
    
    # for each channel
    channels <- vrImageChannelNames(object, name = spat)
    if(!grepl("No Channels", channels)){
      for(ch in channels){
        
        img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
        img <- modify_seeds(img,
                            function(x) {
                              ImageArray::filepath(x) <- restore_absolute_links_images(ImageArray::filepath(x), dir)
                              x
                            })
        vrImages(object, name = spat, channel = ch) <- img 
      } 
    }
  }
  
  # return
  return(object)
}

#' restore_absolute_links
#'
#' @noRd
restore_absolute_links <- function(x, dir){
  if(inherits(x, "DelayedArray")){
    x@filepath <- basename(x@filepath)
  } else if(inherits(x, "IterableMatrix")){
    x@path <- basename(x@path)
  }

  # check object
  if (!is(x, c("Array")) && !is(x, c("IterableMatrix")))
    stop(message("object is not DelayedArray"))
  
  # get path
  if(inherits(x, "DelayedArray")){
    file_path <- file.path(dir, x@filepath)
  } else if(inherits(x, "IterableMatrix")){
    file_path <- file.path(dir, x@path)
  }

  ## file_path_as_absolute() will fail if the file does
  ## not exist.
  if (!file.exists(file_path))
    stop(message("Object points to an HDF5 file ",
              "that does not exist: ", file_path))
  if(inherits(x, "DelayedArray")){
    x@filepath <- file_path_as_absolute(file_path)
    msg <- validate_absolute_path(x@filepath, paste0("'filepath' slot of Object"))
  } else if(inherits(x, "IterableMatrix")){
    x@path <- file_path_as_absolute(file_path)
    msg <- validate_absolute_path(x@path, paste0("'filepath' slot of Object"))
  }

  # validate
  if (!isTRUE(msg))
    stop(paste0(msg))
  
  # return
  x
}

#' restore_absolute_links_images
#'
#' @noRd
restore_absolute_links_images <- function(file_path, dir){
  file_path <- basename(file_path)
  
  # get path
  file_path <- file.path(dir, file_path)
  
  ## file_path_as_absolute() will fail if the file does
  ## not exist.
  if (!file.exists(file_path))
    stop(message("file_path doesnt exist"))
  file_path <- file_path_as_absolute(file_path)
  
  # validate
  msg <- validate_absolute_path(file_path, paste0("'filepath' slot of object"))
  if (!isTRUE(msg))
    stop(paste0(msg))
  file_path
}



####
# Auxiliary ####
####

#' isTRUEorFALSE
#'
#' @noRd
isTRUEorFALSE <- function (x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

#' isSingleString
#'
#' @noRd
isSingleString <- function (x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

#' create_dir
#'
#' @noRd
create_dir <- function (dir){
  if (file.exists(dir)) 
    stop("\"", dir, "\" already exists and is a file, ", 
         "not a directory")
  if (!suppressWarnings(dir.create(dir))) 
    stop("cannot create directory \"", dir, "\"")
}

#' replace_dir
#'
#' @noRd
replace_dir <- function(dir){
  if (unlink(dir, recursive = TRUE) != 0L) 
    stop("failed to delete directory \"", dir, "\"")
  if (!suppressWarnings(dir.create(dir))) 
    stop("cannot create directory \"", dir, "\"")
}

#' check_and_delete_files
#'
#' @noRd
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

#' stop_if_bad_dir
#'
#' @noRd
stop_if_bad_dir <- function(dir, prefix = "")
{
  .THE_EXPECTED_STUFF <- c(
    "an OnDisk based VoltRon object ",
    "previously saved with saveVoltRon",
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
  stop(paste0(msg))
}

#' file_path_as_absolute
#'
#' @noRd
file_path_as_absolute <- function (x) 
{
  if (length(x) != 1L) 
    stop("'x' must be a single character string")
  if (!file.exists(epath <- path.expand(x))) 
    stop(gettextf("file '%s' does not exist", x), domain = NA)
  normalizePath(epath, "/", TRUE)
}

#' validate_absolute_path
#'
#' @noRd
validate_absolute_path <- function(path, what="'path'")
{
  if (!(isSingleString(path) && nzchar(path)))
    return(paste0(what, " must be a single non-empty string"))
  ## Check that 'path' points to an HDF5 file that is accessible.
  if (!file.exists(path))
    return(paste0(what, " (\"", path, "\") must be the path to ",
                  "an existing HDF5 file"))
  if (dir.exists(path))
    return(paste0(what, " (\"", path, "\") must be the path to ",
                  "an HDF5 file, not a directory"))
  if (path != file_path_as_absolute(path))
    return(paste0(what, " (\"", path, "\") must be the absolute ",
                  "canonical path the HDF5 file"))
  TRUE
}