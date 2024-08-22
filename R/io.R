####
# InDisk Methods ####
####

saveVoltRonInDisk <- function (object, assay = NULL, dir = "my_h5_se", prefix = "", format = "hdf5", replace = FALSE, 
                             chunkdim = NULL, level = NULL, as.sparse = NA, verbose = FALSE) 
{
  if (!is(object, "VoltRon")) 
    stop("'object' must be a VoltRon object")
  if (!isSingleString(dir)) 
    stop("'dir' must be a single string specifying the path ", 
              "to the directory where to save the ", class(object), 
              " object (the directory will be created if needed)")
  if (!isSingleString(prefix)) 
    stop("'prefix' must be a single string")
  if (!isTRUEorFALSE(replace)) 
    stop("'replace' must be TRUE or FALSE")
  if (!dir.exists(dir)) {
    create_dir(dir)
  }
  else if (prefix == "") {
    replace_dir(dir, replace)
  }
  rds_path <- file.path(dir, paste0(prefix, "se.rds"))
  indisk_path <- file.path(dir, paste0(prefix, "assays.", ifelse(format == "hdf5", "h5", "zarr")))
  if (prefix != "") 
    check_and_delete_files(rds_path, indisk_path, replace)
  .write_VoltRonInDisk(object, assay = assay, format = format, rds_path = rds_path, indisk_path = indisk_path, 
                       chunkdim = chunkdim, level = level, as.sparse = as.sparse, 
                       verbose = verbose)
}

.write_VoltRonInDisk <- function(object, assay = NULL, format, rds_path, indisk_path, chunkdim=NULL, level=NULL, as.sparse=NA, verbose=FALSE)
{
  if (!is(object, "VoltRon"))
    stop("'object' must be a VoltRon object")
  
  if (!isSingleString(rds_path) || rds_path == "")
    stop("'rds_path' must be a a non-empty string ",
         "specifying the path to the RDS file ",
         "where to write the ", class(object), " object")
  
  if (!isSingleString(indisk_path) || indisk_path == "")
    stop("'indisk_path' must be a a non-empty string ",
         "specifying the path to the HDF5 file ",
         "where to write the assays of the ", class(object), " object")
  
  if (!isTRUEorFALSE(verbose))
    stop("'verbose' must be TRUE or FALSE")
  
  if(format == "hdf5"){
    object <- write_h5_samples(object, assay = assay, h5_path = indisk_path, chunkdim, level, as.sparse, verbose)
  } else if(format == "zarr"){
    object <- write_zarr_samples(object, assay = assay, zarr_path = indisk_path, chunkdim, level, as.sparse, verbose)
  } else {
    stop("'format' should be either 'hdf5' or 'zarr'")
  }
  # .serialize_HDF5VoltRon(x, rds_path, verbose)
  invisible(object)
}

####
# HDF5 Support ####
####

write_h5_samples <- function(object, assay = NULL, h5_path, chunkdim, level,
                             as.sparse, verbose)
{
  # check packages
  if(!requireNamespace('HDF5Array'))
    stop("Please install HDF5Array package!")
  if(!requireNamespace('rhdf5'))
    stop("Please install rhdf5 package!")
  
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
    a <- vrData(assay_object)
    a <- HDF5Array::writeHDF5Array(a, 
                                   h5_path, 
                                   name = paste0(assy, "/data"),
                                   chunkdim=chunkdim, 
                                   level=level,
                                   as.sparse=as.sparse,
                                   with.dimnames=FALSE,
                                   verbose=verbose)
    assay_object@rawdata <- a  
    assay_object@normdata <- a
    
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

writeHDF5ArrayInImage <- function(object, 
                                  h5_path,
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
      # file.h5$create_group(paste0(name, "/", spat))
      cat(paste0("  Writing '", name, "' image channel '", ch, "' for spatial system '", spat,"' \n"))
      rhdf5::h5createGroup(h5_path, group = paste0(name, "/", spat))
      
      # get image and write to h5
      img <- vrImages(object, name = spat, channel = ch, as.raster = TRUE)
      img <- HDF5Array::writeHDF5Array(img, 
                                       h5_path, 
                                       name = paste0(name, "/", spat, "/", ch),
                                       chunkdim=chunkdim, 
                                       level=level,
                                       as.sparse=as.sparse,
                                       with.dimnames=FALSE,
                                       verbose=verbose)
      vrImages(object, name = spat, channel = ch) <- img 
    } 
  }
  
  return(object)
}

####
# ZARR Support ####
####

write_zarr_samples <- function(object, assay = NULL, zarr_path, chunkdim, level,
                             as.sparse, verbose)
{
  # check packages
  if(!requireNamespace('DelayedArray'))
    stop("Please install DelayedArray package!")
  if(!requireNamespace('pizzarr'))
    stop("Please install pizzarr package!")
  
  # sample metadata
  sample_metadata <- SampleMetadata(object)
  
  # assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # open h5
  cat(paste0("ZARR array: ", zarr_path, "\n"))
  zarr.array <- pizzarr::zarr_open(store = zarr_path)
  
  # iterate over assays
  for (assy in assay_names) {
    
    # get assay object
    assay_object <- object[[assy]]
    
    # create assay group in h5
    zarr.array$create_group(assy)
    
    # get data and write
    cat(paste0("  Writing '", assy, "' data \n"))
    a <- vrData(assay_object)
    # a <- HDF5Array::writeHDF5Array(a, 
    #                                h5_path, 
    #                                name = paste0(assy, "/data"),
    #                                chunkdim=chunkdim, 
    #                                level=level,
    #                                as.sparse=as.sparse,
    #                                with.dimnames=FALSE,
    #                                verbose=verbose)
    zarr.array$create_dataset(name = paste0(assy, "/data"),
                              data = a, 
                              shape=dim(a))
    # assay_object@rawdata <- a  
    # assay_object@normdata <- a
    
    # # get image data and write
    # assay_object <- writeHDF5ArrayInImage(object = assay_object, 
    #                                       h5_path,
    #                                       name = assy,
    #                                       chunkdim=chunkdim, 
    #                                       level=level,
    #                                       as.sparse=as.sparse,
    #                                       with.dimnames=FALSE,
    #                                       verbose=verbose)
    
    # write assay back
    object[[assy]] <- assay_object
    
  }
  object
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

replace_dir <- function(dir, replace){
  if (!replace) 
    stop("Directory \"", dir, "\" already exists. ", 
         "Use 'replace=TRUE' to replace it. ", "Its content will be lost!")
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