#' @importClassesFrom data.table data.table

####
# Objects and Classes ####
####

### $ methods ####

#' @method $ vrMetadata
#'
"$.vrMetadata" <- function(x, i, ...) {
  return(NULL)
}

#' @method $<- vrMetadata
#'
#' @importFrom methods new slot
"$<-.vrMetadata" <- function(x, i, ..., value) {

  # molecule metadata
  mol.metadata <- methods::slot(x, "molecule")
  if(nrow(mol.metadata) > 0)
    mol.metadata[[i]] <- value

  # cell metadata
  cell.metadata <- methods::slot(x, "cell")
  if(nrow(cell.metadata) > 0)
    cell.metadata[[i]] <- value

  # spot metadata
  spot.metadata <- methods::slot(x, "spot")
  if(nrow(spot.metadata) > 0)
    spot.metadata[[i]] <- value

  # ROI metadata
  roi.metadata <- methods::slot(x, "ROI")
  if(nrow(roi.metadata) > 0)
    roi.metadata[[i]] <- value

  # ROI metadata
  tile.metadata <- methods::slot(x, "tile")
  if(nrow(tile.metadata) > 0)
    tile.metadata[[i]] <- value

  return(methods::new("vrMetadata", molecule = mol.metadata, cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata, tile = tile.metadata))
}

#' @method $<- vrMetadata
#'
#' @importFrom methods new slot
#'
"[[<-.vrMetadata" <- function(x, i, ..., value) {

  # molecule metadata
  mol.metadata <- methods::slot(x, "molecule")
  if(nrow(mol.metadata) > 0)
    mol.metadata[[i]] <- value

  # cell metadata
  cell.metadata <- methods::slot(x, "cell")
  if(nrow(cell.metadata) > 0)
    cell.metadata[[i]] <- value

  # spot metadata
  spot.metadata <- methods::slot(x, "spot")
  if(nrow(spot.metadata) > 0)
    spot.metadata[[i]] <- value

  # ROI metadata
  roi.metadata <- methods::slot(x, "ROI")
  if(nrow(roi.metadata) > 0)
    roi.metadata[[i]] <- value

  # ROI metadata
  tile.metadata <- methods::slot(x, "tile")
  if(nrow(tile.metadata) > 0)
    tile.metadata[[i]] <- value

  return(methods::new("vrMetadata", molecule = mol.metadata, cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata, tile = tile.metadata))
}

####
# Methods ####
####

vrSpatialPointsvrMetadata <- function(object, assay = NULL) {
  
  # get spatial points
  points <- unlist(lapply(methods::slotNames(object), function(x) {
    if(x %in% c("cell", "spot", "ROI")){
      mdata <- slot(object, name = x)
      if(nrow(mdata) > 0){
        if(!is.null(rownames(mdata))){
          sp <- rownames(mdata)
        } else {
          sp <- as.vector(mdata$id)
        }
        if(!is.null(assay))
          sp <- sp[grepl(paste(paste0(assay, "$"), collapse = "|"), sp)]
        return(sp)  
      }
    } else {
      mdata <- slot(object, name = x)
      if(nrow(mdata) > 0){
        if(inherits(mdata, "data.table")){
          if(!is.null(assay))
            sp <- subset(mdata, subset = assay_id %in% assay)
          return(sp[["id"]])
        } else {
          sp <- as.vector(mdata$id)
          if(!is.null(assay))
            sp <- sp[grepl(paste(paste0(assay, "$"), collapse = "|"), sp)]
          return(sp)
        }
      }
    }
  }))
  
  # return points
  return(points)
}

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @rdname vrSpatialPoints
#' @order 3
#'
#' @importFrom methods slotNames
#'
#' @export
setMethod("vrSpatialPoints", "vrMetadata", vrSpatialPointsvrMetadata)

subsetvrMetadata <- function(x, subset, samples = NULL, assays = NULL, spatialpoints = NULL) {
  
  # start 
  object <- x
  
  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }
  
  # subset all metadata types
  if(!is.null(samples)){
    if(nrow(object@molecule) > 0){
      mol.metadata <- subset_metadata(object@molecule, samples = samples)
    } else {
      mol.metadata <- data.table::data.table()
    }
    cell.metadata <- subset_metadata(object@cell, samples = samples)
    spot.metadata <- subset_metadata(object@spot, samples = samples)
    roi.metadata <- subset_metadata(object@ROI, samples = samples)
    if(nrow(object@tile) > 0){
      tile.metadata <- subset_metadata(object@tile, samples = samples)
    } else {
      tile.metadata <- data.table::data.table()
    }
  } else if(!is.null(assays)){
    assay_names <- unique(lapply(slotToList(object), function(x) {
      if(inherits(x, "data.table")){
        return(unique(as.vector(x$assay_id)))
      } else {
        if(!is.null(rownames(x))){
          return(unique(stringr::str_extract(rownames(x), "Assay[0-9]+")))
        } else {
          return(unique(stringr::str_extract(as.vector(x$id), "Assay[0-9]+"))) 
        }
      }
    }))
    assay_names <- unique(do.call(c,assay_names))
    if(all(assays %in% assay_names)){
      if(nrow(object@molecule) > 0) {
        mol.metadata <- subset_metadata(object@molecule, assays = assays)
      } else {
        mol.metadata <- data.table::data.table()
      }
      cell.metadata <- subset_metadata(object@cell, assays = assays)
      spot.metadata <- subset_metadata(object@spot, assays = assays)
      roi.metadata <- subset_metadata(object@ROI, assays = assays)
      if(nrow(object@tile) > 0) {
        tile.metadata <- object@tile[assay_id %in% assays, ]
      } else {
        tile.metadata <- data.table::data.table()
      }
    } else {
      if(nrow(object@molecule) > 0) {
        mol.metadata <- subset_metadata(object@molecule, assaytypes = assays)
      } else {
        mol.metadata <- data.table::data.table()
      }
      cell.metadata <- subset_metadata(object@cell, assaytypes = assays)
      spot.metadata <- subset_metadata(object@spot, assaytypes = assays)
      roi.metadata <- subset_metadata(object@ROI, assaytypes = assays)
      if(nrow(object@tile) > 0) {
        tile.metadata <- subset_metadata(object@tile, assaytypes = assays)
      } else {
        tile.metadata <- data.table::data.table()
      }
    }
  } else if(!is.null(spatialpoints)){
    if(nrow(object@molecule) > 0){
      mol.metadata <- subset_metadata(object@molecule, spatialpoints = spatialpoints)
    } else {
      mol.metadata <- data.table::data.table()
    }
    cell.metadata <- subset_metadata(object@cell, spatialpoints = spatialpoints)
    spot.metadata <- subset_metadata(object@spot, spatialpoints = spatialpoints)
    roi.metadata <- subset_metadata(object@ROI, spatialpoints = spatialpoints)
    if(nrow(object@tile) > 0){
      tile.metadata <- subset_metadata(object@tile, spatialpoints = spatialpoints)
    } else {
      tile.metadata <- data.table::data.table()
    }
  } else {
    stop("No assay, sample or spatial points were provided!")
  }
  
  # return new metadata
  methods::new("vrMetadata",
               molecule = mol.metadata,
               cell = cell.metadata,
               spot = spot.metadata,
               ROI = roi.metadata,
               tile = tile.metadata)
}

#' Subsetting vrMetadata objects
#'
#' Given a vrMetadata object, subset the object given one of the attributes
#'
#' @param x a vrMetadata object
#' @param subset the subset statement
#' @param samples the set of samples to subset the object
#' @param assays assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param spatialpoints the set of spatial points to subset the object
#'
#' @method subset vrMetadata
#' @order 3
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_extract
#' @importFrom data.table setkey
setMethod("subset", "vrMetadata", subsetvrMetadata)

#' subset_sampleMetadata
#'
#' Subseting sample metadata
#'
#' @param metadata sample metadata of a VoltRon object
#' @param samples the set of samples to subset the object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#'
#' @noRd
subset_sampleMetadata <- function(metadata, samples = NULL, assays = NULL) {

  # subseting on samples, layers and assays
  if(!is.null(samples)){
    if(all(samples %in% metadata$Sample)){
      metadata <- metadata[metadata$Sample %in% samples,]
    } else {
      stop("Some samples with the names '", paste(samples, collapse = ", "), "' are not found in the object")
    }
  } else if(!is.null(assays)) {
    if(all(assays %in% rownames(metadata))){
      metadata <- metadata[assays,]
    } else if(all(assays %in% metadata$Assay)){
      metadata <- metadata[metadata$Assay %in% assays,]
    } else {
      stop("Some assay with the names or types '", paste(assays, collapse = ", "), "' are not found in the object")
    }
  }
  metadata
}

mergevrMetadata <- function(x, y) {

  # start 
  object <- x
  object_list <- y
  
  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)

  # check if all are VoltRon
  if(!all(lapply(object_list, class) == "vrMetadata"))
    stop("All arguements have to be of vrMetadata class")

  # choose objects
  obj1 <- object_list[[1]]
  obj2 <- object_list[[2]]

  # initial combination
  if(length(object_list) > 2){
    combined.metadata <- mergevrMetadata(obj1, obj2)
    for(i in 3:(length(object_list))){
      combined.metadata <- mergevrMetadata(combined.metadata, object_list[[i]])
    }
  } else {
    updateobjects <- updateMetadataAssay(obj1, obj2)
    obj1 <- updateobjects$object1
    obj2 <- updateobjects$object2
    mol.metadata <- rbind_metadata(methods::slot(obj1, "molecule"), methods::slot(obj2, "molecule"))
    cell.metadata <- rbind_metadata(methods::slot(obj1, "cell"), methods::slot(obj2, "cell"))
    spot.metadata <- rbind_metadata(methods::slot(obj1, "spot"), methods::slot(obj2, "spot"))
    roi.metadata <- rbind_metadata(methods::slot(obj1, "ROI"), methods::slot(obj2, "ROI"))
    tile.metadata <- rbind_metadata(methods::slot(obj1, "tile"), methods::slot(obj2, "tile"))
    combined.metadata <- methods::new("vrMetadata", 
                                      molecule = mol.metadata, 
                                      cell = cell.metadata, 
                                      spot = spot.metadata, 
                                      ROI = roi.metadata, 
                                      tile = tile.metadata)
  }

  # return combined object
  return(combined.metadata)
}

#' Merging vrMetadata objects
#'
#' Given a vrMetadata object, and a list of vrMetadata objects, merge all.
#'
#' @param x a vrMetadata object
#' @param y a single or a list of vrMetadata objects
#'
#' @method merge vrMetadata
#'
#' @importFrom dplyr bind_rows
#' @importFrom methods slot
#' @export
setMethod("merge", "vrMetadata", mergevrMetadata)

#' rbind_metadata
#'
#' @param metadata1 metadata1
#' @param metadata2 metadata2
#'
#' @method merge vrMetadata
#'
#' @importFrom dplyr bind_rows
#' @noRd
rbind_metadata <- function(metadata1, metadata2){
  flag1 <- FALSE
  flag2 <- FALSE
  if(!inherits(metadata1, "DataFrame")){
    flag1 <- TRUE
  }
  if(!inherits(metadata2, "DataFrame")){
    flag2 <- TRUE
  }
  if(flag1 && flag2){
    return(dplyr::bind_rows(metadata1,metadata2))
  } else {
    metadata1 <- .make_DF(metadata1)
    metadata2 <- .make_DF(metadata2)
    return(rbind(metadata1, metadata2))
  }
}

#' cbind_metadata
#'
#' @param ... set of metadata objects
#'
#' @method merge vrMetadata
#'
#' @importFrom dplyr bind_cols
#' @importFrom S4Vectors DataFrame
#' @noRd
cbind_metadata <- function(...){
  # initial combination
  metadata_list <- list(...)
  if(length(metadata_list) > 2){
    combined.metadata <- cbind_metadata(metadata_list[[1]], 
                                        metadata_list[[2]])
    for(i in 3:(length(metadata_list))){
      combined.metadata <- cbind_metadata(combined.metadata, 
                                          metadata_list[[i]])
    }
    return(combined.metadata)
  } else {
    flag1 <- FALSE
    flag2 <- FALSE
    if(!inherits(metadata_list[[1]], "DataFrame")){
      flag1 <- TRUE
    }
    if(!inherits(metadata_list[[2]], "DataFrame")){
      flag2 <- TRUE
    }
    if(flag1 && flag2){
      return(dplyr::bind_cols(metadata_list[[1]],metadata_list[[2]]))
    } else {
      if(flag1)
        metadata_list[[1]] <- S4Vectors::DataFrame(metadata_list[[1]])
      if(flag2)
        metadata_list[[2]] <- S4Vectors::DataFrame(metadata_list[[2]])
      return(cbind(metadata_list[[1]], metadata_list[[2]]))
    } 
  }
}

#' subset_metadata
#'
#' @param metadata metadata
#' @param samples the set of samples to subset the object
#' @param assays assay name (exp: Assay1), see \code{SampleMetadata(object)}
#' @param assaytypes assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param spatialpoints the set of spatial points to subset the object
#'
#' @noRd
subset_metadata <- function(metadata, assays = NULL, assaytypes = NULL, samples = NULL, spatialpoints = NULL){
  
  if(inherits(metadata, "data.table")){
    if(nrow(metadata) > 0){
      if(!is.null(assays)){
        metadata <- subset(metadata, subset = assay_id %in% assays)
      } else if(!is.null(assaytypes)){
        metadata <- subset(metadata, subset = Assay %in% assaytypes)
      } else if(!is.null(samples)){
        metadata <- subset(metadata, subset = Sample %in% samples)
      } else if(!is.null(spatialpoints)){
        metadata <- subset(metadata, subset = id %in% spatialpoints)
      } else {
        stop("No assay, sample or spatial points were provided!")
      }  
    } else {
      metadata <- data.table::data.table()
    }
  } else if(inherits(metadata, "DataFrame")){
    if(!is.null(assays)){
      if("assay_id" %in% colnames(metadata)){
        cur_column <- as.vector(metadata$assay_id)
        metadata <- metadata[cur_column %in% assays,]
      } else {
        cur_column <- as.vector(metadata$id)
        metadata <- metadata[stringr::str_extract(cur_column, "Assay[0-9]+") %in% assays, ]
      }
    } else if(!is.null(assaytypes)){
      cur_column <- as.vector(metadata$Assay)
      metadata <- metadata[cur_column %in% assaytypes,]
    } else if(!is.null(samples)){
      cur_column <- as.vector(metadata$Sample)
      metadata <- metadata[cur_column %in% samples,]
    } else if(!is.null(spatialpoints)){
      cur_column <- as.vector(metadata$id)
      metadata <- metadata[cur_column %in% spatialpoints,]
    } else {
      stop("No assay, sample or spatial points were provided!")
    }  
  } else {
    if(nrow(metadata) > 0){
      if(!is.null(assays)){
        if(!is.null(rownames(metadata))){
          metadata <- metadata[stringr::str_extract(rownames(metadata), "Assay[0-9]+") %in% assays, ]
        } else {
          if("assay_id" %in% colnames(metadata)){
            metadata <- subset(metadata, subset = assay_id %in% assays)
          } else {
            metadata <- metadata[stringr::str_extract(metadata$id, "Assay[0-9]+") %in% assays, ]
          }
        }
      } else if(!is.null(assaytypes)){
        metadata <- subset(metadata, subset = Assay %in% assaytypes)
      } else if(!is.null(samples)){
        metadata <- subset(metadata, subset = Sample %in% samples)
      } else if(!is.null(spatialpoints)){
        if(!is.null(rownames(metadata))){
          metadata <- metadata[rownames(metadata) %in% spatialpoints,]
        } else {
          metadata <- metadata[metadata$id %in% spatialpoints,]
        }
      } else {
        stop("No assay, sample or spatial points were provided!")
      }  
    }
  }
  metadata
}

#' #' get_from_metadata
#' #'
#' #' @param metadata metadata
#' #' @param column the column to return
#' #'
#' #' @noRd
#' get_from_metadata <- function(metadata, column){
#'   if(inherits(metadata, "data.table")){
#'     column_vector <- metadata[,get(names(metadata)[which(colnames(metadata) == column)])]
#'   } else {
#'     if("id" %in% colnames(metadata)){
#'       column_vector <- as.vector(metadata[match(rownames(datax), as.vector(metadata$id)),column])
#'     } else{
#'       column_vector <- metadata[rownames(datax),column]
#'     }
#'   }
#'   column_vector
#' }
  
#' merge.sampleMetadata
#'
#' Merging sample.metadata from two VoltRon objects
#'
#' @param metadata_list a list of sample metadata of a VoltRon object
#'
#' @noRd
#'
merge_sampleMetadata <- function(metadata_list) {

  sample_names <- NULL
  sample.metadata <- do.call(rbind, metadata_list)
  rownames(sample.metadata) <- paste0("Assay", seq_len(nrow(sample.metadata)))

  # change sample names if provided
  if(!is.null(sample_names)){

    # check the number sample names
    if(!length(sample_names) %in% c(1,nrow(sample.metadata))){
      stop("Please provide only one sample name or of length of object list!")
    } else {
      sample.metadata$Sample <- sample_names
      section_ids <- rep(NA,nrow(sample.metadata))
      uniq_names <- unique(sample.metadata$Sample)
      for(i in seq_len(length(uniq_names))){
        cur_ind <- which(sample.metadata$Sample == uniq_names[i])
        section_ids[cur_ind] <- seq_len(length(cur_ind))
      }
      sample.metadata$Layer <- paste0("Section", section_ids)
    }
  }
  sample.metadata
}

#' @noRd
.fix_metadata <- function(metadata){
  if(!"id" %in% colnames(metadata))
    if(!is.null(rownames(metadata))){
      metadata$id <- id <- rownames(metadata)
    } else {
      stop("Ill-defined metadata with no rownames, please check the object!")
    }
  # if(!"assay_id" %in% colnames(metadata)){
  #   metadata[["assay_id"]]  <- stringr::str_extract(id, "Assay[0-9]+")
  # }
  suppressMessages(rownames(metadata) <- NULL)
  metadata
}

#' @importFrom S4Vectors DataFrame
#' @noRd
.make_DF <- function(metadata){
  if(!requireNamespace("DelayedArray"))
    stop("Please install DelayedArray package!:", "
         BiocManager::install('DelayedArray')")
  metadata <- as.list(metadata)
  metadata <- lapply(metadata, \(x){
    if(!inherits(x, "DelayedArray")) {
      return(DelayedArray::DelayedArray(array(x)))
    }
    return(x)
  })
  S4Vectors::DataFrame(metadata)
}

### Assay Methods ####

addAssayvrMetadata <- function(object, metadata = NULL, assay, assay_name, sample = "Sample1", layer = "Section1"){

  # get metadata and other info
  assay.type <- vrAssayTypes(assay)
  object_metadata <- methods::slot(object, name = assay.type)
  data <- vrData(assay, norm = FALSE)
  
  # fix object metadata
  object_metadata <- .fix_metadata(object_metadata)

  # add new assay
  assay_ids <- vrAssayNames(object)
  assay_ids <- as.numeric(gsub("Assay", "", assay_ids))
  assay_id <- paste0("Assay", max(assay_ids)+1)
  
  # fix assay metadata
  if(!is.null(metadata))
    metadata <- .fix_metadata(metadata)

  # metadata
  if(inherits(metadata, "data.table")){

    if(!is.null(metadata)){
      count <- if(nrow(data) > 0) Matrix::colSums(data) else NULL
      assay_metadata <- 
        data.table::data.table(
          metadata[, "id", with=FALSE], 
          assay_id = assay_id, 
          Count = count,
          Assay = assay_name, 
          Layer = layer, 
          Sample = sample,
          metadata[, colnames(metadata)[
            !colnames(metadata) %in% 
              c("id", "assay_id", "Count", "Assay", "Layer", "Sample")], 
            with=FALSE
            ])
    }
    
  } else {

    # get original names
    entityID_nopostfix <- stringr::str_replace(vrSpatialPoints(assay), 
                                               pattern = "_Assay[0-9]+", "")
    entityID <- stringr::str_replace(entityID_nopostfix,
                                     pattern = "$", 
                                     paste0("_", assay_id))

    # get metadata id
    if("id" %in% colnames(metadata)){
      metadata_id <- stringr::str_replace(as.vector(metadata$id), 
                                          pattern = "_Assay[0-9]+", "") 
    } else {
      metadata_id <- entityID_nopostfix
    }
    
    # initiate metadata
    if(nrow(data) > 0){
      assay_metadata <- data.frame(id = entityID, 
                                   Count = Matrix::colSums(data), 
                                   row.names = NULL)
    } else {
      assay_metadata <- data.frame(id = entityID, 
                                   row.names = NULL)
    }
    
    # add metadata
    if(!is.null(metadata)){
      if(length(setdiff(metadata_id, entityID_nopostfix)) > 0){
        if(nrow(metadata) != length(entityID_nopostfix)){
          stop("Some spatial points in the metadata 
               does not match with the assay!")
        }
      }
      assay_metadata <- cbind_metadata(
        assay_metadata,
        data.frame(Assay = rep(assay_name, length(entityID)),
                   Layer = rep(layer, length(entityID)),
                   Sample = rep(sample, length(entityID)),
                   assay_id = assay_id),
        metadata[,!colnames(metadata) %in% 
                   c("id", "Count", "assay_id", "Assay", "Layer", "Sample"), 
                 drop = FALSE])
    } else {
      assay_metadata <- cbind_metadata(
        assay_metadata,
        data.frame(Assay = rep(assay_name, length(entityID)),
                   Layer = rep(layer, length(entityID)),
                   Sample = rep(sample, length(entityID)),
                   assay_id = assay_id))
    }
  }

  # add to the main metadata
  object_metadata <- rbind_metadata(object_metadata, assay_metadata)
  methods::slot(object, name = assay.type) <- object_metadata

  # return
  return(object)
}

#' @rdname addAssay
#' @method addAssay vrMetadata
#'
#' @importFrom dplyr bind_rows bind_cols
#' @importFrom methods slot slot<-
#' @importFrom stringr str_replace
#' @importFrom data.table data.table
#' @importFrom Matrix colSums
#'
#' @export
setMethod("addAssay", "vrMetadata", addAssayvrMetadata)

vrAssayNamesvrMetadata <- function(object){
  
  # get assay names from metadata
  assay_names <- NULL
  for(sl in methods::slotNames(object)){
    cur_metadata <- slot(object, name = sl)
    if(sl %in% c("molecule", "tile")){
      cur_names <- cur_metadata$assay_id
    } else {
      if("assay_id" %in% colnames(cur_metadata)){
        cur_names <- as.vector(cur_metadata$assay_id)
      } else if(!is.null(rownames(cur_metadata))){
        cur_names <- stringr::str_extract(rownames(cur_metadata), 
                                          "Assay[0-9]+")
      } else{
        cur_names <- stringr::str_extract(as.vector(cur_metadata$id), 
                                          "Assay[0-9]+")
      }
    }
    assay_names <- c(assay_names, unique(cur_names))
  }
  assay_names
}

#' @rdname vrAssayNames
#' @order 3
#' @importFrom methods slotNames
#' @export
setMethod("vrAssayNames", "vrMetadata", vrAssayNamesvrMetadata)

#' updateMetadataAssay
#'
#' Updating assay names for merge
#'
#' @param object1 vrMetadata object
#' @param object2 vrMetadata object
#'
#' @importFrom stringr str_extract
#' @importFrom methods new
#'
#' @noRd
updateMetadataAssay <- function(object1, object2){

  # get assay types
  object_list <- slotToList(object1)
  assaytype <- unlist(lapply(object_list, function(obj) {
    if(inherits(obj, "data.table")){
      unique(obj$assay_id)
    } else if(inherits(obj, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame"))){
      if("assay_id" %in% colnames(obj)){
        unique(as.vector(obj$assay_id))
      } else {
        unique(stringr::str_extract(as.vector(obj$id), "Assay[0-9]+$"))
      }
    } else {
      unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
    }
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", seq_len(length(assaytype)))
  object1 <- lapply(object_list, function(obj) {
    if(nrow(obj) > 0){
      
      if(inherits(obj, "data.table")){
        
        # change assay id
        temp <- obj$assay_id
        for(i in seq_len(length(assaytype)))
          temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
        obj$assay_id <- temp
        return(obj)
        
      } else if(inherits(obj, c("HDF5DataFrame", 
                                "ZarrDataFrame", 
                                "DataFrame"))){
        
        # change assay id
        if("assay_id" %in% colnames(obj)){
          temp <- as.vector(obj$assay_id)
          for(i in seq_len(length(assaytype)))
            temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
          obj$assay_id <- temp
        }
        
        # change id
        temp <- as.vector(obj$id)
        for(i in seq_len(length(assaytype))){
          temp[grepl(paste0(assaytype[i],"$"), obj$id)] <- 
            gsub(paste0(assaytype[i],"$"), replacement[i],
                 obj$id[grepl(paste0(assaytype[i],"$"),  obj$id)])
        }
        obj$id <- temp
        
        return(obj)
      } else {
        
        # change rownames
        temp <- rownames(obj)
        for(i in seq_len(length(assaytype)))
          temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- 
            gsub(paste0(assaytype[i],"$"), replacement[i],
                 rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
        rownames(obj) <- temp
        
        # change assay id
        if("assay_id" %in% colnames(obj)){
          temp <- obj$assay_id
          for(i in seq_len(length(assaytype)))
            temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
          obj$assay_id <- temp
        }
        return(obj)
      }
    } else {
      return(obj)
    }
  })
  object1 <- methods::new("vrMetadata", 
                          molecule = object1$molecule, 
                          cell = object1$cell, 
                          spot = object1$spot, 
                          ROI = object1$ROI, 
                          tile = object1$tile)

  # get assay types
  object_list <- slotToList(object2)
  assaytype <- unlist(lapply(object_list, function(obj) {
    if(inherits(obj, "data.table")){
      unique(obj$assay_id)
    } else if(inherits(obj, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame"))){
      if("assay_id" %in% colnames(obj)){
        unique(as.vector(obj$assay_id))
      } else {
        unique(stringr::str_extract(as.vector(obj$id), "Assay[0-9]+$"))
      }
    } else {
      unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
    }
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  object2 <- lapply(object_list, function(obj) {
    if(nrow(obj) > 0){
      if(inherits(obj, "data.table")){
        
        # change assay id
        temp <- obj$assay_id
        for(i in seq_len(length(assaytype)))
          temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
        obj$assay_id <- temp
        
        return(obj)
      } else if(inherits(obj, c("HDF5DataFrame", "ZarrDataFrame", "DataFrame"))){
        
        # change assay id
        if("assay_id" %in% colnames(obj)){
          temp <- as.vector(obj$assay_id)
          for(i in seq_len(length(assaytype)))
            temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
          obj$assay_id <- temp
        }
        
        # change id
        temp <- as.vector(obj$id)
        for(i in seq_len(length(assaytype))){
          temp[grepl(paste0(assaytype[i],"$"), obj$id)] <- 
            gsub(paste0(assaytype[i],"$"), replacement[i], 
                 obj$id[grepl(paste0(assaytype[i],"$"), obj$id)])
        }
        obj$id <- temp
        
        return(obj)
      } else {
        
        # change row names
        temp <- rownames(obj)
        for(i in seq_len(length(assaytype)))
          temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- 
            gsub(paste0(assaytype[i],"$"), replacement[i],
                 rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
        rownames(obj) <- temp
        
        # change id
        temp <- obj$id
        for(i in seq_len(length(assaytype))){
          temp[grepl(paste0(assaytype[i],"$"), obj$id)] <- 
            gsub(paste0(assaytype[i],"$"), replacement[i], 
                 obj$id[grepl(paste0(assaytype[i],"$"), obj$id)])
        }
        obj$id <- temp
        
        # change assay id
        if("assay_id" %in% colnames(obj)){
          temp <- obj$assay_id
          for(i in seq_len(length(assaytype)))
            temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
          obj$assay_id <- temp
        }
        obj
      }
    } else {
      return(obj)
    }
  })
  object2 <- methods::new("vrMetadata", 
                          molecule = object2$molecule, 
                          cell = object2$cell, 
                          spot = object2$spot, 
                          ROI = object2$ROI, 
                          tile = object2$tile)

  # return
  return(list(object1 = object1, object2 = object2))
}

changeSampleNamesvrMetadata <- function(object, sample_metadata_table){

  # get old and new samples
  old.samples <- sample_metadata_table$Sample
  new.samples <- sample_metadata_table$NewSample

  # check all types in the vrMetadata object
  new_object <- object
  all_types <- methods::slotNames(object)
  for(type in all_types){
    metadata <- methods::slot(object, name = type)
    new_metadata <-  methods::slot(new_object, name = type)
    if(nrow(new_metadata) > 0){

      # change samples
      for(i in seq_len(length(old.samples)))
        new_metadata$Sample[new_metadata$Sample==old.samples[i]] <- new.samples[i]

      # change layers
      for(i in seq_len(nrow(sample_metadata_table))){
        new_metadata$Layer[grepl(paste0(sample_metadata_table$AssayID[i], "$"), rownames(new_metadata))] <- sample_metadata_table[sample_metadata_table$AssayID[i], "NewLayer"]
      }

      # rewrite metadata type
      methods::slot(new_object, name = type) <- new_metadata
    }
  }

  # return
  return(new_object)
}

#' changeSampleNames.vrMetadata
#'
#' Change the sample names of the vrMetadata object and reorient layers if needed
#' This functions requires the new and old sample and layer names passed from \code{changeSampleNames.VoltRon}
#'
#' @param sample_metadata_table the sample metadata with old and new layers and samples passed from \code{changeSampleNames.VoltRon}
#' 
#' @rdname changeSampleNames
#' @method changeSampleNames vrMetadata
#'
#' @importFrom methods slot slot<- slotNames
#'
#' @noRd
setMethod("changeSampleNames", "vrMetadata", changeSampleNamesvrMetadata)

### Sample Methods ####

vrSampleNamesvrMetadata <- function(object){
  
  # get assay names from metadata
  sample_names <- NULL
  for(sl in methods::slotNames(object)){
    cur_metadata <- slot(object, name = sl)
    sample_names <- c(sample_names, unique(cur_metadata$Sample))
  }
  
  # return
  sample_names
}

#' @rdname vrSampleNames
#' @method vrSampleNames vrMetadata
#'
#' @importFrom methods slotNames
#' @export
setMethod("vrSampleNames", "vrMetadata", vrSampleNamesvrMetadata)

####
# Functions ####
####

#' setVRMetadata
#'
#' @param molecule molecule data frame
#' @param cell cell data frame
#' @param spot spot data frame
#' @param ROI ROI data frame
#' @param tile tile data frame
#'
#' @importFrom methods new
#'
#' @noRd
setVRMetadata <- function(metadata, data, entityID, main.assay, assay.type, sample_name, layer_name, version){
  
  if(is.null(metadata)){

    # set metadata
    vr_metadata <- list(molecule = data.table::data.table(),
                        cell = data.frame(),
                        spot = data.frame(),
                        ROI = data.frame(),
                        tile = data.table::data.table())

    # create entity IDs using Assay index, make it colnames
    entityID <- stringr::str_replace(entityID, pattern = "$", paste0("_Assay1"))

    # create metadata
    # slot(vr_metadata, name = assay.type) <-
    if(version == "v1"){
      vr_metadata[[assay.type]] <- 
        data.frame(Count = Matrix::colSums(data),
                   assay_id = "Assay1",
                   Assay = main.assay,
                   Layer = layer_name,
                   Sample = sample_name,
                   row.names = entityID) 
    } else if (version == "v2"){
      vr_metadata[[assay.type]] <- 
        data.frame(id = entityID,
                   Count = Matrix::colSums(data),
                   assay_id = "Assay1",
                   Assay = main.assay,
                   Layer = layer_name,
                   Sample = sample_name,
                   row.names = entityID) 
    }

  } else {
    if(any(is(metadata) %in% c("data.table", "data.frame", "matrix"))){
      vr_metadata <- list(molecule = data.table::data.table(),
                          cell = data.frame(),
                          spot = data.frame(),
                          ROI = data.frame(),
                          tile = data.table::data.table())

      # if metadata is a data.table
      if(inherits(metadata, "data.table")){

        # if there are no id column, insert entityID
        if(!"id" %in% colnames(metadata)){
          metadata$id <- entityID  
        }
        
        # check ID names
        if(length(setdiff(metadata$id, entityID)) > 0){
          stop("Entity IDs are not matching")
        } else {

          # entity IDs
          metadata <- subset(metadata, subset = entityID %in% id)

          # create entity IDs using Assay index, make it colnames
          set.seed(nrow(metadata$id))
          entityID <- paste0(metadata$id, "_", ids::random_id(bytes = 3, use_openssl = FALSE))

          if(nrow(data) > 0){
            suppressWarnings({
              vr_metadata[[assay.type]] <-
                data.table::data.table(id = entityID,
                                       assay_id = "Assay1",
                                       Count = Matrix::colSums(data),
                                       Assay = main.assay,
                                       Layer = layer_name,
                                       Sample = sample_name,
                                       metadata[,-"id"])
            })
          } else{
            suppressWarnings({
              vr_metadata[[assay.type]] <-
                  data.table::data.table(id = entityID,
                                         assay_id = "Assay1",
                                         Assay = main.assay,
                                         Layer = layer_name,
                                         Sample = sample_name, 
                                         metadata[,-"id"])
            })
          }
        }

      # if metadata is a regular data.frame
      } else if(inherits(metadata, "data.frame")){

        # check row names
        if(length(setdiff(rownames(metadata), entityID)) > 0){
          stop("Entity IDs are not matching")
        } else {

          # entity IDs
          if(version == "v1") {
            metadata <- metadata[entityID,] 
          } else if(version == "v2") {
            
            # if there are no id column, insert entityID
            if(!"id" %in% colnames(metadata)){
              metadata$id <- entityID  
            }
            metadata <- metadata[match(entityID, metadata$id),]
          }

          # create entity IDs using Assay index, make it colnames
          entityID <- stringr::str_replace(entityID, pattern = "$", paste0("_Assay1"))

          # create metadata for version 1
          if(version == "v1"){
            if(nrow(data) > 0){
              vr_metadata[[assay.type]] <-
                data.frame(Count = Matrix::colSums(data),
                           assay_id = "Assay1",
                           Assay = main.assay,
                           Layer = layer_name,
                           Sample = sample_name,
                           metadata, 
                           row.names = entityID)
            } else{
              vr_metadata[[assay.type]] <-
                data.frame(assay_id = "Assay1",
                           Assay = main.assay,
                           Layer = layer_name,
                           Sample = sample_name,
                           metadata,
                           row.names = entityID)
            } 
          
          # create metadata for version 2
          } else if(version == "v2"){
            if(nrow(data) > 0){
              vr_metadata[[assay.type]] <-
                data.frame(id = entityID,
                           Count = Matrix::colSums(data),
                           assay_id = "Assay1",
                           Assay = main.assay,
                           Layer = layer_name,
                           Sample = sample_name,
                           metadata, 
                           row.names = entityID)
            } else{
              vr_metadata[[assay.type]] <-
                data.frame(id = entityID,
                           Assay = main.assay,
                           assay_id = "Assay1",
                           Layer = layer_name,
                           Sample = sample_name,
                           metadata,
                           row.names = entityID)
            }
          }
          
        }
      }
    }
  }
  
  return(
    list(
      entityID = entityID,
      vr_metadata =  methods::new("vrMetadata", 
                                 molecule = vr_metadata$molecule, 
                                 cell = vr_metadata$cell, 
                                 spot = vr_metadata$spot, 
                                 ROI = vr_metadata$ROI, 
                                 tile = vr_metadata$tile)
    )
  )
}

#' setVRSampleMetadata
#'
#' @param samples a list of vrSample object
#'
#' @noRd
setVRSampleMetadata <- function(samples){

  # imput missing sample names
  # sample_name_ind <- sapply(names(samples), is.null)
  sample_name_ind <- vapply(names(samples), is.null, logical(1))
  if(length(sample_name_ind) > 0){
    names_samples <- names(samples)
    if(any(sample_name_ind)){
      null_samples_ind <- which(sample_name_ind)
      names_samples[null_samples_ind] <- paste0("Sample", null_samples_ind)
    }
  } else {
    names_samples <- paste0("Sample", seq_len(length(samples)))
  }

  # get sample metadata
  sample_list <- names(samples)
  sample.metadata <- NULL
  for(i in seq_len(length(sample_list))){
    layer_list <- samples[[sample_list[i]]]@layer
    layer_data <- NULL
    for(j in seq_len(length(layer_list))){
      assay_list <- layer_list[[j]]@assay
      layer_data <- rbind(layer_data, cbind(names(assay_list), names(layer_list)[j]))
    }
    sample.metadata <- rbind(sample.metadata, cbind(layer_data, sample_list[i]))
  }
  sample.metadata <- data.frame(sample.metadata, row.names = paste0("Assay", seq_len(nrow(sample.metadata))))
  colnames(sample.metadata) <- c("Assay", "Layer", "Sample")

  sample.metadata
}
