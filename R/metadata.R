####
# Objects and Classes ####
####

## vrMetadata ####

#' The vrMetadata (VoltRon Metadata) Class
#'
#' @slot samples A list of layers for the this vrSample object
#'
#' @name vrMetadata-class
#' @rdname vrMetadata-class
#' @exportClass vrMetadata
#'
vrMetadata <- setClass(
  Class = 'vrMetadata',
  slots = c(
    cell = 'data.frame',
    spot = 'data.frame',
    ROI = 'data.frame'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrMetadata',
  definition = function(object) {
    cat("VoltRon Metadata Object \n")
    cat("This object includes: \n")
    if(nrow(object@cell) > 0)
      cat("  ", nrow(object@cell), "cells \n")
    if(nrow(object@spot) > 0)
      cat("  ", nrow(object@spot), "spots \n")
    if(nrow(object@ROI) > 0)
      cat("  ", nrow(object@ROI), "ROIs \n")
    return(invisible(x = NULL))
  }
)

### $ methods ####

#' @export
#' @method $ vrMetadata
#'
"$.vrMetadata" <- function(x, i, ...) {
  return(NULL)
}

#' @export
#' @method $<- vrMetadata
#'
"$<-.vrMetadata" <- function(x, i, ..., value) {

  # cell metadata
  cell.metadata <- slot(x, "cell")
  if(nrow(cell.metadata) > 0)
    cell.metadata[[i]] <- value

  # spot metadata
  spot.metadata <- slot(x, "spot")
  if(nrow(spot.metadata) > 0)
    spot.metadata[[i]] <- value

  # ROI metadata
  roi.metadata <- slot(x, "ROI")
  if(nrow(roi.metadata) > 0)
    roi.metadata[[i]] <- value

  return(new("vrMetadata", cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata))
}

#' @export
#' @method $<- vrMetadata
#'
"[[<-.vrMetadata" <- function(x, i, ..., value) {

  # cell metadata
  cell.metadata <- slot(x, "cell")
  if(nrow(cell.metadata) > 0)
    cell.metadata[[i]] <- value

  # spot metadata
  spot.metadata <- slot(x, "spot")
  if(nrow(spot.metadata) > 0)
    spot.metadata[[i]] <- value

  # ROI metadata
  roi.metadata <- slot(x, "ROI")
  if(nrow(roi.metadata) > 0)
    roi.metadata[[i]] <- value

  return(new("vrMetadata", cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata))
}

####
# Methods ####
####

#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrMetadata
#'
#' @export
#'
vrSpatialPoints.vrMetadata <- function(object, ...) {

  # get the combination of cells, spots and ROIs
  points <- c(rownames(object@cell),
                rownames(object@spot),
                rownames(object@ROI))

  return(points)
}

#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#'
#' @method subset vrMetadata
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_extract
#'
subset.vrMetadata <- function(metadata, subset, samples = NULL, assays = NULL, spatialpoints = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subset all metadata types
  if(!is.null(samples)){
    cell.metadata <- metadata@cell[metadata@cell$Sample %in% samples, ]
    spot.metadata <- metadata@spot[metadata@spot$Sample %in% samples, ]
    roi.metadata <- metadata@ROI[metadata@ROI$Sample %in% samples, ]
  } else if(!is.null(assays)){
    assay_names <- unique(lapply(slotToList(metadata), function(x) {
      unique(stringr::str_extract(rownames(x), "Assay[0-9]+"))
    }))
    assay_names <- unique(do.call(c,assay_names))
    if(all(assays %in% assay_names)){
      cell.metadata <- metadata@cell[stringr::str_extract(rownames(metadata@cell), "Assay[0-9]+") %in% assays, ]
      spot.metadata <- metadata@spot[stringr::str_extract(rownames(metadata@spot), "Assay[0-9]+") %in% assays, ]
      roi.metadata <- metadata@ROI[stringr::str_extract(rownames(metadata@ROI), "Assay[0-9]+") %in% assays, ]
    } else {
      cell.metadata <- metadata@cell[metadata@cell$Assay %in% assays, ]
      spot.metadata <- metadata@spot[metadata@spot$Assay %in% assays, ]
      roi.metadata <- metadata@ROI[metadata@ROI$Assay %in% assays, ]
    }
  } else if(!is.null(spatialpoints)){
    if(all(spatialpoints %in% vrSpatialPoints(metadata))){
      cell.metadata <- metadata@cell[rownames(metadata@cell) %in% spatialpoints, ]
      spot.metadata <- metadata@spot[rownames(metadata@spot) %in% spatialpoints, ]
      roi.metadata <- metadata@ROI[rownames(metadata@ROI) %in% spatialpoints, ]
    } else {
      stop("Some spatial points are not found in the metadata and the object")
    }
  } else {
    stop(paste0("No assay or sample name was provided!"))
  }

  # return new metadata
  setVRMetadata(cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata)
}


#' subset.sampleMetadata
#'
#' @param metadata sample metadata of a VoltRon object
#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#'
#' @export
#'
subset.sampleMetadata <- function(metadata, samples = NULL, assays = NULL) {

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

#' @param object_list a list of vrMetadata objects
#'
#' @method merge vrMetadata
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
merge.vrMetadata <- function(object, object_list) {

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
    combined.metadata <- merge(obj1, obj2)
    for(i in 3:(length(object_list))){
      combined.metadata <- merge(combined.metadata, object_list[[i]])
    }
  } else {
    updateobjects <- updateMetadataAssay(obj1, obj2)
    obj1 <- updateobjects$object1
    obj2 <- updateobjects$object2
    cell.metadata <- bind_rows(slot(obj1, "cell"), slot(obj2, "cell"))
    spot.metadata <- bind_rows(slot(obj1, "spot"), slot(obj2, "spot"))
    roi.metadata <- bind_rows(slot(obj1, "ROI"), slot(obj2, "ROI"))
    combined.metadata <- setVRMetadata(cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata)
  }

  # return combined object
  return(combined.metadata)
}

#' merge.sampleMetadata
#'
#' Merging sample.metadata from two VoltRon objects
#'
#' @param metadata_list a list of sample metadata of a VoltRon object
#' @param sample_names new sample names
#'
#' @export
#'
merge_sampleMetadata <- function(metadata_list) {

  sample_names <- NULL
  sample.metadata <- do.call(rbind, metadata_list)
  rownames(sample.metadata) <- paste0("Assay", 1:nrow(sample.metadata))

  # change sample names if provided
  if(!is.null(sample_names)){

    # check the number sample names
    if(!length(sample_names) %in% c(1,nrow(sample.metadata))){
      stop("Please provide only one sample name or of length of object list!")
    } else {
      sample.metadata$Sample <- sample_names
      section_ids <- rep(NA,nrow(sample.metadata))
      uniq_names <- unique(sample.metadata$Sample)
      for(i in 1:length(uniq_names)){
        cur_ind <- which(sample.metadata$Sample == uniq_names[i])
        section_ids[cur_ind] <- 1:length(cur_ind)
      }
      sample.metadata$Layer <- paste0("Section", section_ids)
    }
  }
  sample.metadata
}

### Assay Methods ####

#' @param assay assay
#' @param assay_name assay name
#' @param sample sample name
#' @param layer layer name
#'
#' @rdname addAssay
#' @method addAssay vrMetadata
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
addAssay.vrMetadata <- function(object, assay, assay_name, sample = "Sample1", layer = "Section1"){

  # assay info
  assay.type <- vrAssayTypes(assay)

  # get metadata and other info
  metadata <- slot(object, name = assay.type)
  data <- vrData(assay, norm = FALSE)
  spatialpoints <- vrSpatialPoints(object)

  # add new assay
  assay_ids <- stringr::str_extract(spatialpoints, "Assay[0-9]+")
  assay_ids <- as.numeric(gsub("Assay", "", assay_ids))
  assay_id <- paste0("Assay", max(assay_ids)+1)
  entityID <- gsub("Assay[0-9]+$", assay_id, vrSpatialPoints(assay))

  # metadata
  assay_metadata <- data.frame(Count = colSums(data),
                               Assay = rep(assay_name, length(entityID)),
                               Layer = rep(layer, length(entityID)),
                               Sample = rep(sample, length(entityID)),
                               row.names = entityID)
  metadata <- dplyr::bind_rows(metadata, assay_metadata)

  slot(object, name = assay.type) <- metadata

  # return
  return(object)
}

#' updateMetadataAssay
#'
#' Updating assay names for merge
#'
#' @param object1 vrMetadata object
#' @param object2 vrMetadata object
#'
#' @importFrom stringr str_extract
#' @importFrom stringi stri_replace_all_regex
#'
updateMetadataAssay <- function(object1, object2){

  # get assay types
  object_list <- slotToList(object1)
  assaytype <- unlist(lapply(object_list, function(obj) {
    unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", 1:length(assaytype))
  object1 <- lapply(object_list, function(obj) {
    temp <- rownames(obj)
    for(i in 1:length(assaytype))
      temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                   rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
    rownames(obj) <- temp
    obj
  })
  object1 <- new("vrMetadata", cell = object1$cell, spot = object1$spot, ROI = object1$ROI)

  # get assay types
  object_list <- slotToList(object2)
  assaytype <- unlist(lapply(object_list, function(obj) {
    unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  object2 <- lapply(object_list, function(obj) {
    temp <- rownames(obj)
    for(i in 1:length(assaytype))
      temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                   rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
    rownames(obj) <- temp
    obj
  })
  object2 <- new("vrMetadata", cell = object2$cell, spot = object2$spot, ROI = object2$ROI)

  # return
  return(list(object1 = object1, object2 = object2))
}

#' changeSampleNames.vrMetadata
#'
#' Change the sample names of the vrMetadata object and reorient layers if needed
#' This functions requires the new and old sample and layer names passed from \code{changeSampleNames.VoltRon}
#'
#' @rdname changeSampleNames
#' @method changeSampleNames vrMetadata
#'
#' @param object A VoltRon object
#' @param sample_metadata_table the sample metadata with old and new layers and samples passed from \code{changeSampleNames.VoltRon}
#'
changeSampleNames.vrMetadata <- function(object, sample_metadata_table){

  # get old and new samples
  old.samples <- sample_metadata_table$Sample
  new.samples <- sample_metadata_table$NewSample

  # check all types in the vrMetadata object
  new_object <- object
  all_types <- slotNames(object)
  for(type in all_types){
    metadata <- slot(object, name = type)
    new_metadata <-  slot(new_object, name = type)
    if(nrow(new_metadata) > 0){

      # change samples
      for(i in 1:length(old.samples))
        new_metadata$Sample[new_metadata$Sample==old.samples[i]] <- new.samples[i]

      # change layers
      for(i in 1:nrow(sample_metadata_table)){
        new_metadata$Layer[grepl(paste0(sample_metadata_table$AssayID[i], "$"), rownames(new_metadata))] <- sample_metadata_table[sample_metadata_table$AssayID[i], "NewLayer"]
      }

      # rewrite metadata type
      slot(new_object, name = type) <- new_metadata
    }
  }

  # return
  return(new_object)
}

####
# Functions ####
####

#' setVRMetadata
#'
#' @param cell cell data frame
#' @param spot spot data frame
#' @param ROI ROI data frame
#'
setVRMetadata <- function(cell, spot, ROI){
  new("vrMetadata", cell = cell, spot = spot, ROI = ROI)
}

#' setVRSampleMetadata
#'
#' @param samples a list of vrSample object
#'
setVRSampleMetadata <- function(samples){

  # imput missing sample names
  sample_name_ind <- sapply(names(samples), is.null)
  if(length(sample_name_ind) > 0){
    names_samples <- names(samples)
    if(any(sample_name_ind)){
      null_samples_ind <- which(sample_name_ind)
      names_samples[null_samples_ind] <- paste0("Sample", null_samples_ind)
    }
  } else {
    names_samples <- paste0("Sample", 1:length(samples))
  }

  # get sample metadata
  sample_list <- names(samples)
  sample.metadata <- NULL
  for(i in 1:length(sample_list)){
    layer_list <- samples[[sample_list[i]]]@layer
    layer_data <- NULL
    for(j in 1:length(layer_list)){
      assay_list <- layer_list[[j]]@assay
      layer_data <- rbind(layer_data, cbind(names(assay_list), names(layer_list)[j]))
    }
    sample.metadata <- rbind(sample.metadata, cbind(layer_data, sample_list[i]))
  }
  sample.metadata <- data.frame(sample.metadata, row.names = paste0("Assay", 1:nrow(sample.metadata)))
  colnames(sample.metadata) <- c("Assay", "Layer", "Sample")

  sample.metadata
}
