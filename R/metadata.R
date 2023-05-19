####
# Objects and Classes ####
####

## srMetadata ####

#' The srMetadata (SpaceRover Metadata) Class
#'
#' @slot samples A list of layers for the this srSample object
#'
#' @name srMetadata-class
#' @rdname srMetadata-class
#' @exportClass srMetadata
#'
srMetadata <- setClass(
  Class = 'srMetadata',
  slots = c(
    cell = 'data.frame',
    spot = 'data.frame',
    ROI = 'data.frame'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'srMetadata',
  definition = function(object) {
    cat("SpaceRover Metadata Object \n")
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
#' @method $ srMetadata
#'
"$.srMetadata" <- function(x, i, ...) {
  return(NULL)
}

#' @export
#' @method $<- srMetadata
#'
"$<-.srMetadata" <- function(x, i, ..., value) {

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

  return(new("srMetadata", cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata))
}

#' @export
#' @method $<- srMetadata
#'
"[[<-.srMetadata" <- function(x, i, ..., value) {

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

  return(new("srMetadata", cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata))
}

####
# Methods ####
####

#' @rdname Entities
#' @method Entities srMetadata
#'
#' @export
#'
Entities.srMetadata <- function(object, ...) {

  # get the combination of cells, spots and ROIs
  points <- c(rownames(object@cell),
                rownames(object@spot),
                rownames(object@ROI))

  return(points)
}

#' @method subset srMetadata
#'
#' @aliases subset
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_extract
#'
subset.srMetadata <- function(metadata, subset, samples = NULL, assays = NULL, entities = NULL) {

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
    if(assays %in% assay_names){
      cell.metadata <- metadata@cell[stringr::str_extract(rownames(metadata@cell), "Assay[0-9]+") %in% assays, ]
      spot.metadata <- metadata@spot[stringr::str_extract(rownames(metadata@spot), "Assay[0-9]+") %in% assays, ]
      roi.metadata <- metadata@ROI[stringr::str_extract(rownames(metadata@ROI), "Assay[0-9]+") %in% assays, ]
    } else {
      cell.metadata <- metadata@cell[metadata@cell$Assay %in% assays, ]
      spot.metadata <- metadata@spot[metadata@spot$Assay %in% assays, ]
      roi.metadata <- metadata@ROI[metadata@ROI$Assay %in% assays, ]
    }
  } else if(!is.null(cells)){
    cell.metadata <- metadata@cell[rownames(metadata@cell) %in% entities, ]
    spot.metadata <- metadata@spot[rownames(metadata@spot) %in% entities, ]
    roi.metadata <- metadata@ROI[rownames(metadata@ROI) %in% entities, ]
  } else {
    stop(paste0("No assay or sample name was provided!"))
  }

  # return new metadata
  setSRMetadata(cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata)
}


#' subset.sampleMetadata
#'
#' @param metadata sample metadata of a spaceRover object
#' @param samples samples
#' @param assays assays
#'
#' @export
#'
subset.sampleMetadata <- function(metadata, samples = NULL, assays = NULL) {

  # subseting on samples, layers and assays
  if(!is.null(samples)){
    metadata <- metadata[metadata$Sample %in% samples,]
  } else if(!is.null(assays)) {
    if(assays %in% rownames(metadata)){
      metadata <- metadata[assays,]
    } else if(assays %in% metadata$Assay){
      metadata <- metadata[metadata$Assay %in% assays,]
    } else {
      stop("No assay with the names or types '", paste(assays, collapse = ", "), "' found in the object")
    }
  }
  metadata
}

#' @method merge srMetadata
#'
#' @importFrom dplyr bind_rows
#' @export
#'
merge.srMetadata <- function(object, object_list) {

  # combine all elements
  if(!is.list(object_list))
    object_list <- list(object_list)
  object_list <- c(object, object_list)

  # check if all are spaceRover
  if(!all(lapply(object_list, class) == "srMetadata"))
    stop("All arguements have to be of srMetadata class")

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
    combined.metadata <- setSRMetadata(cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata)
  }

  # return combined object
  return(combined.metadata)
}

#' updateMetadataAssay
#'
#' updating assay names for merge
#'
#' @param object1 srMetadata object
#' @param object2 srMetadata object
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
    rownames(obj) <- stringi::stri_replace_all_regex(rownames(obj),
                                                     pattern=paste0(assaytype,"$"),
                                                     replacement=replacement)
    obj
  })
  object1 <- new("srMetadata", cell = object1$cell, spot = object1$spot, ROI = object1$ROI)

  # get assay types
  object_list <- slotToList(object2)
  assaytype <- unlist(lapply(object_list, function(obj) {
    unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
  }))
  assaytype <- sort(assaytype)

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  object2 <- lapply(object_list, function(obj) {
    rownames(obj) <- stringi::stri_replace_all_regex(rownames(obj),
                                                     pattern=paste0(assaytype,"$"),
                                                     replacement=replacement)
    obj
  })
  object2 <- new("srMetadata", cell = object2$cell, spot = object2$spot, ROI = object2$ROI)

  # return
  return(list(object1 = object1, object2 = object2))
}

#' merge.sampleMetadata
#'
#' @param metadata_list a list of sample metadata of a spaceRover object
#' @param sample_name sample
#'
#' @export
#'
merge.sampleMetadata <- function(metadata_list, sample_name = NULL) {

  sample.metadata <- do.call(rbind, metadata_list)
  rownames(sample.metadata) <- paste0("Assay", 1:nrow(sample.metadata))
  if(!is.null(sample_name)){
    sample.metadata$Sample <- sample_name
    sample.metadata$Layer <- paste0("Section", 1:nrow(sample.metadata))
    unique_assay <- unique(sample.metadata$Assay)
    # if(nrow(sample.metadata) != length(unique_assay)){
    #   for(cur_assay in unique_assay){
    #     cur_assay_ind <- which(sample.metadata$Assay %in% cur_assay)
    #     # sample.metadata$Assay[cur_assay_ind] <- paste0(sample.metadata$Assay[cur_assay_ind], "_", 1:length(cur_assay_ind))
    #   }
    # }
  }
  sample.metadata
}

####
# Functions ####
####

#' setSRMetadata
#'
#' @param cell cell data frame
#' @param spot spot data frame
#' @param ROI ROI data frame
#'
setSRMetadata <- function(cell, spot, ROI){
  new("srMetadata", cell = cell, spot = spot, ROI = ROI)
}

#' setSRSampleMetadata
#'
#' @param samples a list of srSample object
#'
setSRSampleMetadata <- function(samples){

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
