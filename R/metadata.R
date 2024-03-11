####
# Objects and Classes ####
####

## vrMetadata ####

#' The vrMetadata (VoltRon Metadata) Class
#'
#' @slot tile the metadata data table of tiles
#' @slot molecule the metadata data table of molecules
#' @slot cell the metadata data frame of cells
#' @slot spot the metadata data frame of spot
#' @slot ROI the metadata data frame of ROI
#'
#' @name vrMetadata-class
#' @rdname vrMetadata-class
#' @exportClass vrMetadata
#'
vrMetadata <- setClass(
  Class = 'vrMetadata',
  slots = c(
    molecule = 'data.table',
    cell = 'data.frame',
    spot = 'data.frame',
    ROI = 'data.frame',
    tile = 'data.table'
  )
)

### show ####

setMethod(
  f = 'show',
  signature = 'vrMetadata',
  definition = function(object) {
    cat("VoltRon Metadata Object \n")
    cat("This object includes: \n")
    lapply(slotNames(object), function(x){
      if(nrow(slot(object, name = x))){
        cat("  ", nrow(slot(object, name = x)), paste0(x, "s"), "\n")
      }
    })
    return(invisible(x = NULL))
  }
)

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

#' @param assay the assay name or type
#'
#' @rdname vrSpatialPoints
#' @method vrSpatialPoints vrMetadata
#'
vrSpatialPoints.vrMetadata <- function(object, assay = NULL) {

  # points <- c(rownames(object@molecule),
  #               rownames(object@cell),
  #               rownames(object@spot),
  #               rownames(object@ROI))

  # # get assay names if there arent
  # if(!is.null(assay)){
  #   assay_names <- vrAssayNames(object, assay = assay)
  # } else {
  #   assay_names <- NULL
  # }

  # get spatial points
  points <- unlist(lapply(slotNames(object), function(x) {
    if(x %in% c("cell", "spot", "ROI")){
      sp <- rownames(slot(object, name = x))
      if(!is.null(assay))
        sp <- sp[grepl(paste(paste0(assay, "$"), collapse = "|"), sp)]
      return(sp)
    } else {
      mdata <- slot(object, name = x)
      if(nrow(mdata) > 0){
        sp_data <- subset(slot(object, name = x), subset = assay_id %in% assay)
        return(sp_data[["id"]])
      }
    }
  }))

  # return points
  return(points)
}

#' Subsetting vrMetadata objects
#'
#' Given a vrMetadata object, subset the object given one of the attributes
#'
#' @param object a vrMetadata object
#' @param subset the subset statement
#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#' @param spatialpoints the set of spatial points to subset the object
#'
#' @method subset vrMetadata
#'
#' @importFrom rlang enquo
#' @importFrom stringr str_extract
#'
subset.vrMetadata <- function(object, subset, samples = NULL, assays = NULL, spatialpoints = NULL) {

  if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
  }

  # subset all metadata types
  if(!is.null(samples)){
    if(nrow(object@molecule) > 0){
      mol.metadata <- object@molecule[Sample %in% samples, ]
    } else {
      mol.metadata <- data.table::data.table()
    }
    cell.metadata <- object@cell[object@cell$Sample %in% samples, ]
    spot.metadata <- object@spot[object@spot$Sample %in% samples, ]
    roi.metadata <- object@ROI[object@ROI$Sample %in% samples, ]
    # tile.metadata <- object@tile[object@tile$Sample %in% samples, ]
    if(nrow(object@tile) > 0){
      tile.metadata <- object@tile[Sample %in% samples, ]
    } else {
      tile.metadata <- data.table::data.table()
    }
  } else if(!is.null(assays)){
    assay_names <- unique(lapply(slotToList(object), function(x) {
      if(inherits(x, "data.table")){
        unique(x$assay_id)
      } else {
        unique(stringr::str_extract(rownames(x), "Assay[0-9]+"))
      }
    }))
    assay_names <- unique(do.call(c,assay_names))
    if(all(assays %in% assay_names)){
      if(nrow(object@molecule) > 0) {
        mol.metadata <- object@molecule[assay_id %in% assays, ]
      } else {
        mol.metadata <- data.table::data.table()
      }
      cell.metadata <- object@cell[stringr::str_extract(rownames(object@cell), "Assay[0-9]+") %in% assays, ]
      spot.metadata <- object@spot[stringr::str_extract(rownames(object@spot), "Assay[0-9]+") %in% assays, ]
      roi.metadata <- object@ROI[stringr::str_extract(rownames(object@ROI), "Assay[0-9]+") %in% assays, ]
      if(nrow(object@tile) > 0) {
        tile.metadata <- object@tile[assay_id %in% assays, ]
      } else {
        tile.metadata <- data.table::data.table()
      }
    } else {
      if(nrow(object@molecule) > 0) {
        mol.metadata <- object@molecule[Assay %in% assays, ]
      } else {
        mol.metadata <- data.table::data.table()
      }
      cell.metadata <- object@cell[object@cell$Assay %in% assays, ]
      spot.metadata <- object@spot[object@spot$Assay %in% assays, ]
      roi.metadata <- object@ROI[object@ROI$Assay %in% assays, ]
      if(nrow(object@tile) > 0) {
        tile.metadata <- object@tile[Assay %in% assays, ]
      } else {
        tile.metadata <- data.table::data.table()
      }
    }
  } else if(!is.null(spatialpoints)){
    if(all(spatialpoints %in% vrSpatialPoints(object))){
      if(nrow(object@molecule) > 0){
        mol.metadata <- object@molecule[id %in% spatialpoints, ]
      } else {
        mol.metadata <- data.table::data.table()
      }
      cell.metadata <- object@cell[rownames(object@cell) %in% spatialpoints, ]
      spot.metadata <- object@spot[rownames(object@spot) %in% spatialpoints, ]
      roi.metadata <- object@ROI[rownames(object@ROI) %in% spatialpoints, ]
      if(nrow(object@tile) > 0){
        tile.metadata <- object@tile[id %in% spatialpoints, ]
      } else {
        tile.metadata <- data.table::data.table()
      }
    } else {
      stop("Some spatial points are not found in the metadata and the object")
    }
  } else {
    stop(paste0("No assay or sample name was provided!"))
  }

  # return new metadata
  setVRMetadata(molecule = mol.metadata, cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata, tile = tile.metadata)
}


#' subset_sampleMetadata
#'
#' Subseting sample metadata
#'
#' @param metadata sample metadata of a VoltRon object
#' @param samples the set of samples to subset the object
#' @param assays the set of assays to subset the object
#'
#' @noRd
#'
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

#' Merging vrMetadata objects
#'
#' Given a vrMetadata object, and a list of vrMetadata objects, merge all.
#'
#' @param object a vrMetadata object
#' @param object_list a list of vrMetadata objects
#'
#' @method merge vrMetadata
#'
#' @importFrom dplyr bind_rows
#' @importFrom methods slot
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
    mol.metadata <- bind_rows(methods::slot(obj1, "molecule"), methods::slot(obj2, "molecule"))
    cell.metadata <- bind_rows(methods::slot(obj1, "cell"), methods::slot(obj2, "cell"))
    spot.metadata <- bind_rows(methods::slot(obj1, "spot"), methods::slot(obj2, "spot"))
    roi.metadata <- bind_rows(methods::slot(obj1, "ROI"), methods::slot(obj2, "ROI"))
    tile.metadata <- bind_rows(methods::slot(obj1, "tile"), methods::slot(obj2, "tile"))
    combined.metadata <- setVRMetadata(molecule = mol.metadata, cell = cell.metadata, spot = spot.metadata, ROI = roi.metadata, tile = tile.metadata)
  }

  # return combined object
  return(combined.metadata)
}

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
#' @param metadata a predefined metadata
#' @param assay_name assay name
#' @param sample sample name
#' @param layer layer name
#'
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
#'
addAssay.vrMetadata <- function(object, metadata = NULL, assay, assay_name, sample = "Sample1", layer = "Section1"){

  # assay info
  assay.type <- vrAssayTypes(assay)

  # get metadata and other info
  object_metadata <- methods::slot(object, name = assay.type)
  data <- vrData(assay, norm = FALSE)

  # add new assay
  assay_ids <- vrAssayNames(object)
  assay_ids <- as.numeric(gsub("Assay", "", assay_ids))
  assay_id <- paste0("Assay", max(assay_ids)+1)

  # metadata
  if(inherits(metadata, "data.table")){

    if(!is.null(metadata)){

      if(nrow(data) > 0){
        assay_metadata<- data.table::data.table(metadata[, "id", with=FALSE], assay_id = assay_id, Count = Matrix::colSums(data),
                                                metadata[, colnames(metadata)[!colnames(metadata) %in% c("id", "assay_id", "Count", "Assay", "Layer", "Sample")], with=FALSE],
                                                Assay = assay_name, Layer = layer, Sample = sample)
      } else{
        assay_metadata <- data.table::data.table(metadata[, "id", with=FALSE], assay_id = assay_id,
                                                 metadata[, colnames(metadata)[!colnames(metadata) %in% c("id", "assay_id", "Count", "Assay", "Layer", "Sample")], with=FALSE],
                                                 Assay = assay_name, Layer = layer, Sample = sample)
      }

    }
  } else {

    # get original names
    entityID_nopostfix <- stringr::str_replace(vrSpatialPoints(assay), pattern = "_Assay[0-9]+", "")
    entityID <- stringr::str_replace(entityID_nopostfix, pattern = "$", paste0("_", assay_id))

    if(nrow(data) > 0){
      assay_metadata <- data.frame(Count = Matrix::colSums(data))
    } else {
      assay_metadata <- NULL
    }

    if(!is.null(metadata)){
      rownames_metadata <- stringr::str_replace(rownames(metadata), pattern = "_Assay[0-9]+", "")
      if(length(setdiff(rownames_metadata, entityID_nopostfix)) > 0){
        stop("Some spatial points in the metadata does not match with the assay!")
      } else{
        assay_metadata <- dplyr::bind_cols(assay_metadata,
                                           metadata[,!colnames(metadata) %in% c("Count", "Assay", "Layer", "Sample")])
        rownames(assay_metadata) <- entityID
      }
    }

    # complete assay_metadata
    assay_metadata <- dplyr::bind_cols(assay_metadata,
                                       data.frame(Assay = rep(assay_name, length(entityID)),
                                                  Layer = rep(layer, length(entityID)),
                                                  Sample = rep(sample, length(entityID))))
  }

  # add to the main metadata
  object_metadata <- dplyr::bind_rows(object_metadata, assay_metadata)
  methods::slot(object, name = assay.type) <- object_metadata

  # return
  return(object)
}

#' @rdname vrAssayNames
#' @method vrAssayNames vrMetadata
#'
#' @export
#'
vrAssayNames.vrMetadata <- function(object){

  # get assay names from metadata
  assay_names <- NULL
  for(sl in slotNames(object)){
    cur_metadata <- slot(object, name = sl)
    if(sl %in% c("molecule", "tile")){
      cur_names <- cur_metadata$assay_id
    } else {
      cur_names <- stringr::str_extract(rownames(cur_metadata), "Assay[0-9]+")
    }
    assay_names <- c(assay_names, unique(cur_names))
  }
  assay_names
}

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
    } else {
      unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
    }
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", 1:length(assaytype))
  object1 <- lapply(object_list, function(obj) {
    if(inherits(obj, "data.table")){
      temp <- obj$assay_id
      for(i in 1:length(assaytype))
        temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
      obj$assay_id <- temp
      return(obj)
    } else {
      temp <- rownames(obj)
      for(i in 1:length(assaytype))
        temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                     rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
      rownames(obj) <- temp
      return(obj)
    }
  })
  object1 <- methods::new("vrMetadata", molecule = object1$molecule, cell = object1$cell, spot = object1$spot, ROI = object1$ROI, tile = object1$tile)

  # get assay types
  object_list <- slotToList(object2)
  assaytype <- unlist(lapply(object_list, function(obj) {
    if(inherits(obj, "data.table")){
      unique(obj$assay_id)
    } else {
      unique(stringr::str_extract(rownames(obj), "Assay[0-9]+$"))
    }
  }))
  assaytype <- assaytype[order(nchar(assaytype), assaytype)]

  # replace assay names
  replacement <- paste0("Assay", (length(replacement)+1):(length(replacement) + length(assaytype)))
  object2 <- lapply(object_list, function(obj) {
    if(inherits(obj, "data.table")){
      temp <- obj$assay_id
      for(i in 1:length(assaytype))
        temp[grepl(assaytype[i], obj$assay_id)] <- replacement[i]
      obj$assay_id <- temp
      return(obj)
    } else {
      temp <- rownames(obj)
      for(i in 1:length(assaytype))
        temp[grepl(paste0(assaytype[i],"$"), rownames(obj))] <- gsub(paste0(assaytype[i],"$"), replacement[i],
                                                                     rownames(obj)[grepl(paste0(assaytype[i],"$"), rownames(obj))])
      rownames(obj) <- temp
      obj
    }
  })
  object2 <- methods::new("vrMetadata", molecule = object2$molecule, cell = object2$cell, spot = object2$spot, ROI = object2$ROI, tile = object2$tile)

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
#' @importFrom methods slot slot<-
#'
#' @noRd
changeSampleNames.vrMetadata <- function(object, sample_metadata_table){

  # get old and new samples
  old.samples <- sample_metadata_table$Sample
  new.samples <- sample_metadata_table$NewSample

  # check all types in the vrMetadata object
  new_object <- object
  all_types <- slotNames(object)
  for(type in all_types){
    metadata <- methods::slot(object, name = type)
    new_metadata <-  methods::slot(new_object, name = type)
    if(nrow(new_metadata) > 0){

      # change samples
      for(i in 1:length(old.samples))
        new_metadata$Sample[new_metadata$Sample==old.samples[i]] <- new.samples[i]

      # change layers
      for(i in 1:nrow(sample_metadata_table)){
        new_metadata$Layer[grepl(paste0(sample_metadata_table$AssayID[i], "$"), rownames(new_metadata))] <- sample_metadata_table[sample_metadata_table$AssayID[i], "NewLayer"]
      }

      # rewrite metadata type
      methods::slot(new_object, name = type) <- new_metadata
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
#' @param molecule molecule data frame
#' @param cell cell data frame
#' @param spot spot data frame
#' @param ROI ROI data frame
#' @param tile tile data frame
#'
#' @importFrom methods new
#'
#' @noRd
setVRMetadata <- function(molecule, cell, spot, ROI, tile){
  methods::new("vrMetadata", molecule = molecule, cell = cell, spot = spot, ROI = ROI, tile = tile)
}

#' setVRSampleMetadata
#'
#' @param samples a list of vrSample object
#'
#' @noRd
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
