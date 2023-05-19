####
# Spot Deconvolution ####
####

runRCTD <- function(object, sc.object, sc.assay = "RNA", assay = NULL, sc.cluster = "seurat_clusters", sc.nUMI = "nCount_RNA", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # check assays
  if(is.null(assay))
    assay <- object@main.assay

  # get assay names
  if(assay %in% sample.metadata$Assay){
    assay_names <- rownames(sample.metadata)[sample.metadata$Assay %in% assay]
  } else {
    if(assay %in% rownames(sample.metadata)) {
      assay_names <- assay
    } else {
      stop("Assay name or type is not found in the object")
    }
  }

  # run a list of assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    rawdata <- runRCTDSingle(object = cur_assay, sc.object = sc.object, sc.cluster = sc.cluster, sc.nUMI = sc.nUMI, ...)
    rawdata <- new("srAssay", rawdata = rawdata, normdata = rawdata,
                            coords = cur_assay@coords[colnames(rawdata),], image = cur_assay@image, type = cur_assay@type, params = cur_assay@params)
    object <- AddAssay(object, newassay = rawdata, newassay_name = paste(sample.metadata[assy, "Assay"], "decon", sep = "_"),
                       sample = sample.metadata[assy, "Sample"], layer = sample.metadata[assy, "Layer"])
  }

  # return object with new assays
  return(object)
}

runRCTDSingle <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", sc.nUMI = "nCount_RNA", ...){

  # create spatial data
  spatialcounts <- Data(object, type = "raw")
  coords <- as.data.frame(Coordinates(object))
  spatialnUMI <- colSums(spatialcounts)
  spatialdata <- SpatialRNA(coords, spatialcounts, spatialnUMI)

  # create single cell reference
  sccounts <- sc.object[[sc.assay]]@counts
  sccounts <- as.matrix(apply(sccounts,2,as.integer))
  rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
  names(cell_types) <- colnames(sc.object)
  sc.nUMI <- as.integer(sc.object@meta.data[[sc.nUMI]])
  names(sc.nUMI) <- colnames(sc.object)
  reference <- Reference(sccounts, cell_types, sc.nUMI)

  # Run RCTD
  myRCTD <- create.RCTD(spatialdata, reference, ...)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- as.matrix(myRCTD@results$weights)
  norm_weights <- t(sweep(results, 1, rowSums(results), "/"))

  # return
  return(norm_weights)
}
