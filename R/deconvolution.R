####
# Spot Deconvolution ####
####

RunDecon <- function(object, sc.object, sc.assay = "RNA", assay = NULL, sc.cluster = "seurat_clusters", sc.nUMI = "nCount_RNA", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # check assay type
  assay.types <- unique(AssayTypes(object, assay = assay))

  if(length(assay.types) > 1){
    stop("Please make sure that only assays of one assay type (cell/spot/ROI) are being deconvoluted at a time!")
  } else {

    # spots
    if(assay.types == "spot"){
      cat("Running RCTD for spot deconvolution ...\n")
      results <- RunRCTD(object = object, sc.object = sc.object, sc.assay = sc.assay, assay = assay, sc.cluster = sc.cluster, sc.nUMI = sc.nUMI, ...)
    }
  }

  # results
  return(results)
}

RunRCTD <- function(object, sc.object, sc.assay = "RNA", assay = NULL, sc.cluster = "seurat_clusters", sc.nUMI = "nCount_RNA", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # run a list of assays
  for(assy in assay_names){

    cur_assay <- object[[assy]]

    # RCTD
    rawdata <- RunRCTDSingle(object = cur_assay, sc.object = sc.object, sc.cluster = sc.cluster, sc.nUMI = sc.nUMI, ...)

    # Add as new assay
    cat("Adding cell type compositions as new assay:", paste(sample.metadata[assy, "Assay"], "decon", sep = "_"), "...\n")
    rawdata <- new("srAssay",
                   rawdata = rawdata, normdata = rawdata,
                   coords = cur_assay@coords[colnames(rawdata),], image = cur_assay@image,
                   type = cur_assay@type, params = cur_assay@params)
    object <- AddAssay(object,
                       newassay = rawdata,
                       newassay_name = paste(sample.metadata[assy, "Assay"], "decon", sep = "_"),
                       sample = sample.metadata[assy, "Sample"],
                       layer = sample.metadata[assy, "Layer"])
  }

  # return object with new assays
  return(object)
}

RunRCTDSingle <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", sc.nUMI = "nCount_RNA", ...){

  # create spatial data
  cat("Configuring Spatial Assay ...\n")
  spatialcounts <- Data(object, type = "raw")
  coords <- as.data.frame(Coordinates(object))
  spatialnUMI <- colSums(spatialcounts)
  spatialdata <- SpatialRNA(coords, spatialcounts, spatialnUMI)

  # create single cell reference
  cat("Configuring Single Cell Assay (reference) ...\n")
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
  cat("Calculating Cell Type Compositions ...\n")
  myRCTD <- quiet(run.RCTD(myRCTD, doublet_mode = 'full'))
  results <- as.matrix(myRCTD@results$weights)
  norm_weights <- t(sweep(results, 1, rowSums(results), "/"))

  # return
  return(norm_weights)
}
