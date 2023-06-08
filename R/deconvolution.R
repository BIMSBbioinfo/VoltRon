####
# Spot Deconvolution ####
####

#' getDeconvolution
#'
#' Calculate deconvolution of spatial spots and ROIs
#'
#' @param object a VoltRon object
#' @param assay assay
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell data reference
#' @param sc.cluster metadata column variable used for the single cell data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), MuSiC (ROI)
#' @param ... additional parameters passed to \code{}
#'
#' @export
#'
getDeconvolution <- function(object, assay = NULL, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", method = "RCTD", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # check assay type
  assay.types <- unique(vrAssayTypes(object, assay = assay))

  if(length(assay.types) > 1){
    stop("Please make sure that only assays of one assay type (cell/spot/ROI) are being deconvoluted at a time!")
  } else {
    # run a list of assays
    for(assy in assay_names){

      cur_assay <- object[[assy]]

      # RCTD
      rawdata <- getDeconSingle(object = cur_assay, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, method = method, ...)

      # Add as new assay
      cat("Adding cell type compositions as new assay:", paste(sample.metadata[assy, "Assay"], "decon", sep = "_"), "...\n")
      rawdata <- new("vrAssay",
                     rawdata = rawdata, normdata = rawdata,
                     coords = cur_assay@coords[colnames(rawdata),],
                     image = cur_assay@image, segments = cur_assay@segments,
                     type = cur_assay@type, params = cur_assay@params)
      object <- addAssay(object,
                         assay = rawdata,
                         assay_name = paste(sample.metadata[assy, "Assay"], "decon", sep = "_"),
                         sample = sample.metadata[assy, "Sample"],
                         layer = sample.metadata[assy, "Layer"])
    }
  }

  return(object)
}

#' getDeconvolution
#'
#' Calculate deconvolution of spots and ROIs of a single vrAssay object
#'
#' @param object a vrAssay object
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell data reference
#' @param sc.cluster metadata column variable used for the single cell data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), MuSiC (ROI)
#' @param ... additional parameters passed to \code{}
#'
getDeconSingle <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", method = "RCTD", ...){

  # get assay type
  assay.type <- vrAssayTypes(object)

  if(assay.type == "spot"){

    if(method == "RCTD"){
      cat("Running RCTD for spot deconvolution ...\n")
      rawdata <- getRCTD(object = object, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, ...)
    } else if(method == "SPOTlight") {
      cat("Running SPOTlight for spot deconvolution ...\n")
      rawdata <- getSPOTlight(object = object, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, ...)
    } else {
      stop("The selected method is not provided for spot deconvolution. Switching to RCTD")
      rawdata <- getRCTD(object = object, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, ...)
    }

  } else if(assay.type == "ROI"){

    if(method == "MuSiC"){
      cat("Running MuSiC for ROI deconvolution ...\n")
      rawdata <- getMuSiC(object = object, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, ...)
    } else {
      stop("The selected method is not provided for spot deconvolution. Switching to MuSiC")
      rawdata <- getMuSiC(object = object, sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster, ...)
    }

  }

  # return
  return(rawdata)
}

#' getRCTD
#'
#' Calculate RCTD deconvolution for spot transcriptomics
#'
#' @param object a VoltRon object
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell data reference
#' @param sc.cluster metadata column variable used for the single cell data reference
#' @param ... additional parameters passed to \code{create.RCTD} function
#'
getRCTD <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", ...){

  if (!requireNamespace('spacexr'))
    stop("Please install spacexr package to use the RCTD algorithm")

  # create spatial data
  cat("Configuring Spatial Assay ...\n")
  spatialcounts <- vrData(object, norm = FALSE)
  coords <- as.data.frame(vrCoordinates(object))
  spatialnUMI <- colSums(spatialcounts)
  spatialdata <- spacexr::SpatialRNA(coords, spatialcounts, spatialnUMI)

  # create single cell reference
  cat("Configuring Single Cell Assay (reference) ...\n")
  sccounts <- GetAssayData(sc.object[[sc.assay]], slot = "counts")
  sccounts <- as.matrix(apply(sccounts,2,ceiling))
  rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
  names(cell_types) <- colnames(sc.object)
  sc.nUMI <- colSums(sccounts)
  names(sc.nUMI) <- colnames(sc.object)
  reference <- spacexr::Reference(sccounts, cell_types, sc.nUMI)

  # Run RCTD
  myRCTD <- create.RCTD(spatialdata, reference, ...)
  cat("Calculating Cell Type Compositions of spots with RCTD ...\n")
  myRCTD <- quiet(run.RCTD(myRCTD, doublet_mode = 'full'))
  results <- as.matrix(myRCTD@results$weights)
  norm_weights <- t(sweep(results, 1, rowSums(results), "/"))

  # return
  return(norm_weights)
}

getSPOTlight <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", ...){

  if (!requireNamespace('spacexr'))
    stop("Please install spacexr package to use the RCTD algorithm")

  # create spatial data
  cat("Configuring Spatial Assay ...\n")
  spatialcounts <- vrData(object, norm = FALSE)
  coords <- as.data.frame(vrCoordinates(object))
  spatialnUMI <- colSums(spatialcounts)
  spatialdata <- spacexr::SpatialRNA(coords, spatialcounts, spatialnUMI)

  # create single cell reference
  cat("Configuring Single Cell Assay (reference) ...\n")
  sccounts <- GetAssayData(sc.object[[sc.assay]], slot = "counts")
  sccounts <- as.matrix(apply(sccounts,2,ceiling))
  rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
  names(cell_types) <- colnames(sc.object)
  sc.nUMI <- colSums(sccounts)
  names(sc.nUMI) <- colnames(sc.object)
  reference <- spacexr::Reference(sccounts, cell_types, sc.nUMI)

  # Run RCTD
  myRCTD <- create.RCTD(spatialdata, reference, ...)
  cat("Calculating Cell Type Compositions of spots with RCTD ...\n")
  myRCTD <- quiet(run.RCTD(myRCTD, doublet_mode = 'doublet'))
  results <- as.matrix(myRCTD@results$weights)
  norm_weights <- t(sweep(results, 1, rowSums(results), "/"))

  # return
  return(norm_weights)
}

#' getMuSiC
#'
#' @param object A vrAssay object
#' @param sc.object a Seurat object
#' @param sc.assay an assay in Seurat object where single cell count data is pulled
#' @param sc.cluster metadata column in Seurat provides the cell types in single cell data
#' @param sc.Samples metadata column in Seurat that provides the samples in the single cell data
#'
getMuSiC <- function(object, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", sc.Samples = NULL){

  if (!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")
  if (!requireNamespace('MuSiC'))
    stop("Please install MuSiC package for ROI deconvolution")

  if(is.null(sc.Samples))
    stop("Please provide a metadata column for samples for MuSiC algorithm to work")

  # Single cell data
  cat("Configuring Single Cell Assay (reference) ...\n")
  sccounts <- GetAssayData(sc.object[[sc.assay]], slot = "counts")
  sccounts <- as.matrix(apply(sccounts,2,ceiling))
  rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  scRNAseq <- Seurat::as.SingleCellExperiment(CreateSeuratObject(sccounts, meta.data = sc.object@meta.data))

  # deconvolute using
  cat("Calculating Cell Type Compositions of ROIs with MuSiC ...\n")
  results <- music_prop(bulk.mtx = as.matrix(vrData(object)),
                        sc.sce = scRNAseq,
                        clusters = sc.cluster,
                        samples = sc.Samples,
                        verbose = T)
  t(results$Est.prop.weighted)
}
