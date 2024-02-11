####
# Spot Deconvolution ####
####

#' getDeconvolution
#'
#' Calculate deconvolution of spatial spots and ROIs
#'
#' @param object a VoltRon object
#' @param assay assay
#' @param features features
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell data reference
#' @param sc.cluster metadata column variable used for the single cell data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), MuSiC (ROI)
#' @param ... additional parameters passed to \code{getDeconSingle}
#'
#' @export
#'
getDeconvolution <- function(object, assay = NULL, features = NULL, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", method = "RCTD", ...){

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # check assay type
  assay.types <- unique(vrAssayTypes(object, assay = assay))

  if(length(assay.types) > 1){
    stop("Please make sure that only assays of one assay type (cell/spot/ROI) are being deconvoluted at a time!")
  } else {

    # make single cell reference
    reference <- getDeconReference(sc.object = sc.object, sc.assay = sc.assay, sc.cluster = sc.cluster,
                                   method = method, assay.type = assay.types)

    # run a list of assays
    for(assy in assay_names){

      # get assay
      cur_assay <- object[[assy]]

      # RCTD
      # rawdata <- getDeconSingle(object = cur_assay, features = features, sc.object = sc.object, sc.assay = sc.assay,
      #                           sc.cluster = sc.cluster, method = method, ...)
      rawdata <- getDeconSingle(object = cur_assay, features = features, reference = reference, method = method, ...)

      # create new assay
      cat("Adding cell type compositions as new assay:", paste(sample.metadata[assy, "Assay"], "decon", sep = "_"), "...\n")
      spatialpoints <- colnames(rawdata)
      new_assay <- formAssay(data = rawdata,
                             coords = vrCoordinates(cur_assay)[spatialpoints,], segments = vrSegments(cur_assay)[spatialpoints],
                             image = vrImages(cur_assay), type = cur_assay@type, params = cur_assay@params, name = cur_assay@name, main_image = cur_assay@main_image)
      new_assay@image <- cur_assay@image
      new_assay <- subset(new_assay, spatialpoints = spatialpoints)

      # add the new assay
      object <- addAssay(object,
                         assay = new_assay,
                         metadata = Metadata(object, assay = assy)[spatialpoints,],
                         assay_name = paste(sample.metadata[assy, "Assay"], "decon", sep = "_"),
                         sample = sample.metadata[assy, "Sample"],
                         layer = sample.metadata[assy, "Layer"])
    }
  }

  return(object)
}

#' getDeconReference
#'
#' Establish and produce the single cell reference for deconvolution
#'
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell data reference
#' @param sc.cluster metadata column variable used for the single cell data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), MuSiC (ROI)
#' @param assay.type the assay type associated with the single cell deconvolution reference
#'
#' @noRd
getDeconReference <- function(sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", method = "RCTD", assay.type = NULL){

  # Deconvolute for spots
  if(assay.type == "spot"){

    # check method
    if(!method %in% c("RCTD")){
      message("The selected method is not provided for spot deconvolution. Switching to RCTD")
      method <- "RCTD"
    }

    # deconvolution with RCTD
    if(method == "RCTD"){

      if (!requireNamespace('spacexr'))
        stop("Please install spacexr package to use the RCTD algorithm")
      if (!requireNamespace('Seurat'))
        stop("Please install Seurat package for using Seurat objects")

      cat("Configuring Single Cell Assay (reference) ...\n")
      sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
      sccounts <- as.matrix(apply(sccounts,2,ceiling))
      rownames(sccounts) <- rownames(sc.object[[sc.assay]])
      cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
      names(cell_types) <- colnames(sc.object)
      sc.nUMI <- colSums(sccounts)
      names(sc.nUMI) <- colnames(sc.object)
      reference <- spacexr::Reference(sccounts, cell_types, sc.nUMI)
    }

  # Deconvolute for ROIs
  } else if(assay.type == "ROI"){

    # check method
    if(!method %in% c("MuSiC")){
      message("The selected method is not provided for ROI deconvolution. Switching to MuSiC")
      method <- "MuSiC"
    }

    # deconvolution with MuSiC
    if(method == "MuSiC"){

      cat("Configuring Single Cell Assay (reference) ...\n")
      if(inherits(sc.object, "SingleCellExperiment")){
        sc.object$music_decon_clusters <- sc.object[[sc.cluster]]
        reference <- sc.object
      } else if(inherits(sc.object, "Seurat")){
        sc.object$music_decon_clusters <- sc.object@meta.data[[sc.cluster]]
        sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
        sccounts <- as.matrix(apply(sccounts,2,ceiling))
        rownames(sccounts) <- rownames(sc.object[[sc.assay]])
        reference <- Seurat::as.SingleCellExperiment(Seurat::CreateSeuratObject(sccounts, meta.data = sc.object@meta.data))
      } else{
        stop("'sc.object' should either be of a Seurat or SingleCellExperiment class!")
      }

    }
  }

  # return
  return(reference)
}

#' getDeconSingle
#'
#' Calculate deconvolution of spots and ROIs of a single vrAssay object
#'
#' @param object a vrAssay object
#' @param features features
#' @param reference the single cell deconvolution reference, see \code{getDeconReference}
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), MuSiC (ROI)
#' @param ... additional parameters passed to method specific functions
#'
#' @noRd
getDeconSingle <- function(object, features = features, reference, method = "RCTD", ...){

  # get assay type
  assay.type <- vrAssayTypes(object)

  if(assay.type == "spot"){

    # check method
    if(!method %in% c("RCTD")){
      message("The selected method is not provided for spot deconvolution. Switching to RCTD")
      method <- "RCTD"
    }

    if(method == "RCTD"){
      cat("Running RCTD for spot deconvolution ...\n")
      rawdata <- getRCTD(object = object, features = features, reference = reference, ...)
    }

  } else if(assay.type == "ROI"){

    # check method
    if(!method %in% c("MuSiC")){
      message("The selected method is not provided for ROI deconvolution. Switching to MuSiC")
      method <- "MuSiC"
    }

    if(method == "MuSiC"){
      cat("Running MuSiC for ROI deconvolution ...\n")
      rawdata <- getMuSiC(object = object, features = features, reference = reference, ...)
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
#' @param features features
#' @param reference the single cell deconvolution reference, see \code{getDeconReference}
#' @param ... additional parameters passed to \code{create.RCTD} function
#'
#' @noRd
getRCTD <- function(object, features = NULL, reference, ...){

  if (!requireNamespace('spacexr'))
    stop("Please install spacexr package to use the RCTD algorithm")
  if (!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")

  # create spatial data
  cat("Configuring Spatial Assay ...\n")
  spatialcounts <- vrData(object, norm = FALSE)
  coords <- as.data.frame(vrCoordinates(object))
  spatialnUMI <- colSums(spatialcounts)
  spatialdata <- spacexr::SpatialRNA(coords, spatialcounts, spatialnUMI)

  # create single cell reference
  # cat("Configuring Single Cell Assay (reference) ...\n")
  # sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
  # sccounts <- as.matrix(apply(sccounts,2,ceiling))
  # rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  # cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
  # names(cell_types) <- colnames(sc.object)
  # sc.nUMI <- colSums(sccounts)
  # names(sc.nUMI) <- colnames(sc.object)
  # reference <- spacexr::Reference(sccounts, cell_types, sc.nUMI)

  # Run RCTD
  myRCTD <- spacexr::create.RCTD(spatialdata, reference, ...)
  cat("Calculating Cell Type Compositions of spots with RCTD ...\n")
  myRCTD <- quiet(spacexr::run.RCTD(myRCTD, doublet_mode = 'full'))
  results <- as.matrix(myRCTD@results$weights)
  norm_weights <- t(sweep(results, 1, rowSums(results), "/"))

  # return
  return(norm_weights)
}

# getSPOTlight <- function(object, features = NULL, sc.object, sc.assay = "RNA", sc.cluster = "seurat_clusters", ...){
#
#   if (!requireNamespace('spacexr'))
#     stop("Please install spacexr package to use the RCTD algorithm")
#   if (!requireNamespace('Seurat'))
#     stop("Please install Seurat package for using Seurat objects")
#
#   # create spatial data
#   cat("Configuring Spatial Assay ...\n")
#   spatialcounts <- vrData(object, norm = FALSE)
#   coords <- as.data.frame(vrCoordinates(object))
#   spatialnUMI <- colSums(spatialcounts)
#   spatialdata <- spacexr::SpatialRNA(coords, spatialcounts, spatialnUMI)
#
#   # create single cell reference
#   cat("Configuring Single Cell Assay (reference) ...\n")
#   sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
#   sccounts <- as.matrix(apply(sccounts,2,ceiling))
#   rownames(sccounts) <- rownames(sc.object[[sc.assay]])
#   cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
#   names(cell_types) <- colnames(sc.object)
#   sc.nUMI <- colSums(sccounts)
#   names(sc.nUMI) <- colnames(sc.object)
#   reference <- spacexr::Reference(sccounts, cell_types, sc.nUMI)
#
#   # Run RCTD
#   myRCTD <- spacexr::create.RCTD(spatialdata, reference, ...)
#   cat("Calculating Cell Type Compositions of spots with RCTD ...\n")
#   myRCTD <- quiet(spacexr::run.RCTD(myRCTD, doublet_mode = 'doublet'))
#   results <- as.matrix(myRCTD@results$weights)
#   norm_weights <- t(sweep(results, 1, rowSums(results), "/"))
#
#   # return
#   return(norm_weights)
# }

#' getMuSiC
#'
#' Calculate MuSiC deconvolution for ROIs
#'
#' @param object A vrAssay object
#' @param features features
#' @param reference the single cell deconvolution reference, see \code{getDeconReference}
#' @param sc.samples metadata column in Seurat that provides the samples in the single cell data
#'
#' @noRd
getMuSiC <- function(object, features = NULL, reference, sc.samples = NULL){

  if (!requireNamespace('Seurat'))
    stop("Please install Seurat package for using Seurat objects")
  if (!requireNamespace('MuSiC'))
    stop("Please install MuSiC package for ROI deconvolution")
  if (!requireNamespace('SingleCellExperiment'))
    stop("Please install SingleCellExperiment package for ROI deconvolution")

  if(is.null(sc.samples))
    stop("Please provide a metadata column for samples for MuSiC algorithm to work, e.g. sc.samples = Sample")

  if(is.null(features)) {
    features <- vrFeatures(object)
  }

  # Single cell data
  # cat("Configuring Single Cell Assay (reference) ...\n")
  # if(inherits(sc.object, "SingleCellExperiment")){
  #   scRNAseq <- sc.object
  # } else if(inherits(sc.object, "Seurat")){
  #   sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
  #   sccounts <- as.matrix(apply(sccounts,2,ceiling))
  #   rownames(sccounts) <- rownames(sc.object[[sc.assay]])
  #   scRNAseq <- Seurat::as.SingleCellExperiment(Seurat::CreateSeuratObject(sccounts, meta.data = sc.object@meta.data))
  # } else{
  #   stop("'sc.object' should either be of a Seurat or SingleCellExperiment class!")
  # }
  reference <- reference[rownames(reference) %in% features,]

  # data
  datax <- as.matrix(vrData(object))

  # common features
  common_features <- intersect(rownames(reference), rownames(datax))
  common_features <- intersect(common_features, features)
  if(length(common_features) < 5)
    stop("The number of common or selected features are less than 5!")
  reference <- reference[rownames(reference) %in% common_features,]
  datax <- datax[common_features,]

  # deconvolute
  cat("Calculating Cell Type Compositions of ROIs with MuSiC ...\n")
  results <- MuSiC::music_prop(bulk.mtx = datax,
                        sc.sce = reference,
                        # clusters = sc.cluster,
                        clusters = "music_decon_clusters",
                        samples = sc.samples,
                        verbose = T)
  t(results$Est.prop.weighted)
}
