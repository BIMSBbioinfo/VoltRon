####
# Deconvolution ####
####

#' getDeconvolution
#'
#' Calculate deconvolution of spots and ROIs
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium),
#' see \link{SampleMetadata}.
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param features features
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell 
#' data reference
#' @param sc.cluster metadata column variable used for the single cell 
#' data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), 
#' MuSiC (ROI)
#' @param ... additional parameters passed to method specific functions,
#' e.g. RCTD, MuSiC.
#'
#' @export
getDeconvolution <- function(
  object,
  assay = NULL,
  features = NULL,
  sc.object,
  sc.assay = "RNA",
  sc.cluster = "seurat_clusters",
  method = "RCTD",
  ...
) {
  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # check assay type
  assay.types <- unique(vrAssayTypes(object, assay = assay))

  if (length(assay.types) > 1) {
    stop(
      "Please make sure that only assays of one assay type (cell/spot/ROI) are being deconvoluted at a time!"
    )
  } else {
    # make single cell reference
    reference <- getDeconReference(
      sc.object = sc.object,
      sc.assay = sc.assay,
      sc.cluster = sc.cluster,
      method = method,
      assay.type = assay.types
    )

    # run a list of assays
    for (assy in assay_names) {
      # get assay
      cur_assay <- object[[assy]]

      # RCTD
      rawdata <- getDeconSingle(
        object = cur_assay,
        features = features,
        reference = reference,
        method = method,
        sc.cluster = sc.cluster,
        ...
      )

      # add cell type mixtures as new featureset
      object <- addFeature(
        object,
        assay = assy,
        data = rawdata,
        feature_name = "Decon"
      )
    }
  }

  return(object)
}

#' getDeconReference
#'
#' Establish and produce the single cell reference for deconvolution
#'
#' @param sc.object Seurat Object
#' @param sc.assay assay of the Seurat Object used for the single cell 
#' data reference
#' @param sc.cluster metadata column variable used for the single cell 
#' data reference
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), 
#' MuSiC (ROI)
#' @param assay.type the assay type associated with the single cell 
#' deconvolution reference
#'
#' @noRd
getDeconReference <- function(
  sc.object,
  sc.assay = "RNA",
  sc.cluster = "seurat_clusters",
  method = "RCTD",
  assay.type = NULL
) {
  # Deconvolute for spots
  if (assay.type == "spot") {
    reference <- getDeconReferenceSpot(sc.object, sc.assay, sc.cluster, method)

    # Deconvolute for ROIs
  } else if (assay.type == "ROI") {
    # check method
    if (!method %in% c("MuSiC")) {
      message(
        "The selected method is not provided for ROI deconvolution. ", 
        "Switching to MuSiC ..."
      )
      method <- "MuSiC"
    }

    # deconvolution with MuSiC
    if (method == "MuSiC") {
      message("Configuring Single Cell Assay (Reference) ...")
      if (inherits(sc.object, "SingleCellExperiment")) {
        sc.object$music_decon_clusters <- sc.object[[sc.cluster]]
        reference <- sc.object
      } else if (inherits(sc.object, "Seurat")) {
        sc.object$music_decon_clusters <- sc.object@meta.data[[sc.cluster]]
        sccounts <- Seurat::GetAssayData(sc.object[[sc.assay]], slot = "counts")
        sccounts <- as.matrix(apply(sccounts, 2, ceiling))
        rownames(sccounts) <- rownames(sc.object[[sc.assay]])
        reference <- Seurat::as.SingleCellExperiment(Seurat::CreateSeuratObject(
          sccounts,
          meta.data = sc.object@meta.data
        ))
      } else {
        stop(
          "'sc.object' should either be of a Seurat or ", 
          "SingleCellExperiment class!"
        )
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
#' @param reference the single cell deconvolution reference, generated 
#' by \code{getDeconReference}
#' @param method Deconvolution method, RCTD (spot), SPOTlight (spot), 
#' MuSiC (ROI)
#' @param sc.cluster metadata column variable used for the single cell 
#' data reference
#' @param ... additional parameters passed to method specific functions
#'
#' @noRd
getDeconSingle <- function(
  object,
  features = features,
  reference,
  method = "RCTD",
  sc.cluster,
  ...
) {
  # get assay type
  assay.type <- vrAssayTypes(object)

  if (assay.type == "spot") {
    # check method
    if (!method %in% c("RCTD")) {
      message(
        "The selected method is not provided for spot deconvolution. ", 
        "Switching to RCTD"
      )
      method <- "RCTD"
    }

    if (method == "RCTD") {
      message("Running RCTD for spot deconvolution ...")
      rawdata <- getRCTD(
        object = object,
        features = features,
        reference = reference,
        sc.cluster = sc.cluster,
        ...
      )
    }
  } else if (assay.type == "ROI") {
    # check method
    if (!method %in% c("MuSiC")) {
      message(
        "The selected method is not provided for ROI deconvolution. ", 
        "Switching to MuSiC ..."
      )
      method <- "MuSiC"
    }

    if (method == "MuSiC") {
      message("Running MuSiC for ROI deconvolution ...")
      rawdata <- getMuSiC(
        object = object,
        features = features,
        reference = reference,
        sc.cluster,
        ...
      )
    }
  }

  # return
  return(rawdata)
}

####
# Spot Deconvolution ####
####

#' getRCTD
#'
#' Calculate RCTD deconvolution for spot transcriptomics
#'
#' @param object a VoltRon object
#' @param features features
#' @param reference the single cell deconvolution reference, generated 
#' \code{getDeconReference}
#' @param sc.cluster metadata column variable used for the single 
#' cell data reference
#' @param ... additional parameters passed to \code{create.RCTD} function
#'
#' @noRd
getRCTD <- function(object, features = NULL, reference, sc.cluster, ...) {
  if (!requireNamespace('spacexr')) {
    stop(
      "Please install spacexr package to use the RCTD algorithm: ", 
      "devtools::install_github('dmcable/spacexr')"
    )
  }
  if (!requireNamespace('SingleCellExperiment')) {
    stop(
      "Please install SingleCellExperiment package for using SCE objects: 
         Biocmanager::install('SingleCellExperiment')"
    )
  }

  # create spatial data
  message("Configuring Spatial Assay ...")
  spatialcounts <- vrData(object, norm = FALSE)
  coords <- as.matrix(as(vrCoordinates(object), "dgCMatrix"))[, c("x", "y")]
  # spatialnUMI <- colSums(spatialcounts)
  # spatialdata <- spacexr::createSpatialRNA(coords, spatialcounts, spatialnUMI)
  spatialdata <- SpatialExperiment::SpatialExperiment(
    assay = spatialcounts,
    spatialCoords = coords
  )

  # Run RCTD
  myRCTD <- spacexr::createRctd(
    spatialdata,
    reference,
    cell_type_col = sc.cluster
  )
  message("Calculating Cell Type Compositions of spots with RCTD ...")
  myRCTD <- quiet(spacexr::runRctd(myRCTD, rctd_mode = 'full', ...))
  results <- SummarizedExperiment::assay(myRCTD, i = "weights")
  norm_weights <- sweep(results, 2, colSums(results), "/")
  norm_weights <- as.matrix(norm_weights)

  # return
  return(norm_weights)
}

#' @noRd
getDeconReferenceSpot <- function(
  sc.object,
  sc.assay = NULL,
  sc.cluster,
  method
) {
  # check method
  if (!method %in% c("RCTD")) {
    message(
      "The selected method is not provided for spot deconvolution. ", 
      "Switching to RCTD"
    )
    method <- "RCTD"
  }

  # deconvolution with RCTD
  if (method == "RCTD") {
    # check package
    if (!requireNamespace('spacexr')) {
      stop("Please install spacexr package to use the RCTD algorithm")
    }

    message("Configuring Single Cell Assay (reference) ...")
    if (inherits(sc.object, "Seurat")) {
      if (!requireNamespace('Seurat')) {
        stop("Please install Seurat package for using Seurat objects")
      }
      if (is.null(sc.assay)) {
        message(
          "The sc.assay arguement is not provided, using",
          Seurat::DefaultAssay(sc.object),
          "instead!"
        )
      }
      tryCatch(
        {
          sccounts <- sc.object[[sc.assay]]
        },
        error = function(e) {
          stop("The assay", sc.assay, "not found!")
        }
      )
      sccounts <- Seurat::GetAssayData(sccounts, slot = "counts")
      cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])
      rownames(sccounts) <- rownames(sc.object[[sc.assay]])
      # names(cell_types) <- colnames(sc.object)
      # sc.nUMI <- colSums(sccounts)
      # names(sc.nUMI) <- colnames(sc.object)
    } else if (inherits(sc.object, "SingleCellExperiment")) {
      if (!requireNamespace('SingleCellExperiment')) {
        stop(
          "Please install SingleCellExperiment package for using SCE objects"
        )
      }
      if (!is.null(sc.assay)) {
        message(
          "The sc.assay arguement is provided by ignored for SCE objects!"
        )
      }
      sccounts <- SummarizedExperiment::assays(sc.object)[["counts"]]
      meta.data <- as.data.frame(SingleCellExperiment::colData(sc.object))
      cell_types <- as.factor(meta.data[[sc.cluster]])
      # sc.nUMI <- colSums(sccounts)
      # names(sc.nUMI) <- colnames(sc.object)
    } else {
      stop("The reference should be of either Seurat or SingleCellExperiment!")
    }

    # build reference
    sc.nUMI <- colSums(sccounts)
    names(sc.nUMI) <- colnames(sccounts)
    names(cell_types) <- colnames(sccounts)
    cell_types <- data.frame(cell_types)
    colnames(cell_types) <- sc.cluster
    # reference <- spacexr::createReference(sccounts, cell_types, sc.nUMI)
    reference <- SingleCellExperiment::SingleCellExperiment(
      assays = sccounts,
      colData = cell_types
    )
  }

  # return
  reference
}

####
# ROI Deconvolution ####
####

#' getMuSiC
#'
#' Calculate MuSiC deconvolution for ROIs
#'
#' @param object a vrAssay object
#' @param features features
#' @param reference the single cell deconvolution reference, generated 
#' \code{getDeconReference}
#' @param sc.cluster metadata column variable used for the single cell 
#' data reference
#' @param sc.samples metadata column in Seurat that provides the samples 
#' in the single cell data
#'
#' @noRd
getMuSiC <- function(
  object,
  features = NULL,
  reference,
  sc.cluster,
  sc.samples = NULL
) {
  if (!requireNamespace('Seurat')) {
    stop(
      "Please install Seurat package for using Seurat objects: ", 
      "install.packages('Seurat')"
    )
  }
  if (!requireNamespace('MuSiC')) {
    stop(
      "Please install MuSiC package for ROI deconvolution: ", 
      "devtools::install_github('xuranw/MuSiC')"
    )
  }
  if (!requireNamespace('SingleCellExperiment')) {
    stop(
      "Please install SingleCellExperiment package for ROI deconvolution: ", 
      "BiocManager::install('SingleCellExperiment')"
    )
  }

  if (is.null(sc.samples)) {
    stop(
      "Please provide a metadata column for samples for MuSiC algorithm ", 
      "to work, e.g. sc.samples = Sample"
    )
  }

  if (is.null(features)) {
    features <- vrFeatures(object)
  }

  # Single cell reference data
  reference <- reference[rownames(reference) %in% features, ]

  # data
  datax <- as.matrix(vrData(object))

  # common features
  common_features <- intersect(rownames(reference), rownames(datax))
  common_features <- intersect(common_features, features)
  if (length(common_features) < 5) {
    stop("The number of common or selected features are less than 5!")
  }
  reference <- reference[rownames(reference) %in% common_features, ]
  datax <- datax[common_features, ]

  # deconvolute
  message("Calculating Cell Type Compositions of ROIs with MuSiC ...")
  results <- MuSiC::music_prop(
    bulk.mtx = datax,
    sc.sce = reference,
    clusters = sc.cluster,
    samples = sc.samples,
    verbose = T
  )
  t(results$Est.prop.weighted)
}
