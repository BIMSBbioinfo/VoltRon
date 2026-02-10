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
#' @importFrom utils packageVersion
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

    # Validate reference construction
    if (is.null(reference)) {
      stop(
        "Deconvolution reference was not constructed. Check your sc.object and assay parameters."
      )
    }

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
      warning(
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

        # Seurat v5 layer compatibility
        if (
          packageVersion("SeuratObject") >= "5.0.0" &&
            exists("JoinLayers", where = asNamespace("SeuratObject"))
        ) {
          sc.object <- SeuratObject::JoinLayers(sc.object, assay = sc.assay)
          sccounts <- Seurat::GetAssayData(
            sc.object,
            assay = sc.assay,
            layer = "counts"
          )
        } else {
          sccounts <- Seurat::GetAssayData(
            sc.object[[sc.assay]],
            slot = "counts"
          )
        }

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
  } else {
    stop("Unsupported assay type: ", assay.type)
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
  features = NULL,
  reference,
  method = "RCTD",
  sc.cluster,
  ...
) {
  # get assay type
  assay.type <- vrAssayTypes(object)

  # Check that exactly one assay type is returned
  if (length(assay.type) != 1) {
    stop(
      "getDeconSingle requires exactly one assay type. Detected: ",
      paste(assay.type, collapse = ", ")
    )
  }

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
  if (!requireNamespace('spacexr', quietly = TRUE)) {
    stop(
      "Please install spacexr package to use the RCTD algorithm: ",
      "devtools::install_github('dmcable/spacexr')"
    )
  }

  # Prepare target spots for final alignment
  target_spots <- vrSpatialPoints(object)

  message("Configuring Spatial Assay ...")

  # Create spatial data
  spatialcounts <- vrData(object, norm = FALSE)

  # Validate coordinates
  coords_raw <- vrCoordinates(object)

  if (!all(c("x", "y") %in% colnames(coords_raw))) {
    stop("VoltRon coordinates must contain columns named 'x' and 'y'")
  }

  coords <- as.matrix(coords_raw[, c("x", "y")])

  # Enforce Feature Selection
  if (!is.null(features)) {
    common_feats <- intersect(rownames(spatialcounts), features)
    if (length(common_feats) == 0) {
      stop("None of the requested features were found in the spatial assay.")
    }
    spatialcounts <- spatialcounts[common_feats, ]
  }

  # Align counts and coords
  common_spots <- intersect(colnames(spatialcounts), rownames(coords))
  if (length(common_spots) == 0) {
    stop("No matching spots found between count data and coordinates.")
  }

  spatialcounts <- spatialcounts[, common_spots]
  coords <- coords[common_spots, ]
  nUMI <- Matrix::colSums(spatialcounts)

  # --- spacexr >=v2.0 ---
  if (packageVersion("spacexr") >= "2.0.0") {
    message("Using native spacexr v2.0+ structures...")
    coords_df <- as.data.frame(coords)

    spatialdata <- spacexr::SpatialRNA(
      coords = coords_df,
      counts = spatialcounts,
      nUMI = nUMI
    )

    message("Initializing RCTD...")
    myRCTD <- spacexr::create.RCTD(
      spatialRNA = spatialdata,
      reference = reference,
      ...
    )

    message("Calculating Cell Type Compositions of spots with RCTD ...")

    args <- list(...)
    if (!"doublet_mode" %in% names(args)) {
      args$doublet_mode <- 'full'
    }

    myRCTD <- do.call(
      spacexr::run.RCTD,
      c(list(myRCTD), args[names(args) %in% names(formals(spacexr::run.RCTD))])
    )

    results <- myRCTD@results$weights
  } else {
    # Legacy spacexr (< 2.0)
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
      stop("Please install SingleCellExperiment package for using SCE objects")
    }
    if (!requireNamespace('SpatialExperiment', quietly = TRUE)) {
      stop("Please install SpatialExperiment package")
    }

    message("Using legacy spacexr (< 2.0) structures...")
    spatialdata <- SpatialExperiment::SpatialExperiment(
      assay = list(counts = spatialcounts),
      spatialCoords = coords
    )

    if (exists("createRctd", where = asNamespace("spacexr"))) {
      myRCTD <- spacexr::createRctd(
        spatialdata,
        reference,
        cell_type_col = sc.cluster
      )

      message("Calculating Cell Type Compositions of spots with RCTD ...")

      suppressMessages(
        myRCTD <- spacexr::runRctd(myRCTD, rctd_mode = 'full', ...)
      )
    } else {
      myRCTD <- spacexr::create.RCTD(
        spatialdata,
        reference,
        cell_type_col = sc.cluster
      )

      message("Calculating Cell Type Compositions of spots with RCTD ...")

      suppressMessages(
        myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full', ...)
      )
    }

    results <- SummarizedExperiment::assay(myRCTD, i = "weights")
  }

  # --- Spot Alignment & Zero-Filling ---
  # Standardize to [Spots x Features]
  res_rows <- nrow(results)
  res_cols <- ncol(results)
  n_target <- length(target_spots)

  if (res_cols == n_target && res_rows != n_target) {
    results <- t(results)
  }

  # Spot Filling
  if (nrow(results) != n_target) {
    full_data <- matrix(0, nrow = n_target, ncol = ncol(results))
    rownames(full_data) <- target_spots
    colnames(full_data) <- colnames(results)

    common <- intersect(rownames(full_data), rownames(results))

    if (length(common) == 0 && n_target > 0) {
      safe_target <- make.names(target_spots)
      safe_raw <- make.names(rownames(results))
      map_idx <- match(safe_target, safe_raw)
      valid_idx <- which(!is.na(map_idx))

      if (length(valid_idx) > 0) {
        full_data[valid_idx, ] <- results[map_idx[valid_idx], ]
      } else {
        stop(
          "RCTD finished, but spot names do not match the VoltRon object (neither exact nor sanitized matches). Check your sample IDs."
        )
      }
    } else {
      full_data[common, ] <- as.matrix(results[common, ,drop=FALSE])
    }
    results <- full_data
  }

  # Normalize
  row_sums <- rowSums(results)
  row_sums[row_sums == 0] <- 1
  norm_weights <- sweep(results, 1, row_sums, "/")

  # Final Transpose [Features x Spots]
  norm_weights <- t(as.matrix(norm_weights))

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
    warning(
      "The selected method is not provided for spot deconvolution. ",
      "Switching to RCTD"
    )
    method <- "RCTD"
  }

  # deconvolution with RCTD
  if (method == "RCTD") {
    # check package
    if (!requireNamespace('spacexr', quietly = TRUE)) {
      stop("Please install spacexr package to use the RCTD algorithm")
    }

    message("Configuring Single Cell Assay (reference) ...")
    if (inherits(sc.object, "Seurat")) {
      if (!requireNamespace('Seurat', quietly = TRUE)) {
        stop("Please install Seurat package for using Seurat objects")
      }

      if (is.null(sc.assay)) {
        message(
          "The sc.assay argument is not provided, using",
          Seurat::DefaultAssay(sc.object),
          "instead!"
        )
        sc.assay <- Seurat::DefaultAssay(sc.object)
      }

      # Central validation: Check if sc.cluster exists in meta.data
      if (!sc.cluster %in% colnames(sc.object@meta.data)) {
        stop(
          "The cluster column '",
          sc.cluster,
          "' was not found in the Seurat object metadata."
        )
      }

      # Seurat v5 layer compatibility
      if (
        packageVersion("SeuratObject") >= "5.0.0" &&
          exists("JoinLayers", where = asNamespace("SeuratObject"))
      ) {
        sc.object <- SeuratObject::JoinLayers(sc.object, assay = sc.assay)
        sccounts <- Seurat::GetAssayData(
          sc.object,
          assay = sc.assay,
          layer = "counts"
        )
      } else {
        tryCatch(
          {
            sccounts <- sc.object[[sc.assay]]
          },
          error = function(e) {
            stop("The assay ", sc.assay, " not found!")
          }
        )
        sccounts <- Seurat::GetAssayData(sccounts, slot = "counts")
      }

      if (!inherits(sccounts, "sparseMatrix")) {
        sccounts <- as(sccounts, "dgCMatrix")
      }
      cell_types <- as.factor(sc.object@meta.data[[sc.cluster]])

      names(cell_types) <- colnames(sc.object)
      rownames(sccounts) <- rownames(sc.object[[sc.assay]])
    } else if (inherits(sc.object, "SingleCellExperiment")) {
      if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
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
      if (!inherits(sccounts, "sparseMatrix")) {
        sccounts <- as(sccounts, "dgCMatrix")
      }
      names(cell_types) <- colnames(sc.object)
    } else {
      stop("The reference should be of either Seurat or SingleCellExperiment!")
    }

    # build reference
    sc.nUMI <- colSums(sccounts)
    names(sc.nUMI) <- colnames(sccounts)

    if (packageVersion("spacexr") >= "2.0.0") {
      message("Converting to spacexr::Reference object (v2.0+ detected)...")
      reference <- spacexr::Reference(
        counts = sccounts,
        cell_types = cell_types,
        nUMI = sc.nUMI
      )
    } else {
      names(cell_types) <- colnames(sccounts)
      cell_types <- data.frame(cell_types)
      colnames(cell_types) <- sc.cluster

      reference <- SingleCellExperiment::SingleCellExperiment(
        assays = sccounts,
        colData = cell_types
      )
    }
  }

  # return
  return(reference)
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
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop(
      "Please install Seurat package for using Seurat objects: ",
      "install.packages('Seurat')"
    )
  }
  if (!requireNamespace('MuSiC', quietly = TRUE)) {
    stop(
      "Please install MuSiC package for ROI deconvolution: ",
      "devtools::install_github('xuranw/MuSiC')"
    )
  }
  if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
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
