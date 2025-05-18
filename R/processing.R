#' @include allgenerics.R
#'
NULL

####
# Normalization ####
####

normalizeDataVoltRon <- function(object, 
                                 assay = NULL, 
                                 method = "LogNorm", 
                                 desiredQuantile = 0.9, 
                                 scale = 0.2, 
                                 sizefactor = 10000, 
                                 feat_type = NULL) {
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- normalizeData(cur_assay, 
                                    method = method, 
                                    desiredQuantile = desiredQuantile, 
                                    scale = scale, 
                                    sizefactor = sizefactor, 
                                    feat_type = feat_type)
  }
  
  # return
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class 
#  (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param method the normalization method: "LogNorm", 
#' "Q3Norm", "LogQ3Norm" or "CLR"
#' @param desiredQuantile the quantile of the data if "QuanNorm" 
#' or "LogQuanNorm" is selected as \code{method}.
#' @param scale the scale parameter for the hyperbolic arcsine transformation
#' @param sizefactor size factor if \code{method} is selected as \code{LogNorm}
#' @param feat_type the feature set type
#' 
#' @rdname normalizeData
#' @method normalizeData VoltRon
#'
#' @export
setMethod("normalizeData", "VoltRon", normalizeDataVoltRon)

normalizeDatavrAssay <- function(object, 
                                 method = "LogNorm", 
                                 desiredQuantile = 0.9, 
                                 scale = 0.2, 
                                 sizefactor = 10000, 
                                 feat_type = NULL) {
  
  # size factor
  rawdata <- vrData(object, feat_type = feat_type, norm = FALSE)
  
  if(!is.numeric(desiredQuantile)){
    stop("desiredQuantile should be numeric")
  } else {
    if(!findInterval(desiredQuantile, c(0,1)) == 1L){
      stop("desiredQuantile should be between [0,1]")
    }
  }
  
  # normalization method
  if(method == "LogNorm"){
    normdata <- LogNorm(rawdata, colSums(rawdata), sizefactor)
  } else if(method == "Q3Norm") {
    # rawdata[rawdata==0] <- 1
    qs <- getColQuantiles(rawdata, desiredQuantile)
    normdata <- getDivideSweep(rawdata, qs / exp(mean(log(qs))))
  } else if(method == "LogQ3Norm") {
    # rawdata[rawdata==0] <- 1
    qs <- getColQuantiles(rawdata, desiredQuantile)
    normdata <- getDivideSweep(rawdata, qs / exp(mean(log(qs))))
    normdata <- log(normdata + 1)
  } else if(method == "CLR") {
    normdata <- getDivideSweep(rawdata, colSums(rawdata))
    normdata <- apply(normdata, 2, function(x) {
      log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    })
  } else if(method == "hyper.arcsine") {
    normdata <- asinh(rawdata/scale)
  } else {
    stop('Please select one of these methods: "LogNorm",', 
         ' "Q3Norm", "LogQ3Norm" or "CLR"')
  }
  
  # get normalized data
  catch_connect1 <- try(slot(object, name = "data"), silent = TRUE)
  catch_connect2 <- try(slot(object, name = "rawdata"), silent = TRUE)
  if(!is(catch_connect1, 'try-error') && !methods::is(catch_connect1,'error')){
    if(is.null(feat_type))
      feat_type <- vrMainFeatureType(object)
    object@data[[paste0(feat_type, "_norm")]] <- normdata
  } else if(!is(catch_connect2, 'try-error') && !methods::is(catch_connect2,'error')){
    object@normdata <- normdata
  }
  
  # return
  return(object)
}

#' @rdname normalizeData
#' @method normalizeData vrAssay
#'
#' @importFrom stats quantile
#'
#' @export
setMethod("normalizeData", "vrAssay", normalizeDatavrAssay)

#' @rdname normalizeData
#' @method normalizeData vrAssayV2
#'
#' @importFrom stats quantile
#'
#' @export
setMethod("normalizeData", "vrAssayV2", normalizeDatavrAssay)

LogNorm <- function(rawdata, coldepth, sizefactor){
  if(inherits(rawdata, "IterableMatrix")){
    if(!requireNamespace("BPCells"))
      stop("You have to install BPCells!: 
           remotes::install_github('bnprks/BPCells/r')")
    normdata <- BPCells::t(BPCells::t(rawdata)/coldepth)
    normdata <- BPCells::log1p_slow(normdata*sizefactor)
  } else if(inherits(rawdata, "DelayedArray")){
    if(!requireNamespace("DelayedArray"))
      stop("You have to install DelayedArray!: 
           BiocManager::install('DelayedArray')")
    # normdata <- DelayedArray::sweep(rawdata, 2L, coldepth, FUN = "/")
    normdata <- t(t(rawdata)/coldepth)
    normdata <- log(normdata*sizefactor + 1)
  } else {
    normdata <- sweep(rawdata, 2L, coldepth, FUN = "/")
    normdata <- log(normdata*sizefactor + 1)
  }
  return(normdata)
}

getDivideSweep <- function(rawdata, divisor){
  if(inherits(rawdata, "IterableMatrix")){
    if(!requireNamespace("BPCells"))
      stop("You have to install BPCells!: 
           remotes::install_github('bnprks/BPCells/r')")
    return(BPCells::t(BPCells::t(rawdata)/divisor))
  } else if(inherits(rawdata, "DelayedArray")){
    if(!requireNamespace("DelayedArray"))
      stop("You have to install DelayedArray!: 
           BiocManager::install('DelayedArray')")
    return(t(t(rawdata)/divisor))
  } else {
    return(sweep(rawdata, 2L, divisor, FUN = "/"))
  }
  return(rawdata)
}

####
# Features ####
####

getFeaturesVoltRon <- function(object, assay = NULL, max.count = 1, n = 3000){
  
  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)
  
  # get features for all coordinates
  for(assy in assay_names){
    object[[assy]] <- getFeatures(object[[assy]], max.count = max.count, n = n)
  }
  
  # return
  return(object)
}

#' @param assay assay name (exp: Assay1) or assay class 
#' (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param max.count maximum count (across spatial points) for 
#' low count filtering
#' @param n the top number of variable features 
#' 
#' @rdname getFeatures
#'
#' @export
setMethod("getFeatures", "VoltRon", getFeaturesVoltRon)

getFeaturesvrAssay <- function(object, max.count = 1, n = 3000){
  
  # get data and coordinates
  rawdata <- vrData(object, norm = FALSE)
  coords <- vrCoordinates(object)
  features <- vrFeatures(object)
  
  # eliminate genes with low counts
  # keep.genes <- which(apply(rawdata,1,max) > max.count)
  keep.genes <- getMaxCount(rawdata, max.count)
  
  # vst estimation
  vst_data <- getVstData(rawdata)
  loess_data <- vst_data[keep.genes,]
  loess_results <- stats::loess(var~mean, loess_data, span = 0.3)
  vst_data$adj_var <- 0
  vst_data$rank <- 0
  vst_data[keep.genes,]$adj_var <- stats::predict(loess_results)
  vst_data[keep.genes,]$rank <- order(order(vst_data$adj_var[keep.genes], 
                                            decreasing = TRUE))
  
  # set feature data
  vrFeatureData(object) <- vst_data
  
  # return
  return(object)
}

#' @rdname getFeatures
#'
#' @importFrom stats loess predict var
#' @importFrom Matrix rowMeans
#'
#' @export
setMethod("getFeatures", "vrAssay", getFeaturesvrAssay)

#' @rdname getFeatures
#'
#' @importFrom stats loess predict var
#' @importFrom Matrix rowMeans
#'
#' @export
setMethod("getFeatures", "vrAssayV2", getFeaturesvrAssay)

getVstData <- function(rawdata){
  if(inherits(rawdata, "IterableMatrix")){
    if(!requireNamespace("BPCells"))
      stop("You have to install BPCells!: 
           remotes::install_github('bnprks/BPCells/r')")
    mean_data <- BPCells::rowMeans(rawdata)
    # var_data <- BPCells::rowSums(rawdata^2)
    # var_data <- (var_data - mean_data^2/nrow(rawdata))/(nrow(rawdata)-1)
    var_data <- BPCells::rowVars(rawdata)
  } else if(inherits(rawdata, "DelayedArray")){
    if(!requireNamespace("DelayedMatrixStats"))
      stop("You have to install DelayedMatrixStats!: 
           BiocManager::install('DelayedMatrixStats')")
    mean_data <- DelayedMatrixStats::rowMeans2(rawdata)
    var_data <- DelayedMatrixStats::rowVars(rawdata)
  } else {
    mean_data <- Matrix::rowMeans(rawdata)
    var_data <- apply(rawdata, 1, stats::var)
  }
  vst_data <- data.frame(mean = mean_data, var = var_data)
  return(vst_data)
}

getMaxCount <- function(rawdata, max.count){
  if(inherits(rawdata, "IterableMatrix")){
    if(!requireNamespace("BPCells"))
      stop("You have to install BPCells!: 
           remotes::install_github('bnprks/BPCells/r')")
    rawdata <- rawdata > max.count
    keep.genes <- which(BPCells::rowSums(rawdata) > 0)
  } else if(inherits(rawdata, "DelayedArray")){
    if(!requireNamespace("DelayedMatrixStats"))
      stop("You have to install DelayedMatrixStats!: 
           BiocManager::install('DelayedMatrixStats')")
    rawdata <- rawdata > max.count
    keep.genes <- which(DelayedMatrixStats::rowSums2(rawdata) > 0)
  } else {
    keep.genes <- which(apply(rawdata,1,max) > max.count)
  }
  return(keep.genes)
}

#' getVariableFeatures
#'
#' get shared variable features across multiple assays
#'
#' @param object a VoltRon Object
#' @param assay assay name (exp: Assay1) or assay class 
#' (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param n the number of features
#' @param ... additional arguements passed to \link{vrFeatureData}
#'
#' @importFrom dplyr full_join
#' @importFrom utils head
#'
#' @export
getVariableFeatures <- function(object, assay = NULL, n = 3000, ...){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get features for all coordinates
  ranks <- NULL
  for(assy in assay_names){
    feature_data <- vrFeatureData(object[[assy]], ...)
    if(!is.null(feature_data)) {
      if(nrow(feature_data) > 0){
        feature_data$gene <- rownames(feature_data)
      }
    } else {
      feature_data <- data.frame(gene = vrFeatures(object[[assy]]), rank = NA)
    }
    if(is.null(ranks)){
      ranks <- feature_data[,c("gene", "rank")]
    } else {
      ranks <- ranks %>% full_join(feature_data[,c("gene", "rank")], 
                                   by = c("gene" = "gene"))
    }
  }

  # get geometric mean of ranks, i.e. rank product statistic
  rownames_ranks <- ranks$gene
  ranks <- ranks[,!colnames(ranks) %in% "gene", drop = FALSE]
  ranks <- apply(ranks, 1, function(x) exp(mean(log(x))))
  names(ranks) <- rownames_ranks
  ranks <- ranks[ranks != 0]

  # get selected features
  if(length(ranks[!is.na(ranks)]) > 0){
    selected_features <- names(utils::head(sort(ranks, decreasing = FALSE), n))
  } else {
    selected_features <- vrFeatures(object, assay = assay)
  }

  # return
  return(selected_features)
}

####
# vrEmbeddings ####
####

#' getPCA
#'
#' calculate PCA of the VoltRon objects
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class 
#' (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param features the selected features for PCA reduction
#' @param dims the number of dimensions extracted from PCA
#' @param type the key name for the embedding, default: pca
#' @param n.workers the number of cores/workers use for parallelization.
#' @param overwrite Whether the existing embedding with name 'type' 
#' should be overwritten in \link{vrEmbeddings}
#' @param seed seed
#'
#' @importFrom BiocSingular runPCA FastAutoParam
#'
#' @export
getPCA <- function(object, 
                   assay = NULL, 
                   features = NULL, 
                   dims = 30, 
                   type = "pca", 
                   n.workers = 1, 
                   overwrite = FALSE, 
                   seed = 1){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get shared features and subset
  assay_features <- vrFeatures(object, assay = assay)

  # if there are features of a VoltRon object, then get variable features too
  if(length(assay_features) > 0) {
    if(is.null(features))
      features <- getVariableFeatures(object, assay = assay)
    object_subset <- subsetVoltRon(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)

    # adjust extraction features length
    if(dims > length(features)){
      message("Requested more PC dimensions than existing ", 
              "features: dims = length(features) now!")
      dims <- length(features)
    }

  # if there are no features in VoltRon object, return the assay as itself
  } else {
    object_subset <- object
  }

  # get data
  normdata <- vrData(object_subset, assay = assay, norm = TRUE)

  # get PCA embedding
  set.seed(seed)
  if(inherits(normdata, "IterableMatrix")){
    if(!requireNamespace("BPCells"))
      stop("You have to install BPCells!: ", 
           "remotes::install_github('bnprks/BPCells/r')")
    svd <- BPCells::svds(normdata, k=dims, threads = as.integer(n.workers))
    pr.data <- BPCells::multiply_cols(svd$v, svd$d)
  } else {
    if(n.workers > 1){
      if(!requireNamespace("BiocParallel"))
        stop("You have to install BiocParallel!: ", 
             "BiocManager::install('BiocParallel')")
      pr.data <- 
        BiocSingular::runPCA(t(normdata), 
                             rank=dims,
                             scale=TRUE,
                             center=TRUE, 
                             BPPARAM = BiocParallel::MulticoreParam(n.workers), 
                             BSPARAM=BiocSingular::FastAutoParam()) 
    } else {
      pr.data <- 
        BiocSingular::runPCA(t(normdata), 
                             rank=dims,
                             scale=TRUE,
                             center=TRUE, 
                             BSPARAM=BiocSingular::FastAutoParam())
    }
    pr.data <- pr.data$x 
  }
  
  # change colnames
  colnames(pr.data) <- paste0("PC", seq_len(dims))
  rownames(pr.data) <- colnames(normdata)

  # set Embeddings
  vrEmbeddings(object, 
               assay = assay, 
               type = type, 
               overwrite = overwrite) <- pr.data

  # return
  return(object)
}

#' getUMAP
#'
#' calculate UMAP of the VoltRon objects
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class 
#' (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
#' @param data.type the type of data used to calculate UMAP from: 
#' "pca" (default), "raw" or "norm"
#' @param dims the number of dimensions extracted from PCA
#' @param umap.key the name of the umap embedding, default: umap
#' @param overwrite Whether the existing embedding with name 'type' 
#' should be overwritten in \link{vrEmbeddings}
#' @param seed seed
#'
#' @importFrom uwot umap
#' @importFrom Matrix t
#'
#' @export
#'
getUMAP <- function(object, 
                    assay = NULL, 
                    data.type = "pca", 
                    dims = seq_len(30), 
                    umap.key = "umap", 
                    overwrite = FALSE, 
                    seed = 1){

  # get data
  if(data.type %in% c("raw", "norm")){
    data <- vrData(object, assay = assay, norm = (data.type == "norm"))
    data <- as.matrix(as(Matrix::t(data),"dgCMatrix"))
  } else{
    embedding_names <- vrEmbeddingNames(object)
    if(data.type %in% vrEmbeddingNames(object)) {
      data <- vrEmbeddings(object, 
                           assay = assay, 
                           type = data.type, 
                           dims = dims)
    } else {
      stop("Please provide a data type from one of ", 
           "three choices: raw, norm and pca")
    }
  }

  # get umap
  set.seed(seed)
  umap_data <- uwot::umap(data)
  colnames(umap_data) <- c("x", "y")
  vrEmbeddings(object, 
               assay = assay, 
               type = umap.key, 
               overwrite = overwrite) <- umap_data

  # return
  return(object)
}

####
# Image Processing ####
####

#' split_into_tiles
#'
#' split image raster data into tiles
#'
#' @param image_data image raster data
#' @param tile_size tile size
#'
#' @noRd
split_into_tiles <- function(image_data, tile_size = 10) {
  n_rows <- nrow(image_data)
  n_cols <- ncol(image_data)

  # Calculate the number of tiles in rows and columns
  n_row_tiles <- n_rows %/% tile_size
  n_col_tiles <- n_cols %/% tile_size

  # Initialize an empty list to store tiles
  tiles <- list()

  # Loop through the image data matrix to extract tiles
  for (i in seq_len(n_row_tiles)) {
    for (j in seq_len(n_col_tiles)) {
      # Calculate the indices for the current tile
      start_row <- (i - 1) * tile_size + 1
      end_row <- i * tile_size
      start_col <- (j - 1) * tile_size + 1
      end_col <- j * tile_size

      # Extract the current tile from the image data matrix
      tile <- image_data[start_row:end_row, start_col:end_col]

      # Store the tile in the list
      tiles[[length(tiles) + 1]] <- tile
    }
  }

  # Return the list of tiles
  return(tiles)
}
