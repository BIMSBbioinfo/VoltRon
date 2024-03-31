#' @include generics.R
#'
NULL

####
# Normalization ####
####

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param ... additional parameters passed to \code{normalizeData.vrAssay}
#'
#' @rdname normalizeData
#' @method normalizeData VoltRon
#'
#' @export
normalizeData.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- normalizeData(cur_assay, ...)
  }

  # return
  return(object)
}

#' @param method the normalization method: "LogNorm", "Q3Norm", "LogQ3Norm" or "CLR"
#' @param desiredQuantile the quantile of the data if "QuanNorm" or "LogQuanNorm" is selected as \code{method}.
#' @param scale the scale parameter for the hyperbolic arcsine transformation
#' @param sizefactor size factor if \code{method} is selected as \code{LogNorm}
#'
#' @rdname normalizeData
#' @method normalizeData vrAssay
#'
#' @importFrom stats quantile
#'
#' @export
normalizeData.vrAssay <- function(object, method = "LogNorm", desiredQuantile = 0.9, scale = 0.2, sizefactor = 10000) {

  # size factor
  rawdata <- vrData(object)
  coldepth <- colSums(rawdata)

  if(!is.numeric(desiredQuantile)){
    stop("desiredQuantile should be numeric")
  } else {
    if(!findInterval(desiredQuantile, c(0,1)) == 1L){
      stop("desiredQuantile should be between [0,1]")
    }
  }

  # normalization method
  if(method == "LogNorm"){
    depth <- matrix(rep(coldepth, nrow(rawdata)), byrow = T, nrow = nrow(rawdata))
    normdata <- (rawdata/depth)*sizefactor
    normdata <- log(normdata + 1)
  } else if(method == "Q3Norm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
  } else if(method == "LogQ3Norm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
    normdata <- log(normdata + 1)
  } else if(method == "CLR") {
    normdata <- apply(rawdata, 2, function(x) {
      log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    })
  } else if(method == "hyper.arcsine") {
    normdata <- asinh(rawdata/scale)
  } else {
    stop('Please select one of these methods: "LogNorm", "Q3Norm", "LogQ3Norm" or "CLR"')
  }

  # get normalized data
  object@normdata <- normdata

  # return
  return(object)
}

####
# Features ####
####

#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param ... arguements passed to other methods
#' 
#' @rdname getFeatures
#'
#' @export
getFeatures.VoltRon <- function(object, assay = NULL, ...){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get features for all coordinates
  for(assy in assay_names){
    object[[assy]] <- getFeatures(object[[assy]], ...)
  }

  # return
  return(object)
}

#' @param max.count maximum count (across spatial points) for low count filtering
#' @param n the top number of variable features 
#'
#' @rdname getFeatures
#'
#' @importFrom stats loess predict var
#'
#' @export
getFeatures.vrAssay <- function(object, max.count = 1, n = 3000){

  # get data and coordinates
  rawdata <- vrData(object, norm = FALSE)
  coords <- vrCoordinates(object)
  features <- vrFeatures(object)

  # eliminate genes with low counts
  keep.genes <- which(apply(rawdata,1,max) > max.count)
  rawdata_subset <- rawdata[keep.genes,]

  # vst estimation
  vst_data <- data.frame(mean = rowMeans(rawdata), var = apply(rawdata, 1, stats::var))
  loess_data <- vst_data[keep.genes,]
  loess_results <- stats::loess(var~mean, loess_data, span = 0.3)
  vst_data$adj_var <- 0
  vst_data$rank <- 0
  vst_data[keep.genes,]$adj_var <- stats::predict(loess_results)
  vst_data[keep.genes,]$rank <- order(order(vst_data$adj_var[keep.genes], decreasing = TRUE))

  # set feature data
  vrFeatureData(object) <- vst_data

  # return
  return(object)
}


#' getSharedFeatures
#'
#' get shared variable features across multiple assays
#'
#' @param object a Voltron Object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param n the number of features
#' @param ... additional arguements passed to \code{vrFeatureData}
#'
#' @importFrom dplyr full_join
#' @importFrom utils head
#'
getVariableFeatures <- function(object, assay = NULL, n = 3000, ...){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get features for all coordinates
  ranks <- NULL
  for(assy in assay_names){
    feature_data <- vrFeatureData(object[[assy]], ...)
    if(nrow(feature_data) > 0){
      feature_data$gene <- rownames(feature_data)
    } else {
      feature_data <- data.frame(gene = vrFeatures(object[[assy]]), rank = NA)
    }
    if(is.null(ranks)){
      ranks <- feature_data[,c("gene", "rank")]
    } else {
      ranks <- ranks %>% full_join(feature_data[,c("gene", "rank")], by = c("gene" = "gene"))
    }
  }

  # get geometric mean of ranks, i.e. rank product statistic
  ranks <- ranks[,!colnames(ranks) %in% "gene", drop = FALSE]
  ranks <- apply(ranks, 1, function(x) exp(mean(log(x))))
  names(ranks) <- rownames(feature_data)
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
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param features the selected features for PCA reduction
#' @param dims the number of dimensions extracted from PCA
#' @param overwrite Whether the existing embedding with name 'type' should be overwritten in \code{vrEmbeddings}
#' @param seed seed
#' @param ... additional parameters passed to \code{vrEmbeddings}
#'
#' @importFrom irlba irlba
#' @importFrom dplyr left_join
#'
#' @export
#'
getPCA <- function(object, assay = NULL, features = NULL, dims = 30, overwrite = FALSE, seed = 1, ...){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get shared features and subset
  assay_features <- vrFeatures(object, assay = assay)

  # if there are features of a VoltRon object, then get variable features too
  if(length(assay_features) > 0) {
    if(is.null(features))
      features <- getVariableFeatures(object, assay = assay)
    object_subset <- subset(object, features = features)
    vrMainAssay(object_subset) <- vrMainAssay(object)

    # adjust extraction features length
    if(dims > length(features)){
      message("Requested more PC dimensions than existing features: dims = length(features) now!")
      dims <- length(features)
    }

  # if there are no features in VoltRon object, return the assay as itself
  } else {
    object_subset <- object
  }

  # get data
  normdata <- vrData(object_subset, assay = assay, norm = TRUE)

  # scale data before PCA
  scale.data <- apply(normdata, 1, scale)

  # get PCA embedding
  set.seed(seed)
  pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  # loading_matrix <- data.frame(pr.data$rotation, features = features)
  pr.data <- pr.data$x
  colnames(pr.data) <- paste0("PC", 1:dims)
  rownames(pr.data) <- colnames(normdata)
  # rownames(pr.data) <- vrSpatialPoints(object_subset, assay = assay)

  # set Embeddings
  vrEmbeddings(object, type = "pca", overwrite = overwrite, ...) <- pr.data

  # return
  return(object)
}

#' getUMAP
#'
#' calculate UMAP of the VoltRon objects
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \code{SampleMetadata(object)}
#' @param data.type the type of data used to calculate UMAP from: "pca" (default), "raw" or "norm"
#' @param dims the number of dimensions extracted from PCA
#' @param umap.key the name of the umap embedding, default: umap
#' @param overwrite Whether the existing embedding with name 'type' should be overwritten in \code{vrEmbeddings}
#' @param seed seed
#' @param ... additional parameters passed to \code{vrEmbeddings}
#'
#' @importFrom uwot umap
#'
#' @export
#'
getUMAP <- function(object, assay = NULL, data.type = "pca", dims = 1:30, umap.key = "umap", overwrite = FALSE, seed = 1, ...){

  # get data
  if(data.type %in% c("raw", "norm")){
    data <- vrData(object, assay = assay, norm = (data.type == "norm"))
    data <- t(data)
  } else{
    embedding_names <- vrEmbeddingNames(object)
    if(data.type %in% vrEmbeddingNames(object)) {
      data <- vrEmbeddings(object, assay = assay, type = data.type, dims = dims)
    } else {
      stop("Please provide a data type from one of three choices: raw, norm and pca")
    }
  }

  # get umap
  set.seed(seed)
  umap_data <- uwot::umap(data)
  colnames(umap_data) <- c("x", "y")
  vrEmbeddings(object, type = umap.key, overwrite = overwrite, ...) <- umap_data

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
  for (i in 1:n_row_tiles) {
    for (j in 1:n_col_tiles) {
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
