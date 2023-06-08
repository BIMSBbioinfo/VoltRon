#' @include generics.R
#'
NULL

####
# Normalization ####
####

#' @rdname normalizeData
#' @concept preprocessing
#' @method normalizeData VoltRon
#'
#' @export
normalizeData.VoltRon <- function(object, assay = NULL, ...) {

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # normalize assays
  for(assy in assay_names){
    cur_assay <- object[[assy]]
    object[[assy]] <- normalizeData(cur_assay, ...)
  }

  # return
  return(object)
}

#' @rdname normalizeData
#' @concept preprocessing
#' @method normalizeData vrAssay
#'
#' @export
normalizeData.vrAssay <- function(object, method = "LogNorm", desiredQuantile = 0.9) {

  # size factor
  rawdata <- object@rawdata
  sizefactor <- colSums(rawdata)

  # normalization method
  if(method == "LogNorm"){
    sizefactor <- matrix(rep(sizefactor, nrow(rawdata)), byrow = T, nrow = nrow(rawdata))
    normdata <- (rawdata/sizefactor)*10000
    normdata <- log(normdata + 1)
  } else if(method == "QuanNorm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
  } else if(method == "LogQuanNorm") {
    rawdata[rawdata==0] <- 1
    qs <- apply(rawdata, 2, function(x) stats::quantile(x, desiredQuantile))
    normdata <- sweep(rawdata, 2L, qs / exp(mean(log(qs))), FUN = "/")
    normdata <- log(normdata + 1)
  } else if(method == "CLR") {
    normdata <- apply(rawdata, 2, function(x) {
      log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
    })
  }

  # get normalized data
  object@normdata <- normdata

  # return
  return(object)
}

####
# Features ####
####

#' @rdname getSpatialFeatures
#' @concept preprocessing
#' @method getSpatialFeatures VoltRon
#'
#' @export
getSpatialFeatures.VoltRon <- function(object, assay = NULL, ...){

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get features for all coordinates
  for(assy in assay_names){
    object[[assy]] <- getSpatialFeatures(object[[assy]], ...)
  }

  # return
  return(object)
}

#' @rdname getSpatialFeatures
#' @concept preprocessing
#' @method getSpatialFeatures vrAssay
#'
#' @export
getSpatialFeatures.vrAssay <- function(object, max.count = 1, n = 3000){

  # get data and coordinates
  normdata <- Data(object, norm = TRUE)
  rawdata <- Data(object, norm = FALSE)
  coords <- Coordinates(object)
  features <- vrFeatures(object)

  # eliminate genes with low counts
  keep.genes <- which(apply(rawdata,1,max) > max.count)
  rawdata_subset <- rawdata[keep.genes,]

  # vst estimation
  vst_data <- data.frame(mean = rowMeans(rawdata), var = apply(rawdata, 1, var))
  loess_data <- vst_data[keep.genes,]
  loess_results <- loess(var~mean, loess_data, span = 0.3)
  vst_data$adj_var <- 0
  vst_data$rank <- 0
  vst_data[keep.genes,]$adj_var <- predict(loess_results)
  vst_data[keep.genes,]$rank <- order(order(vst_data$adj_var[keep.genes], decreasing = TRUE))

  # set feature data
  FeatureData(object) <- vst_data

  # return
  return(object)
}


#' getSharedFeatures
#'
#' get shared variable features across multiple assays
#'
#' @importFrom dplyr full_join
#'
getSharedFeatures <- function(object, assay = NULL, n = 3000, ...){

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get features for all coordinates
  ranks <- NULL
  for(assy in assay_names){
    feature_data <- FeatureData(object[[assy]], ...)
    feature_data$gene <- rownames(feature_data)
    if(is.null(ranks)){
      ranks <- feature_data[,c("gene", "rank")]
    } else {
      ranks <- ranks %>% full_join(feature_data[,c("gene", "rank")], by = c("gene" = "gene"))
    }
  }

  # get geometric mean of ranks, i.e. rank product statistic
  ranks <- ranks[,!colnames(ranks) %in% "gene", drop = FALSE]
  ranks <- apply(ranks, 1, function(x) exp(mean(log(x))))
  ranks <- ranks[ranks != 0]

  # get selected features
  selected_features <- names(head(sort(ranks, decreasing = FALSE), n))

  # return
  return(selected_features)
}

####
# Embeddings ####
####

#' @rdname PCA
#' @concept embedding
#' @method PCA VoltRon
#'
#' @export
PCA.VoltRon <- function(object, assay = NULL, dims = 30){

  # get assay names
  assay_names <- AssayNames(object, assay = assay)

  # get shared features and subset
  features <- getSharedFeatures(object, assay = assay)
  object <- subset(object, features = features)

  # get data
  normdata <- Data(object, norm = TRUE)

  # scale data before PCA
  scale.data <- apply(normdata, 1, scale)

  # get PCA embedding
  pr.data <- irlba(scale.data, nv=dims, center=colMeans(scale.data))
  pr.data <- pr.data$u
  colnames(pr.data) <- paste0("PC", 1:dims)
  rownames(pr.data) <- colnames(normdata)

  # set Embeddings
  Embeddings(object, type = "pca") <- pr.data

  # return
  return(object)
}

#' @rdname UMAP
#' @concept embedding
#' @method UMAP VoltRon
#'
#' @export
UMAP.VoltRon <- function(object, assay = NULL, dims = 30, seed = 1){

  # set Embeddings
  embedding_data <- Embeddings(object, assay = assay)

  # get umap
  umap_data <- umap::umap(embedding_data, preserve.seed = seed)
  umap_data <- umap_data$layout
  colnames(umap_data) <- c("x", "y")
  Embeddings(object, type = "umap") <- umap_data

  # return
  return(object)
}
