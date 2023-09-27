#' @include generics.R
#'
NULL

####
# Normalization ####
####

#' @param assay assay
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

#' @param assay assay
#' @param method the normalization method: "LogNorm", "Q3Norm", "LogQ3Norm" or "CLR"
#' @param desiredQuantile the quantile of the data if "QuanNorm" or "LogQuanNorm" is selected as \code{method}
#'
#' @rdname normalizeData
#' @method normalizeData vrAssay
#'
#' @importFrom stats quantile
#'
#' @export
#'
normalizeData.vrAssay <- function(object, method = "LogNorm", desiredQuantile = 0.9) {

  # size factor
  rawdata <- vrData(object)
  sizefactor <- colSums(rawdata)

  if(!is.numeric(desiredQuantile)){
    stop("desiredQuantile should be numeric")
  } else {
    if(!findInterval(desiredQuantile, c(0,1)) == 1L){
      stop("desiredQuantile should be between [0,1]")
    }
  }

  # normalization method
  if(method == "LogNorm"){
    sizefactor <- matrix(rep(sizefactor, nrow(rawdata)), byrow = T, nrow = nrow(rawdata))
    normdata <- (rawdata/sizefactor)*10000
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

#' @param assay assay
#'
#' @rdname getFeatures
#' @method getFeatures VoltRon
#'
#' @export
#'
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

#' @param max.count max count (maximum across spatial points) for low count filtering
#' @param n the number of features
#'
#' @rdname getFeatures
#' @method getFeatures vrAssay
#'
#' @importFrom stats loess predict var
#'
#' @export
getFeatures.vrAssay <- function(object, max.count = 1, n = 3000){

  # get data and coordinates
  normdata <- vrData(object, norm = TRUE)
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
#' @param object A Voltron Object
#' @param assay assay
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

#' @param assay assay
#' @param features the selected features for PCA reduction
#' @param dims the number of dimensions extracted from PCA
#' @param seed seed
#'
#' @rdname getPCA
#' @method getPCA VoltRon
#'
#' @importFrom irlba irlba
#' @importFrom dplyr left_join
#'
#' @export
#'
getPCA.VoltRon <- function(object, assay = NULL, features = NULL, dims = 30, seed = 1){

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get shared features and subset
  if(is.null(features)){
    features <- getVariableFeatures(object, assay = assay)
  }
  object_subset <- subset(object, features = features)

  # adjust extraction features length
  if(dims > length(features)){
    message("Requested more PC dimensions than existing features: dims = length(features) now!")
    dims <- length(features)
  }

  # get data
  normdata <- vrData(object_subset, norm = TRUE)

  # scale data before PCA
  scale.data <- apply(normdata, 1, scale)

  # get PCA embedding
  set.seed(seed)
  pr.data <- irlba::prcomp_irlba(scale.data, n=dims, center=colMeans(scale.data))
  loading_matrix <- data.frame(pr.data$rotation, features = features)
  pr.data <- pr.data$x
  colnames(pr.data) <- paste0("PC", 1:dims)
  rownames(pr.data) <- colnames(normdata)

  # update feature matrix
  feature_data <- vrFeatureData(object)
  if(nrow(feature_data) > 0){
    feature_data <- feature_data[, !grepl("PC[0-9]+", colnames(feature_data))]
    feature_data <- data.frame(feature_data, features = vrFeatures(object))
    feature_data <- feature_data %>% left_join(loading_matrix)
    vrFeatureData(object) <- data.frame(feature_data[,colnames(feature_data)[!colnames(feature_data) %in% "features"]],
                                        row.names = feature_data$features)
  }

  # set Embeddings
  vrEmbeddings(object, type = "pca") <- pr.data

  # return
  return(object)
}

#' @param assay assay
#' @param data.type the type of data used to calculate UMAP from: "pca" (default), "raw" or "norm"
#' @param dims the number of dimensions extracted from PCA
#' @param seed seed
#'
#' @rdname getUMAP
#' @method getUMAP VoltRon
#'
#' @importFrom umap umap
#'
#' @export
#'
getUMAP.VoltRon <- function(object, assay = NULL, data.type = "pca", dims = 1:30, seed = 1){

  # get data
  if(data.type %in% c("raw", "norm")){
    data <- vrData(object, assay = assay, norm = (data.type == "norm"))
    data <- t(data)
  } else if(data.type == "pca") {
    data <- vrEmbeddings(object, assay = assay, type = data.type, dims = dims)
  } else {
    stop("Please provide a data type from one of three choices: raw, norm and pca")
  }

  # get umap
  set.seed(seed)
  umap_data <- umap::umap(data, preserve.seed = seed)
  umap_data <- umap_data$layout
  colnames(umap_data) <- c("x", "y")
  vrEmbeddings(object, type = "umap") <- umap_data

  # return
  return(object)
}
