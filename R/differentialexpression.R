####
# ROI DE Analysis ####
####

#' getDiffExp
#'
#' Get differential expression with DESeq2
#'
#' @param object a VoltRon object
#' @param assay assay
#' @param group.by the categorical variable from metadata to get differentially expressed features across
#' @param group.base Optional, the base category in \code{group.by} which is used as control group
#' @param covariates the covariate variable for the design formula
#' @param method the method for DE analysis, e.g. DESeq2
#'
#' @export
#'
getDiffExp <- function(object, assay = NULL, group.by, group.base = NULL, covariates = NULL, method = "DESeq2"){

  # get data and metadata
  data <- vrData(object, assay = assay)
  metadata <- Metadata(object, assay = assay)

  # check groups
  if(!is.null(group.base)){
    if(!group.base %in% unique(metadata[[group.by]]))
      stop("Please specify a group that is included in the group.by column of metadata to define the base group!")
  } else{
    group.base <- levels(factor(metadata[[group.by]]))[1]
  }

  # select a method
  results <-
    switch(method,
           DESeq2 = {
             getDiffExpDESeq2(data, metadata, group.by = group.by, group.base = group.base, covariates = covariates)
           })

  # return
  return(results)
}

#' getDiffExpDESeq2
#'
#' @param data the raw data set
#' @param metadata the metadata
#' @param group.by the categorical variable from metadata to get differentially expressed features across
#' @param group.base Optional, the base category in \code{group.by} which is used as control group
#' @param covariates the covariate variable for the design formula
#'
getDiffExpDESeq2 <- function(data, metadata, group.by, group.base = NULL, covariates){

  if (!requireNamespace('DESeq2'))
    stop("Please install DESeq2 package to find differentially expressed genes with DESeq2 method")

  # experimental design for deseq2

  # make formula
  if(is.null(covariates)){
    group.by.data <- metadata[[group.by]]
    uniq_groups <- unique(group.by.data)
    if(is.null(group.base))
      uniq_groups <- c(group.base, uniq_groups[!uniq_groups %in% group.base])
    group.by.data <- factor(group.by.data, levels = uniq_groups)
    colData <- S4Vectors::DataFrame(group.by.data)
    colnames(colData) <- c(group.by, covariates)
    deseq2.formula <- as.formula(paste0("~", group.by))
  } else {
    design.data <- metadata[,c(group.by, covariates)]
    uniq_groups <- unique(design.data[[group.by]])
    if(is.null(group.base))
      uniq_groups <- c(group.base, uniq_groups[!uniq_groups %in% group.base])
    group.by.data <- factor(design.data[[group.by]], levels = uniq_groups)
    design.data[[group.by]] <- group.by.data
    colData <- S4Vectors::DataFrame(design.data)
    colnames(colData) <- c(group.by, covariates)
    deseq2.formula <- as.formula(paste0("~", group.by, " + ", paste(covariates, collapse = " + ")))
  }

  # run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = deseq2.formula)
  dds <- DESeq(dds)

  all_results <- NULL
  for(i in 1:(length(uniq_groups)-1)){
    for(j in (i+1):length(uniq_groups)){
      comparison <- c(group.by, uniq_groups[i], uniq_groups[j])
      cur_results <- as.data.frame(results(dds, comparison))
      cur_results <- data.frame(cur_results, gene = rownames(cur_results), comparison = paste(comparison, collapse = "_"))
      all_results <- rbind(all_results, cur_results)
    }
  }

  # return
  return(all_results)
}
