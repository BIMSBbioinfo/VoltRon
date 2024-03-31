####
# ROI DE Analysis ####
####

#' getDiffExp
#'
#' Get differential expression with DESeq2
#'
#' @param object a VoltRon object
#' @param assay assay name (exp: Assay1) or assay class (exp: Visium, Xenium), see \link{SampleMetadata}. 
#' if NULL, the default assay will be used, see \link{vrMainAssay}.
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
  if(group.by %in% colnames(metadata)){
    if(!is.null(group.base)){
      if(!group.base %in% unique(metadata[[group.by]]))
        stop("Please specify a group that is included in the group.by column of metadata to define the base group!")
    } else{
      group.base <- levels(factor(metadata[[group.by]]))[1]
    }
  } else {
    stop("Column ", group.by, " cannot be found in metadata!")
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
#' @importFrom stats as.formula
#' @importFrom S4Vectors DataFrame
#'
#' @noRd
getDiffExpDESeq2 <- function(data, metadata, group.by, group.base = NULL, covariates){

  if (!requireNamespace('DESeq2'))
    stop("Please install DESeq2 package to find differentially expressed genes with DESeq2 method")

  # experimental design for deseq2

  # make formula
  if(is.null(covariates)){
    group.by.data <- metadata[[group.by]]
    uniq_groups <- unique(group.by.data)
    if(!is.null(group.base))
      uniq_groups <- c(group.base, uniq_groups[!uniq_groups %in% group.base])
    group.by.data <- factor(group.by.data, levels = uniq_groups)
    colData <- S4Vectors::DataFrame(group.by.data)
    colnames(colData) <- c(group.by, covariates)
    deseq2.formula <- stats::as.formula(paste0("~", group.by))
  } else {
    if(all(covariates %in% colnames(metadata))){
      design.data <- metadata[,c(group.by, covariates)]
      uniq_groups <- unique(design.data[[group.by]])
      if(is.null(group.base))
        uniq_groups <- c(group.base, uniq_groups[!uniq_groups %in% group.base])
      group.by.data <- factor(design.data[[group.by]], levels = uniq_groups)
      design.data[[group.by]] <- group.by.data
      colData <- S4Vectors::DataFrame(design.data)
      colnames(colData) <- c(group.by, covariates)
      deseq2.formula <- stats::as.formula(paste0("~", group.by, " + ", paste(covariates, collapse = " + ")))
    } else {
      stop("Columns ", paste(covariates, collapse = ", "), " cannot be found in metadata!")
    }
  }

  # run DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = colData, design = deseq2.formula)
  dds <- DESeq2::DESeq(dds)

  all_results <- NULL
  for(i in 1:(length(uniq_groups)-1)){
    for(j in (i+1):length(uniq_groups)){
      comparison <- c(group.by, uniq_groups[i], uniq_groups[j])
      cur_results <- as.data.frame(DESeq2::results(dds, comparison))
      cur_results <- data.frame(cur_results, gene = rownames(cur_results), comparison = paste(comparison, collapse = "_"))
      all_results <- rbind(all_results, cur_results)
    }
  }

  # return
  return(all_results)
}
