####
# Nanostring Auxiliary tools ####
####

generate_pkc_lookup <- function (jsons_vec)
{
  lookup_df <- data.frame(RTS_ID = character(), Target = character(),
                          Module = character(), CodeClass = character(), ProbeID = character(),
                          GeneID = character(), SystematicName = character(), stringsAsFactors = FALSE)
  for (curr_idx in seq_len(length(jsons_vec))) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    for (targ in curr_json[["Targets"]]) {
      curr_targ <- targ[["DisplayName"]]
      curr_code_class <- gsub("\\d+$", "", targ[["CodeClass"]])
      for (prb in targ[["Probes"]]) {
        if (curr_json[["AnalyteType"]] == "Protein") {
          curr_RTS_ID <- targ$RTS_ID
        }
        else {
          curr_RTS_ID <- prb$RTS_ID
        }
        curr_probe_ID <- prb$ProbeID
        curr_gene_ID <- paste(prb$GeneID, collapse = ", ")
        if (length(prb$GeneID) < 1) {
          curr_gene_ID <- NA
        }
        curr_syst_name <- paste(prb$SystematicName, collapse = ", ")
        lookup_df[nrow(lookup_df) + 1, ] <- list(curr_RTS_ID,
                                                 curr_targ, curr_module, curr_code_class, curr_probe_ID,
                                                 curr_gene_ID, curr_syst_name)
      }
    }
  }
  return(lookup_df)
}

.dccMetadata <-
  list(schema =
         list("Header" =
                data.frame(labelDescription =
                             c("The version of the file",
                               "The version of the software used to create the file",
                               "The date of the sample"),
                           minVersion = numeric_version(c("0.01", "0.01", "0.01")),
                           row.names =
                             c("FileVersion", "SoftwareVersion", "Date"),
                           stringsAsFactors = FALSE),
              "Scan_Attributes" =
                data.frame(labelDescription =
                             c("The sample ID",
                               "The plate ID",
                               "The well ID"),
                           row.names =
                             c("ID", "Plate_ID", "Well"),
                           minVersion = numeric_version(c(rep("0.01", 3L))),
                           stringsAsFactors = FALSE),
              "NGS_Processing_Attributes" =
                data.frame(labelDescription =
                             c(NA_character_,
                               NA_integer_,
                               NA_integer_,
                               NA_integer_,
                               NA_integer_,
                               NA_real_,
                               NA_real_),
                           minVersion = numeric_version(c(rep("0.01", 7L))),
                           row.names =
                             c("SeqSetId", "Raw", "Trimmed",
                               "Stitched", "Aligned", "umiQ30", "rtsQ30"),
                           stringsAsFactors = FALSE),
              "Code_Summary" =
                data.frame(labelDescription =
                             c(NA_character_, NA_integer_),
                           minVersion = numeric_version(c(rep("0.01", 2L))),
                           row.names = c("RTS_ID", "Count"),
                           stringsAsFactors = FALSE)
         )
  )


.dccMetadata[["protocolData"]] <-
  do.call(rbind,
          unname(head(.dccMetadata[["schema"]], 3L)))[, "labelDescription",
                                                      drop = FALSE]

rownames(.dccMetadata[["protocolData"]])[rownames(.dccMetadata[["protocolData"]]) == "ID"] <- "SampleID"


.codeClassMetadata <-
  c("CodeClass,IsControl,Analyte",
    "Endogenous,FALSE,gx|cnv|fusion",
    "Housekeeping,TRUE,gx|fusion",
    "Positive,TRUE,general",
    "Negative,TRUE,general",
    "Binding,TRUE,general",
    "Purification,TRUE,general",
    "Reserved,TRUE,general",
    "SNV_INPUT_CTL,TRUE,SNV",
    "SNV_NEG,TRUE,SNV",
    "SNV_POS,TRUE,SNV",
    "SNV_UDG_CTL,TRUE,SNV",
    "SNV_PCR_CTL,TRUE,SNV",
    "SNV_REF,FALSE,SNV",
    "SNV_VAR,FALSE,SNV",
    "PROTEIN,FALSE,protein",
    "PROTEIN_NEG,TRUE,protein",
    "PROTEIN_CELL_NORM,TRUE,protein",
    "Restriction Site,TRUE,CNV",
    "Invariant,TRUE,CNV")
.codeClassMetadata <-
  utils::read.csv(textConnection(paste0(.codeClassMetadata, collapse = "\n")),
                  colClasses = c("character", "logical", "character"),
                  stringsAsFactors = FALSE)


.validDccSchema <-
  function(x, fileVersion,
           section = c("Header", "Scan_Attributes", "NGS_Processing_Attributes", "Code_Summary"))
  {
    section <- match.arg(section)
    schema <- .dccMetadata[["schema"]][[section]]
    expectedNames <- row.names(schema)[schema[,"minVersion"] <= fileVersion]
    if (all(expectedNames %in% colnames(x))) {
      TRUE
    } else {
      sprintf("<%s> section must contain %s", section,
              paste0("\"", expectedNames, "\"", collapse = ", "))
    }
  }


.allNA <- function(x) {
  all(is.na(x))
}

.allTRUE <- function(x) {
  is.logical(x) && !anyNA(x) && all(x)
}

.allFALSE <- function(x) {
  is.logical(x) && !anyNA(x) && !any(x)
}

.allZero <- function(x) {
  is.numeric(x) && !anyNA(x) && identical(range(x), c(0, 0))
}

.validNonNegativeInteger <- function(x) {
  is.integer(x) && !anyNA(x) && min(x) >= 0L
}

.validNonNegativeNumber <- function(x) {
  is.numeric(x) && !anyNA(x) && min(x) >= 0
}

.validPositiveNumber <- function(x) {
  is.numeric(x) && !anyNA(x) && min(x) > 0
}

####
# Other Auxiliary tools ####
####

fill.na <- function(x, i = 5) {
  if (is.na(x)[i]) {
    return(round(mean(x, na.rm = TRUE), 0))
  }
  else {
    return(round(x[i], 0))
  }
}

#' slotApply
#'
#' apply to slots
#'
#' @param x object
#' @param FUN function
#' @param ... arguments passed to \code{FUN}
#'
#' @importFrom methods slot
#'
slotApply <- function(x,FUN,...){
  cl <- class(x)
  result <- list()
  for(i in methods::slotNames(cl)){
    result[[i]] <- FUN(methods::slot(x,i),...)
  }
  result
}

#' slotToList
#'
#' slot to list
#'
#' @param x object
#'
#' @importFrom methods slot slotNames
#'
slotToList <- function(x){
  returnlist <- list()
  namesslot <- methods::slotNames(x)
  for(cur_slot in namesslot)
    returnlist[[cur_slot]] <- methods::slot(x, name = cur_slot)
  returnlist
}

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
}

jaccard_similarity <- function(mat) {
  matinv <- 1 - mat
  return((mat %*% mat)/(nrow(mat) - (matinv %*% matinv)))
}

#' make_css
#'
#' make_css from \code{tableHTML} package
#'
#' @param ... css style definitions. Each object you provide must be a list of three elements. The first element will be a vector of the selectors to be styled (e.g. table, th, an id or html class). If the first element is a vector of length greater than one then the selectors will be comma separated in the css. The second element will be a vector of the css definitions and the third element will a vector of the values of those definitions.
#' @param file Character sting. If a file name is provided then the css code will be printed into that file. If the argument is NULL (default) then a string will be returned.
#'
#' @importFrom htmltools HTML
#'
#' @noRd
#'
make_css <- function (..., file = NULL)
{
  css_defs <- list(...)
  for (x in css_defs) {
    if ((!is.list(x)) | (length(x) != 3L)) {
      stop("Each element in ... needs to be a list of three elements")
    }
    if (length(x[[2]]) != length(x[[3]])) {
      stop("The second and third elements of each list need to have the same length")
    }
  }
  all_css <- vapply(css_defs, function(x) {
    css_comp <- paste0(x[[2]], ": ", x[[3]], ";")
    style <- paste(css_comp, collapse = "\n  ")
    to_be_styled <- paste(x[[1]], collapse = ",\n")
    paste0(to_be_styled, " {\n  ", style, "\n}\n")
  }, FUN.VALUE = character(1))
  css_string <- htmltools::HTML(paste(all_css, collapse = "\n"))
  if (is.null(file)) {
    css_string
  }
  else {
    cat(css_string, file = file)
    invisible(NULL)
  }
}
