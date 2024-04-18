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
#' @importFrom methods slot slotNames
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

#' make_css
#'
#' make_css from \code{tableHTML} package
#'
#' @param ... css style definitions. Each object you provide must be a list of three elements. The first element will be a vector of the selectors to be styled (e.g. table, th, an id or html class). If the first element is a vector of length greater than one then the selectors will be comma separated in the css. The second element will be a vector of the css definitions and the third element will a vector of the values of those definitions.
#' @param file Character sting. If a file name is provided then the css code will be printed into that file. If the argument is NULL (default) then a string will be returned.
#'
#' @importFrom shiny HTML
#'
#' @noRd
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
  # css_string <- htmltools::HTML(paste(all_css, collapse = "\n"))
  css_string <- shiny::HTML(paste(all_css, collapse = "\n"))
  if (is.null(file)) {
    css_string
  }
  else {
    cat(css_string, file = file)
    invisible(NULL)
  }
}

#' Fast creation of dummy variables
#'
#' Quickly create dummy (binary) columns from character and
#' factor type columns in the inputted data (and numeric columns if specified.)
#' This function is useful for statistical analysis when you want binary
#' columns rather than character columns. Adapted from the \code{fastDummies} package (https://jacobkap.github.io/fastDummies/)
#'
#' @param .data
#' An object with the data set you want to make dummy columns from.
#' @param select_columns
#' Vector of column names that you want to create dummy variables from.
#' If NULL (default), uses all character and factor columns.
#' @param remove_first_dummy
#' Removes the first dummy of every variable such that only n-1 dummies remain.
#' This avoids multicollinearity issues in models.
#' @param remove_most_frequent_dummy
#' Removes the most frequently observed category such that only n-1 dummies
#' remain. If there is a tie for most frequent, will remove the first
#' (by alphabetical order) category that is tied for most frequent.
#' @param ignore_na
#' If TRUE, ignores any NA values in the column. If FALSE (default), then it
#' will make a dummy column for value_NA and give a 1 in any row which has a
#' NA value.
#' @param split
#' A string to split a column when multiple categories are in the cell. For
#' example, if a variable is Pets and the rows are "cat", "dog", and "turtle",
#' each of these pets would become its own dummy column. If one row is "cat, dog",
#' then a split value of "," this row would have a value of 1 for both the cat
#' and dog dummy columns.
#' @param remove_selected_columns
#' If TRUE (not default), removes the columns used to generate the dummy columns.
#' @param omit_colname_prefix
#' If TRUE (not default) and `length(select_columns) == 1`, omit pre-pending the
#' name of `select_columns` to the names of the newly generated dummy columns
#'
#' @return
#' A data.frame (or tibble or data.table, depending on input data type) with
#' same number of rows as inputted data and original columns plus the newly
#' created dummy columns.
#'
#' @importFrom data.table as.data.table is.data.table chmatch alloc.col set
#' @importFrom stringr str_sort str_order
#'
dummy_cols <- function(.data, select_columns = NULL, remove_first_dummy = FALSE,
          remove_most_frequent_dummy = FALSE, ignore_na = FALSE, split = NULL,
          remove_selected_columns = FALSE, omit_colname_prefix = FALSE)
{
  stopifnot(is.null(select_columns) || is.character(select_columns),
            select_columns != "", is.logical(remove_first_dummy),
            length(remove_first_dummy) == 1, is.logical(remove_selected_columns))
  if (remove_first_dummy == TRUE & remove_most_frequent_dummy ==
      TRUE) {
    stop("Select either 'remove_first_dummy' or 'remove_most_frequent_dummy'\n         to proceed.")
  }
  if (is.vector(.data)) {
    .data <- data.frame(.data = .data, stringsAsFactors = FALSE)
  }
  data_type <- check_type(.data)
  if (!data.table::is.data.table(.data)) {
    .data <- data.table::as.data.table(.data)
  }
  if (!is.null(select_columns)) {
    char_cols <- select_columns
    cols_not_in_data <- char_cols[!char_cols %in% names(.data)]
    char_cols <- char_cols[!char_cols %in% cols_not_in_data]
    if (length(char_cols) == 0) {
      stop("select_columns is/are not in data. Please check data and spelling.")
    }
  }
  else if (ncol(.data) == 1) {
    char_cols <- names(.data)
  }
  else {
    char_cols <- sapply(.data, class)
    char_cols <- char_cols[char_cols %in% c("factor", "character")]
    char_cols <- names(char_cols)
  }
  if (length(char_cols) == 0 && is.null(select_columns)) {
    stop(paste0("No character or factor columns found. ",
                "Please use select_columns to choose columns."))
  }
  if (!is.null(select_columns) && length(cols_not_in_data) >
      0) {
    warning(paste0("NOTE: The following select_columns input(s) ",
                   "is not a column in data.\n"), paste0(names(cols_not_in_data),
                                                         "\t"))
  }
  for (col_name in char_cols) {
    if (is.factor(.data[[col_name]])) {
      unique_vals <- levels(.data[[col_name]])
      if (any(is.na(.data[[col_name]]))) {
        unique_vals <- c(unique_vals, NA)
      }
    }
    else {
      unique_vals <- unique(.data[[col_name]])
      unique_vals <- stringr::str_sort(unique_vals, na_last = TRUE,
                                       locale = "en_US", numeric = TRUE)
    }
    unique_vals <- as.character(unique_vals)
    if (!is.null(split)) {
      unique_vals <- unique(trimws(unlist(strsplit(unique_vals,
                                                   split = split))))
    }
    if (ignore_na) {
      unique_vals <- unique_vals[!is.na(unique_vals)]
    }
    if (remove_most_frequent_dummy) {
      vals <- as.character(.data[[col_name]])
      vals <- data.frame(sort(table(vals), decreasing = TRUE),
                         stringsAsFactors = FALSE)
      top_vals <- vals[vals$Freq %in% max(vals$Freq), ]
      other_vals <- vals$vals[!vals$Freq %in% max(vals$Freq)]
      other_vals <- as.character(other_vals)
      top_vals <- top_vals[stringr::str_order(top_vals$vals,
                                              na_last = TRUE, locale = "en_US", numeric = TRUE),
      ]
      if (nrow(top_vals) == 1) {
        top_vals <- NULL
      }
      else {
        top_vals <- as.character(top_vals$vals[2:nrow(top_vals)])
      }
      unique_vals <- c(top_vals, other_vals)
      unique_vals <- stringr::str_sort(unique_vals, na_last = TRUE,
                                       locale = "en_US", numeric = TRUE)
    }
    if (remove_first_dummy) {
      unique_vals <- unique_vals[-1]
    }
    data.table::alloc.col(.data, ncol(.data) + length(unique_vals))
    .data[, paste0(col_name, "_", unique_vals)] <- 0L
    for (unique_value in unique_vals) {
      data.table::set(.data, i = which(data.table::chmatch(as.character(.data[[col_name]]),
                                                           unique_value, nomatch = 0) == 1L), j = paste0(col_name,
                                                                                                         "_", unique_value), value = 1L)
      if (!is.na(unique_value)) {
        data.table::set(.data, i = which(is.na(.data[[col_name]])),
                        j = paste0(col_name, "_", unique_value), value = NA)
      }
      if (!is.null(split)) {
        max_split_length <- max(sapply(strsplit(as.character(.data[[col_name]]),
                                                split = split), length))
        for (split_length in 1:max_split_length) {
          data.table::set(.data, i = which(data.table::chmatch(as.character(trimws(sapply(strsplit(as.character(.data[[col_name]]),
                                                                                                   split = split), `[`, split_length))), unique_value,
                                                               nomatch = 0) == 1L), j = paste0(col_name,
                                                                                               "_", unique_value), value = 1L)
        }
        if (is.na(unique_value)) {
          .data[[paste0(col_name, "_", unique_value)]][which(!is.na(.data[[col_name]]))] <- 0
        }
      }
    }
  }
  if (remove_selected_columns) {
    .data <- .data[-which(names(.data) %in% char_cols)]
  }
  .data <- fix_data_type(.data, data_type)
  if (omit_colname_prefix) {
    if (length(select_columns) == 1) {
      new_col_index <- as.logical(rowSums(sapply(unique_vals,
                                                 function(x) grepl(paste0(select_columns, "_",
                                                                          x), names(.data)))))
      names(.data)[new_col_index] <- gsub(paste0(select_columns,
                                                 "_"), "", names(.data)[new_col_index])
    }
    else {
      message("Can't omit the colname prefix when recoding more than one column.")
      message("Returning prefixed dummy columns.")
    }
  }
  return(.data)
}

#' @importFrom data.table is.data.table
#' @noRd
check_type <- function(.data) {
  if (data.table::is.data.table(.data)) {
    data_type <- "is_data_table"
  } else if (inherits(.data, "tbl_df")) {
    data_type <- "is_tibble"
  } else {
    data_type <- "is_data_frame"
  }

  return(data_type)
}

#' @importFrom dplyr as_tibble
#' @noRd
fix_data_type <- function(.data, data_type) {
  if (data_type == "is_data_frame") {
    .data <- as.data.frame(.data, stringsAsFactors = FALSE)
  } else if (data_type == "is_tibble") {
    .data <- dplyr::as_tibble(.data)
  }

  return(.data)
}