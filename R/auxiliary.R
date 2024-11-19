####
# Matrix Operations ####
####

getMin <- function(data, ...){
  if(inherits(data, "IterableMatrix")){
    data <- as(data, "dgCMatrix")
  } 
  return(min(data, ...))
}

getMax <- function(data, ...){
  if(inherits(data, "IterableMatrix")){
    data <- as(data, "dgCMatrix")
  } 
  return(max(data, ...))
}

getRange <- function(data, ...){
  return(c(getMin(data, ...), getMax(data, ...)))
}


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
# Basilisk Environment ####
####

#' The Python Basilisk environment
#'
#' Defines a conda environment via Basilisk, which is used to convert R objects to Zarr stores.
#'
#' @importFrom basilisk BasiliskEnvironment
#'
#' @keywords internal
#'
#' @noRd
py_env <- basilisk::BasiliskEnvironment(
  envname="VoltRon_basilisk_env",
  pkgname="VoltRon",
  packages=c(
    "numpy==1.*",
    "pandas==1.*",
    "anndata==0.7.*",
    "h5py==3.*",
    "hdf5==1.*",
    "natsort==7.*",
    "packaging==20.*",
    "scipy==1.*",
    "sqlite==3.*",
    "zarr==2.*",
    "numcodecs==0.*",
    "tifffile==2024.2.12"
  ),
  pip=c(
    "ome-zarr==0.2.1"
  )
)

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

#' CSS string helper
#'
#' Convenience function for building CSS style declarations (i.e. the string
#' that goes into a style attribute, or the parts that go inside curly braces in
#' a full stylesheet).
#'
#' CSS uses `'-'` (minus) as a separator character in property names, but
#' this is an inconvenient character to use in an R function argument name.
#' Instead, you can use `'.'` (period) and/or `'_'` (underscore) as
#' separator characters. For example, `css(font.size = "12px")` yields
#' `"font-size:12px;"`.
#'
#' To mark a property as `!important`, add a `'!'` character to the end
#' of the property name. (Since `'!'` is not normally a character that can be
#' used in an identifier in R, you'll need to put the name in double quotes or
#' backticks.)
#'
#' Argument values will be converted to strings using
#' `paste(collapse = " ")`. Any property with a value of `NULL` or
#' `""` (after paste) will be dropped.
#'
#' @param ... Named style properties, where the name is the property name and
#'   the argument is the property value. See Details for conversion rules.
#' @param collapse_ (Note that the parameter name has a trailing underscore
#'   character.) Character to use to collapse properties into a single string;
#'   likely `""` (the default) for style attributes, and either `"\n"`
#'   or `NULL` for style blocks.
#'
#' @importFrom rlang dots_list
#' 
#' @noRd
css <- function(..., collapse_ = "") {
  props <- rlang::dots_list(...)
  if (length(props) == 0) {
    return(NULL)
  }
  
  if (is.null(names(props)) || any(names(props) == "")) {
    stop("cssList expects all arguments to be named")
  }
  
  # Necessary to make factors show up as level names, not numbers
  props[] <- lapply(props, paste, collapse = " ")
  
  # Drop null args
  props <- props[!sapply(props, empty)]
  if (length(props) == 0) {
    return(NULL)
  }
  
  # Translate camelCase, snake_case, and dot.case to kebab-case
  # For standard CSS properties only, not CSS variables
  is_css_var <- grepl("^--", names(props))
  names(props)[!is_css_var] <- standardize_property_names(names(props)[!is_css_var])
  
  # Create "!important" suffix for each property whose name ends with !, then
  # remove the ! from the property name
  important <- ifelse(grepl("!$", names(props), perl = TRUE), " !important", "")
  names(props) <- sub("!$", "", names(props), perl = TRUE)
  
  paste0(names(props), ":", props, important, ";", collapse = collapse_)
}

empty <- function(x) {
  length(x) == 0 || (is.character(x) && !any(nzchar(x)))
}

standardize_property_names <- function(x) {
  # camelCase to kebab-case
  x <- gsub("([A-Z])", "-\\1", x)
  x <- tolower(x)
  # snake_case and dot.case to kebab-case
  gsub("[._]", "-", x)
}

#' @importFrom grDevices hcl
hue_pal <- function(n, h = c(0, 360) + 15, c = 100, l = 65, h.start = 0, direction = 1) 
{
  if (n == 0) {
    cli::cli_abort("Must request at least one colour from a hue palette.")
  }
  if ((diff(h)%%360) < 1) {
    h[2] <- h[2] - 360/n
  }
  hues <- seq(h[1], h[2], length.out = n)
  hues <- (hues + h.start)%%360
  hcl <- cbind(hues, c, l)
  # pal <- farver::encode_colour(hcl, from = "hcl")
  pal <- apply(hcl, 1, function(x){
    grDevices::hcl(x[1], x[2], x[3])
  })
  if (direction == -1) {
    rev(pal)
  }
  else {
    pal
  }
}

rescale_numeric <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), ...) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

as.raster_array <- function (x, max = 1, ...) 
{
  if (!is.numeric(x)) {
    if (is.raw(x)) {
      storage.mode(x) <- "integer"
      max <- 255L
    }
    else stop("a raster array must be numeric")
  }
  if (length(d <- dim(x)) != 3L) 
    stop("a raster array must have exactly 3 dimensions")
  r <- array(if (d[3L] == 3L) 
    rgb(t(x[, , 1L]), t(x[, , 2L]), t(x[, , 3L]), maxColorValue = max)
    else if (d[3L] == 4L) 
      rgb(t(x[, , 1L]), t(x[, , 2L]), t(x[, , 3L]), t(x[, , 4L]), maxColorValue = max)
    else if (d[3L] == 1L) 
      rgb(t(x[, , 1L]), t(x[, , 1L]), t(x[, , 1L]), maxColorValue = max)
    else stop("a raster array must have exactly 1, 3 or 4 planes"), 
    dim = d[1:2])
  class(r) <- "raster"
  r
}

avgHexColor <- function(colors){
  rgb(t(Reduce(`+`, lapply(colors, col2rgb))/3), maxColorValue=255)
}

####
## ggedit tools ####
## See https://github.com/yonicd/ggedit/tree/master/R
####

#' @title ggplot2 layer proto extraction
#' @description Extract geom, stat and position protos from a ggplot2 layer
#' @param l ggproto
#' @noRd
proto_features <- function(l) {
  a <- sapply(c("position", "geom", "stat"), function(x) {
    class(l[[x]])[1]
  })
  
  data.frame(t(a), stringsAsFactors = FALSE)
}

# forked from https://github.com/yihui/knitr/blob/master/R/defaults.R
#' @importFrom stats setNames
#' @noRd
new_defaults <- function(value = list()) {
  defaults <- value
  
  get <- function(name, default = FALSE, drop = TRUE, regex=FALSE, ...) {
    if (default) defaults <- value # this is only a local version
    if (missing(name)) {
      defaults
    } else {
      if (drop && length(name) == 1) {
        if (regex) {
          name_grep <- grep(name, names(defaults), value = TRUE, ...)
          stats::setNames(defaults[name_grep], name_grep)
        } else {
          defaults[[name]]
        }
      } else {
        stats::setNames(defaults[name], name)
      }
    }
  }
  
  set <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return()
    if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]])) {
      if (length(dots <- dots[[1]]) == 0) {
        return()
      }
    }
    defaults <<- merge(dots)
    invisible(NULL)
  }
  
  # merge <- function(values) merge_list(defaults, values)
  
  restore <- function(target = value) defaults <<- target
  
  append <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return()
    if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]])) {
      if (length(dots <- dots[[1]]) == 0) {
        return()
      }
    }
    dots <- sapply(names(dots), function(x) dots[[x]] <- c(defaults[[x]], dots[[x]]), simplify = FALSE)
    defaults <<- merge(dots)
    invisible(NULL)
  }
  
  list(get = get, set = set, append = append, merge = merge, restore = restore)
}

#' @title Creates an independent copy of a ggplot layer object
#' @description Creates copies of ggplot layers from within ggplot objects that
#' are independent of the parent object.
#' @details ggplot objects are comprimsed of layer objects. Once compiled they
#' are part of the plot object environment and if they are changed internally
#' regardless of where they are in the (ie different environment) it will change
#' the original plot. This function allows to create replicates of the plot layers
#' and edit them independent of the original plot. When setting verbose to TRUE
#' function returns the ggplot2 call as a string to paste in regular ggplot script
#' to generate the layer.
#' @param l ggplot2 object layer
#' @param verbose toggle to control if the output is ggproto object (verbose==FALSE,default) or string of layer call (verbose==TRUE)
#' @param showDefaults toggle to control if the verbose output shows all the input arguments passed to the proto object (if verbose==FALSE then ignored)
#' @return ggproto or string object (conditional on verbose)
#'
#' @importFrom utils capture.output
#' @importFrom rlang sym '!!'
#' @noRd
cloneLayer <- function (l, verbose = FALSE, showDefaults = TRUE) 
{
  geom_opts <- ggedit_opts$get("session_geoms")
  parent.layer <- dplyr::left_join(proto_features(l), dplyr::filter(geom_opts, 
                                                                    !grepl("^stat", !!rlang::sym("fn"))), by = c("position", 
                                                                                                                 "geom", "stat"))
  if (is.na(parent.layer$fn)) 
    parent.layer$fn <- paste0(tolower(strsplit(parent.layer$stat, 
                                               "(?<=Stat)", perl = TRUE)[[1]]), collapse = "_")
  layer.names <- c("mapping", "data", "geom", "position", "stat", 
                   "show.legend", "inherit.aes", "aes_params", "geom_params", 
                   "stat_params")
  x <- sapply(layer.names, function(y) {
    b <- l[[y]]
    if ("waiver" %in% class(b)) 
      b <- NULL
    if (y == "geom") 
      b <- eval(parse(text = parent.layer$geom))
    if (y == "position") 
      b <- gsub(y, "", tolower(class(b)[1]))
    if (y == "stat") 
      b <- eval(parse(text = parent.layer$stat))
    b
  })
  x$params <- append(x$stat_params, x$geom_params)
  x$params <- append(x$params, x$aes_params)
  x$params <- x$params[!duplicated(names(x$params))]
  x$geom_params <- x$aes_params <- x$stat_params <- NULL
  if (verbose) {
    nm <- names(x)
    nm <- nm[!sapply(x, typeof) %in% c("environment", "closure", 
                                       "list")]
    geom_aes <- list(geom = parent.layer$fn, mapping = sapply(names(x$mapping), 
                                                              build_map, y = x$mapping), params = sapply(names(x$params), 
                                                                                                         build_map, y = x$params), layer = sapply(rev(nm), 
                                                                                                                                                  build_map, y = x[rev(nm)]), data = paste0("data = ", 
                                                                                                                                                                                            paste0(capture.output(dput(x$data)), collapse = "\n")))
    strRet <- sprintf("%s(mapping=aes(%s),%s,%s)", paste0(geom_aes$geom, 
                                                          collapse = ","), paste0(geom_aes$mapping, collapse = ","), 
                      paste0(geom_aes$params, collapse = ","), paste0(geom_aes$layer, 
                                                                      collapse = ","))
    if (!showDefaults) {
      geom_proto <- cloneProto(eval(parse(text = paste0(geom_aes$geom, 
                                                        "()"))))
      geom_diff <- sapply(names(geom_aes)[-1], function(x) geom_aes[[x]][!geom_aes[[x]] %in% 
                                                                           geom_proto[[x]]])
      strRet <- sprintf("%s(aes(%s),%s,%s,%s)", paste0(geom_aes$geom, 
                                                       collapse = ","), paste0(geom_diff$mapping, collapse = ","), 
                        paste0(geom_diff$params, collapse = ","), paste0(geom_diff$layer, 
                                                                         collapse = ","), geom_aes$data)
    }
    strRet <- gsub("aes()", "", strRet, fixed = T)
    strRet <- gsub("[,]{2,}", ",", strRet)
    strRet <- gsub("data=NULL", "", strRet)
    strRet <- gsub(",)", ")", strRet)
    strRet <- gsub("\\(,", "(", strRet)
    strRet
  }
  else {
    do.call(layer, x)
  }
}

#' @import dplyr 
#' @importFrom rlang sym '!!'
#' @noRd
cloneProto <- function(l) {
  
  geom_opts <- ggedit_opts$get("session_geoms")
  
  parent.layer <- proto_features(l)  |> 
    dplyr::left_join(
      geom_opts  |>  dplyr::filter(!grepl("^stat", !!rlang::sym('fn'))),
      by = c("position", "geom", "stat")
    )
  
  if (is.na(parent.layer$fn)) {
    parent.layer$fn <- paste0(tolower(strsplit(parent.layer$stat, "(?<=Stat)", perl = TRUE)[[1]]), collapse = "_")
  }
  
  layer.names <- c("mapping", "data", "geom", "position", "stat", "show.legend", "inherit.aes", "aes_params", "geom_params", "stat_params")
  
  x <- sapply(layer.names, function(y) {
    b <- l[[y]]
    
    if ("waiver" %in% class(b)) {
      b = NULL
    }
    
    if (y == "geom") {
      b <- eval(parse(text = parent.layer$geom))
    }
    
    if (y == "position") {
      b <- gsub(y, "", tolower(class(b)[1]))
    }
    
    if (y == "stat") {
      b <- eval(parse(text = parent.layer$stat))
    }
    
    b
  })
  
  x$params <- append(x$stat_params, x$geom_params)
  x$params <- append(x$params, x$aes_params)
  x$params <- x$params[!duplicated(names(x$params))]
  
  x$geom_params <- x$aes_params <- x$stat_params <- NULL
  
  fn <- parent.layer$fn
  
  g <- paste0(fn, "()")
  g <- eval(parse(text = g))
  nm <- names(x)
  
  nm <- nm[!sapply(x, typeof) %in% c("environment", "closure", "list")]
  
  geom_aes <- list(
    geom = fn,
    mapping = sapply(names(x$mapping), build_map,y = x$mapping),
    params = sapply(names(x$params), build_map, y = x$params),
    layer = sapply(rev(nm), build_map, y = x[rev(nm)])
  )
  
  nDF <- cbind(names(g$geom$default_aes), paste(g$geom$default_aes))
  nDF[grep("colour|fill|color", nDF[, 1]), 2] <- paste0("'", col2hcl(nDF[grep("colour|fill|color", nDF[, 1]), 2], alpha = NULL), "'")
  
  geom_aes$default <- paste0(apply(nDF, 1, function(x) paste0(x, collapse = "=")))
  
  geom_aes
}

#' @title Default and current ggedit options
#'
#' @description Options for functions in the ggedit package. When running R code, the object \code{ggedit_opts}
#' (default options) is not modified by chunk headers (local chunk options are
#' merged with default options), whereas \code{ggedit_opts_current} (current options)
#' changes with different chunk headers and it always reflects the options for
#' the current chunk.
#'
#' Normally we set up the global options once in the first code chunk in a
#' document using \code{ggedit_opts$set()}, so that all \emph{latter} chunks will
#' use these options. Note the global options set in one chunk will not affect
#' the options in this chunk itself, and that is why we often need to set global
#' options in a separate chunk.
#' 
#' @note \code{ggedit_opts_current} is read-only in the sense that it does nothing if
#'   you call \code{ggedit_opts_current$set()}; you can only query the options via
#'   \code{ggedit_opts_current$get()}.
#' @rdname ggeditOpts
#' @noRd
ggedit_opts <- new_defaults(list(
  fontDefaults = c(
    "sans",
    "Canonical",
    "mono",
    "Courier",
    "Helvetica",
    "serif",
    "Times",
    "AvantGarde",
    "Bookman",
    "Helvetica-Narrow",
    "NewCenturySchoolbook",
    "Palatino",
    "URWGothic",
    "URWBookman",
    "NimbusMon",
    "URWHelvetica",
    "NimbusSan",
    "NimbusSanCond",
    "CenturySch",
    "URWPalladio",
    "URWTimes",
    "NimbusRom"
  ),
  slideDefaults = list(
    alpha = c(min = 0, max = 1),
    size = c(min = 0, max = 10),
    shape = c(min = 1, max = 25),
    stroke = c(min = 0, max = 10),
    weight = c(min = 0, max = 10),
    linetype = c(min = 1, max = 5),
    width = c(min = 0, max = 1),
    angle = c(min = 0, max = 360),
    hjust = c(min = -10, max = 10),
    vjust = c(min = -10, max = 10),
    stroke = c(min = 0, max = 10),
    lineheight = c(min = 0, max = 10),
    linewidth = c(min = 0, max = 5),
    fontface = c(min = 1, max = 4),
    rel_min_height = c(min = 0, max = 1),
    scale = c(min = 0, max = 100)
  ),
  themeTips = list(
    element_rect = list(
      fill = "fill colour",
      colour = "border colour",
      size = "border size (in pts)",
      linetype = paste0(
        paste(
          seq(0, 6),
          c(
            "blank", "solid", "dashed", "dotted", "dotdash",
            "longdash", "twodash"
          ), sep = ": "
        ),
        collapse = ", "
      )
    ),
    element_line = list(
      colour = "line colour",
      size = "numeric (in pts) or \n relative to global size rel(numeric)",
      linetype = paste0(
        paste(
          seq(0, 6),
          c(
            "blank", "solid", "dashed", "dotted", "dotdash",
            "longdash", "twodash"
          ), sep = ": "
        ),
        collapse = ", "
      ),
      lineend = c("butt(default),round,square")
    ),
    element_text = list(
      family = shiny::HTML('<a href="http://www.cookbook-r.com/Graphs/Fonts/" target="_blank">font family</a>'),
      face = 'font face ("plain", "italic", "bold", "bold.italic")',
      colour = "text colour",
      size = "text size (in pts)",
      hjust = "horizontal justification (in [0, 1])",
      vjust = "vertical justification (in [0, 1])",
      angle = "angle (in [0, 360])",
      lineheight = "numeric line height"
    ),
    justification = list(justification = 'anchor point for positioning legend inside plot <br/> "center" or two-element numeric vector'),
    position = list(position = 'the position of legends. <br/> "left", "right", "bottom", "top", or two-element numeric vector')
  ),
  
  ThemeDefaultClass =
    data.frame(
      item = c("angle", "background", "caption", "colour", "face", "family", "fill", "grid.major", "grid.minor", "hjust", "justification", "key", "key.size", "line", "lineheight", "linetype", "margin", "ontop", "position", "size", "subtitle", "switch.pad.grid", "switch.pad.wrap", "text", "text.x", "text.y", "ticks", "ticks.length", "title", "title.x", "title.y", "vjust", "placement"),
      class = c("numeric", "character", "character", "character", "character", "character", "character", "character", "character", "numeric", "character", "character", "character", "character", "numeric", "numeric", "numeric", "character", "character", "numeric", "character", "character", "character", "character", "character", "character", "numeric", "numeric", "character", "character", "character", "numeric", "character"), stringsAsFactors = FALSE
    ),
  session_geoms = 
    data.frame(
      fn = c("annotation_custom", "annotation_logticks", 
             "annotation_map", "annotation_raster", "geom_abline", "geom_area", 
             "geom_bar", "geom_bin2d", "geom_blank", "geom_boxplot", "geom_col", 
             "geom_contour", "geom_count", "geom_crossbar", "geom_curve", 
             "geom_density", "geom_density_2d", "geom_density2d", "geom_dotplot", 
             "geom_errorbar", "geom_errorbarh", "geom_freqpoly", "geom_hex", 
             "geom_histogram", "geom_hline", "geom_jitter", "geom_label", 
             "geom_line", "geom_linerange", "geom_map", "geom_path", "geom_point", 
             "geom_pointrange", "geom_polygon", "geom_qq", "geom_qq_line", 
             "geom_quantile", "geom_raster", "geom_rect", "geom_ribbon", "geom_rug", 
             "geom_segment", "geom_sf", "geom_smooth", "geom_spoke", "geom_step", 
             "geom_text", "geom_tile", "geom_violin", "geom_vline", "stat_bin", 
             "stat_bin_2d", "stat_bin_hex", "stat_bin2d", "stat_binhex", "stat_boxplot", 
             "stat_contour", "stat_count", "stat_density", "stat_density_2d", 
             "stat_density2d", "stat_ecdf", "stat_ellipse", "stat_function", 
             "stat_identity", "stat_qq", "stat_qq_line", "stat_quantile", 
             "stat_sf", "stat_smooth", "stat_sum", "stat_summary", "stat_summary_2d", 
             "stat_summary_bin", "stat_summary_hex", "stat_unique", "stat_ydensity"),
      geom = c("GeomCustomAnn", "GeomLogticks", "GeomAnnotationMap", 
               "GeomRasterAnn", "GeomAbline", "GeomArea", "GeomBar", "GeomTile", 
               "GeomBlank", "GeomBoxplot", "GeomCol", "GeomContour", "GeomPoint", 
               "GeomCrossbar", "GeomCurve", "GeomDensity", "GeomDensity2d", 
               "GeomDensity2d", "GeomDotplot", "GeomErrorbar", "GeomErrorbarh", 
               "GeomPath", "GeomHex", "GeomBar", "GeomHline", "GeomPoint", "GeomLabel", 
               "GeomLine", "GeomLinerange", "GeomMap", "GeomPath", "GeomPoint", 
               "GeomPointrange", "GeomPolygon", "GeomPoint", "GeomPath", "GeomQuantile", 
               "GeomRaster", "GeomRect", "GeomRibbon", "GeomRug", "GeomSegment", 
               "GeomSf", "GeomSmooth", "GeomSpoke", "GeomStep", "GeomText", 
               "GeomTile", "GeomViolin", "GeomVline", "GeomBar", "GeomTile", 
               "GeomHex", "GeomTile", "GeomHex", "GeomBoxplot", "GeomContour", 
               "GeomBar", "GeomArea", "GeomDensity2d", "GeomDensity2d", "GeomStep", 
               "GeomPath", "GeomPath", "GeomPoint", "GeomPoint", "GeomPath", 
               "GeomQuantile", "GeomRect", "GeomSmooth", "GeomPoint", "GeomPointrange", 
               "GeomTile", "GeomPointrange", "GeomHex", "GeomPoint", "GeomViolin"), 
      position = c("PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionStack", "PositionStack", 
                   "PositionIdentity", "PositionIdentity", "PositionDodge2", "PositionStack", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionStack", "PositionIdentity", "PositionJitter", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionDodge", "PositionIdentity", "PositionStack", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionDodge2", "PositionIdentity", "PositionStack", "PositionStack", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionIdentity", "PositionIdentity", "PositionIdentity", 
                   "PositionIdentity", "PositionDodge"), 
      stat = c("StatIdentity", "StatIdentity", "StatIdentity", "StatIdentity", "StatIdentity", 
               "StatIdentity", "StatCount", "StatBin2d", "StatIdentity", "StatBoxplot", 
               "StatIdentity", "StatContour", "StatSum", "StatIdentity", "StatIdentity", 
               "StatDensity", "StatDensity2d", "StatDensity2d", "StatBindot", 
               "StatIdentity", "StatIdentity", "StatBin", "StatBinhex", "StatBin", 
               "StatIdentity", "StatIdentity", "StatIdentity", "StatIdentity", 
               "StatIdentity", "StatIdentity", "StatIdentity", "StatIdentity", 
               "StatIdentity", "StatIdentity", "StatQq", "StatQqLine", "StatQuantile", 
               "StatIdentity", "StatIdentity", "StatIdentity", "StatIdentity", 
               "StatIdentity", "StatSf", "StatSmooth", "StatIdentity", "StatIdentity", 
               "StatIdentity", "StatIdentity", "StatYdensity", "StatIdentity", 
               "StatBin", "StatBin2d", "StatBinhex", "StatBin2d", "StatBinhex", 
               "StatBoxplot", "StatContour", "StatCount", "StatDensity", "StatDensity2d", 
               "StatDensity2d", "StatEcdf", "StatEllipse", "StatFunction", "StatIdentity", 
               "StatQq", "StatQqLine", "StatQuantile", "StatSf", "StatSmooth", 
               "StatSum", "StatSummary", "StatSummary2d", "StatSummaryBin", 
               "StatSummaryHex", "StatUnique", "StatYdensity"), 
      pkg = rep("ggplot2",  77),stringsAsFactors = FALSE)
))

# Function to convert color to HCL and adjust components
#' @importFrom grDevices col2rgb rgb2hsv hcl rgb
col2hcl <- function(colour, h = NULL, c = NULL, l = NULL, alpha = NULL) {
  # Convert color to RGB
  rgb_col <- grDevices::col2rgb(colour) / 255
  
  # Convert RGB to HSV
  hsv_col <- grDevices::rgb2hsv(rgb_col)
  
  # Convert HSV to HCL
  hue <- hsv_col[1] * 360  # Convert hue to degrees (0-360)
  chroma <- hsv_col[2] * 100  # Chroma is similar to saturation
  luminance <- hsv_col[3] * 100  # Luminance is related to value
  
  # Allow user to override H, C, or L values
  if (!is.null(h)) hue <- h
  if (!is.null(c)) chroma <- c
  if (!is.null(l)) luminance <- l
  
  # Create the HCL color with potentially modified components
  hcl_col <- grDevices::hcl(h = hue, c = chroma, l = luminance)
  
  # Convert HCL back to RGB to apply alpha
  rgb_col_with_alpha <- grDevices::col2rgb(hcl_col) / 255
  
  # Add alpha transparency if provided
  if (!is.null(alpha)) {
    rgba_col <- grDevices::rgb(rgb_col_with_alpha[1], rgb_col_with_alpha[2], rgb_col_with_alpha[3], alpha = alpha)
  } else {
    rgba_col <- grDevices::rgb(rgb_col_with_alpha[1], rgb_col_with_alpha[2], rgb_col_with_alpha[3])
  }
  
  # Return the final RGBA color
  return(rgba_col)
}

#' @importFrom stats as.formula
#' @importFrom rlang quo_name
#' @noRd
build_map <- function(item,y) {
  
  y <- y[[item]]
  
  if (inherits(y,'quosure')){
    return(sprintf('%s = %s',item,rlang::quo_name(y)))
  }
  
  if (inherits(y,'character')){
    return(sprintf("%s = '%s'",item,y))
  }
  
  if (inherits(y, "formula")){
    return(sprintf("formula=stats::as.formula('%s')",
                   paste0(as.character(y)[-1], collapse = "~")))
  }
  
  
  if (inherits(y,'NULL')) {
    return(sprintf('%s = NULL',item))
  }
  
  
  if (inherits(y, c("function", "call", "ggproto"))) {
    return(sprintf("%s = %s",
                   item,
                   paste(capture.output(
                     dput(y)),
                     collapse = "\n")
    ))
  }
  
  if (inherits(y, c("data.frame"))) {
    return(paste0("=", paste(capture.output(dput(y)), collapse = "\n")))
  }
  

  return(sprintf('%s = %s',item, y))
  
}

#' @noRd
capture.output <- function (..., file = NULL, append = FALSE, type = c("output", "message"), split = FALSE) 
{
  type <- match.arg(type)
  rval <- NULL
  closeit <- TRUE
  if (is.null(file)) 
    file <- textConnection("rval", "w", local = TRUE)
  else if (is.character(file)) 
    file <- file(file, if (append) 
      "a"
      else "w")
  else if (inherits(file, "connection")) {
    if (!isOpen(file)) 
      open(file, if (append) 
        "a"
        else "w")
    else closeit <- FALSE
  }
  else stop("'file' must be NULL, a character string or a connection")
  sink(file, type = type, split = split)
  on.exit({
    sink(type = type, split = split)
    if (closeit) close(file)
  })
  for (i in seq_len(...length())) {
    out <- withVisible(...elt(i))
    if (out$visible) 
      print(out$value)
  }
  on.exit()
  sink(type = type, split = split)
  if (closeit) 
    close(file)
  if (is.null(rval)) 
    invisible(NULL)
  else rval
}

avgHexColor <- function(colors, ctrlcolor){
  colors <- lapply(colors, col2rgb)
  rgb(t(Reduce(`+`, colors)/length(colors)), maxColorValue=255)
}

fill_na_with_preceding <- function(x) {
  if (all(is.na(x))) return(x)
  for (i in 2:length(x)) {
    if (is.na(x[i])) {
      x[i] <- x[i - 1]
    }
  }
  return(x)
}