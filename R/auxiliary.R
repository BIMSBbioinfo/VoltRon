fill.na <- function(x, i = 5) {
  if (is.na(x)[i]) {
    return(round(mean(x, na.rm = TRUE), 0))
  }
  else {
    return(round(x[i], 0))
  }
}
