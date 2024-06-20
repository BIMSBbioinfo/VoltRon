#include "Rcpp.h"

// [[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix replaceNaMatrix(Rcpp::NumericMatrix mat, int replace) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (Rcpp::NumericMatrix::is_na(mat(i, j))) {
        mat(i, j) = replace;
      }
    }
  }
  return mat;
}
