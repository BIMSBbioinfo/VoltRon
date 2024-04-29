#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double calculateMoransI(NumericMatrix data, NumericMatrix datadist, double sumW) {

  // Number of rows and columns
  int p = data.nrow();
  int n = data.ncol();

  // Calculate Moran's I
  NumericVector rowMeans(p);
  for (int i = 0; i < p; i++) { // for all features

    // calculate mean
    double rowSum = 0.0;
    for (int mean_ind = 0; mean_ind < n; mean_ind++) {
      rowSum += data(i, mean_ind);
    }
    rowMeans[i] = rowSum / n;

    // calculate I statistics
    double numerator = 0.0;
    double denominator = 0.0;
    double moranI = 0.0;
    for (int jj = 0; jj < n; jj++) {

      // denominator
      denominator += (data(i, jj) - rowMeans[i]) * (data(i, jj) - rowMeans[i]);

      // nominator
      if (jj < n-1){
        for (int kk = jj+1; kk < n; kk++) {
          numerator += datadist(jj,kk) * (data(i, jj) - rowMeans[i]) * (data(i, kk) - rowMeans[i]);
        }
      }
    }
    moranI = (2.0*numerator)/denominator;
    // std::cout << "Pi: " << i << std::endl;
  }

  // return
  return 0.0;
}

