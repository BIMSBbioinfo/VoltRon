#include <Rcpp.h>
#include "my_class.h"
using namespace Rcpp;
using namespace Myclasses;

// [[Rcpp::export]]
void doing_something() {
  Myclass nc;
  nc.do_something();
}
