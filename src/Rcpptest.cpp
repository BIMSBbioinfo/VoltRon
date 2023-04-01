#include <Rcpp.h>
#include "my_class.h"
using namespace Rcpp;
using namespace Newclasses;

// [[Rcpp::export]]
List rcpp_hello_world() {
  CharacterVector x = CharacterVector::create( "foo", "bar" );
  NumericVector y = NumericVector::create( 0.0, 1.0 );
  List z = List::create( x, y );

  return z;
}

// [[Rcpp::export]]
void doing_something() {
  Newclass nc;
  nc.do_something();
}
