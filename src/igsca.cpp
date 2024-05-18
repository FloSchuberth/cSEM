#include <Rcpp.h>
using namespace Rcpp;

// From the documentation of a default rstudio cpp file
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
