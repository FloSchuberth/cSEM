#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat ckronecker(arma::mat X, arma::mat Y) {
  return kron(X, Y);
}