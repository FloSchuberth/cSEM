#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat ckronecker(arma::mat X, arma::mat Y) {
  return arma::kron(X, Y);
}