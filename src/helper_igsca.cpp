#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat kroneckerC(arma::mat &X, arma::mat &Y, arma::uvec &idx) {
  arma::mat prod;
  prod = arma::kron(X, Y);
  
  return prod.cols(idx-1);
}