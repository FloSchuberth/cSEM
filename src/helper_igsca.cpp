#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat kroneckerC(arma::mat &X, arma::mat &Y, arma::uvec &idx) {
  arma::mat prod;
  prod = arma::kron(X, Y);
  
  return prod.cols(idx-1);
}


//[[Rcpp::export]]
arma::sp_mat kroneckerCsp(arma::mat &X, arma::mat &Y, arma::uvec &idx) {
  // specific to kronecker product of Identity matrix with arbitrary matrix
  arma::sp_mat prod;
  prod = arma::kron(X, Y); 
  
  return prod.cols(idx - 1);
}