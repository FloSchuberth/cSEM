// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// SVD on full X: decompose X = UDV', beta = V * D^{-1} * U' * y
// [[Rcpp::export]]
arma::vec arma_svd_X(const arma::mat& X, const arma::vec& y) {
  arma::mat U, V;
  arma::vec d;
  arma::svd_econ(U, d, V, X);
  return V * ((U.t() * y) / d);
}

// SVD on X'X: decompose X'X via SVD, then back-solve
// [[Rcpp::export]]
arma::vec arma_svd_XtX(const arma::mat& X, const arma::vec& y) {
  arma::mat XtX = X.t() * X;
  arma::mat U, V;
  arma::vec d;
  arma::svd(U, d, V, XtX);
  return V * ((U.t() * (X.t() * y)) / d);
}

// Pseudoinverse of X'X (analogous to MASS::ginv)
// [[Rcpp::export]]
arma::vec arma_pinv(const arma::mat& X, const arma::vec& y) {
  arma::mat XtX = X.t() * X;
  return arma::pinv(XtX) * (X.t() * y);
}

// Cholesky-based solve (analogous to base::solve)
// [[Rcpp::export]]
arma::vec arma_solve(const arma::mat& X, const arma::vec& y) {
  arma::mat XtX = X.t() * X;
  arma::vec Xty = X.t() * y;
  return arma::solve(XtX, Xty);
}

// --- IGSCA Step 3: U-matrix update -------------------------------------------

// Full QR approach (2017): forms complete N×N Q, extracts Q₂ — O(N²) memory
// [[Rcpp::export]]
arma::mat arma_update_U_qrQ(const arma::mat& Eta_normed,
                             const arma::mat& Z_normed,
                             const arma::mat& D) {
  arma::uword P = Eta_normed.n_cols;
  arma::mat Q, R;
  arma::qr(Q, R, Eta_normed);                      // Q is N×N
  arma::mat Q2 = Q.cols(P, Q.n_cols - 1);           // N × (N-P)

  arma::mat U_svd, V_svd;
  arma::vec d_svd;
  arma::svd_econ(U_svd, d_svd, V_svd, D * Z_normed.t() * Q2);  // J × (N-P)

  arma::mat Utilde = V_svd * U_svd.t();             // (N-P) × J
  return Q2 * Utilde;                                // N × J
}

// Thin QR + projection approach (Householder-equivalent for Armadillo):
// Uses qr_econ to get Q₁ (N×P), then projects via I - Q₁Q₁' — O(NJ+NP) memory
// [[Rcpp::export]]
arma::mat arma_update_U_qrqty(const arma::mat& Eta_normed,
                               const arma::mat& Z_normed,
                               const arma::mat& D) {
  arma::mat Q1, R;
  arma::qr_econ(Q1, R, Eta_normed);                 // Q1 is N×P (thin)

  // Q₂'Z = Z - Q₁(Q₁'Z) via projection
  arma::mat QtZ_null = Z_normed - Q1 * (Q1.t() * Z_normed);  // N × J

  arma::mat U_svd, V_svd;
  arma::vec d_svd;
  arma::svd_econ(U_svd, d_svd, V_svd, D * QtZ_null.t());     // J × N

  arma::mat Utilde = V_svd * U_svd.t();             // N × J

  // Q₂ Ũ = Ũ - Q₁(Q₁'Ũ) via projection
  return Utilde - Q1 * (Q1.t() * Utilde);           // N × J
}
