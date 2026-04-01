# Drafted using Claude
# Benchmark: svd() vs MASS::ginv() for solving OLS normal equations
#
# OLS normal equations: beta = (X'X)^{-1} X'y
#   R implementations:
#   - svd_X approach:   decompose X = UDV', then beta = V D^{-1} U'y
#   - svd_XtX approach: decompose X'X via svd, then back-solve
#   - ginv approach:    beta = ginv(X'X) %*% X'y  (Moore-Penrose pseudoinverse)
#   - solve baseline:   solve(X'X, X'y)           (Cholesky-based)
#
#   RcppArmadillo implementations (same 4 algorithms, compiled C++):
#   - arma_svd_X:  arma::svd on full X
#   - arma_svd_XtX: arma::svd on X'X
#   - arma_pinv:   arma::pinv(X'X) * X'y
#   - arma_solve:  arma::solve(X'X, X'y)
#
# We vary n (observations) and p (predictors) across orders of magnitude
# to see how time and memory scale.

library(bench)
library(MASS)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(2)
devtools::load_all()  # compile and load RcppArmadillo solvers from src/

# --- Solver implementations ---------------------------------------------------

solve_ols_svd <- function(X, y) {
  s <- svd(X)
  # beta = V %*% diag(1/d) %*% t(U) %*% y
  # Efficient: avoid forming the diagonal matrix
  s$v %*% ((crossprod(s$u, y)) / s$d)
}

solve_ols_ginv <- function(X, y) {
  ginv(crossprod(X)) %*% crossprod(X, y)
}

solve_ols_svd_XtX <- function(X, y) {
  XtX <- crossprod(X)
  s <- svd(XtX)
  # (X'X)^{-1} = V D^{-2} V' for the normal equations
  s$v %*% (crossprod(s$u, crossprod(X, y)) / s$d)
}

solve_ols_solve <- function(X, y) {
  solve(crossprod(X), crossprod(X, y))
}

# --- Benchmark via bench::press -----------------------------------------------
# Grid: n in {1e3, 1e4, 1e5} x p in {5, 50, 500}  (all combos have n > p)

results <- bench::press(
  n = c(1e3, 1e4, 1e5),
  p = c(5L, 50L, 500L),
  {
    set.seed(42)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    y <- X %*% rnorm(p) + rnorm(n, sd = 0.5)

    bench::mark(
      svd_X      = solve_ols_svd(X, y),
      svd_XtX    = solve_ols_svd_XtX(X, y),
      ginv       = solve_ols_ginv(X, y),
      solve      = solve_ols_solve(X, y),
      arma_svd_X   = arma_svd_X(X, y),
      arma_svd_XtX = arma_svd_XtX(X, y),
      arma_pinv    = arma_pinv(X, y),
      arma_solve   = arma_solve(X, y),
      min_iterations = 3,
      check = FALSE
    )
  }
)

print(results)
