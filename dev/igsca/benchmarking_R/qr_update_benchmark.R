# Drafted using Claude
# Benchmark: qr.Q (2017) vs qr.qty (implicit Householder) for IGSCA Step 3
#
# IGSCA Step 3 updates the unique scores matrix U by solving an orthogonal
# Procrustes problem constrained so U is orthogonal to construct scores Eta.
#
#   R implementations:
#   - qr.Q:   forms full N×N Q matrix, extracts Q₂ — O(N²) memory
#   - qr.qty: implicit Householder reflections via qr.qty/qr.qy — O(NJ+NP) memory
#
#   RcppArmadillo implementations:
#   - arma_qrQ:   arma::qr() for full Q, extract Q₂
#   - arma_qrqty: arma::qr_econ() for thin Q₁, project via Q₂ = I - Q₁Q₁'
#
# Grid: N (sample size) × P (constructs) × J (indicators)

library(bench)
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(2)
devtools::load_all()

# --- R implementations --------------------------------------------------------

# Original 2017 approach: full Q matrix via qr.Q(complete = TRUE)
update_U_qrQ <- function(Eta_normed, Z_normed, D) {
  N <- nrow(Eta_normed)
  P <- ncol(Eta_normed)
  qr_obj <- qr(Eta_normed)
  Eta_Q2 <- qr.Q(qr_obj, complete = TRUE)[, (P + 1):N, drop = FALSE]
  svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
  Utilde <- svd_mx$v %*% t(svd_mx$u)
  Eta_Q2 %*% Utilde
}

# New implicit Householder approach: qr.qty / qr.qy
update_U_qrqty <- function(Eta_normed, Z_normed, D) {
  N <- nrow(Eta_normed)
  P <- ncol(Eta_normed)
  J <- ncol(Z_normed)
  qr_eta <- qr(Eta_normed)
  QtZ_null <- qr.qty(qr_eta, Z_normed)[(P + 1):N, , drop = FALSE]
  svd_mx <- svd(D %*% t(QtZ_null))
  Utilde <- svd_mx$v %*% t(svd_mx$u)
  qr.qy(qr_eta, rbind(matrix(0, P, J), Utilde))
}

# --- Benchmark via bench::press -----------------------------------------------
# Grid: N in {1e3, 1e4} x P in {2, 5, 10} x J in {10, 20, 50}
# Note: N = 1e5 with qr.Q would need ~80 GB, so capped at 1e4 for safety.

results <- bench::press(
  N = c(1e3, 1e4),
  P = c(2L, 5L, 10L),
  J = c(10L, 20L, 50L),
  {
    set.seed(42)
    Eta_normed <- scale(matrix(rnorm(N * P), N, P)) / sqrt(N - 1)
    Z_normed <- scale(matrix(rnorm(N * J), N, J)) / sqrt(N - 1)
    D <- diag(runif(J, 0.2, 0.9))

    bench::mark(
      qrQ        = update_U_qrQ(Eta_normed, Z_normed, D),
      qrqty      = update_U_qrqty(Eta_normed, Z_normed, D),
      arma_qrQ   = arma_update_U_qrQ(Eta_normed, Z_normed, D),
      arma_qrqty = arma_update_U_qrqty(Eta_normed, Z_normed, D),
      min_iterations = 3,
      check = FALSE
    )
  }
)

print(results)
View(results)