# Plan: Benchmark qr.qty vs qr.Q for U-matrix Update (R + Armadillo)

## Context
The IGSCA Step 3 (unique scores update) has two R implementations:
1. **Original 2017 (qr.Q):** Forms the full N×N Q matrix via `qr.Q(complete=TRUE)`, extracts Q₂, then does SVD — O(N²) memory, O(N²J) time
2. **New implicit Householder (qr.qty):** Uses `qr.qty()`/`qr.qy()` to apply Householder reflections without forming Q — O(NJ+NP) memory, O(NPJ) time

We want a `bench::press()` comparison of both approaches, plus RcppArmadillo counterparts, varying N and P/J.

## Step 1: New R benchmark file
Create `dev/igsca/benchmarking_R/qr_update_benchmark.R` with:

### R implementations (extracted from `R/helper_igsca.R`):

**`update_U_qrQ(Eta_normed, Z_normed, D)`** — original 2017 approach:
```r
qr_obj <- qr(Eta_normed)
Eta_Q2 <- qr.Q(qr_obj, complete = TRUE)[, (P+1):N, drop = FALSE]
svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
Utilde <- svd_mx$v %*% t(svd_mx$u)
U <- Eta_Q2 %*% Utilde
```

**`update_U_qrqty(Eta_normed, Z_normed, D)`** — new implicit approach:
```r
qr_eta <- qr(Eta_normed)
QtZ_null <- qr.qty(qr_eta, Z_normed)[(P+1):N, , drop = FALSE]
svd_mx <- svd(D %*% t(QtZ_null))
Utilde <- svd_mx$v %*% t(svd_mx$u)
U <- qr.qy(qr_eta, rbind(matrix(0, P, J), Utilde))
```

### RcppArmadillo counterparts in `src/ols_solvers.cpp`:

**`arma_update_U_qrQ(Eta_normed, Z_normed, D)`** — full Q via `arma::qr(Q, R, Eta)`, extract Q₂, same SVD flow

**`arma_update_U_qrqty(Eta_normed, Z_normed, D)`** — use `arma::qr_econ` for thin QR to get Q₁ (N×P), then compute Q₂'Z and Q₂Ũ via projection (Z - Q₁(Q₁'Z)) instead of Householder reflections. Armadillo doesn't expose raw Householder reflectors like R does, so the projection approach is the closest equivalent.

### bench::press grid
```r
bench::press(
  N = c(1e3, 1e4, 1e5),
  P = c(2L, 5L, 10L),
  J = c(10L, 20L, 50L),
  { ... }
)
```
27 combinations. All have N >> P and N >> J.

Setup block generates: `Eta_normed` (N×P), `Z_normed` (N×J), `D` (J×J diagonal).

## Step 2: Add 2 new C++ functions to `src/ols_solvers.cpp`
- `arma_update_U_qrQ`: full QR approach
- `arma_update_U_qrqty`: thin QR + projection approach

Both take `arma::mat Eta_normed`, `arma::mat Z_normed`, `arma::mat D` and return `arma::mat U` (N×J).

## Step 3: Recompile and run
- `devtools::document()` to regenerate exports
- Run the benchmark

## Files to create/modify
- **Create:** `dev/igsca/benchmarking_R/qr_update_benchmark.R`
- **Modify:** `src/ols_solvers.cpp` (add 2 functions)
- **Auto-modified:** `src/RcppExports.cpp`, `R/RcppExports.R` (by `devtools::document()`)

## Verification
- `devtools::load_all()` compiles cleanly
- Quick sanity check: all 4 methods return numerically similar U on a small test case
- Run full benchmark grid
