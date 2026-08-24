# Internal: Delta-method standard errors of the HTMT

Combines the per-construct-pair HTMT gradients with the asymptotic
variance-covariance matrix of the indicator correlations
([`calculateCorVCV()`](https://floschuberth.github.io/cSEM/reference/calculateCorVCV.md))
into the delta-method standard error of each HTMT,
`sqrt(t(g) %*% Sigma %*% g)`. The variance-covariance matrix is formed
once and reused for every pair. Gradient vectors and variance-covariance
matrix share the lower-triangular ordering of the indicator correlation
matrix, so they align by position.

## Usage

``` r
calculateHTMTasymptoticSE(.gradients, .X)
```

## Arguments

- .gradients:

  A named list of gradient vectors, one per construct pair, each aligned
  with the lower triangular of the indicator correlation matrix.

- .X:

  A matrix of processed data (scaled, cleaned and ordered).

## Value

A named numeric vector of standard errors, one per construct pair.

## Details

No lower bound is imposed on the variance: since `Sigma` is positive
semi-definite, `t(g) %*% Sigma %*% g` is non-negative by construction,
so a negative value would signal a genuine problem and is deliberately
allowed to surface as `NaN` (with a warning) rather than being silently
clamped to zero.
