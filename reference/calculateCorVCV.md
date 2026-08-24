# Internal: Distribution-free asymptotic variance-covariance matrix of correlations

Large-sample variance-covariance matrix of the unique sample
correlations of the columns of `.X`, using the estimator of Steiger and
Hakstian (1982) .

## Usage

``` r
calculateCorVCV(.X)
```

## Arguments

- .X:

  A matrix of processed data (scaled, cleaned and ordered).

## Value

The `(L x L)` asymptotic variance-covariance matrix of the correlations,
`L = P(P - 1)/2`, in lower-triangular order.

## Details

The correlations are taken in lower-triangular order of the correlation
matrix of `.X`, the same order used by the HTMT gradient, so that a
gradient vector and this matrix align by position for a delta-method
confidence interval, i.e. `Var(f) = t(g) %*% calculateCorVCV(.X) %*% g`.

## References

Steiger JH, Hakstian AR (1982). “The asymptotic distribution of elements
of a correlation matrix: Theory and application.” *British Journal of
Mathematical and Statistical Psychology*, **35**(2), 208–215.
[doi:10.1111/j.2044-8317.1982.tb00653.x](https://doi.org/10.1111/j.2044-8317.1982.tb00653.x)
.
