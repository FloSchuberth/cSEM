# Internal: Check if matrix is positive semi definite

Check if a matrix is positive semi definite (PSD). A matrix is PSD if
all of its eigenvalues are greater than or equal to zero. One caveat,
the package previously used
[`matrixcalc::is.positive.semi.definite()`](https://rdrr.io/pkg/matrixcalc/man/is.positive.semi.definite.html)
where negative values close to zero `(abs(eigenvalues) < tol)` were not
counted as negative. This behaviour is replicated here as well.

## Usage

``` r
isPositiveSemiDefinite(X, zero.tol = 1e-08)
```

## Arguments

- X:

  Numeric matrix to be checked.

- zero.tol:

  Tolerance for close to zero (negative) values.

## Value

A logical scalar: `TRUE` if `X` is PSD, otherwise `FALSE`.
