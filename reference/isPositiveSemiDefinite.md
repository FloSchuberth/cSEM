# Internal: Check if matrix is positive semi definite

Check if a matrix is positive semi definite (PSD). A matrix is PSD if
all of its eigenvalues are greater than or equal to zero.

## Usage

``` r
isPositiveSemiDefinite(.matrix, .tolerance = 1e-08)
```

## Arguments

- .matrix:

  Numeric matrix to be checked.

- .tolerance:

  Tolerance for close to zero (negative) values.

## Value

A logical scalar: `TRUE` if `X` is PSD, otherwise `FALSE`.
