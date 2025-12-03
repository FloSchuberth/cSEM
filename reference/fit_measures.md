# Model fit measures

Calculate fit measures.

## Usage

``` r
calculateChiSquare(.object, .saturated = FALSE)

calculateChiSquareDf(.object)

calculateCFI(.object)

calculateGFI(.object, .type_gfi = c("ML", "GLS", "ULS"), ...)

calculateCN(.object, .alpha = 0.05, ...)

calculateIFI(.object)

calculateNFI(.object)

calculateNNFI(.object)

calculateRMSEA(.object)

calculateRMSTheta(.object)

calculateSRMR(
  .object = NULL,
  .matrix1 = NULL,
  .matrix2 = NULL,
  .saturated = FALSE,
  ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .type_gfi:

  Character string. Which fitting function should the GFI be based on?
  One of *"ML"* for the maximum likelihood fitting function, *"GLS"* for
  the generalized least squares fitting function or *"ULS"* for the
  unweighted least squares fitting function (same as the squared
  Euclidean distance). Defaults to *"ML"*.

- ...:

  Ignored.

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .matrix1:

  A `matrix` to compare.

- .matrix2:

  A `matrix` to compare.

## Value

A single numeric value.

## Details

See the [Fit
indices](https://floschuberth.github.io/cSEM/articles/Using-assess.html#fit_indices)
section of the [cSEM
website](https://floschuberth.github.io/cSEM/index.html) for details on
the implementation.

## Functions

- `calculateChiSquare()`: The chi square statistic.

- `calculateChiSquareDf()`: The Chi square statistic divided by its
  degrees of freedom.

- `calculateCFI()`: The comparative fit index (CFI).

- `calculateGFI()`: The goodness of fit index (GFI).

- `calculateCN()`: The Hoelter index alias Hoelter's (critical) N (CN).

- `calculateIFI()`: The incremental fit index (IFI).

- `calculateNFI()`: The normed fit index (NFI).

- `calculateNNFI()`: The non-normed fit index (NNFI). Also called the
  Tucker-Lewis index (TLI).

- `calculateRMSEA()`: The root mean square error of approximation
  (RMSEA).

- `calculateRMSTheta()`: The root mean squared residual covariance
  matrix of the outer model residuals (RMS theta).

- `calculateSRMR()`: The standardized root mean square residual (SRMR).
