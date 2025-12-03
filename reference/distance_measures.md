# Calculate difference between S and Sigma_hat

Calculate the difference between the empirical (S) and the model-implied
indicator variance-covariance matrix (Sigma_hat) using different
distance measures.

## Usage

``` r
calculateDG(
  .object = NULL,
  .matrix1 = NULL,
  .matrix2 = NULL,
  .saturated = FALSE,
  ...
)

calculateDL(
  .object = NULL,
  .matrix1 = NULL,
  .matrix2 = NULL,
  .saturated = FALSE,
  ...
)

calculateDML(
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

- .matrix1:

  A `matrix` to compare.

- .matrix2:

  A `matrix` to compare.

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- ...:

  Ignored.

## Value

A single numeric value giving the distance between two matrices.

## Details

The distances may also be computed for any two matrices A and B by
supplying A and B directly via the `.matrix1` and `.matrix2` arguments.
If A and B are supplied `.object` is ignored.

## Functions

- `calculateDG()`: The geodesic distance (dG).

- `calculateDL()`: The squared Euclidean distance

- `calculateDML()`: The distance measure (fit function) used by ML
