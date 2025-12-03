# Degrees of freedom

Calculate the degrees of freedom for a given model from a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object.

## Usage

``` r
calculateDf(
  .object     = NULL,
  .null_model = FALSE,
  ...
  )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .null_model:

  Logical. Should the degrees of freedom for the null model be computed?
  Defaults to `FALSE`.

- ...:

  Ignored.

## Value

A single numeric value.

## Details

Although, composite-based estimators always retrieve parameters of the
postulated models via the estimation of a composite model, the
computation of the degrees of freedom depends on the postulated model.

See: [cSEM
website](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
for details on how the degrees of freedom are calculated.

To compute the degrees of freedom of the null model use
`.null_model = TRUE`. The degrees of freedom of the null model are
identical to the number of non-redundant off-diagonal elements of the
empirical indicator correlation matrix. This implicitly assumes a null
model with model-implied indicator correlation matrix equal to the
identity matrix.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
