# Fornell-Larcker criterion

Computes the Fornell-Larcker matrix.

## Usage

``` r
calculateFLCriterion(
  .object              = NULL,
  .only_common_factors = TRUE,
  ...
  )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .only_common_factors:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- ...:

  Ignored.

## Value

A matrix with the squared construct correlations on the off-diagonal and
the AVEs on the main diagonal.

## Details

The Fornell-Larcker criterion (FL criterion) is a rule suggested by
Fornell and Larcker (1981) to assess discriminant validity. The
Fornell-Larcker criterion is a decision rule based on a comparison
between the squared construct correlations and the average variance
extracted (AVE).

The FL criterion is inherently tied to the common factor model. It is
therefore unclear how to meaningfully interpret the FL criterion in the
context of a model that contains constructs modeled as composites.

## References

Fornell C, Larcker DF (1981). “Evaluating structural equation models
with unobservable variables and measurement error.” *Journal of
Marketing Research*, **XVIII**, 39–50.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
