# Average variance extracted (AVE)

Calculate the average variance extracted (AVE) as proposed by Fornell
and Larcker (1981) . For details see the [cSEM
website](https://floschuberth.github.io/cSEM/articles/Using-assess.html#ave)

## Usage

``` r
calculateAVE(
 .object              = NULL,
 .only_common_factors = TRUE
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

## Value

A named vector of numeric values (the AVEs). If `.object` is a list of
`cSEMResults` objects, a list of AVEs is returned.

## Details

The AVE is inherently tied to the common factor model. It is therefore
unclear how to meaningfully interpret the AVE in the context of a
composite model. It is possible, however, to force computation of the
AVE for constructs modeled as composites by setting
`.only_common_factors = FALSE`.

## References

Fornell C, Larcker DF (1981). “Evaluating structural equation models
with unobservable variables and measurement error.” *Journal of
Marketing Research*, **XVIII**, 39–50.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
