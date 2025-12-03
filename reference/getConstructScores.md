# Get construct scores

**\[stable\]**

## Usage

``` r
getConstructScores(
 .object        = NULL,
 .standardized  = TRUE
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .standardized:

  Logical. Should standardized scores be returned? Defaults to `TRUE`.

## Value

A list of three with elements `Construct_scores`, `W_used`,
`Indicators_used`.

## Details

Get the standardized or unstandardized construct scores.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
