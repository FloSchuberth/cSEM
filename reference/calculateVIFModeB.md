# Calculate variance inflation factors (VIF) for weights obtained by PLS Mode B

Calculate the variance inflation factor (VIF) for weights obtained by
PLS-PM's Mode B.

## Usage

``` r
calculateVIFModeB(.object = NULL)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Value

A named list of vectors containing the VIF values. Each list name is the
name of a construct whose weights were obtained by Mode B. The vectors
contain the VIF values obtained from a regression of each explanatory
variable of a given construct on the remaining explanatory variables of
that construct.

If the weighting approach is not `"PLS-PM"` or for none of the
constructs Mode B is used, the function silently returns `NA`.

## Details

Weight estimates obtained by Mode B can suffer from multicollinearity.
VIF values are commonly used to assess the severity of
multicollinearity.

The function is only applicable to objects of class
`cSEMResults_default`. For other object classes use
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
