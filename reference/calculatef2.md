# Calculate Cohen's f^2

Calculate the effect size for regression analysis (Cohen 1992) known as
Cohen's f^2.

## Usage

``` r
calculatef2(.object = NULL)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Value

A matrix with as many rows as there are structural equations. The number
of columns is equal to the total number of right-hand side variables of
these equations.

## References

Cohen J (1992). “A power primer.” *Psychological Bulletin*, **112**(1),
155–159.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
