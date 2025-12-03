# Internal: Calculate direct, indirect and total effect

The direct effects are equal to the estimated coefficients. The total
effect equals (I-B)^-1 Gamma. The indirect effect equals the difference
between the total effect and the indirect effect. In addition, the
variance accounted for (VAF) is calculated. The VAF is defined as the
ratio of a variables indirect effect to its total effect. Helper for
generic functions
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
and
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).

## Usage

``` r
calculateEffects(
 .object       = NULL,
 .output_type  = c("data.frame", "matrix")
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .output_type:

  Character string. The type of output to return. One of "*complete*" or
  "*structured*". See the Value section for details. Defaults to
  "*complete*".

## Value

A matrix or a data frame of effects.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
