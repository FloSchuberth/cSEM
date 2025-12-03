# Relative Goodness of Fit (relative GoF)

Calculate the Relative Goodness of Fit (GoF) proposed by Vinzi et al.
(2010) . Note that, contrary to what the name suggests, the Relative GoF
is **not** a measure of model fit in the sense of SEM. See e.g. Henseler
and Sarstedt (2012) for a discussion.

## Usage

``` r
calculateRelativeGoF(
 .object              = NULL
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Value

A single numeric value.

## References

Henseler J, Sarstedt M (2012). “Goodness-of-fit Indices for Partial
Least Squares Path Modeling.” *Computational Statistics*, **28**(2),
565–580.
[doi:10.1007/s00180-012-0317-1](https://doi.org/10.1007/s00180-012-0317-1)
.  
  
Vinzi VE, Trinchera L, Amato S (2010). “PLS path modeling: From
foundations to recent developments and open issues for model assessment
and improvement.” In Vinzi VE, Wang H (eds.), *Handbook of Partial Least
Squares*, 47–82. Springer.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
