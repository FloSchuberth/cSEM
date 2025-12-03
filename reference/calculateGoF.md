# Goodness of Fit (GoF)

Calculate the Goodness of Fit (GoF) proposed by Tenenhaus et al. (2004)
. Note that, contrary to what the name suggests, the GoF is **not** a
measure of model fit in the sense of SEM. See e.g. Henseler and Sarstedt
(2012) for a discussion.

## Usage

``` r
calculateGoF(
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

## Details

The GoF is inherently tied to the common factor model. It is therefore
unclear how to meaningfully interpret the GoF in the context of a model
that contains constructs modeled as composites.

## References

Henseler J, Sarstedt M (2012). “Goodness-of-fit Indices for Partial
Least Squares Path Modeling.” *Computational Statistics*, **28**(2),
565–580.
[doi:10.1007/s00180-012-0317-1](https://doi.org/10.1007/s00180-012-0317-1)
.  
  
Tenenhaus M, Amanto S, Vinzi VE (2004). “A Global Goodness-of-Fit Index
for PLS Structural Equation Modelling.” In *Proceedings of the XLII SIS
Scientific Meeting*, 739–742.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
