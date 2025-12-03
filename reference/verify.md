# Verify admissibility

**\[stable\]**

## Usage

``` r
verify(.object)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Value

A logical vector indicating which (if any) problem occurred. A `FALSE`
indicates that the specific problem did not occurred. For models
containing second-order constructs estimated by the two/three-stage
approach, a list of two such vectors (one for the first and one for the
second stage) is returned. Status codes are:

- 1: The algorithm has converged.

- 2: All absolute standardized loading estimates are smaller than or
  equal to 1. A violation implies either a negative variance of the
  measurement error or a correlation larger than 1.

- 3: The construct VCV is positive semi-definite.

- 4: All reliability estimates are smaller than or equal to 1.

- 5: The model-implied indicator VCV is positive semi-definite. This is
  only checked for linear models (including models containing
  second-order constructs).

## Details

Verify admissibility of the results obtained using
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

Results exhibiting one of the following defects are deemed inadmissible:
non-convergence of the algorithm used to obtain weights, loadings and/or
(congeneric) reliabilities larger than 1, a construct
variance-covariance (VCV) and/or model-implied VCV matrix that is not
positive semi-definite.

If `.object` is of class `cSEMResults_2ndorder` (i.e., estimates are
based on a model containing second-order constructs) both the first and
the second stage are checked separately.

Currently, a model-implied indicator VCV matrix for nonlinear model is
not available. `verify()` therefore skips the check for positive
definiteness of the model-implied indicator VCV matrix for nonlinear
models and returns "ok".

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)

## Examples

``` r
### Without higher order constructs --------------------------------------------
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"
  
# Estimate
out <- csem(threecommonfactors, model)
  
# Check admissibility
verify(out) # ok!
#> ________________________________________________________________________________
#> 
#> Verify admissibility:
#> 
#>     admissible
#> 
#> Details:
#> 
#>   Code   Status    Description
#>   1      ok        Convergence achieved                                   
#>   2      ok        All absolute standardized loading estimates <= 1       
#>   3      ok        Construct VCV is positive semi-definite                
#>   4      ok        All reliability estimates <= 1                         
#>   5      ok        Model-implied indicator VCV is positive semi-definite  
#> ________________________________________________________________________________

## Examine the structure of a cSEMVerify object
str(verify(out))
#>  'cSEMVerify' Named logi [1:5] FALSE FALSE FALSE FALSE FALSE
#>  - attr(*, "names")= chr [1:5] "1" "2" "3" "4" ...

### With higher order constructs -----------------------------------------------
# If the model containes higher order constructs both the first and the second-
# stage estimates estimates are checked for admissibility

if (FALSE) { # \dontrun{
require(cSEM.DGP) # download from https://m-e-rademaker.github.io/cSEM.DGP/
  
# Create DGP with 2nd order construct. Loading for indicator y51 is set to 1.1
# to produce a failing first stage model
  
dgp_2ndorder <- "
## Path model / Regressions
eta2 ~ 0.5*eta1
eta3 ~ 0.35*eta1 + 0.4*eta2

## Composite model
eta1 =~ 0.8*y41 + 0.6*y42 + 0.6*y43
eta2 =~ 1.1*y51 + 0.7*y52 + 0.7*y53
c1   =~ 0.8*y11 + 0.4*y12
c2   =~ 0.5*y21 + 0.3*y22

## Higher order composite
eta3 =~ 0.4*c1 + 0.4*c2
"
  
dat <- generateData(dgp_2ndorder) # requires the cSEM.DGP package
out <- csem(dat, .model = dgp_2ndorder)

verify(out) # not ok
} # }
```
