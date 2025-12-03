# Test for overall model fit

**\[maturing\]**

## Usage

``` r
testOMF(
 .object                = NULL, 
 .alpha                 = 0.05,
 .fit_measures          = FALSE,
 .handle_inadmissibles  = c("drop", "ignore", "replace"), 
 .R                     = 499, 
 .saturated             = FALSE,
 .seed                  = NULL,
 ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .fit_measures:

  Logical. (EXPERIMENTAL) Should additional fit measures be included?
  Defaults to `FALSE`.

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*drop*", "*ignore*", or "*replace*". If "*drop*", all
  replications/resamples yielding an inadmissible result will be dropped
  (i.e. the number of results returned will potentially be less than
  `.R`). For "*ignore*" all results are returned even if all or some of
  the replications yielded inadmissible results (i.e. number of results
  returned is equal to `.R`). For "*replace*" resampling continues until
  there are exactly `.R` admissible solutions. Depending on the
  frequency of inadmissible solutions this may significantly increase
  computing time. Defaults to "*drop*".

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- ...:

  Can be used to determine the fitting function used in the calculateGFI
  function.

## Value

A list of class `cSEMTestOMF` containing the following list elements:

- `$Test_statistic`:

  The value of the test statistics.

- `$Critical_value`:

  The corresponding critical values obtained by the bootstrap.

- `$Decision`:

  The test decision. One of: `FALSE` (**Reject**) or `TRUE` (**Do not
  reject**).

- `$Information`:

  The `.R` bootstrap values; The number of admissible results; The seed
  used and the number of total runs.

## Details

Bootstrap-based test for overall model fit originally proposed by Beran
and Srivastava (1985) . See also Dijkstra and Henseler (2015) who first
suggested the test in the context of PLS-PM.

By default, `testOMF()` tests the null hypothesis that the population
indicator correlation matrix equals the population model-implied
indicator correlation matrix. Several discrepancy measures may be used.
By default, `testOMF()` uses four distance measures to assess the
distance between the sample indicator correlation matrix and the
estimated model-implied indicator correlation matrix, namely the
geodesic distance, the squared Euclidean distance, the standardized root
mean square residual (SRMR), and the distance based on the maximum
likelihood fit function. The reference distribution for each test
statistic is obtained by the bootstrap as proposed by Beran and
Srivastava (1985) .

It is possible to perform the bootstrap-based test using fit measures
such as the CFI, RMSEA or the GFI if `.fit_measures = TRUE`. This is
experimental. To the best of our knowledge the applicability and
usefulness of the fit measures for model fit assessment have not been
formally (statistically) assessed yet. Theoretically, the logic of the
test applies to these fit indices as well. Hence, their applicability is
theoretically justified. Only use if you know what you are doing.

If `.saturated = TRUE` the original structural model is ignored and
replaced by a saturated model, i.e., a model in which all constructs are
allowed to correlate freely. This is useful to test misspecification of
the measurement model in isolation.

## References

Beran R, Srivastava MS (1985). “Bootstrap Tests and Confidence Regions
for Functions of a Covariance Matrix.” *The Annals of Statistics*,
**13**(1), 95–115.
[doi:10.1214/aos/1176346579](https://doi.org/10.1214/aos/1176346579) .  
  
Dijkstra TK, Henseler J (2015). “Consistent and Asymptotically Normal
PLS Estimators for Linear Structural Equations.” *Computational
Statistics & Data Analysis*, **81**, 10–23.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
[`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
[`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md),
[`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md),
[`exportToExcel()`](https://floschuberth.github.io/cSEM/reference/exportToExcel.md)

## Examples

``` r
# ===========================================================================
# Basic usage
# ===========================================================================
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Estimate
out <- csem(threecommonfactors, model, .approach_weights = "PLS-PM")

## Test
testOMF(out, .R = 50, .seed = 320)
#> ________________________________________________________________________________
#> --------- Test for overall model fit based on Beran & Srivastava (1985) --------
#> 
#> Null hypothesis:
#> 
#>        ┌──────────────────────────────────────────────────────────────────┐
#>        │                                                                  │
#>        │   H0: The model-implied indicator covariance matrix equals the   │
#>        │   population indicator covariance matrix.                        │
#>        │                                                                  │
#>        └──────────────────────────────────────────────────────────────────┘
#> 
#> Test statistic and critical value: 
#> 
#>                                      Critical value
#>  Distance measure    Test statistic    95%   
#>  dG                      0.0060      0.0192  
#>  SRMR                    0.0158      0.0272  
#>  dL                      0.0112      0.0333  
#>  dML                     0.0320      0.1027  
#>  
#> 
#> Decision: 
#> 
#>                          Significance level
#>  Distance measure             95%        
#>  dG                      Do not reject  
#>  SRMR                    Do not reject  
#>  dL                      Do not reject  
#>  dML                     Do not reject  
#>  
#> Additional information:
#> 
#>  Out of 50 bootstrap replications 50 are admissible.
#>  See ?verify() for what constitutes an inadmissible result.
#> 
#>  The seed used was: 320
#> ________________________________________________________________________________
```
