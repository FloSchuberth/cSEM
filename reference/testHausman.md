# Regression-based Hausman test

**\[experimental\]**

## Usage

``` r
testHausman(
 .object               = NULL,
 .eval_plan            = c("sequential", "multicore", "multisession"),
 .handle_inadmissibles = c("drop", "ignore", "replace"),
 .R                    = 499,
 .resample_method      = c("bootstrap", "jackknife"),
 .seed                 = NULL
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .eval_plan:

  Character string. The evaluation plan to use. One of "*sequential*",
  "*multicore*", or "*multisession*". In the two latter cases all
  available cores will be used. Defaults to "*sequential*".

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

- .resample_method:

  Character string. The resampling method to use. One of: "*none*",
  "*bootstrap*" or "*jackknife*". Defaults to "*none*".

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

## Details

Calculates the regression-based Hausman test to be used to compare OLS
to 2SLS estimates or 2SLS to 3SLS estimates. See e.g., Wooldridge (2010)
(pages 131 f.) for details.

The function is somewhat experimental. Only use if you know what you are
doing.

## References

Wooldridge JM (2010). *Econometric Analysis of Cross Section and Panel
Data*, 2 edition. MIT Press.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)

## Examples

``` r
### Example from Dijkstra & Hensler (2015)
## Prepartion (values are from p. 15-16 of the paper)
Lambda <- t(kronecker(diag(6), c(0.7, 0.7, 0.7)))
Phi <- matrix(c(1.0000, 0.5000, 0.5000, 0.5000, 0.0500, 0.4000, 
                0.5000, 1.0000, 0.5000, 0.5000, 0.5071, 0.6286,
                0.5000, 0.5000, 1.0000, 0.5000, 0.2929, 0.7714,
                0.5000, 0.5000, 0.5000, 1.0000, 0.2571, 0.6286,
                0.0500, 0.5071, 0.2929, 0.2571, 1.0000, sqrt(0.5),
                0.4000, 0.6286, 0.7714, 0.6286, sqrt(0.5), 1.0000), 
              ncol = 6)

## Create population indicator covariance matrix
Sigma <- t(Lambda) %*% Phi %*% Lambda
diag(Sigma) <- 1
dimnames(Sigma) <- list(paste0("x", rep(1:6, each = 3), 1:3),
                        paste0("x", rep(1:6, each = 3), 1:3))

## Generate data
dat <- MASS::mvrnorm(n = 500, mu = rep(0, 18), Sigma = Sigma, empirical = TRUE)
# empirical = TRUE to show that 2SLS is in fact able to recover the true population
# parameters.

## Model to estimate
model <- "
## Structural model (nonrecurisve)
eta5 ~ eta6 + eta1 + eta2
eta6 ~ eta5 + eta3 + eta4

## Measurement model
eta1 =~ x11 + x12 + x13
eta2 =~ x21 + x22 + x23
eta3 =~ x31 + x32 + x33
eta4 =~ x41 + x42 + x43

eta5 =~ x51 + x52 + x53
eta6 =~ x61 + x62 + x63
"

library(cSEM)

## Estimate
res_ols <- csem(dat, .model = model, .approach_paths = "OLS")
sum_res_ols <- summarize(res_ols) 

# Note: For the example the model-implied indicator correlation is irrelevant
#       the warnings can be ignored.

res_2sls <- csem(dat, .model = model, .approach_paths = "2SLS",
                 .instruments = list("eta5" = c('eta1','eta2','eta3','eta4'), 
                                     "eta6" = c('eta1','eta2','eta3','eta4')))
sum_res_2sls <- summarize(res_2sls)
# Note that exogenous constructs are supplied as instruments for themselves!

## Test for endogeneity
test_ha <- testHausman(res_2sls, .R = 200)
test_ha
#> ________________________________________________________________________________
#> ───────────────────────── Regression-based Hausman test ────────────────────────
#> 
#> Null hypothesis:
#> 
#>    ┌──────────────────────────────────────────────────────────────────────────┐
#>    │                                                                          │
#>    │   H0: Variable(s) suspected to be endogenous are uncorrelated with the   │
#>    │   error term (no endogeneity).                                           │
#>    │                                                                          │
#>    └──────────────────────────────────────────────────────────────────────────┘
#> 
#> Regression output: 
#> 
#>  
#>   Dependent construct: 'eta5'
#> 
#>  Independent construct    Estimate  Std. error   t-stat.   p-value
#>  eta1                      -0.3000      0.0994   -3.0186    0.0025
#>  eta2                       0.4999      0.0984    5.0791    0.0000
#>  eta6                       0.2501      0.1136    2.2024    0.0276
#>  Resid_eta6                 0.9784      0.2899    3.3749    0.0007
#> 
#>   Dependent construct: 'eta6'
#> 
#>  Independent construct    Estimate  Std. error   t-stat.   p-value
#>  eta3                       0.4999      0.0511    9.7841    0.0000
#>  eta4                       0.2501      0.0575    4.3511    0.0000
#>  eta5                       0.5001      0.1066    4.6919    0.0000
#>  Resid_eta5                -0.0056      0.1336   -0.0418    0.9667
#> ________________________________________________________________________________
```
