
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![Build
Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)

<!-- WARNING: THIS IS WORK IN PROGRESS. BREAKING CHANGES TO THE API ARE VERY LIKELY.  -->

<!--          Use the package with caution and please report bugs to [the package developers](mailto:manuel.rademaker@uni-wuerzburg.de;f.schuberth@utwente.nl).  -->

<!--          The first stable relase will be version 0.0.1, most likely towards the end -->

<!--          of 2019. -->

## Purpose

Estimate, analyse, test, and study linear, nonlinear, hierachical and
multi-group structural equation models using composite-based approaches
and procedures, including estimation techniques such as partial least
squares path modeling (PLS-PM) and its derivatives (PLSc, OrdPLSc,
robustPLSc), generalized structured component analysis (GSCA),
generalized structured component analysis with uniqueness terms (GSCAm),
generalized canonical correlation analysis (GCCA), principal component
analysis (PCA), factor score regression (FSR) using sum score,
regression or bartlett scores (including bias correction using Croon’s
approach), as well as several tests and typical postestimation
procedures (e.g., verify admissibility of the estimates, assess the
model fit, test the model fit, compute confidence intervals, compare
groups, etc.).

## Installation

The package is available on [CRAN](https://cran.r-project.org/):

``` r
install.packages("cSEM")
```

To install the development version use:

``` r
# install.packages("devtools")
devtools::install_github("M-E-Rademaker/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://m-e-rademaker.github.io/cSEM/).

<!-- ## Philosophy -->

<!-- - First and foremost: `cSEM` has a user-centered design!. "User-centered" mainly  -->

<!--   boils down to: `cSEM` is easy, i.e. intuitive to use by non-R experts!  -->

<!-- - Modern in a sense that the package integrates modern developments within  -->

<!--   the R community. This mainly includes ideas/recommendations/design choices that -->

<!--   fead into the packages of the [tidyverse](https://github.com/tidyverse/tidyverse). -->

<!-- - State of the art in a sense that we seek to quickly implement recent methodological -->

<!--   developments in composite-based SEM.  -->

## Basic usage

The basic usage is illustrated below.

<img src="man/figures/api.png" width="80%" style="display: block; margin: auto;" />

Usully, using `cSEM` is the same 3 step procedure:

> 1.  Pick a dataset and specify a model using [lavaan
>     syntax](http://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use `csem()`
> 3.  Apply one of the postestimation functions listed below on the
>     resulting object.

## Postestimation functions

There are five major postestimation verbs, four test family functions
and two do-family of function:

  - `assess()` : assess the model using common quality criteria
  - `infer()` : calculate common inferencial quantities (e.g., standard
    errors, confidence intervals)
  - `predict()` : predict endogenous indicator values
  - `summarize()` : summarize the results
  - `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently the
following tests are implemented:

  - `testOMF()` : performs a test for overall model fit
  - `testMICOM()` : performs a test for composite measurement invariance
  - `testMGD()` : performs several tests to assess multi-group
    differences
  - `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently two do functions are implemented:

  - `doFloodlightAnalysis()`: performs a floodlight analysis
  - `doRedundancyAnalysis()`: performs a redundancy analysis

All functions require a `cSEMResults` object.

## Example

Models are defined using [lavaan
syntax](http://lavaan.ugent.be/tutorial/syntax1.html) with some slight
modifications (see the [Specifying a
model](https://m-e-rademaker.github.io/cSEM/articles/cSEM.html#using-csem)
section on the [cSEM-website](https://m-e-rademaker.github.io/cSEM/)).
For illustration we use the build-in and well-known `satisfaction`
dataset.

``` r
require(cSEM)
    
## Note: The operator "<~" tells cSEM that the construct to its left is modelled
##       as a composite.
##       The operator "=~" tells cSEM that the construct to its left is modelled
##       as a common factor.
##       The operator "~" tells cSEM which are the dependent (left-hand side) and
##       independent variables (right-hand side).
    
model <- "
# Structural model
EXPE ~ IMAG
QUAL ~ EXPE
VAL  ~ EXPE + QUAL
SAT  ~ IMAG + EXPE + QUAL + VAL 
LOY  ~ IMAG + SAT

# Composite model
IMAG <~ imag1 + imag2 + imag3
EXPE <~ expe1 + expe2 + expe3 
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
VAL  <~ val1  + val2  + val3

# Reflective measurement model
SAT  =~ sat1  + sat2  + sat3  + sat4
LOY  =~ loy1  + loy2  + loy3  + loy4
"
```

The estimation is conducted using the `csem()` function.

``` r
# Estimate using defaults
res <- csem(.data = satisfaction, .model = model)
res
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ## Estimation was successful.
    ## 
    ## The result is a list of class cSEMResults with list elements:
    ## 
    ##  - Estimates
    ##  - Information
    ## 
    ## To get an overview or help type:
    ## 
    ##  - ?cSEMResults
    ##  - str(<object-name>)
    ##  - listviewer::jsondedit(<object-name>, mode = 'view')
    ## 
    ## If you wish to access the list elements directly type e.g. 
    ## 
    ##  - <object-name>$Estimates
    ## 
    ## Available postestimation commands:
    ## 
    ##  - assess(<object-name>)
    ##  - infer(<object-name)
    ##  - predict(<object-name>)
    ##  - summarize(<object-name>)
    ##  - verify(<object-name>)
    ## ________________________________________________________________________________

This is equal to:

``` r
csem(
   .data                        = satisfaction,
   .model                       = model,
   .approach_cor_robust         = "none",
   .approach_nl                 = "sequential",
   .approach_paths              = "OLS",
   .approach_weights            = "PLS-PM",
   .conv_criterion              = "diff_absolute",
   .disattenuate                = TRUE,
   .dominant_indicators         = NULL,
   .estimate_structural         = TRUE,
   .id                          = NULL,
   .iter_max                    = 100,
   .normality                   = FALSE,
   .PLS_approach_cf             = "dist_squared_euclid",
   .PLS_ignore_structural_model = FALSE,
   .PLS_modes                   = NULL,
   .PLS_weight_scheme_inner     = "path",
   .reliabilities               = NULL,
   .starting_values             = NULL,
   .tolerance                   = 1e-05,
   .resample_method             = "none", 
   .resample_method2            = "none",
   .R                           = 499,
   .R2                          = 199,
   .handle_inadmissibles        = "drop",
   .user_funs                   = NULL,
   .eval_plan                   = "sequential",
   .seed                        = NULL,
   .sign_change_option          = "none"
    )
```

The result is always a named list of class `cSEMResults`.

To access list elements use `$`:

``` r
res$Estimates$Loading_estimates 
res$Information$Model
```

A useful tool to examine a list is the [listviewer
package](https://github.com/timelyportfolio/listviewer). If you are new
to `cSEM` this might be a good way to familiarize yourself with the
structure of a `cSEMResults` object.

``` r
listviewer::jsonedit(res, mode = "view") # requires the listviewer package.
```

Apply postestimation functions:

``` r
## Get a summary
summarize(res) 
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ##  General information:
    ##  ------------------------
    ##  Estimation status                = Ok
    ##  Number of observations           = 250
    ##  Weight estimator                 = PLS-PM
    ##  Inner weighting scheme           = path
    ##  Type of indicator correlation    = Pearson
    ##  Path model estimator             = OLS
    ##  Second order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode 
    ## 
    ##  IMAG  Composite      First order   modeB
    ##  EXPE  Composite      First order   modeB
    ##  QUAL  Composite      First order   modeB
    ##  VAL   Composite      First order   modeB
    ##  SAT   Common factor  First order   modeA
    ##  LOY   Common factor  First order   modeA
    ## 
    ## ----------------------------------- Estimates ----------------------------------
    ## 
    ## Estimated path coefficients:
    ## ============================
    ##   Path           Estimate  Std. error   t-stat.   p-value
    ##   EXPE ~ IMAG      0.4714          NA        NA        NA
    ##   QUAL ~ EXPE      0.8344          NA        NA        NA
    ##   VAL ~ EXPE       0.0457          NA        NA        NA
    ##   VAL ~ QUAL       0.7013          NA        NA        NA
    ##   SAT ~ IMAG       0.2450          NA        NA        NA
    ##   SAT ~ EXPE      -0.0172          NA        NA        NA
    ##   SAT ~ QUAL       0.2215          NA        NA        NA
    ##   SAT ~ VAL        0.5270          NA        NA        NA
    ##   LOY ~ IMAG       0.1819          NA        NA        NA
    ##   LOY ~ SAT        0.6283          NA        NA        NA
    ## 
    ## Estimated loadings:
    ## ===================
    ##   Loading          Estimate  Std. error   t-stat.   p-value
    ##   IMAG =~ imag1      0.6306          NA        NA        NA
    ##   IMAG =~ imag2      0.9246          NA        NA        NA
    ##   IMAG =~ imag3      0.9577          NA        NA        NA
    ##   EXPE =~ expe1      0.7525          NA        NA        NA
    ##   EXPE =~ expe2      0.9348          NA        NA        NA
    ##   EXPE =~ expe3      0.7295          NA        NA        NA
    ##   QUAL =~ qual1      0.7861          NA        NA        NA
    ##   QUAL =~ qual2      0.9244          NA        NA        NA
    ##   QUAL =~ qual3      0.7560          NA        NA        NA
    ##   QUAL =~ qual4      0.7632          NA        NA        NA
    ##   QUAL =~ qual5      0.7834          NA        NA        NA
    ##   VAL =~ val1        0.9518          NA        NA        NA
    ##   VAL =~ val2        0.8056          NA        NA        NA
    ##   VAL =~ val3        0.6763          NA        NA        NA
    ##   SAT =~ sat1        0.9243          NA        NA        NA
    ##   SAT =~ sat2        0.8813          NA        NA        NA
    ##   SAT =~ sat3        0.7127          NA        NA        NA
    ##   SAT =~ sat4        0.7756          NA        NA        NA
    ##   LOY =~ loy1        0.9097          NA        NA        NA
    ##   LOY =~ loy2        0.5775          NA        NA        NA
    ##   LOY =~ loy3        0.9043          NA        NA        NA
    ##   LOY =~ loy4        0.4917          NA        NA        NA
    ## 
    ## Estimated weights:
    ## ==================
    ##   Weights          Estimate  Std. error   t-stat.   p-value
    ##   IMAG <~ imag1      0.0156          NA        NA        NA
    ##   IMAG <~ imag2      0.4473          NA        NA        NA
    ##   IMAG <~ imag3      0.6020          NA        NA        NA
    ##   EXPE <~ expe1      0.2946          NA        NA        NA
    ##   EXPE <~ expe2      0.6473          NA        NA        NA
    ##   EXPE <~ expe3      0.2374          NA        NA        NA
    ##   QUAL <~ qual1      0.2370          NA        NA        NA
    ##   QUAL <~ qual2      0.4712          NA        NA        NA
    ##   QUAL <~ qual3      0.1831          NA        NA        NA
    ##   QUAL <~ qual4      0.1037          NA        NA        NA
    ##   QUAL <~ qual5      0.2049          NA        NA        NA
    ##   VAL <~ val1        0.7163          NA        NA        NA
    ##   VAL <~ val2        0.2202          NA        NA        NA
    ##   VAL <~ val3        0.2082          NA        NA        NA
    ##   SAT <~ sat1        0.3209          NA        NA        NA
    ##   SAT <~ sat2        0.3059          NA        NA        NA
    ##   SAT <~ sat3        0.2474          NA        NA        NA
    ##   SAT <~ sat4        0.2692          NA        NA        NA
    ##   LOY <~ loy1        0.3834          NA        NA        NA
    ##   LOY <~ loy2        0.2434          NA        NA        NA
    ##   LOY <~ loy3        0.3812          NA        NA        NA
    ##   LOY <~ loy4        0.2073          NA        NA        NA
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##   Correlation       Estimate  Std. error   t-stat.   p-value
    ##   imag1 ~~ imag2      0.6437          NA        NA        NA
    ##   imag1 ~~ imag3      0.5433          NA        NA        NA
    ##   imag2 ~~ imag3      0.7761          NA        NA        NA
    ##   expe1 ~~ expe2      0.5353          NA        NA        NA
    ##   expe1 ~~ expe3      0.4694          NA        NA        NA
    ##   expe2 ~~ expe3      0.5467          NA        NA        NA
    ##   qual1 ~~ qual2      0.6053          NA        NA        NA
    ##   qual1 ~~ qual3      0.5406          NA        NA        NA
    ##   qual1 ~~ qual4      0.5662          NA        NA        NA
    ##   qual1 ~~ qual5      0.5180          NA        NA        NA
    ##   qual2 ~~ qual3      0.6187          NA        NA        NA
    ##   qual2 ~~ qual4      0.6517          NA        NA        NA
    ##   qual2 ~~ qual5      0.6291          NA        NA        NA
    ##   qual3 ~~ qual4      0.4752          NA        NA        NA
    ##   qual3 ~~ qual5      0.5074          NA        NA        NA
    ##   qual4 ~~ qual5      0.6402          NA        NA        NA
    ##   val1 ~~ val2        0.6344          NA        NA        NA
    ##   val1 ~~ val3        0.4602          NA        NA        NA
    ##   val2 ~~ val3        0.6288          NA        NA        NA
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##   Total effect    Estimate  Std. error   t-stat.   p-value
    ##   EXPE ~ IMAG       0.4714          NA        NA        NA
    ##   QUAL ~ IMAG       0.3933          NA        NA        NA
    ##   QUAL ~ EXPE       0.8344          NA        NA        NA
    ##   VAL ~ IMAG        0.2974          NA        NA        NA
    ##   VAL ~ EXPE        0.6309          NA        NA        NA
    ##   VAL ~ QUAL        0.7013          NA        NA        NA
    ##   SAT ~ IMAG        0.4807          NA        NA        NA
    ##   SAT ~ EXPE        0.5001          NA        NA        NA
    ##   SAT ~ QUAL        0.5911          NA        NA        NA
    ##   SAT ~ VAL         0.5270          NA        NA        NA
    ##   LOY ~ IMAG        0.4840          NA        NA        NA
    ##   LOY ~ EXPE        0.3142          NA        NA        NA
    ##   LOY ~ QUAL        0.3714          NA        NA        NA
    ##   LOY ~ VAL         0.3311          NA        NA        NA
    ##   LOY ~ SAT         0.6283          NA        NA        NA
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value
    ##   QUAL ~ IMAG          0.3933          NA        NA        NA
    ##   VAL ~ IMAG           0.2974          NA        NA        NA
    ##   VAL ~ EXPE           0.5852          NA        NA        NA
    ##   SAT ~ IMAG           0.2357          NA        NA        NA
    ##   SAT ~ EXPE           0.5173          NA        NA        NA
    ##   SAT ~ QUAL           0.3696          NA        NA        NA
    ##   LOY ~ IMAG           0.3020          NA        NA        NA
    ##   LOY ~ EXPE           0.3142          NA        NA        NA
    ##   LOY ~ QUAL           0.3714          NA        NA        NA
    ##   LOY ~ VAL            0.3311          NA        NA        NA
    ## ________________________________________________________________________________

``` r
## Verify admissibility of the results
verify(res) 
```

    ## ________________________________________________________________________________
    ## 
    ## Verify admissibility:
    ## 
    ##   admissible
    ## 
    ## Details:
    ## 
    ##   Code   Status    Description
    ##   1      ok        Convergence achieved                                   
    ##   2      ok        All absolute standardized loading estimates <= 1       
    ##   3      ok        Construct VCV is positive semi-definite                
    ##   4      ok        All reliability estimates <= 1                         
    ##   5      ok        Model-implied indicator VCV is positive semi-definite  
    ## ________________________________________________________________________________

``` r
## Test overall model fit
testOMF(res, .verbose = FALSE)
```

    ## ________________________________________________________________________________
    ## --------- Test for overall model fit based on Beran & Srivastava (1985) --------
    ## 
    ## Null hypothesis:
    ## 
    ##                                       +------------------------------------------------------------+
    ##                                       |                                                            |
    ##                                       |   H0: Population indicator covariance matrix is equal to   |
    ##                                       |   model-implied indicator covariance matrix.               |
    ##                                       |                                                            |
    ##                                       +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3316  
    ##  SRMR                    0.0940      0.0519  
    ##  dL                      2.2340      0.6819  
    ##  
    ## 
    ## Decision: 
    ## 
    ##                          Significance level
    ##  Distance measure          95%   
    ##  dG                      reject  
    ##  SRMR                    reject  
    ##  dL                      reject  
    ##  
    ## Additonal information:
    ## 
    ##  Out of 499 bootstrap replications 473 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -74457539
    ## ________________________________________________________________________________

``` r
## Assess the model
assess(res)
```

    ## ________________________________________________________________________________
    ## 
    ##  Construct        AVE           R2          R2_adj    
    ##  SAT            0.6851        0.7624        0.7585    
    ##  LOY            0.5552        0.5868        0.5834    
    ## 
    ## -------------- Common (internal consistency) reliability estimates -------------
    ## 
    ##  Construct Cronbachs_alpha   Joereskogs_rho   Dijkstra-Henselers_rho_A 
    ##  SAT        0.8940           0.8960                0.9051          
    ##  LOY        0.8194           0.8237                0.8761          
    ## 
    ## ----------- Alternative (internal consistency) reliability estimates -----------
    ## 
    ##  Construct       RhoC         RhoC_mm    RhoC_weighted
    ##  SAT            0.8938        0.8960        0.9051    
    ##  LOY            0.8011        0.8237        0.8761    
    ## 
    ##  Construct  RhoC_weighted_mm     RhoT      RhoT_weighted
    ##  SAT            0.9051        0.8940        0.8869    
    ##  LOY            0.8761        0.8194        0.7850    
    ## 
    ## --------------------------- Distance and fit measures --------------------------
    ## 
    ##  Geodesic distance           = 0.6493432
    ##  Squared Euclidian distance  = 2.23402
    ##  ML distance                 = 2.921932
    ## 
    ##  CFI          = 0.8573048
    ##  GFI          = 0.9642375
    ##  IFI          = 0.8593711
    ##  NFI          = 0.8229918
    ##  NNFI         = 0.8105598
    ##  RMSEA        = 0.1130338
    ##  RMS_theta    = 0.05069299
    ##  SRMR         = 0.09396871
    ## 
    ##  Degrees of freedom    = 174
    ## 
    ## ----------------------- Variance inflation factors (VIFs) ----------------------
    ## 
    ##   Dependent construct: 'VAL'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      3.2928   
    ##  QUAL                      3.2928   
    ##  IMAG                      0.0000   
    ##  VAL                       0.0000   
    ##  SAT                       0.0000   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      3.2985   
    ##  QUAL                      4.4151   
    ##  IMAG                      1.7280   
    ##  VAL                       2.6726   
    ##  SAT                       0.0000   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      0.0000   
    ##  QUAL                      0.0000   
    ##  IMAG                      1.9345   
    ##  VAL                       0.0000   
    ##  SAT                       1.9345   
    ## 
    ## -------------------------- Effect sizes (Cohen's f^2) --------------------------
    ## 
    ##   Dependent construct: 'EXPE'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      0.2856   
    ## 
    ##   Dependent construct: 'QUAL'
    ## 
    ##  Independent construct   Effect size
    ##  EXPE                      2.2928   
    ## 
    ##   Dependent construct: 'VAL'
    ## 
    ##  Independent construct   Effect size
    ##  EXPE                      1.2097   
    ##  QUAL                      1.2097   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      3.2086   
    ##  EXPE                      3.2086   
    ##  QUAL                      3.2086   
    ##  VAL                       3.2086   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      1.4199   
    ##  SAT                       1.4199   
    ## 
    ## ------------------------------ Validity assessment -----------------------------
    ## 
    ##  Heterotrait-montrait ratio of correlation matrix (HTMT matrix)
    ## 
    ##           SAT LOY
    ## SAT 0.0000000   0
    ## LOY 0.7432489   0
    ## 
    ## 
    ##  Fornell-Larcker matrix
    ## 
    ##           SAT       LOY
    ## SAT 0.6851491 0.5696460
    ## LOY 0.5696460 0.5551718
    ## 
    ## ________________________________________________________________________________

``` r
## Predict indicator scores of endogenous constructs
predict(res)
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ##  Number of obs. training          = 225
    ##  Number of obs. test              = 25
    ##  Number of cv folds               = 10
    ##  Number of repetitions            = 10
    ##  Handle inadmissibles             = stop
    ##  Target                           = 'PLS-PM'
    ##  Benchmark                        = 'lm'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name     MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   sat1         1.3523         1.2334       1.7893         1.6205       0.2208
    ##   sat2         1.3079         1.1980       1.7682         1.6292       0.2000
    ##   sat3         1.4093         1.2810       1.7518         1.7276       0.1346
    ##   sat4         1.4157         1.2649       1.7827         1.6386       0.1745
    ##   loy1         1.7740         1.6638       2.2973         2.2325       0.2268
    ##   loy2         1.5159         1.4735       1.9362         1.9786       0.1091
    ##   loy3         1.7819         1.6702       2.3455         2.2736       0.2280
    ##   loy4         1.7114         1.6781       2.1973         2.3119       0.0702
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential quantities are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways to compute resamples:

1.  Setting `.resample_method` in `csem()` to `"jackkinfe"` or
    `"bootstrap"` and subsequently using postestimation functions
    `summarize()` or `infer()`.
2.  The same result is achieved by passing a `cSEMResults` object to
    `resamplecSEMResults()` and subsequently using postestimation
    functions `summarize()` or `infer()`.

<!-- end list -->

``` r
# Setting `.resample_method`
b1 <- csem(.data = satisfaction, .model = model, .resample_method = "bootstrap")
# Using resamplecSEMResults()
b2 <- resamplecSEMResults(res)
```

Now `summarize()` shows inferencial quantities as well:

``` r
summarize(b1)
```

    ## ________________________________________________________________________________
    ## ----------------------------------- Overview -----------------------------------
    ## 
    ##  General information:
    ##  ------------------------
    ##  Estimation status                = Ok
    ##  Number of observations           = 250
    ##  Weight estimator                 = PLS-PM
    ##  Inner weighting scheme           = path
    ##  Type of indicator correlation    = Pearson
    ##  Path model estimator             = OLS
    ##  Second order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample methode                 = bootstrap
    ##  Number of resamples              = 499
    ##  Number of admissible results     = 484
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 1117995718
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode 
    ## 
    ##  IMAG  Composite      First order   modeB
    ##  EXPE  Composite      First order   modeB
    ##  QUAL  Composite      First order   modeB
    ##  VAL   Composite      First order   modeB
    ##  SAT   Common factor  First order   modeA
    ##  LOY   Common factor  First order   modeA
    ## 
    ## ----------------------------------- Estimates ----------------------------------
    ## 
    ## Estimated path coefficients:
    ## ============================
    ##                                                              CI_percentile   
    ##   Path           Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG      0.4714      0.0658    7.1623    0.0000 [ 0.3375; 0.5968 ] 
    ##   QUAL ~ EXPE      0.8344      0.0236   35.3226    0.0000 [ 0.7848; 0.8750 ] 
    ##   VAL ~ EXPE       0.0457      0.0832    0.5495    0.5826 [-0.1040; 0.2203 ] 
    ##   VAL ~ QUAL       0.7013      0.0816    8.5898    0.0000 [ 0.5295; 0.8544 ] 
    ##   SAT ~ IMAG       0.2450      0.0535    4.5767    0.0000 [ 0.1431; 0.3502 ] 
    ##   SAT ~ EXPE      -0.0172      0.0699   -0.2467    0.8052 [-0.1679; 0.1049 ] 
    ##   SAT ~ QUAL       0.2215      0.1056    2.0970    0.0360 [ 0.0302; 0.4244 ] 
    ##   SAT ~ VAL        0.5270      0.0860    6.1240    0.0000 [ 0.3643; 0.6847 ] 
    ##   LOY ~ IMAG       0.1819      0.0829    2.1935    0.0283 [ 0.0142; 0.3280 ] 
    ##   LOY ~ SAT        0.6283      0.0827    7.5931    0.0000 [ 0.4788; 0.7998 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0967    6.5186    0.0000 [ 0.4325; 0.8018 ] 
    ##   IMAG =~ imag2      0.9246      0.0369   25.0679    0.0000 [ 0.8330; 0.9750 ] 
    ##   IMAG =~ imag3      0.9577      0.0274   34.9078    0.0000 [ 0.8886; 0.9905 ] 
    ##   EXPE =~ expe1      0.7525      0.0775    9.7127    0.0000 [ 0.5734; 0.8787 ] 
    ##   EXPE =~ expe2      0.9348      0.0280   33.4261    0.0000 [ 0.8677; 0.9737 ] 
    ##   EXPE =~ expe3      0.7295      0.0727   10.0394    0.0000 [ 0.5704; 0.8514 ] 
    ##   QUAL =~ qual1      0.7861      0.0660   11.9161    0.0000 [ 0.6314; 0.8876 ] 
    ##   QUAL =~ qual2      0.9244      0.0240   38.5465    0.0000 [ 0.8761; 0.9598 ] 
    ##   QUAL =~ qual3      0.7560      0.0607   12.4643    0.0000 [ 0.6131; 0.8533 ] 
    ##   QUAL =~ qual4      0.7632      0.0513   14.8814    0.0000 [ 0.6518; 0.8560 ] 
    ##   QUAL =~ qual5      0.7834      0.0474   16.5310    0.0000 [ 0.6842; 0.8741 ] 
    ##   VAL =~ val1        0.9518      0.0215   44.2078    0.0000 [ 0.8991; 0.9838 ] 
    ##   VAL =~ val2        0.8056      0.0645   12.4851    0.0000 [ 0.6656; 0.8979 ] 
    ##   VAL =~ val3        0.6763      0.0700    9.6578    0.0000 [ 0.5283; 0.8000 ] 
    ##   SAT =~ sat1        0.9243      0.0233   39.7037    0.0000 [ 0.8734; 0.9639 ] 
    ##   SAT =~ sat2        0.8813      0.0284   31.0861    0.0000 [ 0.8208; 0.9281 ] 
    ##   SAT =~ sat3        0.7127      0.0544   13.0998    0.0000 [ 0.5928; 0.8089 ] 
    ##   SAT =~ sat4        0.7756      0.0498   15.5895    0.0000 [ 0.6707; 0.8618 ] 
    ##   LOY =~ loy1        0.9097      0.0499   18.2357    0.0000 [ 0.8004; 0.9885 ] 
    ##   LOY =~ loy2        0.5775      0.0821    7.0313    0.0000 [ 0.4058; 0.7239 ] 
    ##   LOY =~ loy3        0.9043      0.0425   21.2704    0.0000 [ 0.8084; 0.9770 ] 
    ##   LOY =~ loy4        0.4917      0.0906    5.4269    0.0000 [ 0.3409; 0.6702 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1111    0.1408    0.8880 [-0.2023; 0.2386 ] 
    ##   IMAG <~ imag2      0.4473      0.1405    3.1842    0.0015 [ 0.1750; 0.6994 ] 
    ##   IMAG <~ imag3      0.6020      0.1312    4.5872    0.0000 [ 0.3093; 0.8346 ] 
    ##   EXPE <~ expe1      0.2946      0.1141    2.5809    0.0099 [ 0.0651; 0.5137 ] 
    ##   EXPE <~ expe2      0.6473      0.0885    7.3109    0.0000 [ 0.4601; 0.7998 ] 
    ##   EXPE <~ expe3      0.2374      0.0935    2.5392    0.0111 [ 0.0690; 0.4299 ] 
    ##   QUAL <~ qual1      0.2370      0.0857    2.7658    0.0057 [ 0.0726; 0.4133 ] 
    ##   QUAL <~ qual2      0.4712      0.0839    5.6158    0.0000 [ 0.3050; 0.6304 ] 
    ##   QUAL <~ qual3      0.1831      0.0777    2.3561    0.0185 [ 0.0212; 0.3285 ] 
    ##   QUAL <~ qual4      0.1037      0.0623    1.6644    0.0960 [ 0.0039; 0.2465 ] 
    ##   QUAL <~ qual5      0.2049      0.0630    3.2496    0.0012 [ 0.0736; 0.3228 ] 
    ##   VAL <~ val1        0.7163      0.0900    7.9557    0.0000 [ 0.5387; 0.8693 ] 
    ##   VAL <~ val2        0.2202      0.0865    2.5472    0.0109 [ 0.0545; 0.3997 ] 
    ##   VAL <~ val3        0.2082      0.0582    3.5783    0.0003 [ 0.1029; 0.3245 ] 
    ##   SAT <~ sat1        0.3209      0.0154   20.7924    0.0000 [ 0.2951; 0.3551 ] 
    ##   SAT <~ sat2        0.3059      0.0139   22.0595    0.0000 [ 0.2829; 0.3360 ] 
    ##   SAT <~ sat3        0.2474      0.0115   21.6021    0.0000 [ 0.2246; 0.2685 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.1505    0.0000 [ 0.2466; 0.2941 ] 
    ##   LOY <~ loy1        0.3834      0.0242   15.8347    0.0000 [ 0.3365; 0.4284 ] 
    ##   LOY <~ loy2        0.2434      0.0298    8.1550    0.0000 [ 0.1846; 0.2942 ] 
    ##   LOY <~ loy3        0.3812      0.0255   14.9227    0.0000 [ 0.3309; 0.4311 ] 
    ##   LOY <~ loy4        0.2073      0.0340    6.0892    0.0000 [ 0.1448; 0.2744 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0670    9.6062    0.0000 [ 0.4934; 0.7579 ] 
    ##   imag1 ~~ imag3      0.5433      0.0711    7.6382    0.0000 [ 0.3957; 0.6858 ] 
    ##   imag2 ~~ imag3      0.7761      0.0386   20.1151    0.0000 [ 0.6992; 0.8438 ] 
    ##   expe1 ~~ expe2      0.5353      0.0586    9.1266    0.0000 [ 0.4211; 0.6335 ] 
    ##   expe1 ~~ expe3      0.4694      0.0640    7.3378    0.0000 [ 0.3380; 0.5931 ] 
    ##   expe2 ~~ expe3      0.5467      0.0603    9.0690    0.0000 [ 0.4228; 0.6567 ] 
    ##   qual1 ~~ qual2      0.6053      0.0555   10.8997    0.0000 [ 0.5015; 0.7026 ] 
    ##   qual1 ~~ qual3      0.5406      0.0633    8.5404    0.0000 [ 0.4112; 0.6579 ] 
    ##   qual1 ~~ qual4      0.5662      0.0647    8.7524    0.0000 [ 0.4339; 0.6837 ] 
    ##   qual1 ~~ qual5      0.5180      0.0687    7.5374    0.0000 [ 0.3816; 0.6425 ] 
    ##   qual2 ~~ qual3      0.6187      0.0557   11.1039    0.0000 [ 0.4992; 0.7155 ] 
    ##   qual2 ~~ qual4      0.6517      0.0599   10.8854    0.0000 [ 0.5292; 0.7650 ] 
    ##   qual2 ~~ qual5      0.6291      0.0587   10.7235    0.0000 [ 0.5063; 0.7396 ] 
    ##   qual3 ~~ qual4      0.4752      0.0635    7.4792    0.0000 [ 0.3419; 0.5900 ] 
    ##   qual3 ~~ qual5      0.5074      0.0616    8.2430    0.0000 [ 0.3860; 0.6262 ] 
    ##   qual4 ~~ qual5      0.6402      0.0545   11.7448    0.0000 [ 0.5327; 0.7339 ] 
    ##   val1 ~~ val2        0.6344      0.0595   10.6662    0.0000 [ 0.5166; 0.7418 ] 
    ##   val1 ~~ val3        0.4602      0.0696    6.6141    0.0000 [ 0.3262; 0.5909 ] 
    ##   val2 ~~ val3        0.6288      0.0637    9.8680    0.0000 [ 0.5141; 0.7478 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0658    7.1623    0.0000 [ 0.3375; 0.5968 ] 
    ##   QUAL ~ IMAG       0.3933      0.0605    6.4991    0.0000 [ 0.2728; 0.5086 ] 
    ##   QUAL ~ EXPE       0.8344      0.0236   35.3226    0.0000 [ 0.7848; 0.8750 ] 
    ##   VAL ~ IMAG        0.2974      0.0593    5.0185    0.0000 [ 0.1872; 0.4165 ] 
    ##   VAL ~ EXPE        0.6309      0.0490   12.8673    0.0000 [ 0.5336; 0.7240 ] 
    ##   VAL ~ QUAL        0.7013      0.0816    8.5898    0.0000 [ 0.5295; 0.8544 ] 
    ##   SAT ~ IMAG        0.4807      0.0639    7.5252    0.0000 [ 0.3652; 0.6018 ] 
    ##   SAT ~ EXPE        0.5001      0.0566    8.8391    0.0000 [ 0.3934; 0.6061 ] 
    ##   SAT ~ QUAL        0.5911      0.0937    6.3087    0.0000 [ 0.4149; 0.7619 ] 
    ##   SAT ~ VAL         0.5270      0.0860    6.1240    0.0000 [ 0.3643; 0.6847 ] 
    ##   LOY ~ IMAG        0.4840      0.0709    6.8239    0.0000 [ 0.3424; 0.6215 ] 
    ##   LOY ~ EXPE        0.3142      0.0544    5.7733    0.0000 [ 0.2131; 0.4232 ] 
    ##   LOY ~ QUAL        0.3714      0.0848    4.3789    0.0000 [ 0.2315; 0.5571 ] 
    ##   LOY ~ VAL         0.3311      0.0758    4.3664    0.0000 [ 0.1997; 0.4988 ] 
    ##   LOY ~ SAT         0.6283      0.0827    7.5931    0.0000 [ 0.4788; 0.7998 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0605    6.4991    0.0000 [ 0.2728; 0.5086 ] 
    ##   VAL ~ IMAG           0.2974      0.0593    5.0185    0.0000 [ 0.1872; 0.4165 ] 
    ##   VAL ~ EXPE           0.5852      0.0705    8.3001    0.0000 [ 0.4398; 0.7158 ] 
    ##   SAT ~ IMAG           0.2357      0.0476    4.9472    0.0000 [ 0.1454; 0.3314 ] 
    ##   SAT ~ EXPE           0.5173      0.0672    7.7011    0.0000 [ 0.3937; 0.6630 ] 
    ##   SAT ~ QUAL           0.3696      0.0636    5.8110    0.0000 [ 0.2472; 0.4898 ] 
    ##   LOY ~ IMAG           0.3020      0.0554    5.4520    0.0000 [ 0.2091; 0.4181 ] 
    ##   LOY ~ EXPE           0.3142      0.0544    5.7733    0.0000 [ 0.2131; 0.4232 ] 
    ##   LOY ~ QUAL           0.3714      0.0848    4.3789    0.0000 [ 0.2315; 0.5571 ] 
    ##   LOY ~ VAL            0.3311      0.0758    4.3664    0.0000 [ 0.1997; 0.4988 ] 
    ## ________________________________________________________________________________

Several resample-based confidence intervals are implemented, see
`?infer()`:

``` r
infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet
```

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as setting random seeds via the [future
framework](https://github.com/HenrikBengtsson/future). For
multiprocessing simply set `.eval_plan = "multiprocess"` in which case
the maximum number of available cores is used if not on Windows. On
Windows as many separate R instances are opened in the backround as
there are cores available instead. Note that this naturally has some
overhead so for a small number of resamples multiprocessing will not
always be faster compared to sequential (single core) processing (the
default). Seeds are set via the `.seed` argument.

``` r
b <- csem(
  .data            = satisfaction,
  .model           = model, 
  .resample_method = "bootstrap",
  .R               = 999,
  .seed            = 98234,
  .eval_plan       = "multiprocess")
```
