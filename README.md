
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![Build
Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)

WARNING: THIS IS WORK IN PROGRESS. BREAKING CHANGES TO THE API ARE VERY
LIKELY. Use the package with caution and please report bugs to [the
package
developers](mailto:manuel.rademaker@uni-wuerzburg.de;f.schuberth@utwente.nl).
The first stable relase will be version 0.0.1, most likely towards the
end of 2019.

## Purpose

Estimate, analyse, test, and study linear, nonlinear, hierachical
structural equation models using composite-based approaches and
procedures, including estimation techniques such as partial least
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

``` r
# Currently only a development version from GitHub is available:
# install.packages("devtools")
devtools::install_github("M-E-Rademaker/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://m-e-rademaker.github.io/cSEM/) (work in
progress).

## Philosophy

  - First and foremost: `cSEM` has a user-centered design\!
  - “User-centered” mainly boils down to: `cSEM` is easy, i.e. intuitive
    to use by non-R experts\!
    <!--  - There is one central function called `csem()` that provides default choices -->
    <!--    for most of its arguments (similarity to the `sem()` and `cfa()` functions of the [lavaan](http://lavaan.ugent.be/)  -->
    <!--    package is intended). --> <!-- -  -->
    <!--  - cSEM is Well documented (vignettes, HTML output, a website, (eventually) intro course(s) and cheatsheets) -->
    <!--  - Structured output/results  that aims to be "easy"" in a sense that it is -->
    <!--      - ... descriptive/verbose -->
    <!--      - ... (eventually) easy to export to other environments such as MS Word, Latex files etc. (exportability) -->
    <!--      - ... (eventually) easy to migrate from/to/between other PLS/VB/CB-based systems (lavaan, semPLS, ADANCO, SmartPLS) -->
  - The package is designed to be flexible/modular enough so that
    researchers developing new methods can take specific function
    provided by the package and alter them according to their need
    without working their way through a chain of other functions
    (naturally this will not always be possible). Modularity is largly
    inspired by the `matrixpls` package.
  - Modern in a sense that the package integrates modern developments
    within the R community. This mainly includes
    ideas/recommendations/design choices that fead into the packages of
    the [tidyverse](https://github.com/tidyverse/tidyverse).
  - State of the art in a sense that we seek to quickly implement recent
    methodological developments in composite-based SEM.

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

Currently we have five major postestimation verbs, four test family
functions and two do-family of function:

  - `assess()` : assess the model using common quality criteria
  - `infer()` : calculate common inferencial quantities (e.g, standard
    errors)
  - `predict()` : predict endogenous indicator values
  - `summarize()` : summarize the results
  - `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently the
following tests are implemented:

  - `testOMF()` : performs a test for overall model fit
  - `testMICOM()` : performs a test for composite measurement invariance
  - `testMGD` : performs several test to assess multi-group differences
  - `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity.

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

A usefule tool to examine a list is the [listviewer
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
    ##                                                            +------------------------------------------------------------+
    ##                                                            |                                                            |
    ##                                                            |   H0: Population indicator covariance matrix is equal to   |
    ##                                                            |   model-implied indicator covariance matrix.               |
    ##                                                            |                                                            |
    ##                                                            +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3123  
    ##  SRMR                    0.0940      0.0524  
    ##  dL                      2.2340      0.6935  
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
    ##  Out of 499 bootstrap replications 478 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 58379673
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
    ##  EXPE                      0.0014   
    ##  QUAL                      0.3301   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      0.1462   
    ##  EXPE                      0.0004   
    ##  QUAL                      0.0468   
    ##  VAL                       0.4373   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      0.0414   
    ##  SAT                       0.4938   
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
    ##   sat1         1.3516         1.2364       1.7892         1.6230       0.2211
    ##   sat2         1.3070         1.1979       1.7662         1.6306       0.2004
    ##   sat3         1.4084         1.2785       1.7514         1.7249       0.1340
    ##   sat4         1.4159         1.2642       1.7830         1.6371       0.1753
    ##   loy1         1.7744         1.6596       2.2998         2.2254       0.2261
    ##   loy2         1.5168         1.4765       1.9397         1.9838       0.1081
    ##   loy3         1.7847         1.6677       2.3499         2.2680       0.2273
    ##   loy4         1.7098         1.6692       2.1939         2.2970       0.0741
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
    ##  Number of admissible results     = 485
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 853658205
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
    ##   EXPE ~ IMAG      0.4714      0.0663    7.1092    0.0000 [ 0.3453; 0.6059 ] 
    ##   QUAL ~ EXPE      0.8344      0.0232   35.9715    0.0000 [ 0.7850; 0.8715 ] 
    ##   VAL ~ EXPE       0.0457      0.0794    0.5755    0.5650 [-0.1012; 0.2217 ] 
    ##   VAL ~ QUAL       0.7013      0.0789    8.8836    0.0000 [ 0.5401; 0.8605 ] 
    ##   SAT ~ IMAG       0.2450      0.0558    4.3893    0.0000 [ 0.1384; 0.3500 ] 
    ##   SAT ~ EXPE      -0.0172      0.0702   -0.2456    0.8060 [-0.1560; 0.0989 ] 
    ##   SAT ~ QUAL       0.2215      0.1021    2.1695    0.0300 [ 0.0259; 0.4401 ] 
    ##   SAT ~ VAL        0.5270      0.0839    6.2827    0.0000 [ 0.3501; 0.6977 ] 
    ##   LOY ~ IMAG       0.1819      0.0785    2.3180    0.0205 [ 0.0451; 0.3360 ] 
    ##   LOY ~ SAT        0.6283      0.0793    7.9190    0.0000 [ 0.4849; 0.7769 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0968    6.5118    0.0000 [ 0.4220; 0.7825 ] 
    ##   IMAG =~ imag2      0.9246      0.0402   23.0008    0.0000 [ 0.8242; 0.9781 ] 
    ##   IMAG =~ imag3      0.9577      0.0284   33.6931    0.0000 [ 0.8903; 0.9911 ] 
    ##   EXPE =~ expe1      0.7525      0.0805    9.3455    0.0000 [ 0.5605; 0.8712 ] 
    ##   EXPE =~ expe2      0.9348      0.0288   32.5081    0.0000 [ 0.8634; 0.9737 ] 
    ##   EXPE =~ expe3      0.7295      0.0742    9.8310    0.0000 [ 0.5661; 0.8507 ] 
    ##   QUAL =~ qual1      0.7861      0.0708   11.1096    0.0000 [ 0.6148; 0.8881 ] 
    ##   QUAL =~ qual2      0.9244      0.0240   38.5477    0.0000 [ 0.8659; 0.9575 ] 
    ##   QUAL =~ qual3      0.7560      0.0633   11.9351    0.0000 [ 0.6134; 0.8683 ] 
    ##   QUAL =~ qual4      0.7632      0.0532   14.3490    0.0000 [ 0.6386; 0.8518 ] 
    ##   QUAL =~ qual5      0.7834      0.0475   16.4975    0.0000 [ 0.6769; 0.8569 ] 
    ##   VAL =~ val1        0.9518      0.0241   39.4270    0.0000 [ 0.8957; 0.9868 ] 
    ##   VAL =~ val2        0.8056      0.0654   12.3104    0.0000 [ 0.6665; 0.9138 ] 
    ##   VAL =~ val3        0.6763      0.0708    9.5558    0.0000 [ 0.5194; 0.8026 ] 
    ##   SAT =~ sat1        0.9243      0.0217   42.6398    0.0000 [ 0.8786; 0.9623 ] 
    ##   SAT =~ sat2        0.8813      0.0299   29.4698    0.0000 [ 0.8222; 0.9325 ] 
    ##   SAT =~ sat3        0.7127      0.0508   14.0173    0.0000 [ 0.6088; 0.8082 ] 
    ##   SAT =~ sat4        0.7756      0.0509   15.2327    0.0000 [ 0.6647; 0.8679 ] 
    ##   LOY =~ loy1        0.9097      0.0516   17.6290    0.0000 [ 0.7905; 0.9872 ] 
    ##   LOY =~ loy2        0.5775      0.0851    6.7889    0.0000 [ 0.4047; 0.7388 ] 
    ##   LOY =~ loy3        0.9043      0.0423   21.3562    0.0000 [ 0.7971; 0.9708 ] 
    ##   LOY =~ loy4        0.4917      0.0970    5.0710    0.0000 [ 0.3254; 0.6853 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1163    0.1345    0.8930 [-0.2171; 0.2432 ] 
    ##   IMAG <~ imag2      0.4473      0.1492    2.9975    0.0027 [ 0.1613; 0.7182 ] 
    ##   IMAG <~ imag3      0.6020      0.1399    4.3032    0.0000 [ 0.3309; 0.8471 ] 
    ##   EXPE <~ expe1      0.2946      0.1186    2.4848    0.0130 [ 0.0435; 0.5130 ] 
    ##   EXPE <~ expe2      0.6473      0.0872    7.4226    0.0000 [ 0.4490; 0.7977 ] 
    ##   EXPE <~ expe3      0.2374      0.0956    2.4824    0.0130 [ 0.0480; 0.4247 ] 
    ##   QUAL <~ qual1      0.2370      0.0888    2.6705    0.0076 [ 0.0790; 0.4208 ] 
    ##   QUAL <~ qual2      0.4712      0.0819    5.7547    0.0000 [ 0.2922; 0.6148 ] 
    ##   QUAL <~ qual3      0.1831      0.0843    2.1717    0.0299 [ 0.0256; 0.3528 ] 
    ##   QUAL <~ qual4      0.1037      0.0612    1.6959    0.0899 [-0.0076; 0.2301 ] 
    ##   QUAL <~ qual5      0.2049      0.0631    3.2451    0.0012 [ 0.0740; 0.3189 ] 
    ##   VAL <~ val1        0.7163      0.0987    7.2540    0.0000 [ 0.4915; 0.8841 ] 
    ##   VAL <~ val2        0.2202      0.0984    2.2380    0.0252 [ 0.0475; 0.4306 ] 
    ##   VAL <~ val3        0.2082      0.0607    3.4286    0.0006 [ 0.0849; 0.3217 ] 
    ##   SAT <~ sat1        0.3209      0.0154   20.8078    0.0000 [ 0.2928; 0.3542 ] 
    ##   SAT <~ sat2        0.3059      0.0139   22.0880    0.0000 [ 0.2813; 0.3378 ] 
    ##   SAT <~ sat3        0.2474      0.0102   24.2382    0.0000 [ 0.2260; 0.2652 ] 
    ##   SAT <~ sat4        0.2692      0.0120   22.4466    0.0000 [ 0.2451; 0.2939 ] 
    ##   LOY <~ loy1        0.3834      0.0281   13.6358    0.0000 [ 0.3293; 0.4363 ] 
    ##   LOY <~ loy2        0.2434      0.0287    8.4948    0.0000 [ 0.1810; 0.2913 ] 
    ##   LOY <~ loy3        0.3812      0.0278   13.7047    0.0000 [ 0.3270; 0.4381 ] 
    ##   LOY <~ loy4        0.2073      0.0358    5.7868    0.0000 [ 0.1447; 0.2783 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0646    9.9641    0.0000 [ 0.5119; 0.7507 ] 
    ##   imag1 ~~ imag3      0.5433      0.0679    8.0056    0.0000 [ 0.4055; 0.6663 ] 
    ##   imag2 ~~ imag3      0.7761      0.0396   19.5901    0.0000 [ 0.6944; 0.8473 ] 
    ##   expe1 ~~ expe2      0.5353      0.0626    8.5560    0.0000 [ 0.4032; 0.6491 ] 
    ##   expe1 ~~ expe3      0.4694      0.0611    7.6813    0.0000 [ 0.3563; 0.5901 ] 
    ##   expe2 ~~ expe3      0.5467      0.0606    9.0270    0.0000 [ 0.4173; 0.6541 ] 
    ##   qual1 ~~ qual2      0.6053      0.0587   10.3070    0.0000 [ 0.4781; 0.7092 ] 
    ##   qual1 ~~ qual3      0.5406      0.0600    9.0173    0.0000 [ 0.4174; 0.6506 ] 
    ##   qual1 ~~ qual4      0.5662      0.0670    8.4529    0.0000 [ 0.4138; 0.6849 ] 
    ##   qual1 ~~ qual5      0.5180      0.0677    7.6524    0.0000 [ 0.3731; 0.6470 ] 
    ##   qual2 ~~ qual3      0.6187      0.0553   11.1916    0.0000 [ 0.5068; 0.7220 ] 
    ##   qual2 ~~ qual4      0.6517      0.0623   10.4654    0.0000 [ 0.5077; 0.7617 ] 
    ##   qual2 ~~ qual5      0.6291      0.0604   10.4109    0.0000 [ 0.5044; 0.7315 ] 
    ##   qual3 ~~ qual4      0.4752      0.0615    7.7316    0.0000 [ 0.3477; 0.5856 ] 
    ##   qual3 ~~ qual5      0.5074      0.0625    8.1229    0.0000 [ 0.3798; 0.6214 ] 
    ##   qual4 ~~ qual5      0.6402      0.0531   12.0592    0.0000 [ 0.5316; 0.7269 ] 
    ##   val1 ~~ val2        0.6344      0.0545   11.6481    0.0000 [ 0.5299; 0.7380 ] 
    ##   val1 ~~ val3        0.4602      0.0647    7.1146    0.0000 [ 0.3436; 0.5879 ] 
    ##   val2 ~~ val3        0.6288      0.0612   10.2732    0.0000 [ 0.5078; 0.7440 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0663    7.1092    0.0000 [ 0.3453; 0.6059 ] 
    ##   QUAL ~ IMAG       0.3933      0.0616    6.3875    0.0000 [ 0.2771; 0.5267 ] 
    ##   QUAL ~ EXPE       0.8344      0.0232   35.9715    0.0000 [ 0.7850; 0.8715 ] 
    ##   VAL ~ IMAG        0.2974      0.0603    4.9359    0.0000 [ 0.1885; 0.4295 ] 
    ##   VAL ~ EXPE        0.6309      0.0490   12.8811    0.0000 [ 0.5316; 0.7211 ] 
    ##   VAL ~ QUAL        0.7013      0.0789    8.8836    0.0000 [ 0.5401; 0.8605 ] 
    ##   SAT ~ IMAG        0.4807      0.0699    6.8794    0.0000 [ 0.3474; 0.6156 ] 
    ##   SAT ~ EXPE        0.5001      0.0576    8.6895    0.0000 [ 0.3878; 0.6096 ] 
    ##   SAT ~ QUAL        0.5911      0.0904    6.5361    0.0000 [ 0.4081; 0.7756 ] 
    ##   SAT ~ VAL         0.5270      0.0839    6.2827    0.0000 [ 0.3501; 0.6977 ] 
    ##   LOY ~ IMAG        0.4840      0.0677    7.1438    0.0000 [ 0.3613; 0.6244 ] 
    ##   LOY ~ EXPE        0.3142      0.0524    6.0000    0.0000 [ 0.2151; 0.4169 ] 
    ##   LOY ~ QUAL        0.3714      0.0770    4.8234    0.0000 [ 0.2423; 0.5388 ] 
    ##   LOY ~ VAL         0.3311      0.0728    4.5456    0.0000 [ 0.2002; 0.4766 ] 
    ##   LOY ~ SAT         0.6283      0.0793    7.9190    0.0000 [ 0.4849; 0.7769 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0616    6.3875    0.0000 [ 0.2771; 0.5267 ] 
    ##   VAL ~ IMAG           0.2974      0.0603    4.9359    0.0000 [ 0.1885; 0.4295 ] 
    ##   VAL ~ EXPE           0.5852      0.0682    8.5782    0.0000 [ 0.4435; 0.7234 ] 
    ##   SAT ~ IMAG           0.2357      0.0474    4.9766    0.0000 [ 0.1527; 0.3368 ] 
    ##   SAT ~ EXPE           0.5173      0.0646    8.0141    0.0000 [ 0.4015; 0.6543 ] 
    ##   SAT ~ QUAL           0.3696      0.0584    6.3269    0.0000 [ 0.2467; 0.4862 ] 
    ##   LOY ~ IMAG           0.3020      0.0565    5.3459    0.0000 [ 0.2088; 0.4314 ] 
    ##   LOY ~ EXPE           0.3142      0.0524    6.0000    0.0000 [ 0.2151; 0.4169 ] 
    ##   LOY ~ QUAL           0.3714      0.0770    4.8234    0.0000 [ 0.2423; 0.5388 ] 
    ##   LOY ~ VAL            0.3311      0.0728    4.5456    0.0000 [ 0.2002; 0.4766 ] 
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
