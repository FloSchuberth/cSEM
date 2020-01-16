
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
  - `testMGD()` : performs several test to assess multi-group
    differences
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
    ##  dG                      0.6493      0.3257  
    ##  SRMR                    0.0940      0.0516  
    ##  dL                      2.2340      0.6732  
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
    ##  The seed used was: 1549311926
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
    ##   sat1         1.3491         1.2338       1.7848         1.6198       0.2231
    ##   sat2         1.3071         1.1969       1.7633         1.6290       0.2023
    ##   sat3         1.4094         1.2764       1.7528         1.7216       0.1341
    ##   sat4         1.4155         1.2641       1.7838         1.6375       0.1732
    ##   loy1         1.7714         1.6615       2.2939         2.2291       0.2281
    ##   loy2         1.5168         1.4744       1.9384         1.9836       0.1071
    ##   loy3         1.7805         1.6701       2.3443         2.2737       0.2286
    ##   loy4         1.7094         1.6751       2.1962         2.3043       0.0705
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
    ##  Number of admissible results     = 489
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -458190711
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
    ##   EXPE ~ IMAG      0.4714      0.0661    7.1319    0.0000 [ 0.3449; 0.5974 ] 
    ##   QUAL ~ EXPE      0.8344      0.0234   35.7083    0.0000 [ 0.7858; 0.8726 ] 
    ##   VAL ~ EXPE       0.0457      0.0864    0.5293    0.5966 [-0.1022; 0.2341 ] 
    ##   VAL ~ QUAL       0.7013      0.0838    8.3739    0.0000 [ 0.5332; 0.8546 ] 
    ##   SAT ~ IMAG       0.2450      0.0543    4.5081    0.0000 [ 0.1347; 0.3459 ] 
    ##   SAT ~ EXPE      -0.0172      0.0739   -0.2333    0.8155 [-0.1607; 0.1240 ] 
    ##   SAT ~ QUAL       0.2215      0.1009    2.1955    0.0281 [ 0.0453; 0.4323 ] 
    ##   SAT ~ VAL        0.5270      0.0871    6.0468    0.0000 [ 0.3523; 0.6874 ] 
    ##   LOY ~ IMAG       0.1819      0.0828    2.1979    0.0280 [ 0.0353; 0.3464 ] 
    ##   LOY ~ SAT        0.6283      0.0841    7.4682    0.0000 [ 0.4647; 0.7975 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1022    6.1679    0.0000 [ 0.4061; 0.8071 ] 
    ##   IMAG =~ imag2      0.9246      0.0419   22.0798    0.0000 [ 0.8162; 0.9774 ] 
    ##   IMAG =~ imag3      0.9577      0.0306   31.2670    0.0000 [ 0.8701; 0.9926 ] 
    ##   EXPE =~ expe1      0.7525      0.0842    8.9403    0.0000 [ 0.5660; 0.8745 ] 
    ##   EXPE =~ expe2      0.9348      0.0300   31.1274    0.0000 [ 0.8601; 0.9737 ] 
    ##   EXPE =~ expe3      0.7295      0.0712   10.2467    0.0000 [ 0.5732; 0.8431 ] 
    ##   QUAL =~ qual1      0.7861      0.0720   10.9229    0.0000 [ 0.6332; 0.8931 ] 
    ##   QUAL =~ qual2      0.9244      0.0232   39.7888    0.0000 [ 0.8648; 0.9570 ] 
    ##   QUAL =~ qual3      0.7560      0.0608   12.4407    0.0000 [ 0.6254; 0.8654 ] 
    ##   QUAL =~ qual4      0.7632      0.0536   14.2481    0.0000 [ 0.6476; 0.8595 ] 
    ##   QUAL =~ qual5      0.7834      0.0486   16.1246    0.0000 [ 0.6734; 0.8586 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.6263    0.0000 [ 0.8927; 0.9845 ] 
    ##   VAL =~ val2        0.8056      0.0620   12.9973    0.0000 [ 0.6723; 0.9015 ] 
    ##   VAL =~ val3        0.6763      0.0735    9.2040    0.0000 [ 0.5173; 0.8044 ] 
    ##   SAT =~ sat1        0.9243      0.0218   42.3359    0.0000 [ 0.8770; 0.9598 ] 
    ##   SAT =~ sat2        0.8813      0.0279   31.5835    0.0000 [ 0.8125; 0.9289 ] 
    ##   SAT =~ sat3        0.7127      0.0505   14.1084    0.0000 [ 0.6100; 0.8089 ] 
    ##   SAT =~ sat4        0.7756      0.0487   15.9265    0.0000 [ 0.6725; 0.8578 ] 
    ##   LOY =~ loy1        0.9097      0.0472   19.2752    0.0000 [ 0.8161; 0.9851 ] 
    ##   LOY =~ loy2        0.5775      0.0879    6.5671    0.0000 [ 0.4017; 0.7307 ] 
    ##   LOY =~ loy3        0.9043      0.0408   22.1411    0.0000 [ 0.8077; 0.9716 ] 
    ##   LOY =~ loy4        0.4917      0.1003    4.9006    0.0000 [ 0.2956; 0.7008 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1233    0.1269    0.8990 [-0.2260; 0.2539 ] 
    ##   IMAG <~ imag2      0.4473      0.1566    2.8564    0.0043 [ 0.1055; 0.7466 ] 
    ##   IMAG <~ imag3      0.6020      0.1451    4.1500    0.0000 [ 0.3136; 0.8616 ] 
    ##   EXPE <~ expe1      0.2946      0.1223    2.4081    0.0160 [ 0.0533; 0.5238 ] 
    ##   EXPE <~ expe2      0.6473      0.0895    7.2307    0.0000 [ 0.4633; 0.8005 ] 
    ##   EXPE <~ expe3      0.2374      0.0910    2.6083    0.0091 [ 0.0752; 0.4152 ] 
    ##   QUAL <~ qual1      0.2370      0.0929    2.5511    0.0107 [ 0.0577; 0.4295 ] 
    ##   QUAL <~ qual2      0.4712      0.0814    5.7905    0.0000 [ 0.3008; 0.6161 ] 
    ##   QUAL <~ qual3      0.1831      0.0802    2.2821    0.0225 [ 0.0253; 0.3320 ] 
    ##   QUAL <~ qual4      0.1037      0.0612    1.6939    0.0903 [-0.0052; 0.2381 ] 
    ##   QUAL <~ qual5      0.2049      0.0626    3.2731    0.0011 [ 0.0840; 0.3170 ] 
    ##   VAL <~ val1        0.7163      0.0950    7.5412    0.0000 [ 0.5131; 0.8733 ] 
    ##   VAL <~ val2        0.2202      0.0919    2.3951    0.0166 [ 0.0537; 0.4001 ] 
    ##   VAL <~ val3        0.2082      0.0619    3.3633    0.0008 [ 0.0906; 0.3302 ] 
    ##   SAT <~ sat1        0.3209      0.0145   22.1656    0.0000 [ 0.2950; 0.3524 ] 
    ##   SAT <~ sat2        0.3059      0.0131   23.3051    0.0000 [ 0.2839; 0.3338 ] 
    ##   SAT <~ sat3        0.2474      0.0108   22.9725    0.0000 [ 0.2270; 0.2679 ] 
    ##   SAT <~ sat4        0.2692      0.0123   21.8304    0.0000 [ 0.2449; 0.2939 ] 
    ##   LOY <~ loy1        0.3834      0.0262   14.6307    0.0000 [ 0.3339; 0.4349 ] 
    ##   LOY <~ loy2        0.2434      0.0311    7.8265    0.0000 [ 0.1765; 0.2910 ] 
    ##   LOY <~ loy3        0.3812      0.0274   13.9153    0.0000 [ 0.3293; 0.4316 ] 
    ##   LOY <~ loy4        0.2073      0.0370    5.5954    0.0000 [ 0.1374; 0.2831 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0631   10.1997    0.0000 [ 0.5184; 0.7585 ] 
    ##   imag1 ~~ imag3      0.5433      0.0709    7.6617    0.0000 [ 0.4124; 0.6761 ] 
    ##   imag2 ~~ imag3      0.7761      0.0401   19.3714    0.0000 [ 0.6913; 0.8419 ] 
    ##   expe1 ~~ expe2      0.5353      0.0611    8.7657    0.0000 [ 0.4019; 0.6370 ] 
    ##   expe1 ~~ expe3      0.4694      0.0642    7.3079    0.0000 [ 0.3233; 0.5791 ] 
    ##   expe2 ~~ expe3      0.5467      0.0616    8.8756    0.0000 [ 0.4132; 0.6516 ] 
    ##   qual1 ~~ qual2      0.6053      0.0593   10.2067    0.0000 [ 0.4897; 0.7063 ] 
    ##   qual1 ~~ qual3      0.5406      0.0613    8.8258    0.0000 [ 0.4198; 0.6597 ] 
    ##   qual1 ~~ qual4      0.5662      0.0730    7.7606    0.0000 [ 0.4198; 0.6932 ] 
    ##   qual1 ~~ qual5      0.5180      0.0741    6.9890    0.0000 [ 0.3631; 0.6460 ] 
    ##   qual2 ~~ qual3      0.6187      0.0562   11.0029    0.0000 [ 0.5024; 0.7247 ] 
    ##   qual2 ~~ qual4      0.6517      0.0615   10.5891    0.0000 [ 0.5298; 0.7638 ] 
    ##   qual2 ~~ qual5      0.6291      0.0597   10.5385    0.0000 [ 0.5001; 0.7320 ] 
    ##   qual3 ~~ qual4      0.4752      0.0666    7.1385    0.0000 [ 0.3445; 0.5964 ] 
    ##   qual3 ~~ qual5      0.5074      0.0637    7.9702    0.0000 [ 0.3772; 0.6194 ] 
    ##   qual4 ~~ qual5      0.6402      0.0517   12.3885    0.0000 [ 0.5341; 0.7289 ] 
    ##   val1 ~~ val2        0.6344      0.0527   12.0301    0.0000 [ 0.5315; 0.7307 ] 
    ##   val1 ~~ val3        0.4602      0.0675    6.8197    0.0000 [ 0.3348; 0.6016 ] 
    ##   val2 ~~ val3        0.6288      0.0629   10.0037    0.0000 [ 0.5055; 0.7494 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0661    7.1319    0.0000 [ 0.3449; 0.5974 ] 
    ##   QUAL ~ IMAG       0.3933      0.0611    6.4340    0.0000 [ 0.2787; 0.5095 ] 
    ##   QUAL ~ EXPE       0.8344      0.0234   35.7083    0.0000 [ 0.7858; 0.8726 ] 
    ##   VAL ~ IMAG        0.2974      0.0603    4.9359    0.0000 [ 0.1923; 0.4264 ] 
    ##   VAL ~ EXPE        0.6309      0.0492   12.8128    0.0000 [ 0.5306; 0.7259 ] 
    ##   VAL ~ QUAL        0.7013      0.0838    8.3739    0.0000 [ 0.5332; 0.8546 ] 
    ##   SAT ~ IMAG        0.4807      0.0676    7.1090    0.0000 [ 0.3513; 0.6087 ] 
    ##   SAT ~ EXPE        0.5001      0.0578    8.6462    0.0000 [ 0.3893; 0.6118 ] 
    ##   SAT ~ QUAL        0.5911      0.0941    6.2800    0.0000 [ 0.3969; 0.7631 ] 
    ##   SAT ~ VAL         0.5270      0.0871    6.0468    0.0000 [ 0.3523; 0.6874 ] 
    ##   LOY ~ IMAG        0.4840      0.0665    7.2780    0.0000 [ 0.3570; 0.6158 ] 
    ##   LOY ~ EXPE        0.3142      0.0566    5.5530    0.0000 [ 0.2107; 0.4336 ] 
    ##   LOY ~ QUAL        0.3714      0.0833    4.4585    0.0000 [ 0.2149; 0.5453 ] 
    ##   LOY ~ VAL         0.3311      0.0753    4.3991    0.0000 [ 0.1985; 0.4869 ] 
    ##   LOY ~ SAT         0.6283      0.0841    7.4682    0.0000 [ 0.4647; 0.7975 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0611    6.4340    0.0000 [ 0.2787; 0.5095 ] 
    ##   VAL ~ IMAG           0.2974      0.0603    4.9359    0.0000 [ 0.1923; 0.4264 ] 
    ##   VAL ~ EXPE           0.5852      0.0715    8.1873    0.0000 [ 0.4361; 0.7201 ] 
    ##   SAT ~ IMAG           0.2357      0.0480    4.9124    0.0000 [ 0.1546; 0.3354 ] 
    ##   SAT ~ EXPE           0.5173      0.0659    7.8554    0.0000 [ 0.4016; 0.6501 ] 
    ##   SAT ~ QUAL           0.3696      0.0605    6.1047    0.0000 [ 0.2483; 0.4786 ] 
    ##   LOY ~ IMAG           0.3020      0.0587    5.1426    0.0000 [ 0.2012; 0.4276 ] 
    ##   LOY ~ EXPE           0.3142      0.0566    5.5530    0.0000 [ 0.2107; 0.4336 ] 
    ##   LOY ~ QUAL           0.3714      0.0833    4.4585    0.0000 [ 0.2149; 0.5453 ] 
    ##   LOY ~ VAL            0.3311      0.0753    4.3991    0.0000 [ 0.1985; 0.4869 ] 
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
