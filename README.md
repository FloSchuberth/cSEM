
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

Currently we have five major postestimation verbs:

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
    ##                                      +------------------------------------------------------------+
    ##                                      |                                                            |
    ##                                      |   H0: Population indicator covariance matrix is equal to   |
    ##                                      |   model-implied indicator covariance matrix.               |
    ##                                      |                                                            |
    ##                                      +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3127  
    ##  SRMR                    0.0940      0.0528  
    ##  dL                      2.2340      0.7061  
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
    ##  Out of 499 bootstrap replications 468 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 1120327556
    ## ________________________________________________________________________________

``` r
## Assess the model
assess(res)
```

    ## ________________________________________________________________________________
    ## 
    ##  Construct        AVE          RhoC      RhoC_weighted      R2      
    ##  SAT            0.6851        0.8938        0.9051        0.7624    
    ##  LOY            0.5552        0.8011        0.8761        0.5868    
    ## 
    ##  Construct      R2_adj         RhoT      RhoT_weighted
    ##  SAT            0.7585        0.8940        0.8869    
    ##  LOY            0.5834        0.8194        0.7850    
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
    ## --------------------------- Effect sizes (f_squared) ---------------------------
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
    ## 
    ##  Redundancy analysis
    ## 
    ##  Construct       Value    
    ##  IMAG           0.9750    
    ##  EXPE           0.9873    
    ##  QUAL           0.9909    
    ##  VAL            0.9744    
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
    ##  Benchmark                        = 'unit'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name     MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   sat1         1.3612         1.3788       1.8018         1.8333       0.2102
    ##   sat2         1.3100         1.3175       1.7692         1.7963       0.1983
    ##   sat3         1.3990         1.4096       1.7404         1.7581       0.1461
    ##   sat4         1.4167         1.4334       1.7841         1.8056       0.1733
    ##   loy1         1.7629         1.7808       2.2855         2.3041       0.2346
    ##   loy2         1.4895         1.4741       1.9147         1.8967       0.1295
    ##   loy3         1.7760         1.8001       2.3393         2.3635       0.2337
    ##   loy4         1.6857         1.6766       2.1740         2.1571       0.0922
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential quantities are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Some closed form standard error are implemented, however, this
feature is still rather preliminary. It is therefore recommoned to use
resampling instead. `cSEM` mostly relies on the `bootstrap` procedure
(although `jackknife` is implemented as well) to estimate standard
errors, test statistics, and critical quantiles.

`cSEM` offers two ways to compute resamples:

1.  Setting `.resample_method` to `"jackkinfe"` or `"bootstrap"` and
    subsequently using postestimation functions `summarize()` or
    `infer()`.
2.  The same result is achieved by passing a `cSEMResults` object to
    `resamplecSEMResults()` and subsequently using postestimation
    functions `summarize()` or `infer()`.

<!-- end list -->

``` r
# Setting `.resample_method`
b1 <- csem(.data = satisfaction, .model = model, .resample_method = "bootstrap")
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
    ##  Number of admissible results     = 483
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 2141785283
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
    ##   EXPE ~ IMAG      0.4714      0.0672    7.0183    0.0000 [ 0.3360; 0.6000 ] 
    ##   QUAL ~ EXPE      0.8344      0.0239   34.9824    0.0000 [ 0.7843; 0.8752 ] 
    ##   VAL ~ EXPE       0.0457      0.0959    0.4766    0.6337 [-0.1176; 0.2464 ] 
    ##   VAL ~ QUAL       0.7013      0.0901    7.7841    0.0000 [ 0.5123; 0.8555 ] 
    ##   SAT ~ IMAG       0.2450      0.0530    4.6235    0.0000 [ 0.1369; 0.3399 ] 
    ##   SAT ~ EXPE      -0.0172      0.0726   -0.2375    0.8122 [-0.1524; 0.1247 ] 
    ##   SAT ~ QUAL       0.2215      0.1012    2.1901    0.0285 [ 0.0345; 0.4349 ] 
    ##   SAT ~ VAL        0.5270      0.0866    6.0883    0.0000 [ 0.3514; 0.6868 ] 
    ##   LOY ~ IMAG       0.1819      0.0803    2.2661    0.0234 [ 0.0272; 0.3410 ] 
    ##   LOY ~ SAT        0.6283      0.0802    7.8293    0.0000 [ 0.4837; 0.7882 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0989    6.3794    0.0000 [ 0.4209; 0.7930 ] 
    ##   IMAG =~ imag2      0.9246      0.0401   23.0517    0.0000 [ 0.8287; 0.9800 ] 
    ##   IMAG =~ imag3      0.9577      0.0278   34.4895    0.0000 [ 0.8876; 0.9913 ] 
    ##   EXPE =~ expe1      0.7525      0.0810    9.2927    0.0000 [ 0.5620; 0.8751 ] 
    ##   EXPE =~ expe2      0.9348      0.0300   31.1614    0.0000 [ 0.8531; 0.9712 ] 
    ##   EXPE =~ expe3      0.7295      0.0704   10.3641    0.0000 [ 0.5777; 0.8478 ] 
    ##   QUAL =~ qual1      0.7861      0.0700   11.2297    0.0000 [ 0.6256; 0.8965 ] 
    ##   QUAL =~ qual2      0.9244      0.0230   40.2214    0.0000 [ 0.8699; 0.9611 ] 
    ##   QUAL =~ qual3      0.7560      0.0621   12.1822    0.0000 [ 0.6218; 0.8650 ] 
    ##   QUAL =~ qual4      0.7632      0.0553   13.8037    0.0000 [ 0.6407; 0.8561 ] 
    ##   QUAL =~ qual5      0.7834      0.0464   16.8690    0.0000 [ 0.6791; 0.8534 ] 
    ##   VAL =~ val1        0.9518      0.0257   37.0357    0.0000 [ 0.8957; 0.9855 ] 
    ##   VAL =~ val2        0.8056      0.0697   11.5620    0.0000 [ 0.6556; 0.9124 ] 
    ##   VAL =~ val3        0.6763      0.0776    8.7113    0.0000 [ 0.5066; 0.8108 ] 
    ##   SAT =~ sat1        0.9243      0.0222   41.5593    0.0000 [ 0.8766; 0.9637 ] 
    ##   SAT =~ sat2        0.8813      0.0285   30.9148    0.0000 [ 0.8212; 0.9253 ] 
    ##   SAT =~ sat3        0.7127      0.0497   14.3283    0.0000 [ 0.6087; 0.8056 ] 
    ##   SAT =~ sat4        0.7756      0.0482   16.0773    0.0000 [ 0.6693; 0.8620 ] 
    ##   LOY =~ loy1        0.9097      0.0487   18.6657    0.0000 [ 0.7881; 0.9845 ] 
    ##   LOY =~ loy2        0.5775      0.0850    6.7964    0.0000 [ 0.3866; 0.7265 ] 
    ##   LOY =~ loy3        0.9043      0.0425   21.2777    0.0000 [ 0.8115; 0.9749 ] 
    ##   LOY =~ loy4        0.4917      0.1035    4.7495    0.0000 [ 0.2965; 0.6927 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1139    0.1374    0.8907 [-0.2143; 0.2224 ] 
    ##   IMAG <~ imag2      0.4473      0.1523    2.9375    0.0033 [ 0.1548; 0.7407 ] 
    ##   IMAG <~ imag3      0.6020      0.1408    4.2755    0.0000 [ 0.3103; 0.8519 ] 
    ##   EXPE <~ expe1      0.2946      0.1187    2.4817    0.0131 [ 0.0734; 0.5285 ] 
    ##   EXPE <~ expe2      0.6473      0.0853    7.5896    0.0000 [ 0.4594; 0.7929 ] 
    ##   EXPE <~ expe3      0.2374      0.0918    2.5861    0.0097 [ 0.0560; 0.4242 ] 
    ##   QUAL <~ qual1      0.2370      0.0959    2.4728    0.0134 [ 0.0698; 0.4400 ] 
    ##   QUAL <~ qual2      0.4712      0.0791    5.9565    0.0000 [ 0.2996; 0.6326 ] 
    ##   QUAL <~ qual3      0.1831      0.0817    2.2412    0.0250 [ 0.0190; 0.3377 ] 
    ##   QUAL <~ qual4      0.1037      0.0623    1.6643    0.0961 [-0.0134; 0.2292 ] 
    ##   QUAL <~ qual5      0.2049      0.0638    3.2110    0.0013 [ 0.0644; 0.3142 ] 
    ##   VAL <~ val1        0.7163      0.1031    6.9449    0.0000 [ 0.5017; 0.8936 ] 
    ##   VAL <~ val2        0.2202      0.0993    2.2180    0.0266 [ 0.0417; 0.4214 ] 
    ##   VAL <~ val3        0.2082      0.0637    3.2662    0.0011 [ 0.0919; 0.3364 ] 
    ##   SAT <~ sat1        0.3209      0.0145   22.1521    0.0000 [ 0.2955; 0.3531 ] 
    ##   SAT <~ sat2        0.3059      0.0129   23.7075    0.0000 [ 0.2847; 0.3354 ] 
    ##   SAT <~ sat3        0.2474      0.0106   23.2508    0.0000 [ 0.2234; 0.2682 ] 
    ##   SAT <~ sat4        0.2692      0.0120   22.4307    0.0000 [ 0.2463; 0.2938 ] 
    ##   LOY <~ loy1        0.3834      0.0259   14.8249    0.0000 [ 0.3309; 0.4340 ] 
    ##   LOY <~ loy2        0.2434      0.0299    8.1363    0.0000 [ 0.1693; 0.2902 ] 
    ##   LOY <~ loy3        0.3812      0.0286   13.3138    0.0000 [ 0.3268; 0.4363 ] 
    ##   LOY <~ loy4        0.2073      0.0382    5.4204    0.0000 [ 0.1342; 0.2775 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0678    9.4996    0.0000 [ 0.4981; 0.7603 ] 
    ##   imag1 ~~ imag3      0.5433      0.0725    7.4981    0.0000 [ 0.3957; 0.6681 ] 
    ##   imag2 ~~ imag3      0.7761      0.0387   20.0571    0.0000 [ 0.6909; 0.8446 ] 
    ##   expe1 ~~ expe2      0.5353      0.0610    8.7740    0.0000 [ 0.4011; 0.6434 ] 
    ##   expe1 ~~ expe3      0.4694      0.0622    7.5428    0.0000 [ 0.3474; 0.5838 ] 
    ##   expe2 ~~ expe3      0.5467      0.0593    9.2239    0.0000 [ 0.4286; 0.6599 ] 
    ##   qual1 ~~ qual2      0.6053      0.0581   10.4240    0.0000 [ 0.4835; 0.7128 ] 
    ##   qual1 ~~ qual3      0.5406      0.0613    8.8198    0.0000 [ 0.4177; 0.6575 ] 
    ##   qual1 ~~ qual4      0.5662      0.0671    8.4389    0.0000 [ 0.4302; 0.6843 ] 
    ##   qual1 ~~ qual5      0.5180      0.0656    7.9020    0.0000 [ 0.3872; 0.6360 ] 
    ##   qual2 ~~ qual3      0.6187      0.0566   10.9331    0.0000 [ 0.4900; 0.7172 ] 
    ##   qual2 ~~ qual4      0.6517      0.0620   10.5097    0.0000 [ 0.5273; 0.7603 ] 
    ##   qual2 ~~ qual5      0.6291      0.0547   11.4988    0.0000 [ 0.5111; 0.7311 ] 
    ##   qual3 ~~ qual4      0.4752      0.0653    7.2791    0.0000 [ 0.3430; 0.5950 ] 
    ##   qual3 ~~ qual5      0.5074      0.0627    8.0927    0.0000 [ 0.3822; 0.6245 ] 
    ##   qual4 ~~ qual5      0.6402      0.0554   11.5621    0.0000 [ 0.5197; 0.7362 ] 
    ##   val1 ~~ val2        0.6344      0.0576   11.0212    0.0000 [ 0.5078; 0.7256 ] 
    ##   val1 ~~ val3        0.4602      0.0719    6.3999    0.0000 [ 0.3083; 0.5912 ] 
    ##   val2 ~~ val3        0.6288      0.0667    9.4292    0.0000 [ 0.4934; 0.7528 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0672    7.0183    0.0000 [ 0.3360; 0.6000 ] 
    ##   QUAL ~ IMAG       0.3933      0.0622    6.3235    0.0000 [ 0.2656; 0.5175 ] 
    ##   QUAL ~ EXPE       0.8344      0.0239   34.9824    0.0000 [ 0.7843; 0.8752 ] 
    ##   VAL ~ IMAG        0.2974      0.0620    4.7960    0.0000 [ 0.1779; 0.4245 ] 
    ##   VAL ~ EXPE        0.6309      0.0520   12.1323    0.0000 [ 0.5095; 0.7246 ] 
    ##   VAL ~ QUAL        0.7013      0.0901    7.7841    0.0000 [ 0.5123; 0.8555 ] 
    ##   SAT ~ IMAG        0.4807      0.0656    7.3220    0.0000 [ 0.3498; 0.5997 ] 
    ##   SAT ~ EXPE        0.5001      0.0577    8.6687    0.0000 [ 0.3881; 0.6157 ] 
    ##   SAT ~ QUAL        0.5911      0.0973    6.0736    0.0000 [ 0.4048; 0.7635 ] 
    ##   SAT ~ VAL         0.5270      0.0866    6.0883    0.0000 [ 0.3514; 0.6868 ] 
    ##   LOY ~ IMAG        0.4840      0.0673    7.1926    0.0000 [ 0.3521; 0.6163 ] 
    ##   LOY ~ EXPE        0.3142      0.0552    5.6873    0.0000 [ 0.2217; 0.4367 ] 
    ##   LOY ~ QUAL        0.3714      0.0854    4.3480    0.0000 [ 0.2217; 0.5469 ] 
    ##   LOY ~ VAL         0.3311      0.0712    4.6510    0.0000 [ 0.1887; 0.4737 ] 
    ##   LOY ~ SAT         0.6283      0.0802    7.8293    0.0000 [ 0.4837; 0.7882 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0622    6.3235    0.0000 [ 0.2656; 0.5175 ] 
    ##   VAL ~ IMAG           0.2974      0.0620    4.7960    0.0000 [ 0.1779; 0.4245 ] 
    ##   VAL ~ EXPE           0.5852      0.0760    7.6985    0.0000 [ 0.4299; 0.7215 ] 
    ##   SAT ~ IMAG           0.2357      0.0497    4.7432    0.0000 [ 0.1463; 0.3429 ] 
    ##   SAT ~ EXPE           0.5173      0.0650    7.9561    0.0000 [ 0.3903; 0.6439 ] 
    ##   SAT ~ QUAL           0.3696      0.0625    5.9084    0.0000 [ 0.2467; 0.4880 ] 
    ##   LOY ~ IMAG           0.3020      0.0550    5.4919    0.0000 [ 0.2109; 0.4205 ] 
    ##   LOY ~ EXPE           0.3142      0.0552    5.6873    0.0000 [ 0.2217; 0.4367 ] 
    ##   LOY ~ QUAL           0.3714      0.0854    4.3480    0.0000 [ 0.2217; 0.5469 ] 
    ##   LOY ~ VAL            0.3311      0.0712    4.6510    0.0000 [ 0.1887; 0.4737 ] 
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
