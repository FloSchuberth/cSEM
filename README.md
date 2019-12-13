
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
    ##  dG                      0.6493      0.3345  
    ##  SRMR                    0.0940      0.0534  
    ##  dL                      2.2340      0.7220  
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
    ##  Out of 499 bootstrap replications 472 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -178853543
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
    ##  Benchmark                        = 'lm'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name     MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   sat1         1.3499         1.2315       1.7855         1.6165       0.2231
    ##   sat2         1.3074         1.1953       1.7634         1.6249       0.2022
    ##   sat3         1.4085         1.2770       1.7503         1.7218       0.1357
    ##   sat4         1.4136         1.2610       1.7811         1.6335       0.1746
    ##   loy1         1.7711         1.6598       2.2949         2.2234       0.2275
    ##   loy2         1.5168         1.4751       1.9382         1.9804       0.1091
    ##   loy3         1.7825         1.6709       2.3474         2.2715       0.2280
    ##   loy4         1.7095         1.6715       2.1965         2.3045       0.0713
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
    ##  Number of admissible results     = 488
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 1358427687
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
    ##   EXPE ~ IMAG      0.4714      0.0634    7.4382    0.0000 [ 0.3498; 0.5919 ] 
    ##   QUAL ~ EXPE      0.8344      0.0226   36.9157    0.0000 [ 0.7881; 0.8762 ] 
    ##   VAL ~ EXPE       0.0457      0.0880    0.5193    0.6035 [-0.1323; 0.2125 ] 
    ##   VAL ~ QUAL       0.7013      0.0861    8.1430    0.0000 [ 0.5447; 0.8685 ] 
    ##   SAT ~ IMAG       0.2450      0.0549    4.4646    0.0000 [ 0.1434; 0.3654 ] 
    ##   SAT ~ EXPE      -0.0172      0.0731   -0.2356    0.8137 [-0.1704; 0.1150 ] 
    ##   SAT ~ QUAL       0.2215      0.1035    2.1412    0.0323 [ 0.0517; 0.4389 ] 
    ##   SAT ~ VAL        0.5270      0.0893    5.9033    0.0000 [ 0.3406; 0.6856 ] 
    ##   LOY ~ IMAG       0.1819      0.0782    2.3254    0.0201 [ 0.0425; 0.3334 ] 
    ##   LOY ~ SAT        0.6283      0.0820    7.6650    0.0000 [ 0.4763; 0.7900 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0956    6.5994    0.0000 [ 0.4339; 0.7943 ] 
    ##   IMAG =~ imag2      0.9246      0.0370   24.9583    0.0000 [ 0.8371; 0.9756 ] 
    ##   IMAG =~ imag3      0.9577      0.0272   35.2441    0.0000 [ 0.8901; 0.9897 ] 
    ##   EXPE =~ expe1      0.7525      0.0722   10.4266    0.0000 [ 0.5910; 0.8739 ] 
    ##   EXPE =~ expe2      0.9348      0.0282   33.1845    0.0000 [ 0.8620; 0.9733 ] 
    ##   EXPE =~ expe3      0.7295      0.0689   10.5942    0.0000 [ 0.5898; 0.8440 ] 
    ##   QUAL =~ qual1      0.7861      0.0623   12.6226    0.0000 [ 0.6399; 0.8916 ] 
    ##   QUAL =~ qual2      0.9244      0.0225   41.0778    0.0000 [ 0.8699; 0.9543 ] 
    ##   QUAL =~ qual3      0.7560      0.0567   13.3281    0.0000 [ 0.6375; 0.8536 ] 
    ##   QUAL =~ qual4      0.7632      0.0528   14.4488    0.0000 [ 0.6448; 0.8538 ] 
    ##   QUAL =~ qual5      0.7834      0.0457   17.1338    0.0000 [ 0.6839; 0.8616 ] 
    ##   VAL =~ val1        0.9518      0.0218   43.6083    0.0000 [ 0.8990; 0.9837 ] 
    ##   VAL =~ val2        0.8056      0.0614   13.1290    0.0000 [ 0.6721; 0.9023 ] 
    ##   VAL =~ val3        0.6763      0.0719    9.4060    0.0000 [ 0.5252; 0.8154 ] 
    ##   SAT =~ sat1        0.9243      0.0234   39.5280    0.0000 [ 0.8695; 0.9634 ] 
    ##   SAT =~ sat2        0.8813      0.0302   29.1973    0.0000 [ 0.8164; 0.9305 ] 
    ##   SAT =~ sat3        0.7127      0.0526   13.5518    0.0000 [ 0.6006; 0.7995 ] 
    ##   SAT =~ sat4        0.7756      0.0501   15.4666    0.0000 [ 0.6733; 0.8679 ] 
    ##   LOY =~ loy1        0.9097      0.0479   18.9958    0.0000 [ 0.8008; 0.9866 ] 
    ##   LOY =~ loy2        0.5775      0.0832    6.9382    0.0000 [ 0.4064; 0.7246 ] 
    ##   LOY =~ loy3        0.9043      0.0411   21.9785    0.0000 [ 0.8147; 0.9809 ] 
    ##   LOY =~ loy4        0.4917      0.0952    5.1636    0.0000 [ 0.3036; 0.6839 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1139    0.1373    0.8908 [-0.2156; 0.2565 ] 
    ##   IMAG <~ imag2      0.4473      0.1413    3.1649    0.0016 [ 0.1846; 0.7120 ] 
    ##   IMAG <~ imag3      0.6020      0.1312    4.5873    0.0000 [ 0.3296; 0.8242 ] 
    ##   EXPE <~ expe1      0.2946      0.1093    2.6955    0.0070 [ 0.0849; 0.5122 ] 
    ##   EXPE <~ expe2      0.6473      0.0834    7.7575    0.0000 [ 0.4685; 0.8001 ] 
    ##   EXPE <~ expe3      0.2374      0.0850    2.7930    0.0052 [ 0.0659; 0.4037 ] 
    ##   QUAL <~ qual1      0.2370      0.0845    2.8039    0.0050 [ 0.0919; 0.4150 ] 
    ##   QUAL <~ qual2      0.4712      0.0759    6.2072    0.0000 [ 0.3094; 0.6099 ] 
    ##   QUAL <~ qual3      0.1831      0.0744    2.4607    0.0139 [ 0.0302; 0.3234 ] 
    ##   QUAL <~ qual4      0.1037      0.0623    1.6651    0.0959 [-0.0159; 0.2245 ] 
    ##   QUAL <~ qual5      0.2049      0.0641    3.1975    0.0014 [ 0.0741; 0.3188 ] 
    ##   VAL <~ val1        0.7163      0.0909    7.8830    0.0000 [ 0.5299; 0.8688 ] 
    ##   VAL <~ val2        0.2202      0.0893    2.4667    0.0136 [ 0.0427; 0.3998 ] 
    ##   VAL <~ val3        0.2082      0.0629    3.3113    0.0009 [ 0.0926; 0.3416 ] 
    ##   SAT <~ sat1        0.3209      0.0151   21.2997    0.0000 [ 0.2965; 0.3542 ] 
    ##   SAT <~ sat2        0.3059      0.0145   21.0598    0.0000 [ 0.2824; 0.3378 ] 
    ##   SAT <~ sat3        0.2474      0.0110   22.5871    0.0000 [ 0.2238; 0.2669 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.0808    0.0000 [ 0.2440; 0.2933 ] 
    ##   LOY <~ loy1        0.3834      0.0258   14.8725    0.0000 [ 0.3308; 0.4270 ] 
    ##   LOY <~ loy2        0.2434      0.0289    8.4290    0.0000 [ 0.1800; 0.2914 ] 
    ##   LOY <~ loy3        0.3812      0.0256   14.9169    0.0000 [ 0.3297; 0.4262 ] 
    ##   LOY <~ loy4        0.2073      0.0356    5.8167    0.0000 [ 0.1385; 0.2707 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0617   10.4254    0.0000 [ 0.5062; 0.7525 ] 
    ##   imag1 ~~ imag3      0.5433      0.0667    8.1425    0.0000 [ 0.4069; 0.6611 ] 
    ##   imag2 ~~ imag3      0.7761      0.0377   20.5951    0.0000 [ 0.6985; 0.8503 ] 
    ##   expe1 ~~ expe2      0.5353      0.0569    9.3988    0.0000 [ 0.4321; 0.6416 ] 
    ##   expe1 ~~ expe3      0.4694      0.0617    7.6111    0.0000 [ 0.3430; 0.5889 ] 
    ##   expe2 ~~ expe3      0.5467      0.0619    8.8306    0.0000 [ 0.4059; 0.6545 ] 
    ##   qual1 ~~ qual2      0.6053      0.0547   11.0689    0.0000 [ 0.4819; 0.7009 ] 
    ##   qual1 ~~ qual3      0.5406      0.0561    9.6403    0.0000 [ 0.4343; 0.6460 ] 
    ##   qual1 ~~ qual4      0.5662      0.0678    8.3561    0.0000 [ 0.4291; 0.6875 ] 
    ##   qual1 ~~ qual5      0.5180      0.0685    7.5609    0.0000 [ 0.3894; 0.6424 ] 
    ##   qual2 ~~ qual3      0.6187      0.0570   10.8532    0.0000 [ 0.4981; 0.7294 ] 
    ##   qual2 ~~ qual4      0.6517      0.0608   10.7148    0.0000 [ 0.5208; 0.7506 ] 
    ##   qual2 ~~ qual5      0.6291      0.0571   11.0138    0.0000 [ 0.5065; 0.7343 ] 
    ##   qual3 ~~ qual4      0.4752      0.0658    7.2178    0.0000 [ 0.3381; 0.5916 ] 
    ##   qual3 ~~ qual5      0.5074      0.0618    8.2113    0.0000 [ 0.3802; 0.6244 ] 
    ##   qual4 ~~ qual5      0.6402      0.0568   11.2696    0.0000 [ 0.5185; 0.7421 ] 
    ##   val1 ~~ val2        0.6344      0.0539   11.7769    0.0000 [ 0.5211; 0.7316 ] 
    ##   val1 ~~ val3        0.4602      0.0699    6.5845    0.0000 [ 0.3193; 0.5979 ] 
    ##   val2 ~~ val3        0.6288      0.0602   10.4389    0.0000 [ 0.5078; 0.7445 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0634    7.4382    0.0000 [ 0.3498; 0.5919 ] 
    ##   QUAL ~ IMAG       0.3933      0.0591    6.6566    0.0000 [ 0.2843; 0.5085 ] 
    ##   QUAL ~ EXPE       0.8344      0.0226   36.9157    0.0000 [ 0.7881; 0.8762 ] 
    ##   VAL ~ IMAG        0.2974      0.0592    5.0217    0.0000 [ 0.1977; 0.4173 ] 
    ##   VAL ~ EXPE        0.6309      0.0514   12.2773    0.0000 [ 0.5325; 0.7296 ] 
    ##   VAL ~ QUAL        0.7013      0.0861    8.1430    0.0000 [ 0.5447; 0.8685 ] 
    ##   SAT ~ IMAG        0.4807      0.0646    7.4407    0.0000 [ 0.3634; 0.6119 ] 
    ##   SAT ~ EXPE        0.5001      0.0570    8.7712    0.0000 [ 0.3842; 0.6060 ] 
    ##   SAT ~ QUAL        0.5911      0.0967    6.1144    0.0000 [ 0.4231; 0.7941 ] 
    ##   SAT ~ VAL         0.5270      0.0893    5.9033    0.0000 [ 0.3406; 0.6856 ] 
    ##   LOY ~ IMAG        0.4840      0.0655    7.3834    0.0000 [ 0.3550; 0.6196 ] 
    ##   LOY ~ EXPE        0.3142      0.0558    5.6291    0.0000 [ 0.2032; 0.4224 ] 
    ##   LOY ~ QUAL        0.3714      0.0833    4.4565    0.0000 [ 0.2383; 0.5569 ] 
    ##   LOY ~ VAL         0.3311      0.0773    4.2808    0.0000 [ 0.1876; 0.4891 ] 
    ##   LOY ~ SAT         0.6283      0.0820    7.6650    0.0000 [ 0.4763; 0.7900 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0591    6.6566    0.0000 [ 0.2843; 0.5085 ] 
    ##   VAL ~ IMAG           0.2974      0.0592    5.0217    0.0000 [ 0.1977; 0.4173 ] 
    ##   VAL ~ EXPE           0.5852      0.0727    8.0460    0.0000 [ 0.4562; 0.7377 ] 
    ##   SAT ~ IMAG           0.2357      0.0467    5.0509    0.0000 [ 0.1502; 0.3260 ] 
    ##   SAT ~ EXPE           0.5173      0.0665    7.7757    0.0000 [ 0.3949; 0.6533 ] 
    ##   SAT ~ QUAL           0.3696      0.0635    5.8200    0.0000 [ 0.2447; 0.4955 ] 
    ##   LOY ~ IMAG           0.3020      0.0571    5.2914    0.0000 [ 0.2031; 0.4295 ] 
    ##   LOY ~ EXPE           0.3142      0.0558    5.6291    0.0000 [ 0.2032; 0.4224 ] 
    ##   LOY ~ QUAL           0.3714      0.0833    4.4565    0.0000 [ 0.2383; 0.5569 ] 
    ##   LOY ~ VAL            0.3311      0.0773    4.2808    0.0000 [ 0.1876; 0.4891 ] 
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
