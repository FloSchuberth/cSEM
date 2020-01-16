
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
    ##                                                   +------------------------------------------------------------+
    ##                                                   |                                                            |
    ##                                                   |   H0: Population indicator covariance matrix is equal to   |
    ##                                                   |   model-implied indicator covariance matrix.               |
    ##                                                   |                                                            |
    ##                                                   +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3083  
    ##  SRMR                    0.0940      0.0522  
    ##  dL                      2.2340      0.6898  
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
    ##  Out of 499 bootstrap replications 481 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -1641784535
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
    ##   sat1         1.3519         1.2324       1.7866         1.6165       0.2222
    ##   sat2         1.3066         1.1937       1.7641         1.6246       0.2017
    ##   sat3         1.4085         1.2750       1.7514         1.7199       0.1351
    ##   sat4         1.4177         1.2618       1.7850         1.6352       0.1740
    ##   loy1         1.7724         1.6578       2.2969         2.2241       0.2273
    ##   loy2         1.5148         1.4756       1.9373         1.9820       0.1071
    ##   loy3         1.7854         1.6682       2.3490         2.2648       0.2254
    ##   loy4         1.7080         1.6673       2.1936         2.2989       0.0724
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
    ##  Number of admissible results     = 490
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -99579465
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
    ##   EXPE ~ IMAG      0.4714      0.0659    7.1523    0.0000 [ 0.3403; 0.6083 ] 
    ##   QUAL ~ EXPE      0.8344      0.0238   35.0465    0.0000 [ 0.7840; 0.8726 ] 
    ##   VAL ~ EXPE       0.0457      0.0895    0.5110    0.6094 [-0.1174; 0.2308 ] 
    ##   VAL ~ QUAL       0.7013      0.0804    8.7185    0.0000 [ 0.5496; 0.8601 ] 
    ##   SAT ~ IMAG       0.2450      0.0575    4.2601    0.0000 [ 0.1340; 0.3615 ] 
    ##   SAT ~ EXPE      -0.0172      0.0728   -0.2368    0.8128 [-0.1626; 0.1242 ] 
    ##   SAT ~ QUAL       0.2215      0.0997    2.2218    0.0263 [ 0.0378; 0.4295 ] 
    ##   SAT ~ VAL        0.5270      0.0833    6.3272    0.0000 [ 0.3587; 0.6751 ] 
    ##   LOY ~ IMAG       0.1819      0.0777    2.3411    0.0192 [ 0.0344; 0.3423 ] 
    ##   LOY ~ SAT        0.6283      0.0772    8.1427    0.0000 [ 0.4865; 0.7739 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0993    6.3491    0.0000 [ 0.4247; 0.8058 ] 
    ##   IMAG =~ imag2      0.9246      0.0395   23.4188    0.0000 [ 0.8217; 0.9799 ] 
    ##   IMAG =~ imag3      0.9577      0.0299   32.0024    0.0000 [ 0.8798; 0.9909 ] 
    ##   EXPE =~ expe1      0.7525      0.0752   10.0096    0.0000 [ 0.5897; 0.8742 ] 
    ##   EXPE =~ expe2      0.9348      0.0289   32.3835    0.0000 [ 0.8624; 0.9740 ] 
    ##   EXPE =~ expe3      0.7295      0.0709   10.2855    0.0000 [ 0.5541; 0.8391 ] 
    ##   QUAL =~ qual1      0.7861      0.0649   12.1163    0.0000 [ 0.6384; 0.8848 ] 
    ##   QUAL =~ qual2      0.9244      0.0237   39.0428    0.0000 [ 0.8693; 0.9576 ] 
    ##   QUAL =~ qual3      0.7560      0.0581   13.0032    0.0000 [ 0.6308; 0.8542 ] 
    ##   QUAL =~ qual4      0.7632      0.0508   15.0117    0.0000 [ 0.6631; 0.8508 ] 
    ##   QUAL =~ qual5      0.7834      0.0461   16.9946    0.0000 [ 0.6772; 0.8611 ] 
    ##   VAL =~ val1        0.9518      0.0235   40.5632    0.0000 [ 0.8995; 0.9847 ] 
    ##   VAL =~ val2        0.8056      0.0649   12.4110    0.0000 [ 0.6555; 0.9060 ] 
    ##   VAL =~ val3        0.6763      0.0680    9.9399    0.0000 [ 0.5410; 0.8025 ] 
    ##   SAT =~ sat1        0.9243      0.0219   42.2090    0.0000 [ 0.8738; 0.9621 ] 
    ##   SAT =~ sat2        0.8813      0.0308   28.6140    0.0000 [ 0.8060; 0.9282 ] 
    ##   SAT =~ sat3        0.7127      0.0518   13.7656    0.0000 [ 0.6037; 0.8022 ] 
    ##   SAT =~ sat4        0.7756      0.0481   16.1404    0.0000 [ 0.6669; 0.8528 ] 
    ##   LOY =~ loy1        0.9097      0.0502   18.1212    0.0000 [ 0.7798; 0.9835 ] 
    ##   LOY =~ loy2        0.5775      0.0877    6.5885    0.0000 [ 0.3912; 0.7410 ] 
    ##   LOY =~ loy3        0.9043      0.0405   22.3162    0.0000 [ 0.8187; 0.9747 ] 
    ##   LOY =~ loy4        0.4917      0.0950    5.1765    0.0000 [ 0.3058; 0.6895 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1144    0.1368    0.8912 [-0.1869; 0.2705 ] 
    ##   IMAG <~ imag2      0.4473      0.1462    3.0600    0.0022 [ 0.1658; 0.7219 ] 
    ##   IMAG <~ imag3      0.6020      0.1424    4.2273    0.0000 [ 0.3031; 0.8446 ] 
    ##   EXPE <~ expe1      0.2946      0.1151    2.5586    0.0105 [ 0.0654; 0.5006 ] 
    ##   EXPE <~ expe2      0.6473      0.0866    7.4757    0.0000 [ 0.4584; 0.8010 ] 
    ##   EXPE <~ expe3      0.2374      0.0914    2.5971    0.0094 [ 0.0520; 0.4135 ] 
    ##   QUAL <~ qual1      0.2370      0.0868    2.7300    0.0063 [ 0.0883; 0.4261 ] 
    ##   QUAL <~ qual2      0.4712      0.0780    6.0421    0.0000 [ 0.3040; 0.6047 ] 
    ##   QUAL <~ qual3      0.1831      0.0761    2.4065    0.0161 [ 0.0262; 0.3192 ] 
    ##   QUAL <~ qual4      0.1037      0.0577    1.7979    0.0722 [-0.0011; 0.2134 ] 
    ##   QUAL <~ qual5      0.2049      0.0628    3.2616    0.0011 [ 0.0827; 0.3241 ] 
    ##   VAL <~ val1        0.7163      0.0945    7.5805    0.0000 [ 0.5227; 0.8784 ] 
    ##   VAL <~ val2        0.2202      0.0936    2.3525    0.0186 [ 0.0533; 0.4125 ] 
    ##   VAL <~ val3        0.2082      0.0589    3.5330    0.0004 [ 0.0919; 0.3139 ] 
    ##   SAT <~ sat1        0.3209      0.0154   20.7841    0.0000 [ 0.2937; 0.3548 ] 
    ##   SAT <~ sat2        0.3059      0.0130   23.4557    0.0000 [ 0.2843; 0.3351 ] 
    ##   SAT <~ sat3        0.2474      0.0106   23.2614    0.0000 [ 0.2276; 0.2675 ] 
    ##   SAT <~ sat4        0.2692      0.0121   22.2409    0.0000 [ 0.2470; 0.2930 ] 
    ##   LOY <~ loy1        0.3834      0.0270   14.2012    0.0000 [ 0.3330; 0.4344 ] 
    ##   LOY <~ loy2        0.2434      0.0302    8.0605    0.0000 [ 0.1763; 0.2962 ] 
    ##   LOY <~ loy3        0.3812      0.0274   13.9251    0.0000 [ 0.3324; 0.4437 ] 
    ##   LOY <~ loy4        0.2073      0.0349    5.9392    0.0000 [ 0.1368; 0.2739 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0668    9.6372    0.0000 [ 0.4967; 0.7507 ] 
    ##   imag1 ~~ imag3      0.5433      0.0727    7.4746    0.0000 [ 0.3917; 0.6861 ] 
    ##   imag2 ~~ imag3      0.7761      0.0393   19.7679    0.0000 [ 0.6917; 0.8457 ] 
    ##   expe1 ~~ expe2      0.5353      0.0606    8.8262    0.0000 [ 0.4197; 0.6393 ] 
    ##   expe1 ~~ expe3      0.4694      0.0621    7.5568    0.0000 [ 0.3463; 0.5812 ] 
    ##   expe2 ~~ expe3      0.5467      0.0607    9.0007    0.0000 [ 0.4142; 0.6619 ] 
    ##   qual1 ~~ qual2      0.6053      0.0573   10.5728    0.0000 [ 0.4808; 0.7080 ] 
    ##   qual1 ~~ qual3      0.5406      0.0587    9.2071    0.0000 [ 0.4150; 0.6451 ] 
    ##   qual1 ~~ qual4      0.5662      0.0666    8.5022    0.0000 [ 0.4245; 0.6797 ] 
    ##   qual1 ~~ qual5      0.5180      0.0702    7.3766    0.0000 [ 0.3838; 0.6425 ] 
    ##   qual2 ~~ qual3      0.6187      0.0550   11.2390    0.0000 [ 0.5027; 0.7143 ] 
    ##   qual2 ~~ qual4      0.6517      0.0596   10.9352    0.0000 [ 0.5289; 0.7632 ] 
    ##   qual2 ~~ qual5      0.6291      0.0569   11.0465    0.0000 [ 0.5220; 0.7329 ] 
    ##   qual3 ~~ qual4      0.4752      0.0645    7.3622    0.0000 [ 0.3313; 0.5867 ] 
    ##   qual3 ~~ qual5      0.5074      0.0618    8.2134    0.0000 [ 0.3784; 0.6272 ] 
    ##   qual4 ~~ qual5      0.6402      0.0558   11.4634    0.0000 [ 0.5216; 0.7415 ] 
    ##   val1 ~~ val2        0.6344      0.0533   11.9031    0.0000 [ 0.5242; 0.7331 ] 
    ##   val1 ~~ val3        0.4602      0.0642    7.1635    0.0000 [ 0.3359; 0.5824 ] 
    ##   val2 ~~ val3        0.6288      0.0630    9.9777    0.0000 [ 0.5013; 0.7386 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0659    7.1523    0.0000 [ 0.3403; 0.6083 ] 
    ##   QUAL ~ IMAG       0.3933      0.0622    6.3241    0.0000 [ 0.2724; 0.5218 ] 
    ##   QUAL ~ EXPE       0.8344      0.0238   35.0465    0.0000 [ 0.7840; 0.8726 ] 
    ##   VAL ~ IMAG        0.2974      0.0623    4.7748    0.0000 [ 0.1851; 0.4340 ] 
    ##   VAL ~ EXPE        0.6309      0.0530   11.8979    0.0000 [ 0.5223; 0.7257 ] 
    ##   VAL ~ QUAL        0.7013      0.0804    8.7185    0.0000 [ 0.5496; 0.8601 ] 
    ##   SAT ~ IMAG        0.4807      0.0671    7.1659    0.0000 [ 0.3482; 0.6068 ] 
    ##   SAT ~ EXPE        0.5001      0.0579    8.6358    0.0000 [ 0.3910; 0.6107 ] 
    ##   SAT ~ QUAL        0.5911      0.0915    6.4599    0.0000 [ 0.4237; 0.7835 ] 
    ##   SAT ~ VAL         0.5270      0.0833    6.3272    0.0000 [ 0.3587; 0.6751 ] 
    ##   LOY ~ IMAG        0.4840      0.0674    7.1795    0.0000 [ 0.3670; 0.6197 ] 
    ##   LOY ~ EXPE        0.3142      0.0538    5.8440    0.0000 [ 0.2168; 0.4228 ] 
    ##   LOY ~ QUAL        0.3714      0.0803    4.6249    0.0000 [ 0.2346; 0.5464 ] 
    ##   LOY ~ VAL         0.3311      0.0715    4.6325    0.0000 [ 0.1971; 0.4697 ] 
    ##   LOY ~ SAT         0.6283      0.0772    8.1427    0.0000 [ 0.4865; 0.7739 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0622    6.3241    0.0000 [ 0.2724; 0.5218 ] 
    ##   VAL ~ IMAG           0.2974      0.0623    4.7748    0.0000 [ 0.1851; 0.4340 ] 
    ##   VAL ~ EXPE           0.5852      0.0690    8.4754    0.0000 [ 0.4657; 0.7241 ] 
    ##   SAT ~ IMAG           0.2357      0.0493    4.7778    0.0000 [ 0.1495; 0.3432 ] 
    ##   SAT ~ EXPE           0.5173      0.0686    7.5376    0.0000 [ 0.3870; 0.6568 ] 
    ##   SAT ~ QUAL           0.3696      0.0621    5.9463    0.0000 [ 0.2474; 0.4939 ] 
    ##   LOY ~ IMAG           0.3020      0.0535    5.6401    0.0000 [ 0.2158; 0.4246 ] 
    ##   LOY ~ EXPE           0.3142      0.0538    5.8440    0.0000 [ 0.2168; 0.4228 ] 
    ##   LOY ~ QUAL           0.3714      0.0803    4.6249    0.0000 [ 0.2346; 0.5464 ] 
    ##   LOY ~ VAL            0.3311      0.0715    4.6325    0.0000 [ 0.1971; 0.4697 ] 
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
