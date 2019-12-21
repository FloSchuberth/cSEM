
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

NOTE: the package will be submitted to CRAN as version 0.1.0 once the
CRAN team returns from vacation (Jan 6, 2020). Expect the package to be
on CRAN some days later.

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
    ##                                   +------------------------------------------------------------+
    ##                                   |                                                            |
    ##                                   |   H0: Population indicator covariance matrix is equal to   |
    ##                                   |   model-implied indicator covariance matrix.               |
    ##                                   |                                                            |
    ##                                   +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3369  
    ##  SRMR                    0.0940      0.0529  
    ##  dL                      2.2340      0.7067  
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
    ##  Out of 499 bootstrap replications 483 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -1888861271
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
    ##   sat1         1.3505         1.2359       1.7863         1.6210       0.2219
    ##   sat2         1.3070         1.1997       1.7642         1.6286       0.2014
    ##   sat3         1.4067         1.2739       1.7487         1.7185       0.1360
    ##   sat4         1.4158         1.2645       1.7824         1.6390       0.1745
    ##   loy1         1.7712         1.6587       2.2947         2.2250       0.2274
    ##   loy2         1.5146         1.4735       1.9359         1.9822       0.1089
    ##   loy3         1.7817         1.6682       2.3449         2.2694       0.2285
    ##   loy4         1.7132         1.6761       2.1995         2.3101       0.0700
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
    ##  Random seed                      = -978869286
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
    ##   EXPE ~ IMAG      0.4714      0.0643    7.3336    0.0000 [ 0.3484; 0.5993 ] 
    ##   QUAL ~ EXPE      0.8344      0.0226   36.9133    0.0000 [ 0.7851; 0.8750 ] 
    ##   VAL ~ EXPE       0.0457      0.0836    0.5468    0.5845 [-0.1068; 0.2138 ] 
    ##   VAL ~ QUAL       0.7013      0.0799    8.7793    0.0000 [ 0.5390; 0.8465 ] 
    ##   SAT ~ IMAG       0.2450      0.0528    4.6357    0.0000 [ 0.1438; 0.3439 ] 
    ##   SAT ~ EXPE      -0.0172      0.0700   -0.2461    0.8056 [-0.1504; 0.1216 ] 
    ##   SAT ~ QUAL       0.2215      0.0941    2.3537    0.0186 [ 0.0426; 0.4066 ] 
    ##   SAT ~ VAL        0.5270      0.0834    6.3160    0.0000 [ 0.3577; 0.6815 ] 
    ##   LOY ~ IMAG       0.1819      0.0785    2.3164    0.0205 [ 0.0372; 0.3347 ] 
    ##   LOY ~ SAT        0.6283      0.0783    8.0217    0.0000 [ 0.4865; 0.7824 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1018    6.1923    0.0000 [ 0.4094; 0.8237 ] 
    ##   IMAG =~ imag2      0.9246      0.0418   22.1154    0.0000 [ 0.8216; 0.9784 ] 
    ##   IMAG =~ imag3      0.9577      0.0297   32.2836    0.0000 [ 0.8795; 0.9924 ] 
    ##   EXPE =~ expe1      0.7525      0.0808    9.3102    0.0000 [ 0.5684; 0.8774 ] 
    ##   EXPE =~ expe2      0.9348      0.0288   32.4590    0.0000 [ 0.8582; 0.9717 ] 
    ##   EXPE =~ expe3      0.7295      0.0700   10.4229    0.0000 [ 0.5780; 0.8414 ] 
    ##   QUAL =~ qual1      0.7861      0.0679   11.5755    0.0000 [ 0.6196; 0.8892 ] 
    ##   QUAL =~ qual2      0.9244      0.0220   42.1023    0.0000 [ 0.8690; 0.9536 ] 
    ##   QUAL =~ qual3      0.7560      0.0601   12.5688    0.0000 [ 0.6230; 0.8581 ] 
    ##   QUAL =~ qual4      0.7632      0.0514   14.8340    0.0000 [ 0.6520; 0.8458 ] 
    ##   QUAL =~ qual5      0.7834      0.0455   17.2011    0.0000 [ 0.6830; 0.8562 ] 
    ##   VAL =~ val1        0.9518      0.0236   40.3614    0.0000 [ 0.8942; 0.9839 ] 
    ##   VAL =~ val2        0.8056      0.0669   12.0458    0.0000 [ 0.6546; 0.9142 ] 
    ##   VAL =~ val3        0.6763      0.0764    8.8543    0.0000 [ 0.5182; 0.8142 ] 
    ##   SAT =~ sat1        0.9243      0.0219   42.1849    0.0000 [ 0.8784; 0.9618 ] 
    ##   SAT =~ sat2        0.8813      0.0284   31.0445    0.0000 [ 0.8223; 0.9266 ] 
    ##   SAT =~ sat3        0.7127      0.0524   13.6095    0.0000 [ 0.5952; 0.8011 ] 
    ##   SAT =~ sat4        0.7756      0.0497   15.6169    0.0000 [ 0.6711; 0.8633 ] 
    ##   LOY =~ loy1        0.9097      0.0529   17.2116    0.0000 [ 0.7828; 0.9883 ] 
    ##   LOY =~ loy2        0.5775      0.0850    6.7915    0.0000 [ 0.3994; 0.7277 ] 
    ##   LOY =~ loy3        0.9043      0.0430   21.0436    0.0000 [ 0.8142; 0.9735 ] 
    ##   LOY =~ loy4        0.4917      0.1031    4.7705    0.0000 [ 0.2762; 0.6905 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1228    0.1274    0.8986 [-0.1914; 0.2753 ] 
    ##   IMAG <~ imag2      0.4473      0.1616    2.7680    0.0056 [ 0.1205; 0.7777 ] 
    ##   IMAG <~ imag3      0.6020      0.1448    4.1568    0.0000 [ 0.3091; 0.8603 ] 
    ##   EXPE <~ expe1      0.2946      0.1209    2.4376    0.0148 [ 0.0493; 0.5229 ] 
    ##   EXPE <~ expe2      0.6473      0.0870    7.4389    0.0000 [ 0.4587; 0.7863 ] 
    ##   EXPE <~ expe3      0.2374      0.0945    2.5111    0.0120 [ 0.0481; 0.4108 ] 
    ##   QUAL <~ qual1      0.2370      0.0906    2.6172    0.0089 [ 0.0832; 0.4242 ] 
    ##   QUAL <~ qual2      0.4712      0.0784    6.0081    0.0000 [ 0.3028; 0.6067 ] 
    ##   QUAL <~ qual3      0.1831      0.0810    2.2603    0.0238 [ 0.0183; 0.3369 ] 
    ##   QUAL <~ qual4      0.1037      0.0603    1.7191    0.0856 [-0.0058; 0.2228 ] 
    ##   QUAL <~ qual5      0.2049      0.0620    3.3048    0.0010 [ 0.0697; 0.3084 ] 
    ##   VAL <~ val1        0.7163      0.0986    7.2642    0.0000 [ 0.5064; 0.8772 ] 
    ##   VAL <~ val2        0.2202      0.0923    2.3849    0.0171 [ 0.0412; 0.4084 ] 
    ##   VAL <~ val3        0.2082      0.0630    3.3027    0.0010 [ 0.0843; 0.3328 ] 
    ##   SAT <~ sat1        0.3209      0.0149   21.5387    0.0000 [ 0.2967; 0.3535 ] 
    ##   SAT <~ sat2        0.3059      0.0138   22.2139    0.0000 [ 0.2815; 0.3372 ] 
    ##   SAT <~ sat3        0.2474      0.0116   21.3898    0.0000 [ 0.2233; 0.2682 ] 
    ##   SAT <~ sat4        0.2692      0.0119   22.6134    0.0000 [ 0.2480; 0.2936 ] 
    ##   LOY <~ loy1        0.3834      0.0262   14.6342    0.0000 [ 0.3349; 0.4317 ] 
    ##   LOY <~ loy2        0.2434      0.0304    8.0024    0.0000 [ 0.1814; 0.2961 ] 
    ##   LOY <~ loy3        0.3812      0.0281   13.5659    0.0000 [ 0.3289; 0.4414 ] 
    ##   LOY <~ loy4        0.2073      0.0388    5.3419    0.0000 [ 0.1267; 0.2838 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0598   10.7702    0.0000 [ 0.5254; 0.7569 ] 
    ##   imag1 ~~ imag3      0.5433      0.0675    8.0480    0.0000 [ 0.4090; 0.6779 ] 
    ##   imag2 ~~ imag3      0.7761      0.0355   21.8927    0.0000 [ 0.6971; 0.8360 ] 
    ##   expe1 ~~ expe2      0.5353      0.0616    8.6905    0.0000 [ 0.4215; 0.6427 ] 
    ##   expe1 ~~ expe3      0.4694      0.0597    7.8638    0.0000 [ 0.3503; 0.5790 ] 
    ##   expe2 ~~ expe3      0.5467      0.0595    9.1919    0.0000 [ 0.4212; 0.6554 ] 
    ##   qual1 ~~ qual2      0.6053      0.0567   10.6680    0.0000 [ 0.4890; 0.7094 ] 
    ##   qual1 ~~ qual3      0.5406      0.0610    8.8652    0.0000 [ 0.4107; 0.6477 ] 
    ##   qual1 ~~ qual4      0.5662      0.0679    8.3446    0.0000 [ 0.4268; 0.6960 ] 
    ##   qual1 ~~ qual5      0.5180      0.0662    7.8261    0.0000 [ 0.3891; 0.6449 ] 
    ##   qual2 ~~ qual3      0.6187      0.0535   11.5712    0.0000 [ 0.5103; 0.7114 ] 
    ##   qual2 ~~ qual4      0.6517      0.0582   11.1971    0.0000 [ 0.5355; 0.7581 ] 
    ##   qual2 ~~ qual5      0.6291      0.0528   11.9090    0.0000 [ 0.5231; 0.7268 ] 
    ##   qual3 ~~ qual4      0.4752      0.0611    7.7713    0.0000 [ 0.3475; 0.5877 ] 
    ##   qual3 ~~ qual5      0.5074      0.0616    8.2408    0.0000 [ 0.3876; 0.6157 ] 
    ##   qual4 ~~ qual5      0.6402      0.0535   11.9592    0.0000 [ 0.5303; 0.7327 ] 
    ##   val1 ~~ val2        0.6344      0.0562   11.2853    0.0000 [ 0.5211; 0.7348 ] 
    ##   val1 ~~ val3        0.4602      0.0674    6.8267    0.0000 [ 0.3294; 0.5853 ] 
    ##   val2 ~~ val3        0.6288      0.0650    9.6814    0.0000 [ 0.5024; 0.7453 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0643    7.3336    0.0000 [ 0.3484; 0.5993 ] 
    ##   QUAL ~ IMAG       0.3933      0.0597    6.5853    0.0000 [ 0.2820; 0.5114 ] 
    ##   QUAL ~ EXPE       0.8344      0.0226   36.9133    0.0000 [ 0.7851; 0.8750 ] 
    ##   VAL ~ IMAG        0.2974      0.0596    4.9931    0.0000 [ 0.1906; 0.4214 ] 
    ##   VAL ~ EXPE        0.6309      0.0500   12.6248    0.0000 [ 0.5399; 0.7220 ] 
    ##   VAL ~ QUAL        0.7013      0.0799    8.7793    0.0000 [ 0.5390; 0.8465 ] 
    ##   SAT ~ IMAG        0.4807      0.0646    7.4412    0.0000 [ 0.3452; 0.6007 ] 
    ##   SAT ~ EXPE        0.5001      0.0582    8.5897    0.0000 [ 0.3924; 0.6172 ] 
    ##   SAT ~ QUAL        0.5911      0.0901    6.5571    0.0000 [ 0.4060; 0.7486 ] 
    ##   SAT ~ VAL         0.5270      0.0834    6.3160    0.0000 [ 0.3577; 0.6815 ] 
    ##   LOY ~ IMAG        0.4840      0.0686    7.0590    0.0000 [ 0.3587; 0.6166 ] 
    ##   LOY ~ EXPE        0.3142      0.0514    6.1177    0.0000 [ 0.2288; 0.4248 ] 
    ##   LOY ~ QUAL        0.3714      0.0809    4.5888    0.0000 [ 0.2228; 0.5371 ] 
    ##   LOY ~ VAL         0.3311      0.0717    4.6148    0.0000 [ 0.2028; 0.4791 ] 
    ##   LOY ~ SAT         0.6283      0.0783    8.0217    0.0000 [ 0.4865; 0.7824 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0597    6.5853    0.0000 [ 0.2820; 0.5114 ] 
    ##   VAL ~ IMAG           0.2974      0.0596    4.9931    0.0000 [ 0.1906; 0.4214 ] 
    ##   VAL ~ EXPE           0.5852      0.0670    8.7303    0.0000 [ 0.4456; 0.7087 ] 
    ##   SAT ~ IMAG           0.2357      0.0486    4.8504    0.0000 [ 0.1562; 0.3375 ] 
    ##   SAT ~ EXPE           0.5173      0.0622    8.3233    0.0000 [ 0.4020; 0.6406 ] 
    ##   SAT ~ QUAL           0.3696      0.0630    5.8642    0.0000 [ 0.2483; 0.4973 ] 
    ##   LOY ~ IMAG           0.3020      0.0529    5.7087    0.0000 [ 0.2157; 0.4232 ] 
    ##   LOY ~ EXPE           0.3142      0.0514    6.1177    0.0000 [ 0.2288; 0.4248 ] 
    ##   LOY ~ QUAL           0.3714      0.0809    4.5888    0.0000 [ 0.2228; 0.5371 ] 
    ##   LOY ~ VAL            0.3311      0.0717    4.6148    0.0000 [ 0.2028; 0.4791 ] 
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
