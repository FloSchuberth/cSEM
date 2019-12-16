
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
    ##   +------------------------------------------------------------+
    ##   |                                                            |
    ##   |   H0: Population indicator covariance matrix is equal to   |
    ##   |   model-implied indicator covariance matrix.               |
    ##   |                                                            |
    ##   +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3291  
    ##  SRMR                    0.0940      0.0530  
    ##  dL                      2.2340      0.7095  
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
    ##  Out of 499 bootstrap replications 471 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -617035516
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
    ##   sat1         1.3518         1.2346       1.7866         1.6199       0.2217
    ##   sat2         1.3059         1.1946       1.7632         1.6246       0.2019
    ##   sat3         1.4062         1.2761       1.7486         1.7201       0.1357
    ##   sat4         1.4141         1.2610       1.7819         1.6334       0.1740
    ##   loy1         1.7724         1.6617       2.2955         2.2282       0.2273
    ##   loy2         1.5157         1.4732       1.9373         1.9798       0.1080
    ##   loy3         1.7836         1.6702       2.3472         2.2706       0.2264
    ##   loy4         1.7108         1.6674       2.1976         2.2987       0.0711
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
    ##  Number of admissible results     = 490
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 142891711
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
    ##   EXPE ~ IMAG      0.4714      0.0675    6.9820    0.0000 [ 0.3353; 0.5952 ] 
    ##   QUAL ~ EXPE      0.8344      0.0237   35.1734    0.0000 [ 0.7802; 0.8736 ] 
    ##   VAL ~ EXPE       0.0457      0.0875    0.5223    0.6015 [-0.1105; 0.2162 ] 
    ##   VAL ~ QUAL       0.7013      0.0838    8.3696    0.0000 [ 0.5195; 0.8537 ] 
    ##   SAT ~ IMAG       0.2450      0.0557    4.3988    0.0000 [ 0.1439; 0.3527 ] 
    ##   SAT ~ EXPE      -0.0172      0.0708   -0.2435    0.8076 [-0.1686; 0.1116 ] 
    ##   SAT ~ QUAL       0.2215      0.1033    2.1447    0.0320 [ 0.0488; 0.4429 ] 
    ##   SAT ~ VAL        0.5270      0.0875    6.0235    0.0000 [ 0.3505; 0.6820 ] 
    ##   LOY ~ IMAG       0.1819      0.0806    2.2571    0.0240 [ 0.0307; 0.3479 ] 
    ##   LOY ~ SAT        0.6283      0.0856    7.3437    0.0000 [ 0.4394; 0.7805 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0972    6.4846    0.0000 [ 0.4208; 0.7998 ] 
    ##   IMAG =~ imag2      0.9246      0.0427   21.6616    0.0000 [ 0.8231; 0.9779 ] 
    ##   IMAG =~ imag3      0.9577      0.0291   32.9653    0.0000 [ 0.8877; 0.9919 ] 
    ##   EXPE =~ expe1      0.7525      0.0803    9.3745    0.0000 [ 0.5844; 0.8862 ] 
    ##   EXPE =~ expe2      0.9348      0.0282   33.1816    0.0000 [ 0.8595; 0.9720 ] 
    ##   EXPE =~ expe3      0.7295      0.0771    9.4646    0.0000 [ 0.5500; 0.8615 ] 
    ##   QUAL =~ qual1      0.7861      0.0693   11.3467    0.0000 [ 0.6358; 0.8936 ] 
    ##   QUAL =~ qual2      0.9244      0.0236   39.1099    0.0000 [ 0.8661; 0.9598 ] 
    ##   QUAL =~ qual3      0.7560      0.0625   12.0987    0.0000 [ 0.6098; 0.8546 ] 
    ##   QUAL =~ qual4      0.7632      0.0541   14.0954    0.0000 [ 0.6502; 0.8512 ] 
    ##   QUAL =~ qual5      0.7834      0.0468   16.7491    0.0000 [ 0.6779; 0.8576 ] 
    ##   VAL =~ val1        0.9518      0.0226   42.1467    0.0000 [ 0.8988; 0.9848 ] 
    ##   VAL =~ val2        0.8056      0.0665   12.1218    0.0000 [ 0.6512; 0.9075 ] 
    ##   VAL =~ val3        0.6763      0.0765    8.8402    0.0000 [ 0.5086; 0.7962 ] 
    ##   SAT =~ sat1        0.9243      0.0233   39.6140    0.0000 [ 0.8731; 0.9631 ] 
    ##   SAT =~ sat2        0.8813      0.0302   29.1501    0.0000 [ 0.8201; 0.9294 ] 
    ##   SAT =~ sat3        0.7127      0.0515   13.8392    0.0000 [ 0.5953; 0.8049 ] 
    ##   SAT =~ sat4        0.7756      0.0533   14.5480    0.0000 [ 0.6601; 0.8691 ] 
    ##   LOY =~ loy1        0.9097      0.0499   18.2272    0.0000 [ 0.7975; 0.9887 ] 
    ##   LOY =~ loy2        0.5775      0.0864    6.6856    0.0000 [ 0.4031; 0.7384 ] 
    ##   LOY =~ loy3        0.9043      0.0424   21.3043    0.0000 [ 0.8116; 0.9839 ] 
    ##   LOY =~ loy4        0.4917      0.1006    4.8899    0.0000 [ 0.2968; 0.6853 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1123    0.1393    0.8892 [-0.1894; 0.2370 ] 
    ##   IMAG <~ imag2      0.4473      0.1556    2.8755    0.0040 [ 0.1274; 0.7352 ] 
    ##   IMAG <~ imag3      0.6020      0.1448    4.1587    0.0000 [ 0.3022; 0.8651 ] 
    ##   EXPE <~ expe1      0.2946      0.1221    2.4130    0.0158 [ 0.0666; 0.5326 ] 
    ##   EXPE <~ expe2      0.6473      0.0842    7.6893    0.0000 [ 0.4721; 0.8029 ] 
    ##   EXPE <~ expe3      0.2374      0.0982    2.4165    0.0157 [ 0.0342; 0.4164 ] 
    ##   QUAL <~ qual1      0.2370      0.0913    2.5955    0.0094 [ 0.0753; 0.4231 ] 
    ##   QUAL <~ qual2      0.4712      0.0777    6.0612    0.0000 [ 0.3172; 0.6291 ] 
    ##   QUAL <~ qual3      0.1831      0.0819    2.2344    0.0255 [ 0.0043; 0.3326 ] 
    ##   QUAL <~ qual4      0.1037      0.0622    1.6680    0.0953 [-0.0071; 0.2317 ] 
    ##   QUAL <~ qual5      0.2049      0.0619    3.3083    0.0009 [ 0.0743; 0.3135 ] 
    ##   VAL <~ val1        0.7163      0.0960    7.4609    0.0000 [ 0.5233; 0.8805 ] 
    ##   VAL <~ val2        0.2202      0.0940    2.3436    0.0191 [ 0.0357; 0.4031 ] 
    ##   VAL <~ val3        0.2082      0.0588    3.5412    0.0004 [ 0.0938; 0.3208 ] 
    ##   SAT <~ sat1        0.3209      0.0162   19.7495    0.0000 [ 0.2938; 0.3601 ] 
    ##   SAT <~ sat2        0.3059      0.0137   22.2560    0.0000 [ 0.2824; 0.3336 ] 
    ##   SAT <~ sat3        0.2474      0.0109   22.6103    0.0000 [ 0.2264; 0.2689 ] 
    ##   SAT <~ sat4        0.2692      0.0122   21.9863    0.0000 [ 0.2449; 0.2935 ] 
    ##   LOY <~ loy1        0.3834      0.0273   14.0697    0.0000 [ 0.3311; 0.4378 ] 
    ##   LOY <~ loy2        0.2434      0.0300    8.1078    0.0000 [ 0.1809; 0.2964 ] 
    ##   LOY <~ loy3        0.3812      0.0285   13.3679    0.0000 [ 0.3318; 0.4388 ] 
    ##   LOY <~ loy4        0.2073      0.0368    5.6352    0.0000 [ 0.1323; 0.2756 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0629   10.2330    0.0000 [ 0.5143; 0.7576 ] 
    ##   imag1 ~~ imag3      0.5433      0.0671    8.0927    0.0000 [ 0.4022; 0.6717 ] 
    ##   imag2 ~~ imag3      0.7761      0.0386   20.1112    0.0000 [ 0.7020; 0.8473 ] 
    ##   expe1 ~~ expe2      0.5353      0.0632    8.4634    0.0000 [ 0.4144; 0.6565 ] 
    ##   expe1 ~~ expe3      0.4694      0.0642    7.3159    0.0000 [ 0.3373; 0.5899 ] 
    ##   expe2 ~~ expe3      0.5467      0.0641    8.5335    0.0000 [ 0.4077; 0.6607 ] 
    ##   qual1 ~~ qual2      0.6053      0.0598   10.1179    0.0000 [ 0.4814; 0.7104 ] 
    ##   qual1 ~~ qual3      0.5406      0.0622    8.6915    0.0000 [ 0.4041; 0.6601 ] 
    ##   qual1 ~~ qual4      0.5662      0.0666    8.5015    0.0000 [ 0.4360; 0.6834 ] 
    ##   qual1 ~~ qual5      0.5180      0.0672    7.7100    0.0000 [ 0.3740; 0.6383 ] 
    ##   qual2 ~~ qual3      0.6187      0.0548   11.2919    0.0000 [ 0.4992; 0.7152 ] 
    ##   qual2 ~~ qual4      0.6517      0.0600   10.8659    0.0000 [ 0.5363; 0.7630 ] 
    ##   qual2 ~~ qual5      0.6291      0.0574   10.9608    0.0000 [ 0.5128; 0.7335 ] 
    ##   qual3 ~~ qual4      0.4752      0.0635    7.4792    0.0000 [ 0.3396; 0.5901 ] 
    ##   qual3 ~~ qual5      0.5074      0.0636    7.9803    0.0000 [ 0.3776; 0.6278 ] 
    ##   qual4 ~~ qual5      0.6402      0.0561   11.4053    0.0000 [ 0.5275; 0.7411 ] 
    ##   val1 ~~ val2        0.6344      0.0562   11.2934    0.0000 [ 0.5216; 0.7359 ] 
    ##   val1 ~~ val3        0.4602      0.0725    6.3511    0.0000 [ 0.3084; 0.5947 ] 
    ##   val2 ~~ val3        0.6288      0.0646    9.7369    0.0000 [ 0.4976; 0.7403 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0675    6.9820    0.0000 [ 0.3353; 0.5952 ] 
    ##   QUAL ~ IMAG       0.3933      0.0628    6.2601    0.0000 [ 0.2663; 0.5157 ] 
    ##   QUAL ~ EXPE       0.8344      0.0237   35.1734    0.0000 [ 0.7802; 0.8736 ] 
    ##   VAL ~ IMAG        0.2974      0.0608    4.8885    0.0000 [ 0.1854; 0.4301 ] 
    ##   VAL ~ EXPE        0.6309      0.0497   12.7065    0.0000 [ 0.5283; 0.7273 ] 
    ##   VAL ~ QUAL        0.7013      0.0838    8.3696    0.0000 [ 0.5195; 0.8537 ] 
    ##   SAT ~ IMAG        0.4807      0.0682    7.0489    0.0000 [ 0.3279; 0.6053 ] 
    ##   SAT ~ EXPE        0.5001      0.0569    8.7821    0.0000 [ 0.3905; 0.6089 ] 
    ##   SAT ~ QUAL        0.5911      0.0970    6.0928    0.0000 [ 0.3893; 0.7755 ] 
    ##   SAT ~ VAL         0.5270      0.0875    6.0235    0.0000 [ 0.3505; 0.6820 ] 
    ##   LOY ~ IMAG        0.4840      0.0675    7.1716    0.0000 [ 0.3524; 0.6121 ] 
    ##   LOY ~ EXPE        0.3142      0.0560    5.6061    0.0000 [ 0.2030; 0.4229 ] 
    ##   LOY ~ QUAL        0.3714      0.0875    4.2447    0.0000 [ 0.2001; 0.5499 ] 
    ##   LOY ~ VAL         0.3311      0.0759    4.3608    0.0000 [ 0.1876; 0.4955 ] 
    ##   LOY ~ SAT         0.6283      0.0856    7.3437    0.0000 [ 0.4394; 0.7805 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0628    6.2601    0.0000 [ 0.2663; 0.5157 ] 
    ##   VAL ~ IMAG           0.2974      0.0608    4.8885    0.0000 [ 0.1854; 0.4301 ] 
    ##   VAL ~ EXPE           0.5852      0.0707    8.2762    0.0000 [ 0.4351; 0.7151 ] 
    ##   SAT ~ IMAG           0.2357      0.0488    4.8283    0.0000 [ 0.1480; 0.3393 ] 
    ##   SAT ~ EXPE           0.5173      0.0697    7.4269    0.0000 [ 0.3916; 0.6631 ] 
    ##   SAT ~ QUAL           0.3696      0.0621    5.9498    0.0000 [ 0.2481; 0.4949 ] 
    ##   LOY ~ IMAG           0.3020      0.0578    5.2212    0.0000 [ 0.1790; 0.4205 ] 
    ##   LOY ~ EXPE           0.3142      0.0560    5.6061    0.0000 [ 0.2030; 0.4229 ] 
    ##   LOY ~ QUAL           0.3714      0.0875    4.2447    0.0000 [ 0.2001; 0.5499 ] 
    ##   LOY ~ VAL            0.3311      0.0759    4.3608    0.0000 [ 0.1876; 0.4955 ] 
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
