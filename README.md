
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
    ##  dG                      0.6493      0.3314  
    ##  SRMR                    0.0940      0.0526  
    ##  dL                      2.2340      0.7013  
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
    ##  Out of 499 bootstrap replications 479 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -2103018466
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
    ##   sat1         1.3527         1.2311       1.7877         1.6177       0.2213
    ##   sat2         1.3062         1.1956       1.7618         1.6270       0.2025
    ##   sat3         1.4076         1.2738       1.7516         1.7188       0.1339
    ##   sat4         1.4155         1.2609       1.7829         1.6328       0.1734
    ##   loy1         1.7725         1.6572       2.2963         2.2222       0.2257
    ##   loy2         1.5173         1.4723       1.9402         1.9806       0.1051
    ##   loy3         1.7859         1.6699       2.3514         2.2690       0.2250
    ##   loy4         1.7106         1.6661       2.1957         2.2964       0.0738
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
    ##  Number of admissible results     = 486
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -1748932873
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
    ##   EXPE ~ IMAG      0.4714      0.0633    7.4490    0.0000 [ 0.3511; 0.5875 ] 
    ##   QUAL ~ EXPE      0.8344      0.0246   33.9290    0.0000 [ 0.7766; 0.8737 ] 
    ##   VAL ~ EXPE       0.0457      0.0814    0.5617    0.5743 [-0.1045; 0.2178 ] 
    ##   VAL ~ QUAL       0.7013      0.0778    9.0162    0.0000 [ 0.5538; 0.8444 ] 
    ##   SAT ~ IMAG       0.2450      0.0559    4.3785    0.0000 [ 0.1314; 0.3582 ] 
    ##   SAT ~ EXPE      -0.0172      0.0686   -0.2511    0.8017 [-0.1584; 0.1221 ] 
    ##   SAT ~ QUAL       0.2215      0.0965    2.2953    0.0217 [ 0.0486; 0.4345 ] 
    ##   SAT ~ VAL        0.5270      0.0840    6.2699    0.0000 [ 0.3685; 0.6925 ] 
    ##   LOY ~ IMAG       0.1819      0.0790    2.3040    0.0212 [ 0.0331; 0.3423 ] 
    ##   LOY ~ SAT        0.6283      0.0807    7.7836    0.0000 [ 0.4625; 0.7802 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1038    6.0741    0.0000 [ 0.4086; 0.8318 ] 
    ##   IMAG =~ imag2      0.9246      0.0410   22.5692    0.0000 [ 0.8207; 0.9814 ] 
    ##   IMAG =~ imag3      0.9577      0.0299   31.9847    0.0000 [ 0.8759; 0.9912 ] 
    ##   EXPE =~ expe1      0.7525      0.0839    8.9667    0.0000 [ 0.5703; 0.8841 ] 
    ##   EXPE =~ expe2      0.9348      0.0295   31.7140    0.0000 [ 0.8628; 0.9720 ] 
    ##   EXPE =~ expe3      0.7295      0.0747    9.7625    0.0000 [ 0.5604; 0.8560 ] 
    ##   QUAL =~ qual1      0.7861      0.0713   11.0314    0.0000 [ 0.6186; 0.8936 ] 
    ##   QUAL =~ qual2      0.9244      0.0243   38.0555    0.0000 [ 0.8684; 0.9606 ] 
    ##   QUAL =~ qual3      0.7560      0.0648   11.6712    0.0000 [ 0.6085; 0.8561 ] 
    ##   QUAL =~ qual4      0.7632      0.0507   15.0540    0.0000 [ 0.6519; 0.8420 ] 
    ##   QUAL =~ qual5      0.7834      0.0470   16.6830    0.0000 [ 0.6778; 0.8593 ] 
    ##   VAL =~ val1        0.9518      0.0225   42.2948    0.0000 [ 0.9030; 0.9846 ] 
    ##   VAL =~ val2        0.8056      0.0661   12.1867    0.0000 [ 0.6625; 0.9100 ] 
    ##   VAL =~ val3        0.6763      0.0725    9.3301    0.0000 [ 0.5170; 0.7959 ] 
    ##   SAT =~ sat1        0.9243      0.0226   40.9744    0.0000 [ 0.8760; 0.9620 ] 
    ##   SAT =~ sat2        0.8813      0.0287   30.6734    0.0000 [ 0.8154; 0.9274 ] 
    ##   SAT =~ sat3        0.7127      0.0522   13.6423    0.0000 [ 0.5939; 0.8087 ] 
    ##   SAT =~ sat4        0.7756      0.0490   15.8170    0.0000 [ 0.6735; 0.8661 ] 
    ##   LOY =~ loy1        0.9097      0.0498   18.2830    0.0000 [ 0.7857; 0.9868 ] 
    ##   LOY =~ loy2        0.5775      0.0911    6.3406    0.0000 [ 0.4139; 0.7439 ] 
    ##   LOY =~ loy3        0.9043      0.0427   21.1993    0.0000 [ 0.8120; 0.9825 ] 
    ##   LOY =~ loy4        0.4917      0.1033    4.7598    0.0000 [ 0.2832; 0.6886 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1232    0.1269    0.8990 [-0.1881; 0.2941 ] 
    ##   IMAG <~ imag2      0.4473      0.1530    2.9239    0.0035 [ 0.1325; 0.7433 ] 
    ##   IMAG <~ imag3      0.6020      0.1456    4.1352    0.0000 [ 0.2886; 0.8481 ] 
    ##   EXPE <~ expe1      0.2946      0.1208    2.4395    0.0147 [ 0.0785; 0.5386 ] 
    ##   EXPE <~ expe2      0.6473      0.0895    7.2327    0.0000 [ 0.4545; 0.7953 ] 
    ##   EXPE <~ expe3      0.2374      0.0978    2.4261    0.0153 [ 0.0431; 0.4129 ] 
    ##   QUAL <~ qual1      0.2370      0.0958    2.4750    0.0133 [ 0.0569; 0.4521 ] 
    ##   QUAL <~ qual2      0.4712      0.0809    5.8262    0.0000 [ 0.3191; 0.6375 ] 
    ##   QUAL <~ qual3      0.1831      0.0839    2.1816    0.0291 [ 0.0087; 0.3354 ] 
    ##   QUAL <~ qual4      0.1037      0.0591    1.7556    0.0792 [-0.0125; 0.2221 ] 
    ##   QUAL <~ qual5      0.2049      0.0650    3.1493    0.0016 [ 0.0581; 0.3226 ] 
    ##   VAL <~ val1        0.7163      0.0930    7.7039    0.0000 [ 0.5373; 0.8760 ] 
    ##   VAL <~ val2        0.2202      0.0938    2.3489    0.0188 [ 0.0541; 0.4171 ] 
    ##   VAL <~ val3        0.2082      0.0589    3.5316    0.0004 [ 0.0842; 0.3160 ] 
    ##   SAT <~ sat1        0.3209      0.0145   22.1118    0.0000 [ 0.2946; 0.3494 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5108    0.0000 [ 0.2799; 0.3350 ] 
    ##   SAT <~ sat3        0.2474      0.0109   22.6287    0.0000 [ 0.2243; 0.2698 ] 
    ##   SAT <~ sat4        0.2692      0.0124   21.6631    0.0000 [ 0.2469; 0.2961 ] 
    ##   LOY <~ loy1        0.3834      0.0273   14.0273    0.0000 [ 0.3288; 0.4340 ] 
    ##   LOY <~ loy2        0.2434      0.0320    7.6185    0.0000 [ 0.1828; 0.2995 ] 
    ##   LOY <~ loy3        0.3812      0.0286   13.3206    0.0000 [ 0.3247; 0.4381 ] 
    ##   LOY <~ loy4        0.2073      0.0386    5.3754    0.0000 [ 0.1288; 0.2813 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0631   10.1963    0.0000 [ 0.5090; 0.7515 ] 
    ##   imag1 ~~ imag3      0.5433      0.0693    7.8404    0.0000 [ 0.3929; 0.6724 ] 
    ##   imag2 ~~ imag3      0.7761      0.0368   21.0658    0.0000 [ 0.6999; 0.8406 ] 
    ##   expe1 ~~ expe2      0.5353      0.0676    7.9208    0.0000 [ 0.3900; 0.6478 ] 
    ##   expe1 ~~ expe3      0.4694      0.0600    7.8272    0.0000 [ 0.3491; 0.5838 ] 
    ##   expe2 ~~ expe3      0.5467      0.0629    8.6962    0.0000 [ 0.4062; 0.6544 ] 
    ##   qual1 ~~ qual2      0.6053      0.0588   10.2970    0.0000 [ 0.4792; 0.7099 ] 
    ##   qual1 ~~ qual3      0.5406      0.0621    8.7026    0.0000 [ 0.4140; 0.6539 ] 
    ##   qual1 ~~ qual4      0.5662      0.0648    8.7382    0.0000 [ 0.4190; 0.6866 ] 
    ##   qual1 ~~ qual5      0.5180      0.0663    7.8168    0.0000 [ 0.3765; 0.6273 ] 
    ##   qual2 ~~ qual3      0.6187      0.0572   10.8080    0.0000 [ 0.4971; 0.7121 ] 
    ##   qual2 ~~ qual4      0.6517      0.0600   10.8618    0.0000 [ 0.5209; 0.7525 ] 
    ##   qual2 ~~ qual5      0.6291      0.0572   10.9984    0.0000 [ 0.5108; 0.7322 ] 
    ##   qual3 ~~ qual4      0.4752      0.0666    7.1348    0.0000 [ 0.3320; 0.5917 ] 
    ##   qual3 ~~ qual5      0.5074      0.0648    7.8320    0.0000 [ 0.3627; 0.6324 ] 
    ##   qual4 ~~ qual5      0.6402      0.0551   11.6186    0.0000 [ 0.5282; 0.7350 ] 
    ##   val1 ~~ val2        0.6344      0.0568   11.1642    0.0000 [ 0.5170; 0.7321 ] 
    ##   val1 ~~ val3        0.4602      0.0710    6.4793    0.0000 [ 0.3180; 0.5957 ] 
    ##   val2 ~~ val3        0.6288      0.0614   10.2341    0.0000 [ 0.5071; 0.7439 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0633    7.4490    0.0000 [ 0.3511; 0.5875 ] 
    ##   QUAL ~ IMAG       0.3933      0.0592    6.6431    0.0000 [ 0.2757; 0.5032 ] 
    ##   QUAL ~ EXPE       0.8344      0.0246   33.9290    0.0000 [ 0.7766; 0.8737 ] 
    ##   VAL ~ IMAG        0.2974      0.0589    5.0452    0.0000 [ 0.1899; 0.4184 ] 
    ##   VAL ~ EXPE        0.6309      0.0514   12.2759    0.0000 [ 0.5314; 0.7323 ] 
    ##   VAL ~ QUAL        0.7013      0.0778    9.0162    0.0000 [ 0.5538; 0.8444 ] 
    ##   SAT ~ IMAG        0.4807      0.0657    7.3120    0.0000 [ 0.3573; 0.6075 ] 
    ##   SAT ~ EXPE        0.5001      0.0588    8.5103    0.0000 [ 0.3756; 0.6142 ] 
    ##   SAT ~ QUAL        0.5911      0.0889    6.6458    0.0000 [ 0.4182; 0.7738 ] 
    ##   SAT ~ VAL         0.5270      0.0840    6.2699    0.0000 [ 0.3685; 0.6925 ] 
    ##   LOY ~ IMAG        0.4840      0.0645    7.4992    0.0000 [ 0.3587; 0.6023 ] 
    ##   LOY ~ EXPE        0.3142      0.0540    5.8147    0.0000 [ 0.2164; 0.4210 ] 
    ##   LOY ~ QUAL        0.3714      0.0783    4.7461    0.0000 [ 0.2316; 0.5586 ] 
    ##   LOY ~ VAL         0.3311      0.0742    4.4624    0.0000 [ 0.1989; 0.4827 ] 
    ##   LOY ~ SAT         0.6283      0.0807    7.7836    0.0000 [ 0.4625; 0.7802 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0592    6.6431    0.0000 [ 0.2757; 0.5032 ] 
    ##   VAL ~ IMAG           0.2974      0.0589    5.0452    0.0000 [ 0.1899; 0.4184 ] 
    ##   VAL ~ EXPE           0.5852      0.0668    8.7638    0.0000 [ 0.4491; 0.7083 ] 
    ##   SAT ~ IMAG           0.2357      0.0475    4.9607    0.0000 [ 0.1551; 0.3415 ] 
    ##   SAT ~ EXPE           0.5173      0.0639    8.0908    0.0000 [ 0.3984; 0.6438 ] 
    ##   SAT ~ QUAL           0.3696      0.0585    6.3176    0.0000 [ 0.2482; 0.4896 ] 
    ##   LOY ~ IMAG           0.3020      0.0550    5.4944    0.0000 [ 0.2070; 0.4143 ] 
    ##   LOY ~ EXPE           0.3142      0.0540    5.8147    0.0000 [ 0.2164; 0.4210 ] 
    ##   LOY ~ QUAL           0.3714      0.0783    4.7461    0.0000 [ 0.2316; 0.5586 ] 
    ##   LOY ~ VAL            0.3311      0.0742    4.4624    0.0000 [ 0.1989; 0.4827 ] 
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
