
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
    ##                                                          +------------------------------------------------------------+
    ##                                                          |                                                            |
    ##                                                          |   H0: Population indicator covariance matrix is equal to   |
    ##                                                          |   model-implied indicator covariance matrix.               |
    ##                                                          |                                                            |
    ##                                                          +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3287  
    ##  SRMR                    0.0940      0.0529  
    ##  dL                      2.2340      0.7092  
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
    ##  Out of 499 bootstrap replications 477 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -693085397
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
    ##  Number of obs. total             = 250
    ##  Number of obs. training          = 225
    ##  Number of obs. test              = 25
    ##  Number of cv folds               = 10
    ##  Number of repetitions            = 10
    ##  Handle inadmissibles             = stop
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name     MAE_target  RMSE_target    MAE_lm   RMSE_lm   Q2_predict
    ##   sat1         1.3493       1.7847    1.2303    1.6163      -0.2192
    ##   sat2         1.3056       1.7611    1.1977    1.6278      -0.1706
    ##   sat3         1.4078       1.7513    1.2784    1.7232      -0.0328
    ##   sat4         1.4146       1.7824    1.2606    1.6334      -0.1907
    ##   loy1         1.7698       2.2937    1.6595    2.2248      -0.0629
    ##   loy2         1.5144       1.9367    1.4731    1.9788       0.0420
    ##   loy3         1.7829       2.3461    1.6688    2.2697      -0.0685
    ##   loy4         1.7084       2.1926    1.6712    2.3019       0.0927
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
    ##  Number of admissible results     = 485
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -458200941
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
    ##   EXPE ~ IMAG      0.4714      0.0639    7.3770    0.0000 [ 0.3603; 0.6105 ] 
    ##   QUAL ~ EXPE      0.8344      0.0237   35.2256    0.0000 [ 0.7848; 0.8766 ] 
    ##   VAL ~ EXPE       0.0457      0.0867    0.5274    0.5979 [-0.1016; 0.2350 ] 
    ##   VAL ~ QUAL       0.7013      0.0823    8.5256    0.0000 [ 0.5349; 0.8445 ] 
    ##   SAT ~ IMAG       0.2450      0.0546    4.4826    0.0000 [ 0.1454; 0.3607 ] 
    ##   SAT ~ EXPE      -0.0172      0.0711   -0.2422    0.8086 [-0.1537; 0.1161 ] 
    ##   SAT ~ QUAL       0.2215      0.1071    2.0678    0.0387 [ 0.0512; 0.4615 ] 
    ##   SAT ~ VAL        0.5270      0.0867    6.0809    0.0000 [ 0.3528; 0.6842 ] 
    ##   LOY ~ IMAG       0.1819      0.0813    2.2373    0.0253 [ 0.0406; 0.3601 ] 
    ##   LOY ~ SAT        0.6283      0.0800    7.8554    0.0000 [ 0.4601; 0.7746 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0992    6.3582    0.0000 [ 0.4323; 0.8048 ] 
    ##   IMAG =~ imag2      0.9246      0.0406   22.7582    0.0000 [ 0.8258; 0.9786 ] 
    ##   IMAG =~ imag3      0.9577      0.0311   30.8158    0.0000 [ 0.8745; 0.9919 ] 
    ##   EXPE =~ expe1      0.7525      0.0786    9.5691    0.0000 [ 0.5550; 0.8712 ] 
    ##   EXPE =~ expe2      0.9348      0.0270   34.6345    0.0000 [ 0.8719; 0.9739 ] 
    ##   EXPE =~ expe3      0.7295      0.0700   10.4183    0.0000 [ 0.5678; 0.8470 ] 
    ##   QUAL =~ qual1      0.7861      0.0706   11.1321    0.0000 [ 0.6153; 0.8846 ] 
    ##   QUAL =~ qual2      0.9244      0.0225   41.0016    0.0000 [ 0.8688; 0.9562 ] 
    ##   QUAL =~ qual3      0.7560      0.0608   12.4258    0.0000 [ 0.6268; 0.8623 ] 
    ##   QUAL =~ qual4      0.7632      0.0552   13.8344    0.0000 [ 0.6331; 0.8592 ] 
    ##   QUAL =~ qual5      0.7834      0.0471   16.6371    0.0000 [ 0.6820; 0.8611 ] 
    ##   VAL =~ val1        0.9518      0.0235   40.5873    0.0000 [ 0.8966; 0.9854 ] 
    ##   VAL =~ val2        0.8056      0.0664   12.1348    0.0000 [ 0.6526; 0.9194 ] 
    ##   VAL =~ val3        0.6763      0.0707    9.5593    0.0000 [ 0.5214; 0.8051 ] 
    ##   SAT =~ sat1        0.9243      0.0228   40.5356    0.0000 [ 0.8733; 0.9647 ] 
    ##   SAT =~ sat2        0.8813      0.0297   29.6257    0.0000 [ 0.8172; 0.9258 ] 
    ##   SAT =~ sat3        0.7127      0.0477   14.9430    0.0000 [ 0.6126; 0.7951 ] 
    ##   SAT =~ sat4        0.7756      0.0522   14.8491    0.0000 [ 0.6621; 0.8667 ] 
    ##   LOY =~ loy1        0.9097      0.0475   19.1561    0.0000 [ 0.8043; 0.9857 ] 
    ##   LOY =~ loy2        0.5775      0.0877    6.5888    0.0000 [ 0.3924; 0.7263 ] 
    ##   LOY =~ loy3        0.9043      0.0438   20.6282    0.0000 [ 0.8090; 0.9766 ] 
    ##   LOY =~ loy4        0.4917      0.0978    5.0291    0.0000 [ 0.3125; 0.6938 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1109    0.1410    0.8879 [-0.2103; 0.2350 ] 
    ##   IMAG <~ imag2      0.4473      0.1550    2.8866    0.0039 [ 0.1337; 0.7459 ] 
    ##   IMAG <~ imag3      0.6020      0.1467    4.1038    0.0000 [ 0.2869; 0.8648 ] 
    ##   EXPE <~ expe1      0.2946      0.1195    2.4653    0.0137 [ 0.0463; 0.4905 ] 
    ##   EXPE <~ expe2      0.6473      0.0809    7.9970    0.0000 [ 0.4796; 0.7926 ] 
    ##   EXPE <~ expe3      0.2374      0.0912    2.6029    0.0092 [ 0.0791; 0.4146 ] 
    ##   QUAL <~ qual1      0.2370      0.0894    2.6507    0.0080 [ 0.0655; 0.4124 ] 
    ##   QUAL <~ qual2      0.4712      0.0751    6.2762    0.0000 [ 0.2927; 0.6108 ] 
    ##   QUAL <~ qual3      0.1831      0.0794    2.3050    0.0212 [ 0.0422; 0.3480 ] 
    ##   QUAL <~ qual4      0.1037      0.0614    1.6897    0.0911 [-0.0117; 0.2236 ] 
    ##   QUAL <~ qual5      0.2049      0.0655    3.1255    0.0018 [ 0.0579; 0.3178 ] 
    ##   VAL <~ val1        0.7163      0.0977    7.3282    0.0000 [ 0.5104; 0.8947 ] 
    ##   VAL <~ val2        0.2202      0.0965    2.2827    0.0224 [ 0.0484; 0.4247 ] 
    ##   VAL <~ val3        0.2082      0.0591    3.5235    0.0004 [ 0.0952; 0.3299 ] 
    ##   SAT <~ sat1        0.3209      0.0152   21.0707    0.0000 [ 0.2948; 0.3552 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5149    0.0000 [ 0.2836; 0.3385 ] 
    ##   SAT <~ sat3        0.2474      0.0105   23.5632    0.0000 [ 0.2260; 0.2660 ] 
    ##   SAT <~ sat4        0.2692      0.0119   22.6449    0.0000 [ 0.2460; 0.2929 ] 
    ##   LOY <~ loy1        0.3834      0.0265   14.4483    0.0000 [ 0.3338; 0.4375 ] 
    ##   LOY <~ loy2        0.2434      0.0300    8.1206    0.0000 [ 0.1796; 0.2937 ] 
    ##   LOY <~ loy3        0.3812      0.0277   13.7841    0.0000 [ 0.3249; 0.4327 ] 
    ##   LOY <~ loy4        0.2073      0.0364    5.6970    0.0000 [ 0.1367; 0.2809 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0621   10.3654    0.0000 [ 0.5218; 0.7624 ] 
    ##   imag1 ~~ imag3      0.5433      0.0693    7.8407    0.0000 [ 0.4189; 0.6868 ] 
    ##   imag2 ~~ imag3      0.7761      0.0395   19.6288    0.0000 [ 0.6982; 0.8520 ] 
    ##   expe1 ~~ expe2      0.5353      0.0609    8.7932    0.0000 [ 0.4270; 0.6589 ] 
    ##   expe1 ~~ expe3      0.4694      0.0610    7.6965    0.0000 [ 0.3522; 0.5853 ] 
    ##   expe2 ~~ expe3      0.5467      0.0615    8.8862    0.0000 [ 0.4264; 0.6617 ] 
    ##   qual1 ~~ qual2      0.6053      0.0585   10.3511    0.0000 [ 0.4904; 0.7144 ] 
    ##   qual1 ~~ qual3      0.5406      0.0626    8.6357    0.0000 [ 0.4155; 0.6576 ] 
    ##   qual1 ~~ qual4      0.5662      0.0712    7.9497    0.0000 [ 0.4048; 0.6906 ] 
    ##   qual1 ~~ qual5      0.5180      0.0716    7.2359    0.0000 [ 0.3649; 0.6429 ] 
    ##   qual2 ~~ qual3      0.6187      0.0572   10.8217    0.0000 [ 0.5012; 0.7143 ] 
    ##   qual2 ~~ qual4      0.6517      0.0627   10.3941    0.0000 [ 0.5264; 0.7674 ] 
    ##   qual2 ~~ qual5      0.6291      0.0578   10.8923    0.0000 [ 0.5145; 0.7302 ] 
    ##   qual3 ~~ qual4      0.4752      0.0651    7.2975    0.0000 [ 0.3399; 0.5845 ] 
    ##   qual3 ~~ qual5      0.5074      0.0626    8.1003    0.0000 [ 0.3766; 0.6167 ] 
    ##   qual4 ~~ qual5      0.6402      0.0569   11.2470    0.0000 [ 0.5276; 0.7382 ] 
    ##   val1 ~~ val2        0.6344      0.0540   11.7425    0.0000 [ 0.5224; 0.7303 ] 
    ##   val1 ~~ val3        0.4602      0.0671    6.8620    0.0000 [ 0.3249; 0.5924 ] 
    ##   val2 ~~ val3        0.6288      0.0599   10.5020    0.0000 [ 0.5034; 0.7433 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0639    7.3770    0.0000 [ 0.3603; 0.6105 ] 
    ##   QUAL ~ IMAG       0.3933      0.0598    6.5766    0.0000 [ 0.2944; 0.5248 ] 
    ##   QUAL ~ EXPE       0.8344      0.0237   35.2256    0.0000 [ 0.7848; 0.8766 ] 
    ##   VAL ~ IMAG        0.2974      0.0592    5.0273    0.0000 [ 0.2018; 0.4351 ] 
    ##   VAL ~ EXPE        0.6309      0.0493   12.7903    0.0000 [ 0.5330; 0.7280 ] 
    ##   VAL ~ QUAL        0.7013      0.0823    8.5256    0.0000 [ 0.5349; 0.8445 ] 
    ##   SAT ~ IMAG        0.4807      0.0653    7.3621    0.0000 [ 0.3576; 0.6141 ] 
    ##   SAT ~ EXPE        0.5001      0.0572    8.7393    0.0000 [ 0.3899; 0.6105 ] 
    ##   SAT ~ QUAL        0.5911      0.0967    6.1152    0.0000 [ 0.3938; 0.7827 ] 
    ##   SAT ~ VAL         0.5270      0.0867    6.0809    0.0000 [ 0.3528; 0.6842 ] 
    ##   LOY ~ IMAG        0.4840      0.0648    7.4691    0.0000 [ 0.3635; 0.6200 ] 
    ##   LOY ~ EXPE        0.3142      0.0531    5.9138    0.0000 [ 0.2160; 0.4155 ] 
    ##   LOY ~ QUAL        0.3714      0.0808    4.5962    0.0000 [ 0.2344; 0.5491 ] 
    ##   LOY ~ VAL         0.3311      0.0720    4.5961    0.0000 [ 0.1867; 0.4819 ] 
    ##   LOY ~ SAT         0.6283      0.0800    7.8554    0.0000 [ 0.4601; 0.7746 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0598    6.5766    0.0000 [ 0.2944; 0.5248 ] 
    ##   VAL ~ IMAG           0.2974      0.0592    5.0273    0.0000 [ 0.2018; 0.4351 ] 
    ##   VAL ~ EXPE           0.5852      0.0704    8.3100    0.0000 [ 0.4438; 0.7111 ] 
    ##   SAT ~ IMAG           0.2357      0.0469    5.0211    0.0000 [ 0.1519; 0.3376 ] 
    ##   SAT ~ EXPE           0.5173      0.0707    7.3188    0.0000 [ 0.4005; 0.6649 ] 
    ##   SAT ~ QUAL           0.3696      0.0594    6.2167    0.0000 [ 0.2492; 0.4801 ] 
    ##   LOY ~ IMAG           0.3020      0.0564    5.3575    0.0000 [ 0.2052; 0.4189 ] 
    ##   LOY ~ EXPE           0.3142      0.0531    5.9138    0.0000 [ 0.2160; 0.4155 ] 
    ##   LOY ~ QUAL           0.3714      0.0808    4.5962    0.0000 [ 0.2344; 0.5491 ] 
    ##   LOY ~ VAL            0.3311      0.0720    4.5961    0.0000 [ 0.1867; 0.4819 ] 
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
