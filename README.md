
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
  - `predict()` : predict indicator values (not yet implemented)
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
    ##  Type of indicator correlation    = Bravais-Pearson
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
# ## Test overall model fit
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
    ##  dG                      0.6493      0.3202  
    ##  SRMR                    0.0940      0.0515  
    ##  dL                      2.2340      0.6718  
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
    ##  Out of 499 bootstrap replications 475 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -758394284
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
    ##  Type of indicator correlation    = Bravais-Pearson
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
    ##  Random seed                      = -1476110302
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
    ##   EXPE ~ IMAG      0.4714      0.0607    7.7601    0.0000 [ 0.3485; 0.5900 ] 
    ##   QUAL ~ EXPE      0.8344      0.0220   37.9301    0.0000 [ 0.7881; 0.8731 ] 
    ##   VAL ~ EXPE       0.0457      0.0811    0.5637    0.5730 [-0.1138; 0.2098 ] 
    ##   VAL ~ QUAL       0.7013      0.0780    8.9881    0.0000 [ 0.5478; 0.8348 ] 
    ##   SAT ~ IMAG       0.2450      0.0551    4.4498    0.0000 [ 0.1353; 0.3568 ] 
    ##   SAT ~ EXPE      -0.0172      0.0773   -0.2230    0.8236 [-0.1727; 0.1287 ] 
    ##   SAT ~ QUAL       0.2215      0.1085    2.0418    0.0412 [ 0.0378; 0.4536 ] 
    ##   SAT ~ VAL        0.5270      0.0894    5.8968    0.0000 [ 0.3473; 0.6909 ] 
    ##   LOY ~ IMAG       0.1819      0.0779    2.3351    0.0195 [ 0.0347; 0.3345 ] 
    ##   LOY ~ SAT        0.6283      0.0797    7.8785    0.0000 [ 0.4663; 0.7880 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0991    6.3650    0.0000 [ 0.4274; 0.7989 ] 
    ##   IMAG =~ imag2      0.9246      0.0396   23.3448    0.0000 [ 0.8307; 0.9813 ] 
    ##   IMAG =~ imag3      0.9577      0.0292   32.8379    0.0000 [ 0.8851; 0.9919 ] 
    ##   EXPE =~ expe1      0.7525      0.0758    9.9285    0.0000 [ 0.5845; 0.8854 ] 
    ##   EXPE =~ expe2      0.9348      0.0280   33.4242    0.0000 [ 0.8604; 0.9715 ] 
    ##   EXPE =~ expe3      0.7295      0.0737    9.9017    0.0000 [ 0.5643; 0.8482 ] 
    ##   QUAL =~ qual1      0.7861      0.0632   12.4305    0.0000 [ 0.6510; 0.8942 ] 
    ##   QUAL =~ qual2      0.9244      0.0231   40.0011    0.0000 [ 0.8683; 0.9589 ] 
    ##   QUAL =~ qual3      0.7560      0.0626   12.0854    0.0000 [ 0.6050; 0.8478 ] 
    ##   QUAL =~ qual4      0.7632      0.0522   14.6221    0.0000 [ 0.6547; 0.8523 ] 
    ##   QUAL =~ qual5      0.7834      0.0458   17.1020    0.0000 [ 0.6806; 0.8644 ] 
    ##   VAL =~ val1        0.9518      0.0225   42.3964    0.0000 [ 0.8966; 0.9851 ] 
    ##   VAL =~ val2        0.8056      0.0648   12.4405    0.0000 [ 0.6561; 0.9076 ] 
    ##   VAL =~ val3        0.6763      0.0762    8.8695    0.0000 [ 0.5099; 0.8064 ] 
    ##   SAT =~ sat1        0.9243      0.0224   41.2588    0.0000 [ 0.8758; 0.9626 ] 
    ##   SAT =~ sat2        0.8813      0.0292   30.1861    0.0000 [ 0.8186; 0.9288 ] 
    ##   SAT =~ sat3        0.7127      0.0542   13.1430    0.0000 [ 0.6096; 0.8057 ] 
    ##   SAT =~ sat4        0.7756      0.0502   15.4361    0.0000 [ 0.6618; 0.8598 ] 
    ##   LOY =~ loy1        0.9097      0.0534   17.0268    0.0000 [ 0.7814; 0.9857 ] 
    ##   LOY =~ loy2        0.5775      0.0879    6.5710    0.0000 [ 0.3750; 0.7278 ] 
    ##   LOY =~ loy3        0.9043      0.0450   20.1058    0.0000 [ 0.7933; 0.9693 ] 
    ##   LOY =~ loy4        0.4917      0.0959    5.1301    0.0000 [ 0.3259; 0.6784 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1133    0.1381    0.8902 [-0.2009; 0.2441 ] 
    ##   IMAG <~ imag2      0.4473      0.1506    2.9699    0.0030 [ 0.1724; 0.7465 ] 
    ##   IMAG <~ imag3      0.6020      0.1402    4.2954    0.0000 [ 0.2941; 0.8290 ] 
    ##   EXPE <~ expe1      0.2946      0.1146    2.5711    0.0101 [ 0.0714; 0.5187 ] 
    ##   EXPE <~ expe2      0.6473      0.0864    7.4898    0.0000 [ 0.4479; 0.7962 ] 
    ##   EXPE <~ expe3      0.2374      0.0924    2.5683    0.0102 [ 0.0525; 0.4016 ] 
    ##   QUAL <~ qual1      0.2370      0.0878    2.7001    0.0069 [ 0.0837; 0.4210 ] 
    ##   QUAL <~ qual2      0.4712      0.0807    5.8428    0.0000 [ 0.3073; 0.6229 ] 
    ##   QUAL <~ qual3      0.1831      0.0805    2.2753    0.0229 [ 0.0128; 0.3219 ] 
    ##   QUAL <~ qual4      0.1037      0.0584    1.7753    0.0759 [-0.0163; 0.2231 ] 
    ##   QUAL <~ qual5      0.2049      0.0653    3.1366    0.0017 [ 0.0663; 0.3200 ] 
    ##   VAL <~ val1        0.7163      0.0960    7.4654    0.0000 [ 0.5175; 0.8778 ] 
    ##   VAL <~ val2        0.2202      0.0932    2.3632    0.0181 [ 0.0616; 0.4116 ] 
    ##   VAL <~ val3        0.2082      0.0570    3.6499    0.0003 [ 0.1020; 0.3286 ] 
    ##   SAT <~ sat1        0.3209      0.0159   20.2190    0.0000 [ 0.2964; 0.3537 ] 
    ##   SAT <~ sat2        0.3059      0.0142   21.4725    0.0000 [ 0.2842; 0.3396 ] 
    ##   SAT <~ sat3        0.2474      0.0110   22.4270    0.0000 [ 0.2257; 0.2674 ] 
    ##   SAT <~ sat4        0.2692      0.0121   22.2020    0.0000 [ 0.2470; 0.2921 ] 
    ##   LOY <~ loy1        0.3834      0.0270   14.1988    0.0000 [ 0.3308; 0.4386 ] 
    ##   LOY <~ loy2        0.2434      0.0315    7.7216    0.0000 [ 0.1697; 0.2943 ] 
    ##   LOY <~ loy3        0.3812      0.0278   13.6928    0.0000 [ 0.3263; 0.4352 ] 
    ##   LOY <~ loy4        0.2073      0.0354    5.8602    0.0000 [ 0.1430; 0.2758 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0665    9.6842    0.0000 [ 0.5028; 0.7552 ] 
    ##   imag1 ~~ imag3      0.5433      0.0693    7.8419    0.0000 [ 0.4081; 0.6755 ] 
    ##   imag2 ~~ imag3      0.7761      0.0398   19.4943    0.0000 [ 0.6867; 0.8428 ] 
    ##   expe1 ~~ expe2      0.5353      0.0597    8.9638    0.0000 [ 0.4153; 0.6483 ] 
    ##   expe1 ~~ expe3      0.4694      0.0590    7.9547    0.0000 [ 0.3600; 0.5871 ] 
    ##   expe2 ~~ expe3      0.5467      0.0603    9.0630    0.0000 [ 0.4237; 0.6608 ] 
    ##   qual1 ~~ qual2      0.6053      0.0547   11.0726    0.0000 [ 0.5038; 0.7063 ] 
    ##   qual1 ~~ qual3      0.5406      0.0577    9.3766    0.0000 [ 0.4212; 0.6483 ] 
    ##   qual1 ~~ qual4      0.5662      0.0637    8.8934    0.0000 [ 0.4389; 0.6822 ] 
    ##   qual1 ~~ qual5      0.5180      0.0654    7.9173    0.0000 [ 0.3984; 0.6454 ] 
    ##   qual2 ~~ qual3      0.6187      0.0559   11.0663    0.0000 [ 0.4944; 0.7066 ] 
    ##   qual2 ~~ qual4      0.6517      0.0617   10.5605    0.0000 [ 0.5222; 0.7639 ] 
    ##   qual2 ~~ qual5      0.6291      0.0565   11.1275    0.0000 [ 0.5160; 0.7372 ] 
    ##   qual3 ~~ qual4      0.4752      0.0650    7.3052    0.0000 [ 0.3383; 0.5815 ] 
    ##   qual3 ~~ qual5      0.5074      0.0612    8.2965    0.0000 [ 0.3872; 0.6249 ] 
    ##   qual4 ~~ qual5      0.6402      0.0567   11.2851    0.0000 [ 0.5178; 0.7344 ] 
    ##   val1 ~~ val2        0.6344      0.0542   11.7012    0.0000 [ 0.5274; 0.7331 ] 
    ##   val1 ~~ val3        0.4602      0.0738    6.2398    0.0000 [ 0.3113; 0.6053 ] 
    ##   val2 ~~ val3        0.6288      0.0642    9.7880    0.0000 [ 0.4912; 0.7417 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0607    7.7601    0.0000 [ 0.3485; 0.5900 ] 
    ##   QUAL ~ IMAG       0.3933      0.0564    6.9697    0.0000 [ 0.2857; 0.5052 ] 
    ##   QUAL ~ EXPE       0.8344      0.0220   37.9301    0.0000 [ 0.7881; 0.8731 ] 
    ##   VAL ~ IMAG        0.2974      0.0569    5.2290    0.0000 [ 0.1910; 0.4176 ] 
    ##   VAL ~ EXPE        0.6309      0.0493   12.8002    0.0000 [ 0.5300; 0.7209 ] 
    ##   VAL ~ QUAL        0.7013      0.0780    8.9881    0.0000 [ 0.5478; 0.8348 ] 
    ##   SAT ~ IMAG        0.4807      0.0634    7.5787    0.0000 [ 0.3515; 0.6089 ] 
    ##   SAT ~ EXPE        0.5001      0.0576    8.6898    0.0000 [ 0.3954; 0.6061 ] 
    ##   SAT ~ QUAL        0.5911      0.0986    5.9959    0.0000 [ 0.4014; 0.7841 ] 
    ##   SAT ~ VAL         0.5270      0.0894    5.8968    0.0000 [ 0.3473; 0.6909 ] 
    ##   LOY ~ IMAG        0.4840      0.0656    7.3719    0.0000 [ 0.3462; 0.6153 ] 
    ##   LOY ~ EXPE        0.3142      0.0551    5.7018    0.0000 [ 0.2217; 0.4308 ] 
    ##   LOY ~ QUAL        0.3714      0.0849    4.3753    0.0000 [ 0.2327; 0.5663 ] 
    ##   LOY ~ VAL         0.3311      0.0768    4.3091    0.0000 [ 0.1921; 0.4920 ] 
    ##   LOY ~ SAT         0.6283      0.0797    7.8785    0.0000 [ 0.4663; 0.7880 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0564    6.9697    0.0000 [ 0.2857; 0.5052 ] 
    ##   VAL ~ IMAG           0.2974      0.0569    5.2290    0.0000 [ 0.1910; 0.4176 ] 
    ##   VAL ~ EXPE           0.5852      0.0671    8.7242    0.0000 [ 0.4541; 0.7008 ] 
    ##   SAT ~ IMAG           0.2357      0.0454    5.1939    0.0000 [ 0.1603; 0.3279 ] 
    ##   SAT ~ EXPE           0.5173      0.0724    7.1407    0.0000 [ 0.3912; 0.6591 ] 
    ##   SAT ~ QUAL           0.3696      0.0644    5.7423    0.0000 [ 0.2442; 0.4879 ] 
    ##   LOY ~ IMAG           0.3020      0.0514    5.8737    0.0000 [ 0.2111; 0.4118 ] 
    ##   LOY ~ EXPE           0.3142      0.0551    5.7018    0.0000 [ 0.2217; 0.4308 ] 
    ##   LOY ~ QUAL           0.3714      0.0849    4.3753    0.0000 [ 0.2327; 0.5663 ] 
    ##   LOY ~ VAL            0.3311      0.0768    4.3091    0.0000 [ 0.1921; 0.4920 ] 
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
