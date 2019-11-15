
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
# ## Test overall model fit
testOMF(res, .verbose = FALSE)
```

    ## ________________________________________________________________________________
    ## --------- Test for overall model fit based on Beran & Srivastava (1985) --------
    ## 
    ## Null hypothesis:
    ## 
    ##                   +------------------------------------------------------------+
    ##                   |                                                            |
    ##                   |   H0: Population indicator covariance matrix is equal to   |
    ##                   |   model-implied indicator covariance matrix.               |
    ##                   |                                                            |
    ##                   +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3059  
    ##  SRMR                    0.0940      0.0517  
    ##  dL                      2.2340      0.6768  
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
    ##  The seed used was: 620575944
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
    ##  Number of admissible results     = 487
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -1629688704
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
    ##   EXPE ~ IMAG      0.4714      0.0637    7.4019    0.0000 [ 0.3562; 0.5986 ] 
    ##   QUAL ~ EXPE      0.8344      0.0223   37.4962    0.0000 [ 0.7872; 0.8750 ] 
    ##   VAL ~ EXPE       0.0457      0.0863    0.5296    0.5964 [-0.1157; 0.2408 ] 
    ##   VAL ~ QUAL       0.7013      0.0801    8.7547    0.0000 [ 0.5263; 0.8444 ] 
    ##   SAT ~ IMAG       0.2450      0.0557    4.3958    0.0000 [ 0.1404; 0.3557 ] 
    ##   SAT ~ EXPE      -0.0172      0.0734   -0.2347    0.8145 [-0.1642; 0.1264 ] 
    ##   SAT ~ QUAL       0.2215      0.1019    2.1736    0.0297 [ 0.0512; 0.4421 ] 
    ##   SAT ~ VAL        0.5270      0.0896    5.8796    0.0000 [ 0.3370; 0.7036 ] 
    ##   LOY ~ IMAG       0.1819      0.0830    2.1921    0.0284 [ 0.0191; 0.3476 ] 
    ##   LOY ~ SAT        0.6283      0.0849    7.3967    0.0000 [ 0.4643; 0.8014 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0986    6.3966    0.0000 [ 0.4159; 0.8073 ] 
    ##   IMAG =~ imag2      0.9246      0.0407   22.7206    0.0000 [ 0.8302; 0.9788 ] 
    ##   IMAG =~ imag3      0.9577      0.0274   34.9315    0.0000 [ 0.8824; 0.9905 ] 
    ##   EXPE =~ expe1      0.7525      0.0788    9.5456    0.0000 [ 0.5668; 0.8684 ] 
    ##   EXPE =~ expe2      0.9348      0.0284   32.9174    0.0000 [ 0.8641; 0.9699 ] 
    ##   EXPE =~ expe3      0.7295      0.0728   10.0184    0.0000 [ 0.5580; 0.8404 ] 
    ##   QUAL =~ qual1      0.7861      0.0674   11.6623    0.0000 [ 0.6287; 0.8872 ] 
    ##   QUAL =~ qual2      0.9244      0.0223   41.3774    0.0000 [ 0.8669; 0.9558 ] 
    ##   QUAL =~ qual3      0.7560      0.0612   12.3527    0.0000 [ 0.6279; 0.8562 ] 
    ##   QUAL =~ qual4      0.7632      0.0525   14.5383    0.0000 [ 0.6497; 0.8504 ] 
    ##   QUAL =~ qual5      0.7834      0.0426   18.3743    0.0000 [ 0.6888; 0.8535 ] 
    ##   VAL =~ val1        0.9518      0.0224   42.5708    0.0000 [ 0.9030; 0.9837 ] 
    ##   VAL =~ val2        0.8056      0.0613   13.1466    0.0000 [ 0.6733; 0.9096 ] 
    ##   VAL =~ val3        0.6763      0.0731    9.2498    0.0000 [ 0.5282; 0.8039 ] 
    ##   SAT =~ sat1        0.9243      0.0232   39.8801    0.0000 [ 0.8718; 0.9619 ] 
    ##   SAT =~ sat2        0.8813      0.0283   31.1650    0.0000 [ 0.8211; 0.9285 ] 
    ##   SAT =~ sat3        0.7127      0.0532   13.4025    0.0000 [ 0.5983; 0.8116 ] 
    ##   SAT =~ sat4        0.7756      0.0497   15.6082    0.0000 [ 0.6765; 0.8665 ] 
    ##   LOY =~ loy1        0.9097      0.0493   18.4378    0.0000 [ 0.7928; 0.9829 ] 
    ##   LOY =~ loy2        0.5775      0.0831    6.9519    0.0000 [ 0.3988; 0.7123 ] 
    ##   LOY =~ loy3        0.9043      0.0436   20.7389    0.0000 [ 0.8067; 0.9769 ] 
    ##   LOY =~ loy4        0.4917      0.0954    5.1556    0.0000 [ 0.3178; 0.6932 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1168    0.1340    0.8934 [-0.2015; 0.2293 ] 
    ##   IMAG <~ imag2      0.4473      0.1525    2.9340    0.0033 [ 0.1345; 0.7365 ] 
    ##   IMAG <~ imag3      0.6020      0.1378    4.3680    0.0000 [ 0.3215; 0.8457 ] 
    ##   EXPE <~ expe1      0.2946      0.1127    2.6139    0.0090 [ 0.0772; 0.5027 ] 
    ##   EXPE <~ expe2      0.6473      0.0831    7.7901    0.0000 [ 0.4559; 0.7925 ] 
    ##   EXPE <~ expe3      0.2374      0.0896    2.6487    0.0081 [ 0.0610; 0.4056 ] 
    ##   QUAL <~ qual1      0.2370      0.0869    2.7271    0.0064 [ 0.0818; 0.4261 ] 
    ##   QUAL <~ qual2      0.4712      0.0778    6.0592    0.0000 [ 0.3151; 0.6318 ] 
    ##   QUAL <~ qual3      0.1831      0.0764    2.3967    0.0165 [ 0.0309; 0.3289 ] 
    ##   QUAL <~ qual4      0.1037      0.0634    1.6357    0.1019 [-0.0232; 0.2300 ] 
    ##   QUAL <~ qual5      0.2049      0.0575    3.5636    0.0004 [ 0.0806; 0.3062 ] 
    ##   VAL <~ val1        0.7163      0.0914    7.8330    0.0000 [ 0.5190; 0.8760 ] 
    ##   VAL <~ val2        0.2202      0.0907    2.4278    0.0152 [ 0.0685; 0.4216 ] 
    ##   VAL <~ val3        0.2082      0.0607    3.4309    0.0006 [ 0.0834; 0.3218 ] 
    ##   SAT <~ sat1        0.3209      0.0152   21.1336    0.0000 [ 0.2928; 0.3544 ] 
    ##   SAT <~ sat2        0.3059      0.0142   21.5025    0.0000 [ 0.2804; 0.3353 ] 
    ##   SAT <~ sat3        0.2474      0.0112   22.1718    0.0000 [ 0.2256; 0.2686 ] 
    ##   SAT <~ sat4        0.2692      0.0115   23.3400    0.0000 [ 0.2488; 0.2939 ] 
    ##   LOY <~ loy1        0.3834      0.0260   14.7629    0.0000 [ 0.3314; 0.4312 ] 
    ##   LOY <~ loy2        0.2434      0.0293    8.3208    0.0000 [ 0.1809; 0.2916 ] 
    ##   LOY <~ loy3        0.3812      0.0269   14.1888    0.0000 [ 0.3306; 0.4356 ] 
    ##   LOY <~ loy4        0.2073      0.0352    5.8850    0.0000 [ 0.1427; 0.2797 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0618   10.4074    0.0000 [ 0.5202; 0.7504 ] 
    ##   imag1 ~~ imag3      0.5433      0.0706    7.6960    0.0000 [ 0.4033; 0.6781 ] 
    ##   imag2 ~~ imag3      0.7761      0.0390   19.8958    0.0000 [ 0.6933; 0.8423 ] 
    ##   expe1 ~~ expe2      0.5353      0.0623    8.5929    0.0000 [ 0.3997; 0.6468 ] 
    ##   expe1 ~~ expe3      0.4694      0.0626    7.5036    0.0000 [ 0.3414; 0.5853 ] 
    ##   expe2 ~~ expe3      0.5467      0.0610    8.9545    0.0000 [ 0.4134; 0.6588 ] 
    ##   qual1 ~~ qual2      0.6053      0.0575   10.5269    0.0000 [ 0.4880; 0.7069 ] 
    ##   qual1 ~~ qual3      0.5406      0.0593    9.1248    0.0000 [ 0.4341; 0.6624 ] 
    ##   qual1 ~~ qual4      0.5662      0.0676    8.3724    0.0000 [ 0.4290; 0.6953 ] 
    ##   qual1 ~~ qual5      0.5180      0.0693    7.4799    0.0000 [ 0.3722; 0.6419 ] 
    ##   qual2 ~~ qual3      0.6187      0.0550   11.2522    0.0000 [ 0.4982; 0.7153 ] 
    ##   qual2 ~~ qual4      0.6517      0.0589   11.0596    0.0000 [ 0.5338; 0.7676 ] 
    ##   qual2 ~~ qual5      0.6291      0.0553   11.3849    0.0000 [ 0.5059; 0.7207 ] 
    ##   qual3 ~~ qual4      0.4752      0.0644    7.3735    0.0000 [ 0.3400; 0.5832 ] 
    ##   qual3 ~~ qual5      0.5074      0.0611    8.3050    0.0000 [ 0.3859; 0.6230 ] 
    ##   qual4 ~~ qual5      0.6402      0.0563   11.3798    0.0000 [ 0.5132; 0.7391 ] 
    ##   val1 ~~ val2        0.6344      0.0536   11.8261    0.0000 [ 0.5191; 0.7282 ] 
    ##   val1 ~~ val3        0.4602      0.0678    6.7858    0.0000 [ 0.3372; 0.6004 ] 
    ##   val2 ~~ val3        0.6288      0.0653    9.6370    0.0000 [ 0.4957; 0.7486 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0637    7.4019    0.0000 [ 0.3562; 0.5986 ] 
    ##   QUAL ~ IMAG       0.3933      0.0589    6.6801    0.0000 [ 0.2881; 0.5111 ] 
    ##   QUAL ~ EXPE       0.8344      0.0223   37.4962    0.0000 [ 0.7872; 0.8750 ] 
    ##   VAL ~ IMAG        0.2974      0.0592    5.0252    0.0000 [ 0.1949; 0.4223 ] 
    ##   VAL ~ EXPE        0.6309      0.0509   12.3916    0.0000 [ 0.5215; 0.7250 ] 
    ##   VAL ~ QUAL        0.7013      0.0801    8.7547    0.0000 [ 0.5263; 0.8444 ] 
    ##   SAT ~ IMAG        0.4807      0.0654    7.3516    0.0000 [ 0.3448; 0.6080 ] 
    ##   SAT ~ EXPE        0.5001      0.0595    8.4117    0.0000 [ 0.3953; 0.6133 ] 
    ##   SAT ~ QUAL        0.5911      0.0926    6.3849    0.0000 [ 0.4095; 0.7707 ] 
    ##   SAT ~ VAL         0.5270      0.0896    5.8796    0.0000 [ 0.3370; 0.7036 ] 
    ##   LOY ~ IMAG        0.4840      0.0682    7.0911    0.0000 [ 0.3545; 0.6200 ] 
    ##   LOY ~ EXPE        0.3142      0.0587    5.3547    0.0000 [ 0.1997; 0.4298 ] 
    ##   LOY ~ QUAL        0.3714      0.0849    4.3766    0.0000 [ 0.2144; 0.5478 ] 
    ##   LOY ~ VAL         0.3311      0.0762    4.3422    0.0000 [ 0.1882; 0.4747 ] 
    ##   LOY ~ SAT         0.6283      0.0849    7.3967    0.0000 [ 0.4643; 0.8014 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0589    6.6801    0.0000 [ 0.2881; 0.5111 ] 
    ##   VAL ~ IMAG           0.2974      0.0592    5.0252    0.0000 [ 0.1949; 0.4223 ] 
    ##   VAL ~ EXPE           0.5852      0.0687    8.5244    0.0000 [ 0.4438; 0.7060 ] 
    ##   SAT ~ IMAG           0.2357      0.0492    4.7905    0.0000 [ 0.1565; 0.3482 ] 
    ##   SAT ~ EXPE           0.5173      0.0687    7.5262    0.0000 [ 0.3962; 0.6565 ] 
    ##   SAT ~ QUAL           0.3696      0.0643    5.7489    0.0000 [ 0.2397; 0.4821 ] 
    ##   LOY ~ IMAG           0.3020      0.0547    5.5178    0.0000 [ 0.2111; 0.4147 ] 
    ##   LOY ~ EXPE           0.3142      0.0587    5.3547    0.0000 [ 0.1997; 0.4298 ] 
    ##   LOY ~ QUAL           0.3714      0.0849    4.3766    0.0000 [ 0.2144; 0.5478 ] 
    ##   LOY ~ VAL            0.3311      0.0762    4.3422    0.0000 [ 0.1882; 0.4747 ] 
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
