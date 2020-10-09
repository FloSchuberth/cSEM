
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![R-CMD-check](https://github.com/M-E-Rademaker/cSEM/workflows/R-CMD-check/badge.svg)](https://github.com/M-E-Rademaker/cSEM/actions)
<!-- [![Build Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM) -->
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/cSEM)](https://cran.r-project.org/package=cSEM)
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

<!-- ## Philosophy -->

<!-- - First and foremost: `cSEM` has a user-centered design!. "User-centered" mainly  -->

<!--   boils down to: `cSEM` is easy, i.e. intuitive to use by non-R experts!  -->

<!-- - Modern in a sense that the package integrates modern developments within  -->

<!--   the R community. This mainly includes ideas/recommendations/design choices that -->

<!--   fead into the packages of the [tidyverse](https://github.com/tidyverse/tidyverse). -->

<!-- - State of the art in a sense that we seek to quickly implement recent methodological -->

<!--   developments in composite-based SEM.  -->

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

There are five major postestimation verbs, three test family functions
and four do-family of function:

  - `assess()` : assess the model using common quality criteria
  - `infer()` : calculate common inferencial quantities (e.g., standard
    errors, confidence intervals)
  - `predict()` : predict endogenous indicator values
  - `summarize()` : summarize the results
  - `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

  - `testOMF()` : performs a test for overall model fit
  - `testMICOM()` : performs a test for composite measurement invariance
  - `testMGD()` : performs several tests to assess multi-group
    differences
  - `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently, three do functions are implemented:

  - `doIPMA()`: performs an importance-performance matrix analysis
  - `doNonlinearEffectsAnalysis()`: performs a nonlinear effects
    analysis such as floodlight and surface analysis
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

A useful tool to examine a list is the [listviewer
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
    ##  Inner weighting scheme           = "path"
    ##  Type of indicator correlation    = Pearson
    ##  Path model estimator             = OLS
    ##  Second-order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode      
    ## 
    ##  IMAG  Composite      First order   "modeB"   
    ##  EXPE  Composite      First order   "modeB"   
    ##  QUAL  Composite      First order   "modeB"   
    ##  VAL   Composite      First order   "modeB"   
    ##  SAT   Common factor  First order   "modeA"   
    ##  LOY   Common factor  First order   "modeA"   
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
    ##   Weight           Estimate  Std. error   t-stat.   p-value
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
testOMF(res)
```

    ## ________________________________________________________________________________
    ## --------- Test for overall model fit based on Beran & Srivastava (1985) --------
    ## 
    ## Null hypothesis:
    ## 
    ##                                           +------------------------------------------------------------------+
    ##                                           |                                                                  |
    ##                                           |   H0: The model-implied indicator covariance matrix equals the   |
    ##                                           |   population indicator covariance matrix.                        |
    ##                                           |                                                                  |
    ##                                           +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3320  
    ##  SRMR                    0.0940      0.0536  
    ##  dL                      2.2340      0.7260  
    ##  dML                     2.9219      1.6613  
    ##  
    ## 
    ## Decision: 
    ## 
    ##                          Significance level
    ##  Distance measure          95%   
    ##  dG                      reject  
    ##  SRMR                    reject  
    ##  dL                      reject  
    ##  dML                     reject  
    ##  
    ## Additional information:
    ## 
    ##  Out of 499 bootstrap replications 478 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -6378864
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
    ##  EXPE             NA          0.2222        0.2190    
    ##  QUAL             NA          0.6963        0.6951    
    ##  VAL              NA          0.5474        0.5438    
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
    ##  SAT            0.8960        0.8938        0.9051    
    ##  LOY            0.8237        0.8011        0.8761    
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
    ##  Chi_square     = 727.5611
    ##  Chi_square_df  = 3.954137
    ##  CFI            = 0.8598825
    ##  CN             = 75.14588
    ##  GFI            = 0.7280612
    ##  IFI            = 0.8615598
    ##  NFI            = 0.8229918
    ##  NNFI           = 0.8240917
    ##  RMSEA          = 0.108922
    ##  RMS_theta      = 0.05069299
    ##  SRMR           = 0.09396871
    ## 
    ##  Degrees of freedom    = 184
    ## 
    ## --------------------------- Model selection criteria ---------------------------
    ## 
    ##  Construct        AIC          AICc          AICu     
    ##  EXPE          -59.8152      192.2824      -57.8072   
    ##  QUAL          -294.9343     -42.8367      -292.9263  
    ##  VAL           -193.2127      58.9506      -190.1945  
    ##  SAT           -350.2874     -97.9418      -345.2368  
    ##  LOY           -215.9322      36.2311      -212.9141  
    ## 
    ##  Construct        BIC           FPE           GM      
    ##  EXPE          -52.7723       0.7872       259.8087   
    ##  QUAL          -287.8914      0.3074       271.8568   
    ##  VAL           -182.6483      0.4617       312.7010   
    ##  SAT           -332.6801      0.2463       278.2973   
    ##  LOY           -205.3678      0.4216       291.0665   
    ## 
    ##  Construct        HQ            HQc       Mallows_Cp  
    ##  EXPE          -56.9806      -56.8695       2.7658    
    ##  QUAL          -292.0997     -291.9886      14.8139   
    ##  VAL           -188.9608     -188.7516      52.1366   
    ##  SAT           -343.2010     -342.7088      10.6900   
    ##  LOY           -211.6804     -211.4711      30.5022   
    ## 
    ## ----------------------- Variance inflation factors (VIFs) ----------------------
    ## 
    ##   Dependent construct: 'VAL'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      3.2928   
    ##  QUAL                      3.2928   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      3.2985   
    ##  QUAL                      4.4151   
    ##  IMAG                      1.7280   
    ##  VAL                       2.6726   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct    VIF value 
    ##  IMAG                      1.9345   
    ##  SAT                       1.9345   
    ## 
    ## ------------ Variance inflation factors (VIFs) for modeB constructs ------------
    ## 
    ##   Construct: 'IMAG'
    ## 
    ##  Weight    VIF value 
    ##  imag1      2.8359   
    ##  imag2      5.5535   
    ##  imag3      4.5088   
    ## 
    ##   Construct: 'EXPE'
    ## 
    ##  Weight    VIF value 
    ##  expe1      2.3551   
    ##  expe2      2.7116   
    ##  expe3      2.4116   
    ## 
    ##   Construct: 'QUAL'
    ## 
    ##  Weight    VIF value 
    ##  qual1      3.0835   
    ##  qual2      4.4376   
    ##  qual3      2.9575   
    ##  qual4      3.7341   
    ##  qual5      3.4566   
    ## 
    ##   Construct: 'VAL'
    ## 
    ##  Weight    VIF value 
    ##  val1       2.7725   
    ##  val2       3.8349   
    ##  val3       2.7307   
    ## 
    ## -------------------------- Effect sizes (Cohen's f^2) --------------------------
    ## 
    ##   Dependent construct: 'EXPE'
    ## 
    ##  Independent construct       f^2    
    ##  IMAG                      0.2856   
    ## 
    ##   Dependent construct: 'QUAL'
    ## 
    ##  Independent construct       f^2    
    ##  EXPE                      2.2928   
    ## 
    ##   Dependent construct: 'VAL'
    ## 
    ##  Independent construct       f^2    
    ##  EXPE                      0.0014   
    ##  QUAL                      0.3301   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct       f^2    
    ##  IMAG                      0.1462   
    ##  EXPE                      0.0004   
    ##  QUAL                      0.0468   
    ##  VAL                       0.4373   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct       f^2    
    ##  IMAG                      0.0414   
    ##  SAT                       0.4938   
    ## 
    ## ------------------------------ Validity assessment -----------------------------
    ## 
    ##  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
    ## 
    ##           SAT LOY
    ## SAT 1.0000000   0
    ## LOY 0.7432489   1
    ## 
    ## 
    ##  Fornell-Larcker matrix
    ## 
    ##           SAT       LOY
    ## SAT 0.6851491 0.5696460
    ## LOY 0.5696460 0.5551718
    ## 
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
    ##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   expe1         1.4541         1.5699       1.9046         2.0963       0.0575
    ##   expe2         1.4111         1.4787       1.9311         2.0298       0.2016
    ##   expe3         1.6333         1.7291       2.1274         2.2233       0.1237
    ##   qual1         1.4749         1.5464       1.9269         2.0642       0.1175
    ##   qual2         1.5790         1.5414       2.0373         2.0625       0.2189
    ##   qual3         1.7307         1.7317       2.2215         2.2863       0.1196
    ##   qual4         1.2331         1.1955       1.5940         1.6279       0.2353
    ##   qual5         1.5056         1.5058       1.9353         1.9613       0.1973
    ##   val1          1.4435         1.3614       1.8674         1.7648       0.2503
    ##   val2          1.2230         1.2073       1.6443         1.7136       0.1762
    ##   val3          1.4804         1.3812       1.9675         1.9365       0.1495
    ##   sat1          1.2471         1.2347       1.6455         1.6197       0.3392
    ##   sat2          1.2308         1.1989       1.6404         1.6300       0.3089
    ##   sat3          1.3403         1.2760       1.6725         1.7210       0.2110
    ##   sat4          1.3177         1.2602       1.6673         1.6326       0.2761
    ##   loy1          1.6896         1.6570       2.2298         2.2264       0.2697
    ##   loy2          1.4859         1.4707       1.9120         1.9776       0.1306
    ##   loy3          1.7014         1.6687       2.2785         2.2734       0.2705
    ##   loy4          1.6868         1.6677       2.1750         2.2994       0.0870
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
    ##  Inner weighting scheme           = "path"
    ##  Type of indicator correlation    = Pearson
    ##  Path model estimator             = OLS
    ##  Second-order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample method                  = "bootstrap"
    ##  Number of resamples              = 499
    ##  Number of admissible results     = 481
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = -2137598003
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode      
    ## 
    ##  IMAG  Composite      First order   "modeB"   
    ##  EXPE  Composite      First order   "modeB"   
    ##  QUAL  Composite      First order   "modeB"   
    ##  VAL   Composite      First order   "modeB"   
    ##  SAT   Common factor  First order   "modeA"   
    ##  LOY   Common factor  First order   "modeA"   
    ## 
    ## ----------------------------------- Estimates ----------------------------------
    ## 
    ## Estimated path coefficients:
    ## ============================
    ##                                                              CI_percentile   
    ##   Path           Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG      0.4714      0.0657    7.1795    0.0000 [ 0.3414; 0.5963 ] 
    ##   QUAL ~ EXPE      0.8344      0.0219   38.0309    0.0000 [ 0.7892; 0.8725 ] 
    ##   VAL ~ EXPE       0.0457      0.0854    0.5352    0.5925 [-0.1109; 0.2232 ] 
    ##   VAL ~ QUAL       0.7013      0.0837    8.3837    0.0000 [ 0.5351; 0.8523 ] 
    ##   SAT ~ IMAG       0.2450      0.0531    4.6125    0.0000 [ 0.1426; 0.3407 ] 
    ##   SAT ~ EXPE      -0.0172      0.0756   -0.2278    0.8198 [-0.1492; 0.1310 ] 
    ##   SAT ~ QUAL       0.2215      0.1074    2.0623    0.0392 [ 0.0254; 0.4443 ] 
    ##   SAT ~ VAL        0.5270      0.0860    6.1252    0.0000 [ 0.3518; 0.6828 ] 
    ##   LOY ~ IMAG       0.1819      0.0800    2.2733    0.0230 [ 0.0282; 0.3496 ] 
    ##   LOY ~ SAT        0.6283      0.0820    7.6647    0.0000 [ 0.4710; 0.7971 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0988    6.3813    0.0000 [ 0.4362; 0.8219 ] 
    ##   IMAG =~ imag2      0.9246      0.0393   23.4984    0.0000 [ 0.8215; 0.9739 ] 
    ##   IMAG =~ imag3      0.9577      0.0294   32.5363    0.0000 [ 0.8822; 0.9927 ] 
    ##   EXPE =~ expe1      0.7525      0.0734   10.2506    0.0000 [ 0.6055; 0.8759 ] 
    ##   EXPE =~ expe2      0.9348      0.0285   32.7604    0.0000 [ 0.8605; 0.9723 ] 
    ##   EXPE =~ expe3      0.7295      0.0712   10.2463    0.0000 [ 0.5746; 0.8432 ] 
    ##   QUAL =~ qual1      0.7861      0.0630   12.4741    0.0000 [ 0.6455; 0.8897 ] 
    ##   QUAL =~ qual2      0.9244      0.0236   39.1187    0.0000 [ 0.8682; 0.9579 ] 
    ##   QUAL =~ qual3      0.7560      0.0612   12.3534    0.0000 [ 0.6207; 0.8555 ] 
    ##   QUAL =~ qual4      0.7632      0.0551   13.8438    0.0000 [ 0.6383; 0.8569 ] 
    ##   QUAL =~ qual5      0.7834      0.0432   18.1246    0.0000 [ 0.6904; 0.8553 ] 
    ##   VAL =~ val1        0.9518      0.0233   40.9183    0.0000 [ 0.8968; 0.9822 ] 
    ##   VAL =~ val2        0.8056      0.0602   13.3868    0.0000 [ 0.6718; 0.9080 ] 
    ##   VAL =~ val3        0.6763      0.0703    9.6147    0.0000 [ 0.5250; 0.7973 ] 
    ##   SAT =~ sat1        0.9243      0.0235   39.3340    0.0000 [ 0.8772; 0.9630 ] 
    ##   SAT =~ sat2        0.8813      0.0279   31.5950    0.0000 [ 0.8231; 0.9312 ] 
    ##   SAT =~ sat3        0.7127      0.0511   13.9539    0.0000 [ 0.6180; 0.8169 ] 
    ##   SAT =~ sat4        0.7756      0.0504   15.3912    0.0000 [ 0.6674; 0.8664 ] 
    ##   LOY =~ loy1        0.9097      0.0510   17.8265    0.0000 [ 0.7813; 0.9850 ] 
    ##   LOY =~ loy2        0.5775      0.0856    6.7434    0.0000 [ 0.4092; 0.7319 ] 
    ##   LOY =~ loy3        0.9043      0.0417   21.7007    0.0000 [ 0.8102; 0.9764 ] 
    ##   LOY =~ loy4        0.4917      0.0995    4.9407    0.0000 [ 0.2885; 0.6582 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1246    0.1256    0.9001 [-0.1990; 0.2763 ] 
    ##   IMAG <~ imag2      0.4473      0.1508    2.9666    0.0030 [ 0.1293; 0.7219 ] 
    ##   IMAG <~ imag3      0.6020      0.1400    4.2995    0.0000 [ 0.3106; 0.8665 ] 
    ##   EXPE <~ expe1      0.2946      0.1154    2.5529    0.0107 [ 0.0801; 0.5236 ] 
    ##   EXPE <~ expe2      0.6473      0.0874    7.4034    0.0000 [ 0.4540; 0.7837 ] 
    ##   EXPE <~ expe3      0.2374      0.0948    2.5046    0.0123 [ 0.0505; 0.4159 ] 
    ##   QUAL <~ qual1      0.2370      0.0910    2.6042    0.0092 [ 0.0790; 0.4271 ] 
    ##   QUAL <~ qual2      0.4712      0.0814    5.7864    0.0000 [ 0.3004; 0.6094 ] 
    ##   QUAL <~ qual3      0.1831      0.0840    2.1797    0.0293 [ 0.0087; 0.3433 ] 
    ##   QUAL <~ qual4      0.1037      0.0623    1.6656    0.0958 [-0.0145; 0.2259 ] 
    ##   QUAL <~ qual5      0.2049      0.0609    3.3655    0.0008 [ 0.0739; 0.3165 ] 
    ##   VAL <~ val1        0.7163      0.0924    7.7497    0.0000 [ 0.5042; 0.8647 ] 
    ##   VAL <~ val2        0.2202      0.0879    2.5061    0.0122 [ 0.0705; 0.3988 ] 
    ##   VAL <~ val3        0.2082      0.0610    3.4137    0.0006 [ 0.0734; 0.3142 ] 
    ##   SAT <~ sat1        0.3209      0.0151   21.2601    0.0000 [ 0.2923; 0.3500 ] 
    ##   SAT <~ sat2        0.3059      0.0134   22.9075    0.0000 [ 0.2836; 0.3379 ] 
    ##   SAT <~ sat3        0.2474      0.0106   23.3086    0.0000 [ 0.2266; 0.2672 ] 
    ##   SAT <~ sat4        0.2692      0.0123   21.8760    0.0000 [ 0.2472; 0.2939 ] 
    ##   LOY <~ loy1        0.3834      0.0274   13.9728    0.0000 [ 0.3275; 0.4350 ] 
    ##   LOY <~ loy2        0.2434      0.0290    8.4009    0.0000 [ 0.1838; 0.2955 ] 
    ##   LOY <~ loy3        0.3812      0.0272   14.0070    0.0000 [ 0.3331; 0.4314 ] 
    ##   LOY <~ loy4        0.2073      0.0375    5.5243    0.0000 [ 0.1332; 0.2717 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0623   10.3261    0.0000 [ 0.5168; 0.7605 ] 
    ##   imag1 ~~ imag3      0.5433      0.0682    7.9616    0.0000 [ 0.4112; 0.6747 ] 
    ##   imag2 ~~ imag3      0.7761      0.0397   19.5411    0.0000 [ 0.6920; 0.8426 ] 
    ##   expe1 ~~ expe2      0.5353      0.0564    9.4873    0.0000 [ 0.4279; 0.6493 ] 
    ##   expe1 ~~ expe3      0.4694      0.0562    8.3464    0.0000 [ 0.3581; 0.5759 ] 
    ##   expe2 ~~ expe3      0.5467      0.0582    9.3954    0.0000 [ 0.4305; 0.6539 ] 
    ##   qual1 ~~ qual2      0.6053      0.0530   11.4181    0.0000 [ 0.5087; 0.7034 ] 
    ##   qual1 ~~ qual3      0.5406      0.0601    8.9958    0.0000 [ 0.4212; 0.6507 ] 
    ##   qual1 ~~ qual4      0.5662      0.0672    8.4311    0.0000 [ 0.4386; 0.6888 ] 
    ##   qual1 ~~ qual5      0.5180      0.0650    7.9693    0.0000 [ 0.3937; 0.6456 ] 
    ##   qual2 ~~ qual3      0.6187      0.0551   11.2205    0.0000 [ 0.5130; 0.7162 ] 
    ##   qual2 ~~ qual4      0.6517      0.0612   10.6442    0.0000 [ 0.5299; 0.7575 ] 
    ##   qual2 ~~ qual5      0.6291      0.0561   11.2219    0.0000 [ 0.5163; 0.7307 ] 
    ##   qual3 ~~ qual4      0.4752      0.0630    7.5466    0.0000 [ 0.3467; 0.5926 ] 
    ##   qual3 ~~ qual5      0.5074      0.0582    8.7201    0.0000 [ 0.3931; 0.6129 ] 
    ##   qual4 ~~ qual5      0.6402      0.0546   11.7173    0.0000 [ 0.5197; 0.7346 ] 
    ##   val1 ~~ val2        0.6344      0.0551   11.5234    0.0000 [ 0.5265; 0.7311 ] 
    ##   val1 ~~ val3        0.4602      0.0696    6.6095    0.0000 [ 0.3178; 0.5869 ] 
    ##   val2 ~~ val3        0.6288      0.0629    9.9917    0.0000 [ 0.5023; 0.7478 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0657    7.1795    0.0000 [ 0.3414; 0.5963 ] 
    ##   QUAL ~ IMAG       0.3933      0.0605    6.5012    0.0000 [ 0.2781; 0.5139 ] 
    ##   QUAL ~ EXPE       0.8344      0.0219   38.0309    0.0000 [ 0.7892; 0.8725 ] 
    ##   VAL ~ IMAG        0.2974      0.0605    4.9189    0.0000 [ 0.1907; 0.4226 ] 
    ##   VAL ~ EXPE        0.6309      0.0514   12.2840    0.0000 [ 0.5281; 0.7283 ] 
    ##   VAL ~ QUAL        0.7013      0.0837    8.3837    0.0000 [ 0.5351; 0.8523 ] 
    ##   SAT ~ IMAG        0.4807      0.0675    7.1193    0.0000 [ 0.3442; 0.6118 ] 
    ##   SAT ~ EXPE        0.5001      0.0572    8.7385    0.0000 [ 0.3813; 0.6086 ] 
    ##   SAT ~ QUAL        0.5911      0.1007    5.8691    0.0000 [ 0.3900; 0.7777 ] 
    ##   SAT ~ VAL         0.5270      0.0860    6.1252    0.0000 [ 0.3518; 0.6828 ] 
    ##   LOY ~ IMAG        0.4840      0.0677    7.1506    0.0000 [ 0.3625; 0.6170 ] 
    ##   LOY ~ EXPE        0.3142      0.0550    5.7121    0.0000 [ 0.2178; 0.4426 ] 
    ##   LOY ~ QUAL        0.3714      0.0891    4.1684    0.0000 [ 0.2240; 0.5819 ] 
    ##   LOY ~ VAL         0.3311      0.0737    4.4927    0.0000 [ 0.1984; 0.4862 ] 
    ##   LOY ~ SAT         0.6283      0.0820    7.6647    0.0000 [ 0.4710; 0.7971 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0605    6.5012    0.0000 [ 0.2781; 0.5139 ] 
    ##   VAL ~ IMAG           0.2974      0.0605    4.9189    0.0000 [ 0.1907; 0.4226 ] 
    ##   VAL ~ EXPE           0.5852      0.0717    8.1593    0.0000 [ 0.4358; 0.7077 ] 
    ##   SAT ~ IMAG           0.2357      0.0492    4.7932    0.0000 [ 0.1449; 0.3359 ] 
    ##   SAT ~ EXPE           0.5173      0.0704    7.3469    0.0000 [ 0.3911; 0.6698 ] 
    ##   SAT ~ QUAL           0.3696      0.0612    6.0395    0.0000 [ 0.2452; 0.4883 ] 
    ##   LOY ~ IMAG           0.3020      0.0559    5.4058    0.0000 [ 0.2035; 0.4243 ] 
    ##   LOY ~ EXPE           0.3142      0.0550    5.7121    0.0000 [ 0.2178; 0.4426 ] 
    ##   LOY ~ QUAL           0.3714      0.0891    4.1684    0.0000 [ 0.2240; 0.5819 ] 
    ##   LOY ~ VAL            0.3311      0.0737    4.4927    0.0000 [ 0.1984; 0.4862 ] 
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
