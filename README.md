
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![R-CMD-check](https://github.com/M-E-Rademaker/cSEM/workflows/R-CMD-check/badge.svg)](https://github.com/M-E-Rademaker/cSEM/actions)
<!-- [![Build Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM) -->
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)
![Lifecycle
Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
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
    ##                            +------------------------------------------------------------------+
    ##                            |                                                                  |
    ##                            |   H0: The model-implied indicator covariance matrix equals the   |
    ##                            |   population indicator covariance matrix.                        |
    ##                            |                                                                  |
    ##                            +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3333  
    ##  SRMR                    0.0940      0.0523  
    ##  dL                      2.2340      0.6916  
    ##  dML                     2.9219      1.6389  
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
    ##  Out of 499 bootstrap replications 476 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 293914126
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
    ##   expe1         1.4544         1.5666       1.9066         2.0944       0.0554
    ##   expe2         1.4116         1.4816       1.9359         2.0300       0.2001
    ##   expe3         1.6323         1.7288       2.1271         2.2239       0.1229
    ##   qual1         1.4770         1.5458       1.9294         2.0640       0.1165
    ##   qual2         1.5807         1.5379       2.0415         2.0601       0.2169
    ##   qual3         1.7323         1.7315       2.2246         2.2878       0.1181
    ##   qual4         1.2356         1.1972       1.5997         1.6338       0.2326
    ##   qual5         1.5061         1.5028       1.9359         1.9587       0.1979
    ##   val1          1.4489         1.3655       1.8727         1.7684       0.2482
    ##   val2          1.2263         1.2062       1.6489         1.7131       0.1738
    ##   val3          1.4806         1.3786       1.9680         1.9302       0.1484
    ##   sat1          1.2459         1.2365       1.6449         1.6224       0.3401
    ##   sat2          1.2316         1.1968       1.6400         1.6293       0.3091
    ##   sat3          1.3404         1.2792       1.6734         1.7280       0.2102
    ##   sat4          1.3185         1.2636       1.6687         1.6362       0.2770
    ##   loy1          1.6938         1.6595       2.2342         2.2257       0.2682
    ##   loy2          1.4863         1.4727       1.9126         1.9789       0.1322
    ##   loy3          1.7070         1.6723       2.2845         2.2748       0.2683
    ##   loy4          1.6880         1.6690       2.1758         2.2990       0.0877
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
    ##  Number of admissible results     = 482
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 1856154751
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
    ##   EXPE ~ IMAG      0.4714      0.0670    7.0314    0.0000 [ 0.3351; 0.5959 ] 
    ##   QUAL ~ EXPE      0.8344      0.0243   34.4097    0.0000 [ 0.7764; 0.8743 ] 
    ##   VAL ~ EXPE       0.0457      0.0886    0.5158    0.6060 [-0.1102; 0.2362 ] 
    ##   VAL ~ QUAL       0.7013      0.0852    8.2288    0.0000 [ 0.5201; 0.8439 ] 
    ##   SAT ~ IMAG       0.2450      0.0568    4.3154    0.0000 [ 0.1336; 0.3481 ] 
    ##   SAT ~ EXPE      -0.0172      0.0697   -0.2471    0.8048 [-0.1702; 0.1093 ] 
    ##   SAT ~ QUAL       0.2215      0.1012    2.1888    0.0286 [ 0.0539; 0.4415 ] 
    ##   SAT ~ VAL        0.5270      0.0906    5.8149    0.0000 [ 0.3457; 0.7039 ] 
    ##   LOY ~ IMAG       0.1819      0.0809    2.2491    0.0245 [ 0.0178; 0.3231 ] 
    ##   LOY ~ SAT        0.6283      0.0844    7.4423    0.0000 [ 0.4788; 0.7952 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1010    6.2461    0.0000 [ 0.4135; 0.8142 ] 
    ##   IMAG =~ imag2      0.9246      0.0405   22.8443    0.0000 [ 0.8222; 0.9778 ] 
    ##   IMAG =~ imag3      0.9577      0.0295   32.5018    0.0000 [ 0.8785; 0.9906 ] 
    ##   EXPE =~ expe1      0.7525      0.0808    9.3148    0.0000 [ 0.5783; 0.8801 ] 
    ##   EXPE =~ expe2      0.9348      0.0295   31.6560    0.0000 [ 0.8604; 0.9719 ] 
    ##   EXPE =~ expe3      0.7295      0.0742    9.8295    0.0000 [ 0.5525; 0.8487 ] 
    ##   QUAL =~ qual1      0.7861      0.0680   11.5528    0.0000 [ 0.6340; 0.8953 ] 
    ##   QUAL =~ qual2      0.9244      0.0226   40.9651    0.0000 [ 0.8628; 0.9538 ] 
    ##   QUAL =~ qual3      0.7560      0.0629   12.0243    0.0000 [ 0.6034; 0.8579 ] 
    ##   QUAL =~ qual4      0.7632      0.0541   14.1117    0.0000 [ 0.6415; 0.8528 ] 
    ##   QUAL =~ qual5      0.7834      0.0455   17.2028    0.0000 [ 0.6865; 0.8660 ] 
    ##   VAL =~ val1        0.9518      0.0248   38.4300    0.0000 [ 0.8898; 0.9839 ] 
    ##   VAL =~ val2        0.8056      0.0639   12.6121    0.0000 [ 0.6732; 0.9131 ] 
    ##   VAL =~ val3        0.6763      0.0697    9.7006    0.0000 [ 0.5308; 0.8084 ] 
    ##   SAT =~ sat1        0.9243      0.0213   43.4006    0.0000 [ 0.8763; 0.9572 ] 
    ##   SAT =~ sat2        0.8813      0.0300   29.3369    0.0000 [ 0.8110; 0.9276 ] 
    ##   SAT =~ sat3        0.7127      0.0529   13.4662    0.0000 [ 0.6038; 0.8095 ] 
    ##   SAT =~ sat4        0.7756      0.0492   15.7746    0.0000 [ 0.6643; 0.8633 ] 
    ##   LOY =~ loy1        0.9097      0.0494   18.4175    0.0000 [ 0.7918; 0.9855 ] 
    ##   LOY =~ loy2        0.5775      0.0849    6.8061    0.0000 [ 0.4134; 0.7257 ] 
    ##   LOY =~ loy3        0.9043      0.0431   21.0034    0.0000 [ 0.8158; 0.9743 ] 
    ##   LOY =~ loy4        0.4917      0.0983    5.0026    0.0000 [ 0.3078; 0.6902 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1160    0.1349    0.8927 [-0.2179; 0.2412 ] 
    ##   IMAG <~ imag2      0.4473      0.1473    3.0362    0.0024 [ 0.1664; 0.7498 ] 
    ##   IMAG <~ imag3      0.6020      0.1380    4.3634    0.0000 [ 0.3135; 0.8519 ] 
    ##   EXPE <~ expe1      0.2946      0.1179    2.4979    0.0125 [ 0.0733; 0.5285 ] 
    ##   EXPE <~ expe2      0.6473      0.0882    7.3354    0.0000 [ 0.4595; 0.7917 ] 
    ##   EXPE <~ expe3      0.2374      0.0975    2.4340    0.0149 [ 0.0356; 0.4479 ] 
    ##   QUAL <~ qual1      0.2370      0.0911    2.6023    0.0093 [ 0.0935; 0.4393 ] 
    ##   QUAL <~ qual2      0.4712      0.0788    5.9795    0.0000 [ 0.2955; 0.6010 ] 
    ##   QUAL <~ qual3      0.1831      0.0811    2.2563    0.0241 [ 0.0164; 0.3320 ] 
    ##   QUAL <~ qual4      0.1037      0.0605    1.7134    0.0866 [-0.0034; 0.2267 ] 
    ##   QUAL <~ qual5      0.2049      0.0604    3.3914    0.0007 [ 0.0880; 0.3176 ] 
    ##   VAL <~ val1        0.7163      0.0988    7.2496    0.0000 [ 0.4865; 0.8632 ] 
    ##   VAL <~ val2        0.2202      0.0968    2.2738    0.0230 [ 0.0541; 0.4248 ] 
    ##   VAL <~ val3        0.2082      0.0601    3.4645    0.0005 [ 0.0846; 0.3302 ] 
    ##   SAT <~ sat1        0.3209      0.0150   21.3691    0.0000 [ 0.2958; 0.3527 ] 
    ##   SAT <~ sat2        0.3059      0.0135   22.5902    0.0000 [ 0.2842; 0.3345 ] 
    ##   SAT <~ sat3        0.2474      0.0110   22.5776    0.0000 [ 0.2259; 0.2691 ] 
    ##   SAT <~ sat4        0.2692      0.0124   21.7550    0.0000 [ 0.2445; 0.2947 ] 
    ##   LOY <~ loy1        0.3834      0.0269   14.2356    0.0000 [ 0.3292; 0.4343 ] 
    ##   LOY <~ loy2        0.2434      0.0293    8.2965    0.0000 [ 0.1812; 0.2920 ] 
    ##   LOY <~ loy3        0.3812      0.0277   13.7400    0.0000 [ 0.3261; 0.4357 ] 
    ##   LOY <~ loy4        0.2073      0.0359    5.7797    0.0000 [ 0.1398; 0.2803 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0653    9.8643    0.0000 [ 0.5145; 0.7567 ] 
    ##   imag1 ~~ imag3      0.5433      0.0707    7.6831    0.0000 [ 0.4074; 0.6772 ] 
    ##   imag2 ~~ imag3      0.7761      0.0379   20.4816    0.0000 [ 0.6885; 0.8356 ] 
    ##   expe1 ~~ expe2      0.5353      0.0616    8.6873    0.0000 [ 0.4041; 0.6400 ] 
    ##   expe1 ~~ expe3      0.4694      0.0604    7.7775    0.0000 [ 0.3450; 0.5799 ] 
    ##   expe2 ~~ expe3      0.5467      0.0628    8.7102    0.0000 [ 0.4065; 0.6438 ] 
    ##   qual1 ~~ qual2      0.6053      0.0570   10.6112    0.0000 [ 0.4864; 0.7059 ] 
    ##   qual1 ~~ qual3      0.5406      0.0606    8.9276    0.0000 [ 0.4152; 0.6542 ] 
    ##   qual1 ~~ qual4      0.5662      0.0707    8.0130    0.0000 [ 0.4297; 0.6935 ] 
    ##   qual1 ~~ qual5      0.5180      0.0704    7.3629    0.0000 [ 0.3813; 0.6521 ] 
    ##   qual2 ~~ qual3      0.6187      0.0584   10.5928    0.0000 [ 0.4942; 0.7112 ] 
    ##   qual2 ~~ qual4      0.6517      0.0606   10.7619    0.0000 [ 0.5268; 0.7630 ] 
    ##   qual2 ~~ qual5      0.6291      0.0582   10.8008    0.0000 [ 0.5078; 0.7319 ] 
    ##   qual3 ~~ qual4      0.4752      0.0631    7.5325    0.0000 [ 0.3385; 0.5849 ] 
    ##   qual3 ~~ qual5      0.5074      0.0598    8.4908    0.0000 [ 0.3917; 0.6108 ] 
    ##   qual4 ~~ qual5      0.6402      0.0557   11.4928    0.0000 [ 0.5158; 0.7298 ] 
    ##   val1 ~~ val2        0.6344      0.0542   11.6980    0.0000 [ 0.5321; 0.7430 ] 
    ##   val1 ~~ val3        0.4602      0.0679    6.7823    0.0000 [ 0.3250; 0.5876 ] 
    ##   val2 ~~ val3        0.6288      0.0641    9.8137    0.0000 [ 0.5034; 0.7435 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0670    7.0314    0.0000 [ 0.3351; 0.5959 ] 
    ##   QUAL ~ IMAG       0.3933      0.0618    6.3594    0.0000 [ 0.2671; 0.5132 ] 
    ##   QUAL ~ EXPE       0.8344      0.0243   34.4097    0.0000 [ 0.7764; 0.8743 ] 
    ##   VAL ~ IMAG        0.2974      0.0607    4.9018    0.0000 [ 0.1786; 0.4188 ] 
    ##   VAL ~ EXPE        0.6309      0.0495   12.7392    0.0000 [ 0.5296; 0.7151 ] 
    ##   VAL ~ QUAL        0.7013      0.0852    8.2288    0.0000 [ 0.5201; 0.8439 ] 
    ##   SAT ~ IMAG        0.4807      0.0688    6.9854    0.0000 [ 0.3446; 0.5986 ] 
    ##   SAT ~ EXPE        0.5001      0.0552    9.0559    0.0000 [ 0.3902; 0.6042 ] 
    ##   SAT ~ QUAL        0.5911      0.0946    6.2466    0.0000 [ 0.4077; 0.7796 ] 
    ##   SAT ~ VAL         0.5270      0.0906    5.8149    0.0000 [ 0.3457; 0.7039 ] 
    ##   LOY ~ IMAG        0.4840      0.0694    6.9721    0.0000 [ 0.3433; 0.6079 ] 
    ##   LOY ~ EXPE        0.3142      0.0549    5.7260    0.0000 [ 0.2163; 0.4298 ] 
    ##   LOY ~ QUAL        0.3714      0.0890    4.1730    0.0000 [ 0.2166; 0.5616 ] 
    ##   LOY ~ VAL         0.3311      0.0783    4.2300    0.0000 [ 0.1886; 0.4983 ] 
    ##   LOY ~ SAT         0.6283      0.0844    7.4423    0.0000 [ 0.4788; 0.7952 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0618    6.3594    0.0000 [ 0.2671; 0.5132 ] 
    ##   VAL ~ IMAG           0.2974      0.0607    4.9018    0.0000 [ 0.1786; 0.4188 ] 
    ##   VAL ~ EXPE           0.5852      0.0735    7.9590    0.0000 [ 0.4213; 0.7090 ] 
    ##   SAT ~ IMAG           0.2357      0.0480    4.9100    0.0000 [ 0.1409; 0.3336 ] 
    ##   SAT ~ EXPE           0.5173      0.0666    7.7688    0.0000 [ 0.3883; 0.6571 ] 
    ##   SAT ~ QUAL           0.3696      0.0669    5.5270    0.0000 [ 0.2407; 0.4967 ] 
    ##   LOY ~ IMAG           0.3020      0.0570    5.2940    0.0000 [ 0.1969; 0.4264 ] 
    ##   LOY ~ EXPE           0.3142      0.0549    5.7260    0.0000 [ 0.2163; 0.4298 ] 
    ##   LOY ~ QUAL           0.3714      0.0890    4.1730    0.0000 [ 0.2166; 0.5616 ] 
    ##   LOY ~ VAL            0.3311      0.0783    4.2300    0.0000 [ 0.1886; 0.4983 ] 
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
