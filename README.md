
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
>     syntax](https://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use `csem()`
> 3.  Apply one of the postestimation functions listed below on the
>     resulting object.

## Postestimation functions

There are five major postestimation verbs, three test family functions
and three do-family of function:

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
syntax](https://lavaan.ugent.be/tutorial/syntax1.html) with some slight
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
    ##  dG                      0.6493      0.3245  
    ##  SRMR                    0.0940      0.0524  
    ##  dL                      2.2340      0.6952  
    ##  dML                     2.9219      1.6185  
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
    ##  Out of 499 bootstrap replications 477 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 891891725
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
    ##   expe1         1.4543         1.5652       1.9036         2.0918       0.0570
    ##   expe2         1.4129         1.4800       1.9337         2.0261       0.2013
    ##   expe3         1.6326         1.7315       2.1255         2.2254       0.1252
    ##   qual1         1.4733         1.5447       1.9249         2.0615       0.1179
    ##   qual2         1.5792         1.5381       2.0406         2.0616       0.2176
    ##   qual3         1.7330         1.7305       2.2222         2.2855       0.1195
    ##   qual4         1.2344         1.1948       1.5957         1.6257       0.2347
    ##   qual5         1.5047         1.5012       1.9339         1.9524       0.1967
    ##   val1          1.4462         1.3657       1.8705         1.7686       0.2496
    ##   val2          1.2258         1.2054       1.6482         1.7122       0.1729
    ##   val3          1.4806         1.3806       1.9684         1.9330       0.1489
    ##   sat1          1.2467         1.2361       1.6438         1.6221       0.3402
    ##   sat2          1.2292         1.1968       1.6368         1.6273       0.3113
    ##   sat3          1.3388         1.2748       1.6699         1.7201       0.2120
    ##   sat4          1.3145         1.2613       1.6635         1.6342       0.2805
    ##   loy1          1.6905         1.6574       2.2340         2.2220       0.2680
    ##   loy2          1.4842         1.4719       1.9105         1.9775       0.1327
    ##   loy3          1.7028         1.6677       2.2797         2.2685       0.2702
    ##   loy4          1.6874         1.6654       2.1756         2.2940       0.0878
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
    ##  Number of admissible results     = 489
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 842990549
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
    ##   EXPE ~ IMAG      0.4714      0.0623    7.5621    0.0000 [ 0.3470; 0.5977 ] 
    ##   QUAL ~ EXPE      0.8344      0.0234   35.6303    0.0000 [ 0.7892; 0.8718 ] 
    ##   VAL ~ EXPE       0.0457      0.0846    0.5403    0.5890 [-0.1084; 0.2299 ] 
    ##   VAL ~ QUAL       0.7013      0.0784    8.9506    0.0000 [ 0.5368; 0.8502 ] 
    ##   SAT ~ IMAG       0.2450      0.0558    4.3924    0.0000 [ 0.1304; 0.3496 ] 
    ##   SAT ~ EXPE      -0.0172      0.0710   -0.2428    0.8081 [-0.1609; 0.1173 ] 
    ##   SAT ~ QUAL       0.2215      0.0984    2.2518    0.0243 [ 0.0532; 0.4280 ] 
    ##   SAT ~ VAL        0.5270      0.0842    6.2554    0.0000 [ 0.3594; 0.6774 ] 
    ##   LOY ~ IMAG       0.1819      0.0814    2.2342    0.0255 [ 0.0328; 0.3488 ] 
    ##   LOY ~ SAT        0.6283      0.0817    7.6949    0.0000 [ 0.4652; 0.7746 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0968    6.5131    0.0000 [ 0.4110; 0.7978 ] 
    ##   IMAG =~ imag2      0.9246      0.0392   23.5839    0.0000 [ 0.8368; 0.9790 ] 
    ##   IMAG =~ imag3      0.9577      0.0281   34.0459    0.0000 [ 0.8810; 0.9905 ] 
    ##   EXPE =~ expe1      0.7525      0.0739   10.1781    0.0000 [ 0.5710; 0.8703 ] 
    ##   EXPE =~ expe2      0.9348      0.0270   34.6519    0.0000 [ 0.8704; 0.9716 ] 
    ##   EXPE =~ expe3      0.7295      0.0736    9.9161    0.0000 [ 0.5655; 0.8512 ] 
    ##   QUAL =~ qual1      0.7861      0.0640   12.2923    0.0000 [ 0.6435; 0.8826 ] 
    ##   QUAL =~ qual2      0.9244      0.0218   42.3976    0.0000 [ 0.8692; 0.9568 ] 
    ##   QUAL =~ qual3      0.7560      0.0634   11.9210    0.0000 [ 0.6088; 0.8548 ] 
    ##   QUAL =~ qual4      0.7632      0.0520   14.6809    0.0000 [ 0.6513; 0.8531 ] 
    ##   QUAL =~ qual5      0.7834      0.0482   16.2425    0.0000 [ 0.6773; 0.8545 ] 
    ##   VAL =~ val1        0.9518      0.0233   40.7882    0.0000 [ 0.8949; 0.9825 ] 
    ##   VAL =~ val2        0.8056      0.0607   13.2795    0.0000 [ 0.6800; 0.9067 ] 
    ##   VAL =~ val3        0.6763      0.0731    9.2486    0.0000 [ 0.5265; 0.7993 ] 
    ##   SAT =~ sat1        0.9243      0.0228   40.4635    0.0000 [ 0.8716; 0.9597 ] 
    ##   SAT =~ sat2        0.8813      0.0294   30.0008    0.0000 [ 0.8121; 0.9318 ] 
    ##   SAT =~ sat3        0.7127      0.0542   13.1468    0.0000 [ 0.5968; 0.8020 ] 
    ##   SAT =~ sat4        0.7756      0.0488   15.8852    0.0000 [ 0.6626; 0.8546 ] 
    ##   LOY =~ loy1        0.9097      0.0475   19.1362    0.0000 [ 0.8023; 0.9829 ] 
    ##   LOY =~ loy2        0.5775      0.0835    6.9170    0.0000 [ 0.3967; 0.7344 ] 
    ##   LOY =~ loy3        0.9043      0.0421   21.5012    0.0000 [ 0.8078; 0.9712 ] 
    ##   LOY =~ loy4        0.4917      0.0928    5.2964    0.0000 [ 0.3218; 0.6737 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1166    0.1342    0.8932 [-0.2235; 0.2426 ] 
    ##   IMAG <~ imag2      0.4473      0.1463    3.0565    0.0022 [ 0.1685; 0.7181 ] 
    ##   IMAG <~ imag3      0.6020      0.1366    4.4089    0.0000 [ 0.3084; 0.8311 ] 
    ##   EXPE <~ expe1      0.2946      0.1147    2.5692    0.0102 [ 0.0727; 0.5230 ] 
    ##   EXPE <~ expe2      0.6473      0.0840    7.7075    0.0000 [ 0.4600; 0.7901 ] 
    ##   EXPE <~ expe3      0.2374      0.0950    2.4985    0.0125 [ 0.0240; 0.4014 ] 
    ##   QUAL <~ qual1      0.2370      0.0890    2.6633    0.0077 [ 0.0787; 0.4108 ] 
    ##   QUAL <~ qual2      0.4712      0.0768    6.1356    0.0000 [ 0.3004; 0.6111 ] 
    ##   QUAL <~ qual3      0.1831      0.0812    2.2548    0.0241 [ 0.0033; 0.3334 ] 
    ##   QUAL <~ qual4      0.1037      0.0645    1.6085    0.1077 [-0.0159; 0.2317 ] 
    ##   QUAL <~ qual5      0.2049      0.0656    3.1224    0.0018 [ 0.0831; 0.3264 ] 
    ##   VAL <~ val1        0.7163      0.0934    7.6659    0.0000 [ 0.5069; 0.8638 ] 
    ##   VAL <~ val2        0.2202      0.0879    2.5050    0.0122 [ 0.0738; 0.4080 ] 
    ##   VAL <~ val3        0.2082      0.0623    3.3436    0.0008 [ 0.0817; 0.3272 ] 
    ##   SAT <~ sat1        0.3209      0.0155   20.7568    0.0000 [ 0.2959; 0.3541 ] 
    ##   SAT <~ sat2        0.3059      0.0141   21.7558    0.0000 [ 0.2841; 0.3380 ] 
    ##   SAT <~ sat3        0.2474      0.0117   21.0768    0.0000 [ 0.2238; 0.2695 ] 
    ##   SAT <~ sat4        0.2692      0.0116   23.2330    0.0000 [ 0.2461; 0.2919 ] 
    ##   LOY <~ loy1        0.3834      0.0242   15.8698    0.0000 [ 0.3377; 0.4269 ] 
    ##   LOY <~ loy2        0.2434      0.0297    8.1888    0.0000 [ 0.1795; 0.2886 ] 
    ##   LOY <~ loy3        0.3812      0.0266   14.3510    0.0000 [ 0.3274; 0.4315 ] 
    ##   LOY <~ loy4        0.2073      0.0347    5.9685    0.0000 [ 0.1429; 0.2774 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0615   10.4704    0.0000 [ 0.5201; 0.7575 ] 
    ##   imag1 ~~ imag3      0.5433      0.0667    8.1414    0.0000 [ 0.4219; 0.6646 ] 
    ##   imag2 ~~ imag3      0.7761      0.0402   19.2995    0.0000 [ 0.6946; 0.8432 ] 
    ##   expe1 ~~ expe2      0.5353      0.0548    9.7731    0.0000 [ 0.4204; 0.6372 ] 
    ##   expe1 ~~ expe3      0.4694      0.0591    7.9377    0.0000 [ 0.3502; 0.5874 ] 
    ##   expe2 ~~ expe3      0.5467      0.0587    9.3106    0.0000 [ 0.4285; 0.6555 ] 
    ##   qual1 ~~ qual2      0.6053      0.0556   10.8895    0.0000 [ 0.4787; 0.7080 ] 
    ##   qual1 ~~ qual3      0.5406      0.0608    8.8891    0.0000 [ 0.4053; 0.6542 ] 
    ##   qual1 ~~ qual4      0.5662      0.0633    8.9457    0.0000 [ 0.4378; 0.6802 ] 
    ##   qual1 ~~ qual5      0.5180      0.0673    7.6936    0.0000 [ 0.3756; 0.6471 ] 
    ##   qual2 ~~ qual3      0.6187      0.0547   11.3100    0.0000 [ 0.5020; 0.7089 ] 
    ##   qual2 ~~ qual4      0.6517      0.0598   10.8960    0.0000 [ 0.5257; 0.7561 ] 
    ##   qual2 ~~ qual5      0.6291      0.0547   11.4977    0.0000 [ 0.5044; 0.7233 ] 
    ##   qual3 ~~ qual4      0.4752      0.0629    7.5534    0.0000 [ 0.3493; 0.5842 ] 
    ##   qual3 ~~ qual5      0.5074      0.0628    8.0788    0.0000 [ 0.3844; 0.6222 ] 
    ##   qual4 ~~ qual5      0.6402      0.0546   11.7207    0.0000 [ 0.5271; 0.7342 ] 
    ##   val1 ~~ val2        0.6344      0.0549   11.5579    0.0000 [ 0.5234; 0.7387 ] 
    ##   val1 ~~ val3        0.4602      0.0673    6.8340    0.0000 [ 0.3136; 0.5915 ] 
    ##   val2 ~~ val3        0.6288      0.0605   10.3853    0.0000 [ 0.5122; 0.7373 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0623    7.5621    0.0000 [ 0.3470; 0.5977 ] 
    ##   QUAL ~ IMAG       0.3933      0.0582    6.7630    0.0000 [ 0.2851; 0.5149 ] 
    ##   QUAL ~ EXPE       0.8344      0.0234   35.6303    0.0000 [ 0.7892; 0.8718 ] 
    ##   VAL ~ IMAG        0.2974      0.0581    5.1166    0.0000 [ 0.1923; 0.4164 ] 
    ##   VAL ~ EXPE        0.6309      0.0487   12.9538    0.0000 [ 0.5330; 0.7266 ] 
    ##   VAL ~ QUAL        0.7013      0.0784    8.9506    0.0000 [ 0.5368; 0.8502 ] 
    ##   SAT ~ IMAG        0.4807      0.0672    7.1482    0.0000 [ 0.3415; 0.6014 ] 
    ##   SAT ~ EXPE        0.5001      0.0552    9.0654    0.0000 [ 0.3899; 0.6090 ] 
    ##   SAT ~ QUAL        0.5911      0.0926    6.3833    0.0000 [ 0.4205; 0.7679 ] 
    ##   SAT ~ VAL         0.5270      0.0842    6.2554    0.0000 [ 0.3594; 0.6774 ] 
    ##   LOY ~ IMAG        0.4840      0.0710    6.8183    0.0000 [ 0.3480; 0.6148 ] 
    ##   LOY ~ EXPE        0.3142      0.0540    5.8202    0.0000 [ 0.2122; 0.4293 ] 
    ##   LOY ~ QUAL        0.3714      0.0828    4.4839    0.0000 [ 0.2363; 0.5527 ] 
    ##   LOY ~ VAL         0.3311      0.0738    4.4886    0.0000 [ 0.1918; 0.4761 ] 
    ##   LOY ~ SAT         0.6283      0.0817    7.6949    0.0000 [ 0.4652; 0.7746 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0582    6.7630    0.0000 [ 0.2851; 0.5149 ] 
    ##   VAL ~ IMAG           0.2974      0.0581    5.1166    0.0000 [ 0.1923; 0.4164 ] 
    ##   VAL ~ EXPE           0.5852      0.0681    8.5879    0.0000 [ 0.4400; 0.7238 ] 
    ##   SAT ~ IMAG           0.2357      0.0463    5.0907    0.0000 [ 0.1504; 0.3332 ] 
    ##   SAT ~ EXPE           0.5173      0.0667    7.7552    0.0000 [ 0.3951; 0.6567 ] 
    ##   SAT ~ QUAL           0.3696      0.0625    5.9168    0.0000 [ 0.2451; 0.4865 ] 
    ##   LOY ~ IMAG           0.3020      0.0547    5.5220    0.0000 [ 0.2025; 0.4111 ] 
    ##   LOY ~ EXPE           0.3142      0.0540    5.8202    0.0000 [ 0.2122; 0.4293 ] 
    ##   LOY ~ QUAL           0.3714      0.0828    4.4839    0.0000 [ 0.2363; 0.5527 ] 
    ##   LOY ~ VAL            0.3311      0.0738    4.4886    0.0000 [ 0.1918; 0.4761 ] 
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
