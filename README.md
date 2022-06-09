
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

## Purpose

Estimate, analyse, test, and study linear, nonlinear, hierarchical and
multi-group structural equation models using composite-based approaches
and procedures, including estimation techniques such as partial least
squares path modeling (PLS-PM) and its derivatives (PLSc, OrdPLSc,
robustPLSc), generalized structured component analysis (GSCA),
generalized structured component analysis with uniqueness terms (GSCAm),
generalized canonical correlation analysis (GCCA), principal component
analysis (PCA), factor score regression (FSR) using sum score,
regression or Bartlett scores (including bias correction using Croon’s
approach), as well as several tests and typical postestimation
procedures (e.g., verify admissibility of the estimates, assess the
model fit, test the model fit, compute confidence intervals, compare
groups, etc.).

## News (2022-05-09):

-   `predict()` function is now able to predict categorical indicators
    (a procedure known as OrdPLScPredict).

-   Use singular value decomposition in GSCAm to deal with large
    datasets, which allows for large datasets (thanks to Heungsun Hwang)

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

-   `assess()` : assess the model using common quality criteria
-   `infer()` : calculate common inferential quantities (e.g., standard
    errors, confidence intervals)
-   `predict()` : predict endogenous indicator values
-   `summarize()` : summarize the results
-   `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

-   `testOMF()` : performs a test for overall model fit
-   `testMICOM()` : performs a test for composite measurement invariance
-   `testMGD()` : performs several tests to assess multi-group
    differences
-   `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently, three do functions are implemented:

-   `doIPMA()`: performs an importance-performance matrix analysis
-   `doNonlinearEffectsAnalysis()`: performs a nonlinear effects
    analysis such as floodlight and surface analysis
-   `doRedundancyAnalysis()`: performs a redundancy analysis

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
    
## Note: The operator "<~" tells cSEM that the construct to its left is modeled
##       as a composite.
##       The operator "=~" tells cSEM that the construct to its left is modeled
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
    ##        +------------------------------------------------------------------+
    ##        |                                                                  |
    ##        |   H0: The model-implied indicator covariance matrix equals the   |
    ##        |   population indicator covariance matrix.                        |
    ##        |                                                                  |
    ##        +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3239  
    ##  SRMR                    0.0940      0.0519  
    ##  dL                      2.2340      0.6811  
    ##  dML                     2.9219      1.6311  
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
    ##  Out of 499 bootstrap replications 471 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -598052067
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
    ##  Squared Euclidean distance  = 2.23402
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
    ## -------------- Variance inflation factors (VIFs) for modeB weights -------------
    ## 
    ##   Construct: 'IMAG'
    ## 
    ##  Weight    VIF value 
    ##  imag1      1.7215   
    ##  imag2      3.0515   
    ##  imag3      2.5356   
    ## 
    ##   Construct: 'EXPE'
    ## 
    ##  Weight    VIF value 
    ##  expe1      1.4949   
    ##  expe2      1.6623   
    ##  expe3      1.5212   
    ## 
    ##   Construct: 'QUAL'
    ## 
    ##  Weight    VIF value 
    ##  qual1      1.8401   
    ##  qual2      2.5005   
    ##  qual3      1.7796   
    ##  qual4      2.1557   
    ##  qual5      2.0206   
    ## 
    ##   Construct: 'VAL'
    ## 
    ##  Weight    VIF value 
    ##  val1       1.6912   
    ##  val2       2.2049   
    ##  val3       1.6714   
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
    ## ----------------------- Discriminant validity assessment -----------------------
    ## 
    ##  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
    ## 
    ##           SAT LOY
    ## SAT 1.0000000   0
    ## LOY 0.7432489   1
    ## 
    ## 
    ##  Advanced heterotrait-monotrait ratio of correlations matrix (HTMT2 matrix)
    ## 
    ##           SAT LOY
    ## SAT 1.0000000   0
    ## LOY 0.7140046   1
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
    ##  Number of repetitions            = 1
    ##  Handle inadmissibles             = stop
    ##  Estimator target                 = 'PLS-PM'
    ##  Estimator benchmark              = 'lm'
    ##  Disattenuation target            = 'TRUE'
    ##  Disattenuation benchmark         = 'FALSE'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   expe1         1.4598         1.5917       1.9112         2.0997       0.0499
    ##   expe2         1.4118         1.4975       1.9298         2.0268       0.2020
    ##   expe3         1.6254         1.7362       2.1165         2.2141       0.1286
    ##   qual1         1.4713         1.5631       1.9302         2.0760       0.1153
    ##   qual2         1.5688         1.5367       2.0291         2.0508       0.2230
    ##   qual3         1.7244         1.7378       2.2121         2.2819       0.1255
    ##   qual4         1.2306         1.1985       1.5913         1.6291       0.2374
    ##   qual5         1.5000         1.5074       1.9259         1.9441       0.1991
    ##   val1          1.4459         1.3632       1.8660         1.7583       0.2523
    ##   val2          1.2325         1.2286       1.6597         1.7302       0.1701
    ##   val3          1.4780         1.4026       1.9692         1.9487       0.1477
    ##   sat1          1.2425         1.2271       1.6450         1.6148       0.3412
    ##   sat2          1.2306         1.1994       1.6359         1.6239       0.3128
    ##   sat3          1.3403         1.2828       1.6727         1.7142       0.2090
    ##   sat4          1.3183         1.2613       1.6716         1.6393       0.2716
    ##   loy1          1.6896         1.6582       2.2365         2.2312       0.2672
    ##   loy2          1.4821         1.4827       1.9080         1.9790       0.1348
    ##   loy3          1.6983         1.6683       2.2793         2.2712       0.2691
    ##   loy4          1.6837         1.6824       2.1760         2.3054       0.0863
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential statistics are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways for resampling:

1.  Setting `.resample_method` in `csem()` to `"jackknife"` or
    `"bootstrap"` and subsequently using postestimation functions
    `summarize()` or `infer()`.
2.  The same result is achieved by passing a `cSEMResults` object to
    `resamplecSEMResults()` and subsequently using postestimation
    functions `summarize()` or `infer()`.

``` r
# Setting `.resample_method`
b1 <- csem(.data = satisfaction, .model = model, .resample_method = "bootstrap")
# Using resamplecSEMResults()
b2 <- resamplecSEMResults(res)
```

The `summarize()` function reports the inferential statistics:

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
    ##  Number of admissible results     = 484
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 458992545
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
    ##   EXPE ~ IMAG      0.4714      0.0701    6.7239    0.0000 [ 0.3396; 0.5968 ] 
    ##   QUAL ~ EXPE      0.8344      0.0235   35.5787    0.0000 [ 0.7854; 0.8763 ] 
    ##   VAL ~ EXPE       0.0457      0.0880    0.5194    0.6035 [-0.1283; 0.2197 ] 
    ##   VAL ~ QUAL       0.7013      0.0861    8.1440    0.0000 [ 0.5177; 0.8506 ] 
    ##   SAT ~ IMAG       0.2450      0.0520    4.7097    0.0000 [ 0.1529; 0.3459 ] 
    ##   SAT ~ EXPE      -0.0172      0.0675   -0.2553    0.7985 [-0.1506; 0.1139 ] 
    ##   SAT ~ QUAL       0.2215      0.0971    2.2811    0.0225 [ 0.0410; 0.4249 ] 
    ##   SAT ~ VAL        0.5270      0.0859    6.1310    0.0000 [ 0.3454; 0.6930 ] 
    ##   LOY ~ IMAG       0.1819      0.0797    2.2825    0.0225 [ 0.0269; 0.3274 ] 
    ##   LOY ~ SAT        0.6283      0.0825    7.6138    0.0000 [ 0.4793; 0.7904 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0949    6.6458    0.0000 [ 0.4203; 0.8064 ] 
    ##   IMAG =~ imag2      0.9246      0.0408   22.6710    0.0000 [ 0.8203; 0.9758 ] 
    ##   IMAG =~ imag3      0.9577      0.0288   33.2485    0.0000 [ 0.8892; 0.9935 ] 
    ##   EXPE =~ expe1      0.7525      0.0737   10.2134    0.0000 [ 0.5780; 0.8695 ] 
    ##   EXPE =~ expe2      0.9348      0.0288   32.5057    0.0000 [ 0.8624; 0.9715 ] 
    ##   EXPE =~ expe3      0.7295      0.0705   10.3457    0.0000 [ 0.5594; 0.8338 ] 
    ##   QUAL =~ qual1      0.7861      0.0679   11.5784    0.0000 [ 0.6387; 0.8892 ] 
    ##   QUAL =~ qual2      0.9244      0.0230   40.1452    0.0000 [ 0.8689; 0.9566 ] 
    ##   QUAL =~ qual3      0.7560      0.0578   13.0758    0.0000 [ 0.6329; 0.8452 ] 
    ##   QUAL =~ qual4      0.7632      0.0531   14.3591    0.0000 [ 0.6400; 0.8520 ] 
    ##   QUAL =~ qual5      0.7834      0.0480   16.3363    0.0000 [ 0.6829; 0.8686 ] 
    ##   VAL =~ val1        0.9518      0.0245   38.7998    0.0000 [ 0.8899; 0.9840 ] 
    ##   VAL =~ val2        0.8056      0.0620   12.9964    0.0000 [ 0.6711; 0.9100 ] 
    ##   VAL =~ val3        0.6763      0.0730    9.2688    0.0000 [ 0.5341; 0.8152 ] 
    ##   SAT =~ sat1        0.9243      0.0217   42.5938    0.0000 [ 0.8748; 0.9628 ] 
    ##   SAT =~ sat2        0.8813      0.0294   29.9613    0.0000 [ 0.8160; 0.9321 ] 
    ##   SAT =~ sat3        0.7127      0.0534   13.3572    0.0000 [ 0.5974; 0.8088 ] 
    ##   SAT =~ sat4        0.7756      0.0484   16.0162    0.0000 [ 0.6748; 0.8533 ] 
    ##   LOY =~ loy1        0.9097      0.0473   19.2243    0.0000 [ 0.8068; 0.9836 ] 
    ##   LOY =~ loy2        0.5775      0.0861    6.7107    0.0000 [ 0.4047; 0.7367 ] 
    ##   LOY =~ loy3        0.9043      0.0430   21.0296    0.0000 [ 0.8088; 0.9672 ] 
    ##   LOY =~ loy4        0.4917      0.1009    4.8755    0.0000 [ 0.3013; 0.6972 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1097    0.1426    0.8866 [-0.1993; 0.2235 ] 
    ##   IMAG <~ imag2      0.4473      0.1487    3.0076    0.0026 [ 0.1250; 0.6928 ] 
    ##   IMAG <~ imag3      0.6020      0.1399    4.3027    0.0000 [ 0.3283; 0.8596 ] 
    ##   EXPE <~ expe1      0.2946      0.1112    2.6501    0.0080 [ 0.0837; 0.4959 ] 
    ##   EXPE <~ expe2      0.6473      0.0842    7.6840    0.0000 [ 0.4546; 0.7877 ] 
    ##   EXPE <~ expe3      0.2374      0.0903    2.6299    0.0085 [ 0.0704; 0.4090 ] 
    ##   QUAL <~ qual1      0.2370      0.0901    2.6316    0.0085 [ 0.0855; 0.4304 ] 
    ##   QUAL <~ qual2      0.4712      0.0765    6.1624    0.0000 [ 0.2972; 0.6050 ] 
    ##   QUAL <~ qual3      0.1831      0.0783    2.3387    0.0193 [ 0.0391; 0.3303 ] 
    ##   QUAL <~ qual4      0.1037      0.0587    1.7682    0.0770 [-0.0191; 0.2215 ] 
    ##   QUAL <~ qual5      0.2049      0.0634    3.2307    0.0012 [ 0.0820; 0.3255 ] 
    ##   VAL <~ val1        0.7163      0.0976    7.3430    0.0000 [ 0.4969; 0.8728 ] 
    ##   VAL <~ val2        0.2202      0.0915    2.4054    0.0162 [ 0.0557; 0.4167 ] 
    ##   VAL <~ val3        0.2082      0.0589    3.5364    0.0004 [ 0.0949; 0.3193 ] 
    ##   SAT <~ sat1        0.3209      0.0152   21.0759    0.0000 [ 0.2976; 0.3580 ] 
    ##   SAT <~ sat2        0.3059      0.0139   22.0723    0.0000 [ 0.2826; 0.3382 ] 
    ##   SAT <~ sat3        0.2474      0.0111   22.3542    0.0000 [ 0.2246; 0.2683 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.0954    0.0000 [ 0.2468; 0.2954 ] 
    ##   LOY <~ loy1        0.3834      0.0270   14.1967    0.0000 [ 0.3294; 0.4342 ] 
    ##   LOY <~ loy2        0.2434      0.0295    8.2502    0.0000 [ 0.1828; 0.2955 ] 
    ##   LOY <~ loy3        0.3812      0.0282   13.4974    0.0000 [ 0.3236; 0.4349 ] 
    ##   LOY <~ loy4        0.2073      0.0367    5.6469    0.0000 [ 0.1361; 0.2768 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0617   10.4315    0.0000 [ 0.5219; 0.7554 ] 
    ##   imag1 ~~ imag3      0.5433      0.0685    7.9344    0.0000 [ 0.4158; 0.6643 ] 
    ##   imag2 ~~ imag3      0.7761      0.0383   20.2660    0.0000 [ 0.6993; 0.8461 ] 
    ##   expe1 ~~ expe2      0.5353      0.0583    9.1804    0.0000 [ 0.4058; 0.6353 ] 
    ##   expe1 ~~ expe3      0.4694      0.0604    7.7718    0.0000 [ 0.3499; 0.5793 ] 
    ##   expe2 ~~ expe3      0.5467      0.0574    9.5214    0.0000 [ 0.4135; 0.6360 ] 
    ##   qual1 ~~ qual2      0.6053      0.0557   10.8618    0.0000 [ 0.4894; 0.7000 ] 
    ##   qual1 ~~ qual3      0.5406      0.0555    9.7478    0.0000 [ 0.4181; 0.6360 ] 
    ##   qual1 ~~ qual4      0.5662      0.0668    8.4699    0.0000 [ 0.4223; 0.6800 ] 
    ##   qual1 ~~ qual5      0.5180      0.0675    7.6708    0.0000 [ 0.3868; 0.6430 ] 
    ##   qual2 ~~ qual3      0.6187      0.0541   11.4361    0.0000 [ 0.5011; 0.7056 ] 
    ##   qual2 ~~ qual4      0.6517      0.0626   10.4131    0.0000 [ 0.5315; 0.7656 ] 
    ##   qual2 ~~ qual5      0.6291      0.0571   11.0091    0.0000 [ 0.5189; 0.7382 ] 
    ##   qual3 ~~ qual4      0.4752      0.0646    7.3556    0.0000 [ 0.3455; 0.5847 ] 
    ##   qual3 ~~ qual5      0.5074      0.0612    8.2962    0.0000 [ 0.3836; 0.6183 ] 
    ##   qual4 ~~ qual5      0.6402      0.0576   11.1192    0.0000 [ 0.5179; 0.7456 ] 
    ##   val1 ~~ val2        0.6344      0.0522   12.1432    0.0000 [ 0.5198; 0.7292 ] 
    ##   val1 ~~ val3        0.4602      0.0689    6.6755    0.0000 [ 0.3321; 0.6053 ] 
    ##   val2 ~~ val3        0.6288      0.0616   10.2077    0.0000 [ 0.5091; 0.7421 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0701    6.7239    0.0000 [ 0.3396; 0.5968 ] 
    ##   QUAL ~ IMAG       0.3933      0.0650    6.0512    0.0000 [ 0.2756; 0.5093 ] 
    ##   QUAL ~ EXPE       0.8344      0.0235   35.5787    0.0000 [ 0.7854; 0.8763 ] 
    ##   VAL ~ IMAG        0.2974      0.0633    4.6957    0.0000 [ 0.1875; 0.4241 ] 
    ##   VAL ~ EXPE        0.6309      0.0515   12.2580    0.0000 [ 0.5284; 0.7267 ] 
    ##   VAL ~ QUAL        0.7013      0.0861    8.1440    0.0000 [ 0.5177; 0.8506 ] 
    ##   SAT ~ IMAG        0.4807      0.0667    7.2037    0.0000 [ 0.3437; 0.6051 ] 
    ##   SAT ~ EXPE        0.5001      0.0534    9.3638    0.0000 [ 0.3983; 0.6042 ] 
    ##   SAT ~ QUAL        0.5911      0.0922    6.4139    0.0000 [ 0.4048; 0.7652 ] 
    ##   SAT ~ VAL         0.5270      0.0859    6.1310    0.0000 [ 0.3454; 0.6930 ] 
    ##   LOY ~ IMAG        0.4840      0.0677    7.1479    0.0000 [ 0.3535; 0.6133 ] 
    ##   LOY ~ EXPE        0.3142      0.0560    5.6078    0.0000 [ 0.2204; 0.4310 ] 
    ##   LOY ~ QUAL        0.3714      0.0835    4.4486    0.0000 [ 0.2201; 0.5382 ] 
    ##   LOY ~ VAL         0.3311      0.0736    4.4995    0.0000 [ 0.1991; 0.4897 ] 
    ##   LOY ~ SAT         0.6283      0.0825    7.6138    0.0000 [ 0.4793; 0.7904 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0650    6.0512    0.0000 [ 0.2756; 0.5093 ] 
    ##   VAL ~ IMAG           0.2974      0.0633    4.6957    0.0000 [ 0.1875; 0.4241 ] 
    ##   VAL ~ EXPE           0.5852      0.0742    7.8871    0.0000 [ 0.4370; 0.7159 ] 
    ##   SAT ~ IMAG           0.2357      0.0479    4.9247    0.0000 [ 0.1528; 0.3341 ] 
    ##   SAT ~ EXPE           0.5173      0.0670    7.7192    0.0000 [ 0.3930; 0.6511 ] 
    ##   SAT ~ QUAL           0.3696      0.0628    5.8819    0.0000 [ 0.2546; 0.4928 ] 
    ##   LOY ~ IMAG           0.3020      0.0562    5.3702    0.0000 [ 0.2040; 0.4254 ] 
    ##   LOY ~ EXPE           0.3142      0.0560    5.6078    0.0000 [ 0.2204; 0.4310 ] 
    ##   LOY ~ QUAL           0.3714      0.0835    4.4486    0.0000 [ 0.2201; 0.5382 ] 
    ##   LOY ~ VAL            0.3311      0.0736    4.4995    0.0000 [ 0.1991; 0.4897 ] 
    ## ________________________________________________________________________________

Several bootstrap-based confidence intervals are implemented, see
`?infer()`:

``` r
infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet
```

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as setting random seeds via the [future
framework](https://github.com/HenrikBengtsson/future). For
multiprocessing simply set `.eval_plan = "multisession"` in which case
the maximum number of available cores is used if not on Windows. On
Windows as many separate R instances are opened in the background as
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
  .eval_plan       = "multisession")
```
