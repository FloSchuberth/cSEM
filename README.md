
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
    ##  Estimation status                  = Ok
    ##  Number of observations             = 250
    ##  Weight estimator                   = PLS-PM
    ##  Inner weighting scheme             = "path"
    ##  Type of indicator correlation      = Pearson
    ##  Path model estimator               = OLS
    ##  Second-order approach              = NA
    ##  Type of path model                 = Linear
    ##  Disattenuated                      = Yes (PLSc)
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
    ##        ┌──────────────────────────────────────────────────────────────────┐
    ##        │                                                                  │
    ##        │   H0: The model-implied indicator covariance matrix equals the   │
    ##        │   population indicator covariance matrix.                        │
    ##        │                                                                  │
    ##        └──────────────────────────────────────────────────────────────────┘
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3199  
    ##  SRMR                    0.0940      0.0518  
    ##  dL                      2.2340      0.6776  
    ##  dML                     2.9219      1.5557  
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
    ##  Out of 499 bootstrap replications 468 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -464649022
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
    ##  Geodesic distance             = 0.6493432
    ##  Squared Euclidean distance    = 2.23402
    ##  ML distance                   = 2.921932
    ## 
    ##  Chi_square       = 727.5611
    ##  Chi_square_df    = 3.954137
    ##  CFI              = 0.8598825
    ##  CN               = 75.14588
    ##  GFI              = 0.7280612
    ##  IFI              = 0.8615598
    ##  NFI              = 0.8229918
    ##  NNFI             = 0.8240917
    ##  RMSEA            = 0.108922
    ##  RMS_theta        = 0.05069299
    ##  SRMR             = 0.09396871
    ## 
    ##  Degrees of freedom       = 184
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
    ##  Number of obs. training            = 225
    ##  Number of obs. test                = 25
    ##  Number of cv folds                 = 10
    ##  Number of repetitions              = 1
    ##  Handle inadmissibles               = stop
    ##  Estimator target                   = 'PLS-PM'
    ##  Estimator benchmark                = 'lm'
    ##  Disattenuation target              = 'TRUE'
    ##  Disattenuation benchmark           = 'FALSE'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   expe1         1.4636         1.5861       1.9173         2.0928       0.0492
    ##   expe2         1.4183         1.4928       1.9404         2.0247       0.1954
    ##   expe3         1.6283         1.7381       2.1238         2.2202       0.1250
    ##   qual1         1.4801         1.5542       1.9341         2.0550       0.1145
    ##   qual2         1.5918         1.5564       2.0624         2.0859       0.2094
    ##   qual3         1.7337         1.7504       2.2261         2.2952       0.1197
    ##   qual4         1.2387         1.1993       1.6026         1.6525       0.2256
    ##   qual5         1.5191         1.5114       1.9505         1.9670       0.1840
    ##   val1          1.4563         1.3747       1.8868         1.7819       0.2410
    ##   val2          1.2320         1.2238       1.6538         1.7264       0.1702
    ##   val3          1.4847         1.3894       1.9731         1.9326       0.1471
    ##   sat1          1.2529         1.2369       1.6558         1.6316       0.3348
    ##   sat2          1.2406         1.2109       1.6564         1.6407       0.2983
    ##   sat3          1.3463         1.3083       1.6848         1.7557       0.2059
    ##   sat4          1.3176         1.2601       1.6650         1.6322       0.2825
    ##   loy1          1.6933         1.6679       2.2309         2.2317       0.2723
    ##   loy2          1.4846         1.4868       1.9080         1.9831       0.1343
    ##   loy3          1.7210         1.6872       2.3015         2.2909       0.2698
    ##   loy4          1.6908         1.6881       2.1798         2.3134       0.0830
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
    ##  Estimation status                  = Ok
    ##  Number of observations             = 250
    ##  Weight estimator                   = PLS-PM
    ##  Inner weighting scheme             = "path"
    ##  Type of indicator correlation      = Pearson
    ##  Path model estimator               = OLS
    ##  Second-order approach              = NA
    ##  Type of path model                 = Linear
    ##  Disattenuated                      = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample method                    = "bootstrap"
    ##  Number of resamples                = 499
    ##  Number of admissible results       = 482
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 1326272388
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
    ##   EXPE ~ IMAG      0.4714      0.0613    7.6922    0.0000 [ 0.3391; 0.5903 ] 
    ##   QUAL ~ EXPE      0.8344      0.0223   37.3800    0.0000 [ 0.7881; 0.8753 ] 
    ##   VAL ~ EXPE       0.0457      0.0842    0.5432    0.5870 [-0.1063; 0.2152 ] 
    ##   VAL ~ QUAL       0.7013      0.0840    8.3475    0.0000 [ 0.5285; 0.8503 ] 
    ##   SAT ~ IMAG       0.2450      0.0542    4.5235    0.0000 [ 0.1421; 0.3475 ] 
    ##   SAT ~ EXPE      -0.0172      0.0761   -0.2265    0.8208 [-0.1707; 0.1277 ] 
    ##   SAT ~ QUAL       0.2215      0.1085    2.0426    0.0411 [ 0.0390; 0.4510 ] 
    ##   SAT ~ VAL        0.5270      0.0901    5.8505    0.0000 [ 0.3501; 0.6904 ] 
    ##   LOY ~ IMAG       0.1819      0.0802    2.2675    0.0234 [ 0.0330; 0.3350 ] 
    ##   LOY ~ SAT        0.6283      0.0806    7.7916    0.0000 [ 0.4737; 0.7962 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1008    6.2566    0.0000 [ 0.4187; 0.8201 ] 
    ##   IMAG =~ imag2      0.9246      0.0396   23.3305    0.0000 [ 0.8217; 0.9769 ] 
    ##   IMAG =~ imag3      0.9577      0.0290   32.9821    0.0000 [ 0.8775; 0.9926 ] 
    ##   EXPE =~ expe1      0.7525      0.0788    9.5530    0.0000 [ 0.5794; 0.8728 ] 
    ##   EXPE =~ expe2      0.9348      0.0270   34.5866    0.0000 [ 0.8659; 0.9724 ] 
    ##   EXPE =~ expe3      0.7295      0.0713   10.2253    0.0000 [ 0.5887; 0.8494 ] 
    ##   QUAL =~ qual1      0.7861      0.0680   11.5650    0.0000 [ 0.6183; 0.8861 ] 
    ##   QUAL =~ qual2      0.9244      0.0217   42.5752    0.0000 [ 0.8645; 0.9543 ] 
    ##   QUAL =~ qual3      0.7560      0.0581   13.0115    0.0000 [ 0.6186; 0.8535 ] 
    ##   QUAL =~ qual4      0.7632      0.0518   14.7319    0.0000 [ 0.6384; 0.8486 ] 
    ##   QUAL =~ qual5      0.7834      0.0488   16.0554    0.0000 [ 0.6749; 0.8587 ] 
    ##   VAL =~ val1        0.9518      0.0228   41.7463    0.0000 [ 0.8998; 0.9830 ] 
    ##   VAL =~ val2        0.8056      0.0635   12.6818    0.0000 [ 0.6695; 0.9076 ] 
    ##   VAL =~ val3        0.6763      0.0738    9.1579    0.0000 [ 0.5191; 0.8028 ] 
    ##   SAT =~ sat1        0.9243      0.0227   40.7096    0.0000 [ 0.8759; 0.9624 ] 
    ##   SAT =~ sat2        0.8813      0.0280   31.5064    0.0000 [ 0.8192; 0.9281 ] 
    ##   SAT =~ sat3        0.7127      0.0521   13.6786    0.0000 [ 0.5993; 0.8070 ] 
    ##   SAT =~ sat4        0.7756      0.0485   16.0069    0.0000 [ 0.6761; 0.8575 ] 
    ##   LOY =~ loy1        0.9097      0.0505   17.9971    0.0000 [ 0.7932; 0.9796 ] 
    ##   LOY =~ loy2        0.5775      0.0862    6.6979    0.0000 [ 0.3940; 0.7281 ] 
    ##   LOY =~ loy3        0.9043      0.0422   21.4115    0.0000 [ 0.8147; 0.9772 ] 
    ##   LOY =~ loy4        0.4917      0.0947    5.1946    0.0000 [ 0.3163; 0.6824 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1199    0.1305    0.8962 [-0.2377; 0.2731 ] 
    ##   IMAG <~ imag2      0.4473      0.1484    3.0138    0.0026 [ 0.1282; 0.7429 ] 
    ##   IMAG <~ imag3      0.6020      0.1404    4.2873    0.0000 [ 0.3127; 0.8576 ] 
    ##   EXPE <~ expe1      0.2946      0.1139    2.5855    0.0097 [ 0.0645; 0.5052 ] 
    ##   EXPE <~ expe2      0.6473      0.0844    7.6655    0.0000 [ 0.4627; 0.7903 ] 
    ##   EXPE <~ expe3      0.2374      0.0952    2.4928    0.0127 [ 0.0400; 0.4311 ] 
    ##   QUAL <~ qual1      0.2370      0.0900    2.6349    0.0084 [ 0.0713; 0.4184 ] 
    ##   QUAL <~ qual2      0.4712      0.0810    5.8162    0.0000 [ 0.3034; 0.6182 ] 
    ##   QUAL <~ qual3      0.1831      0.0783    2.3393    0.0193 [ 0.0049; 0.3223 ] 
    ##   QUAL <~ qual4      0.1037      0.0604    1.7161    0.0861 [-0.0055; 0.2228 ] 
    ##   QUAL <~ qual5      0.2049      0.0665    3.0806    0.0021 [ 0.0699; 0.3305 ] 
    ##   VAL <~ val1        0.7163      0.0944    7.5851    0.0000 [ 0.5330; 0.8730 ] 
    ##   VAL <~ val2        0.2202      0.0922    2.3884    0.0169 [ 0.0471; 0.4119 ] 
    ##   VAL <~ val3        0.2082      0.0609    3.4176    0.0006 [ 0.0868; 0.3222 ] 
    ##   SAT <~ sat1        0.3209      0.0147   21.7866    0.0000 [ 0.2969; 0.3549 ] 
    ##   SAT <~ sat2        0.3059      0.0141   21.7153    0.0000 [ 0.2825; 0.3380 ] 
    ##   SAT <~ sat3        0.2474      0.0110   22.4943    0.0000 [ 0.2246; 0.2691 ] 
    ##   SAT <~ sat4        0.2692      0.0117   22.9935    0.0000 [ 0.2474; 0.2912 ] 
    ##   LOY <~ loy1        0.3834      0.0266   14.4290    0.0000 [ 0.3278; 0.4325 ] 
    ##   LOY <~ loy2        0.2434      0.0305    7.9881    0.0000 [ 0.1770; 0.2956 ] 
    ##   LOY <~ loy3        0.3812      0.0271   14.0831    0.0000 [ 0.3286; 0.4345 ] 
    ##   LOY <~ loy4        0.2073      0.0348    5.9602    0.0000 [ 0.1406; 0.2759 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0638   10.0911    0.0000 [ 0.5133; 0.7530 ] 
    ##   imag1 ~~ imag3      0.5433      0.0696    7.8079    0.0000 [ 0.3992; 0.6757 ] 
    ##   imag2 ~~ imag3      0.7761      0.0387   20.0364    0.0000 [ 0.6848; 0.8386 ] 
    ##   expe1 ~~ expe2      0.5353      0.0609    8.7881    0.0000 [ 0.4153; 0.6509 ] 
    ##   expe1 ~~ expe3      0.4694      0.0630    7.4477    0.0000 [ 0.3493; 0.5916 ] 
    ##   expe2 ~~ expe3      0.5467      0.0593    9.2256    0.0000 [ 0.4258; 0.6619 ] 
    ##   qual1 ~~ qual2      0.6053      0.0552   10.9619    0.0000 [ 0.4905; 0.6955 ] 
    ##   qual1 ~~ qual3      0.5406      0.0609    8.8784    0.0000 [ 0.4221; 0.6467 ] 
    ##   qual1 ~~ qual4      0.5662      0.0692    8.1779    0.0000 [ 0.4357; 0.6856 ] 
    ##   qual1 ~~ qual5      0.5180      0.0677    7.6494    0.0000 [ 0.3866; 0.6416 ] 
    ##   qual2 ~~ qual3      0.6187      0.0520   11.8872    0.0000 [ 0.5115; 0.7056 ] 
    ##   qual2 ~~ qual4      0.6517      0.0589   11.0576    0.0000 [ 0.5233; 0.7495 ] 
    ##   qual2 ~~ qual5      0.6291      0.0563   11.1649    0.0000 [ 0.5124; 0.7280 ] 
    ##   qual3 ~~ qual4      0.4752      0.0644    7.3785    0.0000 [ 0.3330; 0.5906 ] 
    ##   qual3 ~~ qual5      0.5074      0.0617    8.2210    0.0000 [ 0.3895; 0.6231 ] 
    ##   qual4 ~~ qual5      0.6402      0.0559   11.4491    0.0000 [ 0.5110; 0.7369 ] 
    ##   val1 ~~ val2        0.6344      0.0553   11.4678    0.0000 [ 0.5184; 0.7310 ] 
    ##   val1 ~~ val3        0.4602      0.0710    6.4849    0.0000 [ 0.3241; 0.5978 ] 
    ##   val2 ~~ val3        0.6288      0.0626   10.0421    0.0000 [ 0.4920; 0.7419 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0613    7.6922    0.0000 [ 0.3391; 0.5903 ] 
    ##   QUAL ~ IMAG       0.3933      0.0566    6.9468    0.0000 [ 0.2795; 0.5032 ] 
    ##   QUAL ~ EXPE       0.8344      0.0223   37.3800    0.0000 [ 0.7881; 0.8753 ] 
    ##   VAL ~ IMAG        0.2974      0.0559    5.3174    0.0000 [ 0.1969; 0.4162 ] 
    ##   VAL ~ EXPE        0.6309      0.0479   13.1757    0.0000 [ 0.5289; 0.7152 ] 
    ##   VAL ~ QUAL        0.7013      0.0840    8.3475    0.0000 [ 0.5285; 0.8503 ] 
    ##   SAT ~ IMAG        0.4807      0.0631    7.6191    0.0000 [ 0.3523; 0.5963 ] 
    ##   SAT ~ EXPE        0.5001      0.0582    8.5879    0.0000 [ 0.3922; 0.6053 ] 
    ##   SAT ~ QUAL        0.5911      0.0993    5.9525    0.0000 [ 0.3961; 0.7878 ] 
    ##   SAT ~ VAL         0.5270      0.0901    5.8505    0.0000 [ 0.3501; 0.6904 ] 
    ##   LOY ~ IMAG        0.4840      0.0655    7.3841    0.0000 [ 0.3637; 0.6041 ] 
    ##   LOY ~ EXPE        0.3142      0.0546    5.7573    0.0000 [ 0.2162; 0.4344 ] 
    ##   LOY ~ QUAL        0.3714      0.0862    4.3085    0.0000 [ 0.2247; 0.5653 ] 
    ##   LOY ~ VAL         0.3311      0.0757    4.3717    0.0000 [ 0.2043; 0.4768 ] 
    ##   LOY ~ SAT         0.6283      0.0806    7.7916    0.0000 [ 0.4737; 0.7962 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0566    6.9468    0.0000 [ 0.2795; 0.5032 ] 
    ##   VAL ~ IMAG           0.2974      0.0559    5.3174    0.0000 [ 0.1969; 0.4162 ] 
    ##   VAL ~ EXPE           0.5852      0.0722    8.1051    0.0000 [ 0.4417; 0.7187 ] 
    ##   SAT ~ IMAG           0.2357      0.0464    5.0853    0.0000 [ 0.1541; 0.3343 ] 
    ##   SAT ~ EXPE           0.5173      0.0701    7.3833    0.0000 [ 0.3806; 0.6454 ] 
    ##   SAT ~ QUAL           0.3696      0.0648    5.7004    0.0000 [ 0.2340; 0.4930 ] 
    ##   LOY ~ IMAG           0.3020      0.0559    5.3991    0.0000 [ 0.2092; 0.4286 ] 
    ##   LOY ~ EXPE           0.3142      0.0546    5.7573    0.0000 [ 0.2162; 0.4344 ] 
    ##   LOY ~ QUAL           0.3714      0.0862    4.3085    0.0000 [ 0.2247; 0.5653 ] 
    ##   LOY ~ VAL            0.3311      0.0757    4.3717    0.0000 [ 0.2043; 0.4768 ] 
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
