
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![R-CMD-check](https://github.com/FloSchuberth/cSEM/workflows/R-CMD-check/badge.svg)](https://github.com/FloSchuberth/cSEM/actions)
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
approach), as well as several tests and typical post-estimation
procedures (e.g., verify admissibility of the estimates, assess the
model fit, test the model fit, compute confidence intervals, compare
groups, etc.).

## News (2026-07-27):

- Replace the polycor package by a more efficient implementation to
  calculate polychoric/polyserial correlations. Thanks to Kjell S.
  Slupphaug who contributed this implementation.

- Adjust p-value calculation in testMGD in case of permutation-based
  tests to prevent that p-values can be exactly 0. Thanks to Michael

- Fix bug in BasicCIResample(). Thanks to Michael.

- Implementation of doModelSearch() to perform AGAS-PLS. Thanks to
  Gloria.

- Release of cSEM version 0.6.1

- Release of cSEM Version 0.6.0

- Implementation of a `plot()` function to visualize cSEM models. Thanks
  to Nguyen.

- Enhancement of the `predict()` function

## Installation

The package is available on [CRAN](https://cran.r-project.org/):

``` r
install.packages("cSEM")
```

To install the development version, which is recommended, use:

``` r
# install.packages("pak")
pak::pak("FloSchuberth/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://floschuberth.github.io/cSEM/).

## Basic usage

The basic usage is illustrated below.

<img src="man/figures/api.png" alt="" width="80%" style="display: block; margin: auto;" />

Usually, using `cSEM` is the same 3 step procedure:

> 1.  Pick a dataset and specify a model using [lavaan
>     syntax](https://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use `csem()`
> 3.  Apply one of the post-estimation functions listed below on the
>     resulting object.

## Post-Estimation Functions

There are five major post-estimation verbs, three test family functions
and three do-family of function:

- `assess()` : assess the model using common quality criteria
- `infer()` : calculate common inferential quantities (e.g., standard
  errors, confidence intervals)
- `predict()` : predict endogenous indicator values
- `plot()` : Plot the cSEM model
- `summarize()` : summarize the results
- `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

- `testCVPAT()` performs a cross-validated predictive ability test
- `testOMF()` : performs a test for overall model fit
- `testMICOM()` : performs a test for composite measurement invariance
- `testMGD()` : performs several tests to assess multi-group differences
- `testHausman()` : performs the regression-based Hausman test to test
  for endogeneity

Other miscellaneous post-estimation functions belong do the do-family of
functions. Currently, three do functions are implemented:

- `doIPMA()`: performs an importance-performance matrix analysis
- `doNonlinearEffectsAnalysis()`: performs a nonlinear effects analysis
  such as floodlight and surface analysis
- `doRedundancyAnalysis()`: performs a redundancy analysis

All functions require a `cSEMResults` object.

## Example

Models are defined using [lavaan
syntax](https://lavaan.ugent.be/tutorial/syntax1.html) with some slight
modifications (see the [Specifying a
model](https://floschuberth.github.io/cSEM/articles/cSEM.html#using-csem)
section on the [cSEM-website](https://floschuberth.github.io/cSEM/)).
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
package](https://github.com/timelyportfolio/listviewer/). If you are new
to `cSEM` this might be a good way to familiarize yourself with the
structure of a `cSEMResults` object.

``` r
listviewer::jsonedit(res, mode = "view") # requires the listviewer package.
```

Apply post-estimation functions:

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
    ##  dG                      0.6493      0.3128  
    ##  SRMR                    0.0940      0.0529  
    ##  dL                      2.2340      0.7076  
    ##  dML                     2.9219      1.5456  
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
    ##  Out of 499 bootstrap replications 465 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -604905912
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
    ##  Approach to predict                = 'earliest'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   expe1         1.4579         1.6123       1.9042         2.1228       0.0540
    ##   expe2         1.4099         1.5005       1.9320         2.0326       0.2002
    ##   expe3         1.6230         1.7365       2.1158         2.2171       0.1288
    ##   qual1         1.4725         1.5475       1.9211         2.0582       0.1173
    ##   qual2         1.5772         1.5394       2.0343         2.0547       0.2210
    ##   qual3         1.7292         1.7371       2.2171         2.2784       0.1225
    ##   qual4         1.2370         1.1900       1.5982         1.6232       0.2344
    ##   qual5         1.5005         1.5079       1.9306         1.9488       0.1990
    ##   val1          1.4406         1.3675       1.8643         1.7727       0.2562
    ##   val2          1.2238         1.2168       1.6418         1.7189       0.1803
    ##   val3          1.4787         1.3955       1.9621         1.9403       0.1516
    ##   sat1          1.2445         1.2248       1.6431         1.6201       0.3449
    ##   sat2          1.2266         1.2000       1.6345         1.6272       0.3151
    ##   sat3          1.3393         1.2907       1.6699         1.7223       0.2148
    ##   sat4          1.3129         1.2634       1.6585         1.6353       0.2866
    ##   loy1          1.6856         1.6550       2.2337         2.2251       0.2710
    ##   loy2          1.4879         1.4961       1.9118         1.9904       0.1286
    ##   loy3          1.6974         1.6608       2.2776         2.2561       0.2722
    ##   loy4          1.6898         1.6768       2.1796         2.3013       0.0841
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential statistics are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways for resampling:

1.  Setting `.resample_method` in `csem()` to `"jackknife"` or
    `"bootstrap"` and subsequently using post-estimation functions
    `summarize()` or `infer()`.
2.  The same result is achieved by passing a `cSEMResults` object to
    `resamplecSEMResults()` and subsequently using post-estimation
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
    ##  Random seed                        = -1305138351
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
    ##   EXPE ~ IMAG      0.4714      0.0609    7.7370    0.0000 [ 0.3475; 0.5927 ] 
    ##   QUAL ~ EXPE      0.8344      0.0238   35.1019    0.0000 [ 0.7825; 0.8748 ] 
    ##   VAL ~ EXPE       0.0457      0.0859    0.5323    0.5945 [-0.1201; 0.2177 ] 
    ##   VAL ~ QUAL       0.7013      0.0819    8.5682    0.0000 [ 0.5446; 0.8653 ] 
    ##   SAT ~ IMAG       0.2450      0.0549    4.4632    0.0000 [ 0.1259; 0.3457 ] 
    ##   SAT ~ EXPE      -0.0172      0.0713   -0.2418    0.8089 [-0.1572; 0.1258 ] 
    ##   SAT ~ QUAL       0.2215      0.0977    2.2666    0.0234 [ 0.0612; 0.4331 ] 
    ##   SAT ~ VAL        0.5270      0.0841    6.2674    0.0000 [ 0.3491; 0.6836 ] 
    ##   LOY ~ IMAG       0.1819      0.0765    2.3786    0.0174 [ 0.0443; 0.3424 ] 
    ##   LOY ~ SAT        0.6283      0.0788    7.9717    0.0000 [ 0.4780; 0.7770 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1001    6.3025    0.0000 [ 0.4246; 0.8178 ] 
    ##   IMAG =~ imag2      0.9246      0.0371   24.9062    0.0000 [ 0.8303; 0.9829 ] 
    ##   IMAG =~ imag3      0.9577      0.0281   34.0697    0.0000 [ 0.8834; 0.9913 ] 
    ##   EXPE =~ expe1      0.7525      0.0763    9.8621    0.0000 [ 0.5754; 0.8696 ] 
    ##   EXPE =~ expe2      0.9348      0.0290   32.2431    0.0000 [ 0.8635; 0.9713 ] 
    ##   EXPE =~ expe3      0.7295      0.0695   10.4998    0.0000 [ 0.5780; 0.8368 ] 
    ##   QUAL =~ qual1      0.7861      0.0683   11.5181    0.0000 [ 0.6234; 0.8888 ] 
    ##   QUAL =~ qual2      0.9244      0.0236   39.1678    0.0000 [ 0.8636; 0.9561 ] 
    ##   QUAL =~ qual3      0.7560      0.0605   12.4997    0.0000 [ 0.6139; 0.8535 ] 
    ##   QUAL =~ qual4      0.7632      0.0545   13.9955    0.0000 [ 0.6451; 0.8571 ] 
    ##   QUAL =~ qual5      0.7834      0.0484   16.1937    0.0000 [ 0.6799; 0.8678 ] 
    ##   VAL =~ val1        0.9518      0.0225   42.2177    0.0000 [ 0.9041; 0.9868 ] 
    ##   VAL =~ val2        0.8056      0.0604   13.3400    0.0000 [ 0.6671; 0.9000 ] 
    ##   VAL =~ val3        0.6763      0.0698    9.6868    0.0000 [ 0.5342; 0.7943 ] 
    ##   SAT =~ sat1        0.9243      0.0220   41.9547    0.0000 [ 0.8793; 0.9604 ] 
    ##   SAT =~ sat2        0.8813      0.0309   28.5373    0.0000 [ 0.8173; 0.9362 ] 
    ##   SAT =~ sat3        0.7127      0.0522   13.6640    0.0000 [ 0.5987; 0.8061 ] 
    ##   SAT =~ sat4        0.7756      0.0487   15.9390    0.0000 [ 0.6755; 0.8655 ] 
    ##   LOY =~ loy1        0.9097      0.0492   18.4977    0.0000 [ 0.8004; 0.9865 ] 
    ##   LOY =~ loy2        0.5775      0.0812    7.1140    0.0000 [ 0.4054; 0.7217 ] 
    ##   LOY =~ loy3        0.9043      0.0425   21.2916    0.0000 [ 0.8125; 0.9741 ] 
    ##   LOY =~ loy4        0.4917      0.0859    5.7229    0.0000 [ 0.3349; 0.6567 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1128    0.1387    0.8897 [-0.1990; 0.2417 ] 
    ##   IMAG <~ imag2      0.4473      0.1425    3.1382    0.0017 [ 0.1841; 0.7451 ] 
    ##   IMAG <~ imag3      0.6020      0.1369    4.3980    0.0000 [ 0.2901; 0.8344 ] 
    ##   EXPE <~ expe1      0.2946      0.1153    2.5551    0.0106 [ 0.0605; 0.5055 ] 
    ##   EXPE <~ expe2      0.6473      0.0846    7.6481    0.0000 [ 0.4688; 0.7906 ] 
    ##   EXPE <~ expe3      0.2374      0.0895    2.6524    0.0080 [ 0.0371; 0.4061 ] 
    ##   QUAL <~ qual1      0.2370      0.0913    2.5964    0.0094 [ 0.0819; 0.4238 ] 
    ##   QUAL <~ qual2      0.4712      0.0786    5.9930    0.0000 [ 0.3057; 0.6199 ] 
    ##   QUAL <~ qual3      0.1831      0.0812    2.2538    0.0242 [ 0.0020; 0.3355 ] 
    ##   QUAL <~ qual4      0.1037      0.0618    1.6780    0.0934 [-0.0178; 0.2257 ] 
    ##   QUAL <~ qual5      0.2049      0.0626    3.2745    0.0011 [ 0.0662; 0.3119 ] 
    ##   VAL <~ val1        0.7163      0.0933    7.6746    0.0000 [ 0.5209; 0.8745 ] 
    ##   VAL <~ val2        0.2202      0.0890    2.4748    0.0133 [ 0.0593; 0.4064 ] 
    ##   VAL <~ val3        0.2082      0.0580    3.5903    0.0003 [ 0.1042; 0.3169 ] 
    ##   SAT <~ sat1        0.3209      0.0153   20.9953    0.0000 [ 0.2973; 0.3595 ] 
    ##   SAT <~ sat2        0.3059      0.0133   23.0571    0.0000 [ 0.2832; 0.3331 ] 
    ##   SAT <~ sat3        0.2474      0.0114   21.7955    0.0000 [ 0.2231; 0.2682 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.1405    0.0000 [ 0.2461; 0.2931 ] 
    ##   LOY <~ loy1        0.3834      0.0240   16.0044    0.0000 [ 0.3329; 0.4310 ] 
    ##   LOY <~ loy2        0.2434      0.0286    8.5045    0.0000 [ 0.1807; 0.2921 ] 
    ##   LOY <~ loy3        0.3812      0.0255   14.9544    0.0000 [ 0.3286; 0.4268 ] 
    ##   LOY <~ loy4        0.2073      0.0323    6.4093    0.0000 [ 0.1481; 0.2739 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0657    9.7985    0.0000 [ 0.5103; 0.7521 ] 
    ##   imag1 ~~ imag3      0.5433      0.0712    7.6300    0.0000 [ 0.4060; 0.6708 ] 
    ##   imag2 ~~ imag3      0.7761      0.0394   19.7051    0.0000 [ 0.6979; 0.8458 ] 
    ##   expe1 ~~ expe2      0.5353      0.0593    9.0192    0.0000 [ 0.4142; 0.6438 ] 
    ##   expe1 ~~ expe3      0.4694      0.0590    7.9533    0.0000 [ 0.3605; 0.5895 ] 
    ##   expe2 ~~ expe3      0.5467      0.0582    9.3883    0.0000 [ 0.4377; 0.6679 ] 
    ##   qual1 ~~ qual2      0.6053      0.0571   10.5932    0.0000 [ 0.4959; 0.7114 ] 
    ##   qual1 ~~ qual3      0.5406      0.0595    9.0890    0.0000 [ 0.4201; 0.6472 ] 
    ##   qual1 ~~ qual4      0.5662      0.0678    8.3557    0.0000 [ 0.4286; 0.6917 ] 
    ##   qual1 ~~ qual5      0.5180      0.0666    7.7809    0.0000 [ 0.3940; 0.6523 ] 
    ##   qual2 ~~ qual3      0.6187      0.0532   11.6179    0.0000 [ 0.5053; 0.7101 ] 
    ##   qual2 ~~ qual4      0.6517      0.0620   10.5190    0.0000 [ 0.5292; 0.7601 ] 
    ##   qual2 ~~ qual5      0.6291      0.0577   10.9017    0.0000 [ 0.5084; 0.7349 ] 
    ##   qual3 ~~ qual4      0.4752      0.0647    7.3498    0.0000 [ 0.3463; 0.5902 ] 
    ##   qual3 ~~ qual5      0.5074      0.0594    8.5454    0.0000 [ 0.3809; 0.6124 ] 
    ##   qual4 ~~ qual5      0.6402      0.0554   11.5590    0.0000 [ 0.5171; 0.7374 ] 
    ##   val1 ~~ val2        0.6344      0.0519   12.2277    0.0000 [ 0.5331; 0.7269 ] 
    ##   val1 ~~ val3        0.4602      0.0664    6.9264    0.0000 [ 0.3267; 0.5929 ] 
    ##   val2 ~~ val3        0.6288      0.0627   10.0335    0.0000 [ 0.5094; 0.7392 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0609    7.7370    0.0000 [ 0.3475; 0.5927 ] 
    ##   QUAL ~ IMAG       0.3933      0.0570    6.8965    0.0000 [ 0.2817; 0.5099 ] 
    ##   QUAL ~ EXPE       0.8344      0.0238   35.1019    0.0000 [ 0.7825; 0.8748 ] 
    ##   VAL ~ IMAG        0.2974      0.0570    5.2199    0.0000 [ 0.1924; 0.4230 ] 
    ##   VAL ~ EXPE        0.6309      0.0491   12.8562    0.0000 [ 0.5381; 0.7281 ] 
    ##   VAL ~ QUAL        0.7013      0.0819    8.5682    0.0000 [ 0.5446; 0.8653 ] 
    ##   SAT ~ IMAG        0.4807      0.0620    7.7534    0.0000 [ 0.3594; 0.5974 ] 
    ##   SAT ~ EXPE        0.5001      0.0559    8.9440    0.0000 [ 0.3922; 0.6121 ] 
    ##   SAT ~ QUAL        0.5911      0.0926    6.3842    0.0000 [ 0.4105; 0.7736 ] 
    ##   SAT ~ VAL         0.5270      0.0841    6.2674    0.0000 [ 0.3491; 0.6836 ] 
    ##   LOY ~ IMAG        0.4840      0.0622    7.7770    0.0000 [ 0.3801; 0.6118 ] 
    ##   LOY ~ EXPE        0.3142      0.0540    5.8148    0.0000 [ 0.2172; 0.4214 ] 
    ##   LOY ~ QUAL        0.3714      0.0822    4.5184    0.0000 [ 0.2333; 0.5540 ] 
    ##   LOY ~ VAL         0.3311      0.0722    4.5865    0.0000 [ 0.1961; 0.4841 ] 
    ##   LOY ~ SAT         0.6283      0.0788    7.9717    0.0000 [ 0.4780; 0.7770 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0570    6.8965    0.0000 [ 0.2817; 0.5099 ] 
    ##   VAL ~ IMAG           0.2974      0.0570    5.2199    0.0000 [ 0.1924; 0.4230 ] 
    ##   VAL ~ EXPE           0.5852      0.0702    8.3385    0.0000 [ 0.4534; 0.7171 ] 
    ##   SAT ~ IMAG           0.2357      0.0461    5.1140    0.0000 [ 0.1521; 0.3306 ] 
    ##   SAT ~ EXPE           0.5173      0.0674    7.6786    0.0000 [ 0.4068; 0.6689 ] 
    ##   SAT ~ QUAL           0.3696      0.0616    6.0025    0.0000 [ 0.2486; 0.5060 ] 
    ##   LOY ~ IMAG           0.3020      0.0533    5.6642    0.0000 [ 0.2096; 0.4096 ] 
    ##   LOY ~ EXPE           0.3142      0.0540    5.8148    0.0000 [ 0.2172; 0.4214 ] 
    ##   LOY ~ QUAL           0.3714      0.0822    4.5184    0.0000 [ 0.2333; 0.5540 ] 
    ##   LOY ~ VAL            0.3311      0.0722    4.5865    0.0000 [ 0.1961; 0.4841 ] 
    ## ________________________________________________________________________________

Several bootstrap-based confidence intervals are implemented, see
`?infer()`:

``` r
infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet
```

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as setting random seeds via the [future
framework](https://github.com/futureverse/future/). For multiprocessing
simply set `.eval_plan = "multisession"` in which case the maximum
number of available cores is used if not on Windows. On Windows as many
separate R instances are opened in the background as there are cores
available instead. Note that this naturally has some overhead so for a
small number of resamples multiprocessing will not always be faster
compared to sequential (single core) processing (the default). Seeds are
set via the `.seed` argument.

``` r
b <- csem(
  .data            = satisfaction,
  .model           = model, 
  .resample_method = "bootstrap",
  .R               = 999,
  .seed            = 98234,
  .eval_plan       = "multisession")
```
