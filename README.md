
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

- Fix smaller issue in the print function for `assess()`.

- Replace helper function from the matrixcalc package. Thanks to
  Kjell S. Slupphaug for this contribution.

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
    ##  dG                      0.6493      0.3263  
    ##  SRMR                    0.0940      0.0544  
    ##  dL                      2.2340      0.7501  
    ##  dML                     2.9219      1.5996  
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
    ##  Out of 499 bootstrap replications 480 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 1579570724
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
    ##   expe1         1.4561         1.5963       1.9053         2.1011       0.0579
    ##   expe2         1.4130         1.5075       1.9344         2.0489       0.1983
    ##   expe3         1.6384         1.7620       2.1331         2.2508       0.1231
    ##   qual1         1.4758         1.5659       1.9243         2.0798       0.1168
    ##   qual2         1.5754         1.5568       2.0369         2.0879       0.2164
    ##   qual3         1.7360         1.7434       2.2226         2.2934       0.1171
    ##   qual4         1.2321         1.2024       1.5907         1.6457       0.2349
    ##   qual5         1.5024         1.5159       1.9333         1.9559       0.1971
    ##   val1          1.4479         1.3769       1.8757         1.7851       0.2454
    ##   val2          1.2234         1.2179       1.6446         1.7183       0.1750
    ##   val3          1.4771         1.3897       1.9636         1.9401       0.1481
    ##   sat1          1.2454         1.2340       1.6448         1.6164       0.3383
    ##   sat2          1.2291         1.2008       1.6363         1.6207       0.3110
    ##   sat3          1.3437         1.2928       1.6762         1.7209       0.2112
    ##   sat4          1.3195         1.2598       1.6720         1.6320       0.2707
    ##   loy1          1.6952         1.6706       2.2351         2.2344       0.2646
    ##   loy2          1.4788         1.4857       1.9045         1.9803       0.1356
    ##   loy3          1.6998         1.6654       2.2744         2.2633       0.2713
    ##   loy4          1.6881         1.6963       2.1797         2.3275       0.0846
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
    ##  Number of admissible results       = 491
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 1748902371
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
    ##   EXPE ~ IMAG      0.4714      0.0669    7.0466    0.0000 [ 0.3387; 0.6042 ] 
    ##   QUAL ~ EXPE      0.8344      0.0239   34.9215    0.0000 [ 0.7840; 0.8765 ] 
    ##   VAL ~ EXPE       0.0457      0.0875    0.5225    0.6013 [-0.1177; 0.2313 ] 
    ##   VAL ~ QUAL       0.7013      0.0836    8.3887    0.0000 [ 0.5288; 0.8514 ] 
    ##   SAT ~ IMAG       0.2450      0.0586    4.1780    0.0000 [ 0.1383; 0.3616 ] 
    ##   SAT ~ EXPE      -0.0172      0.0749   -0.2300    0.8181 [-0.1718; 0.1234 ] 
    ##   SAT ~ QUAL       0.2215      0.1061    2.0881    0.0368 [ 0.0288; 0.4497 ] 
    ##   SAT ~ VAL        0.5270      0.0876    6.0170    0.0000 [ 0.3472; 0.6681 ] 
    ##   LOY ~ IMAG       0.1819      0.0823    2.2114    0.0270 [ 0.0177; 0.3345 ] 
    ##   LOY ~ SAT        0.6283      0.0837    7.5063    0.0000 [ 0.4694; 0.7843 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1032    6.1098    0.0000 [ 0.4233; 0.8064 ] 
    ##   IMAG =~ imag2      0.9246      0.0380   24.3628    0.0000 [ 0.8312; 0.9763 ] 
    ##   IMAG =~ imag3      0.9577      0.0294   32.6041    0.0000 [ 0.8787; 0.9919 ] 
    ##   EXPE =~ expe1      0.7525      0.0751   10.0212    0.0000 [ 0.5720; 0.8665 ] 
    ##   EXPE =~ expe2      0.9348      0.0282   33.1731    0.0000 [ 0.8605; 0.9731 ] 
    ##   EXPE =~ expe3      0.7295      0.0678   10.7624    0.0000 [ 0.5748; 0.8369 ] 
    ##   QUAL =~ qual1      0.7861      0.0663   11.8494    0.0000 [ 0.6335; 0.8918 ] 
    ##   QUAL =~ qual2      0.9244      0.0226   40.9167    0.0000 [ 0.8720; 0.9584 ] 
    ##   QUAL =~ qual3      0.7560      0.0611   12.3782    0.0000 [ 0.6188; 0.8516 ] 
    ##   QUAL =~ qual4      0.7632      0.0504   15.1487    0.0000 [ 0.6550; 0.8510 ] 
    ##   QUAL =~ qual5      0.7834      0.0447   17.5265    0.0000 [ 0.6981; 0.8642 ] 
    ##   VAL =~ val1        0.9518      0.0225   42.3646    0.0000 [ 0.8915; 0.9810 ] 
    ##   VAL =~ val2        0.8056      0.0637   12.6538    0.0000 [ 0.6643; 0.9113 ] 
    ##   VAL =~ val3        0.6763      0.0711    9.5102    0.0000 [ 0.5309; 0.8003 ] 
    ##   SAT =~ sat1        0.9243      0.0222   41.7135    0.0000 [ 0.8775; 0.9595 ] 
    ##   SAT =~ sat2        0.8813      0.0283   31.1366    0.0000 [ 0.8204; 0.9309 ] 
    ##   SAT =~ sat3        0.7127      0.0562   12.6872    0.0000 [ 0.5880; 0.8066 ] 
    ##   SAT =~ sat4        0.7756      0.0495   15.6614    0.0000 [ 0.6733; 0.8596 ] 
    ##   LOY =~ loy1        0.9097      0.0507   17.9439    0.0000 [ 0.8002; 0.9858 ] 
    ##   LOY =~ loy2        0.5775      0.0852    6.7822    0.0000 [ 0.3826; 0.7394 ] 
    ##   LOY =~ loy3        0.9043      0.0418   21.6458    0.0000 [ 0.8187; 0.9747 ] 
    ##   LOY =~ loy4        0.4917      0.0986    4.9854    0.0000 [ 0.3008; 0.6676 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1225    0.1278    0.8983 [-0.2093; 0.2499 ] 
    ##   IMAG <~ imag2      0.4473      0.1440    3.1055    0.0019 [ 0.1768; 0.7104 ] 
    ##   IMAG <~ imag3      0.6020      0.1387    4.3396    0.0000 [ 0.3179; 0.8563 ] 
    ##   EXPE <~ expe1      0.2946      0.1130    2.6071    0.0091 [ 0.0739; 0.5029 ] 
    ##   EXPE <~ expe2      0.6473      0.0837    7.7309    0.0000 [ 0.4545; 0.7967 ] 
    ##   EXPE <~ expe3      0.2374      0.0865    2.7432    0.0061 [ 0.0569; 0.3921 ] 
    ##   QUAL <~ qual1      0.2370      0.0879    2.6957    0.0070 [ 0.0731; 0.4197 ] 
    ##   QUAL <~ qual2      0.4712      0.0788    5.9836    0.0000 [ 0.3122; 0.6149 ] 
    ##   QUAL <~ qual3      0.1831      0.0748    2.4470    0.0144 [ 0.0325; 0.3146 ] 
    ##   QUAL <~ qual4      0.1037      0.0607    1.7098    0.0873 [-0.0125; 0.2245 ] 
    ##   QUAL <~ qual5      0.2049      0.0624    3.2843    0.0010 [ 0.0821; 0.3242 ] 
    ##   VAL <~ val1        0.7163      0.0912    7.8539    0.0000 [ 0.5131; 0.8568 ] 
    ##   VAL <~ val2        0.2202      0.0936    2.3532    0.0186 [ 0.0602; 0.4437 ] 
    ##   VAL <~ val3        0.2082      0.0589    3.5327    0.0004 [ 0.0857; 0.3124 ] 
    ##   SAT <~ sat1        0.3209      0.0158   20.2791    0.0000 [ 0.2938; 0.3563 ] 
    ##   SAT <~ sat2        0.3059      0.0141   21.7414    0.0000 [ 0.2813; 0.3351 ] 
    ##   SAT <~ sat3        0.2474      0.0114   21.7585    0.0000 [ 0.2238; 0.2672 ] 
    ##   SAT <~ sat4        0.2692      0.0116   23.2225    0.0000 [ 0.2471; 0.2911 ] 
    ##   LOY <~ loy1        0.3834      0.0263   14.5627    0.0000 [ 0.3289; 0.4345 ] 
    ##   LOY <~ loy2        0.2434      0.0302    8.0644    0.0000 [ 0.1746; 0.2974 ] 
    ##   LOY <~ loy3        0.3812      0.0269   14.1883    0.0000 [ 0.3311; 0.4322 ] 
    ##   LOY <~ loy4        0.2073      0.0368    5.6278    0.0000 [ 0.1287; 0.2751 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0696    9.2443    0.0000 [ 0.5023; 0.7633 ] 
    ##   imag1 ~~ imag3      0.5433      0.0705    7.7102    0.0000 [ 0.4003; 0.6768 ] 
    ##   imag2 ~~ imag3      0.7761      0.0380   20.4236    0.0000 [ 0.6881; 0.8426 ] 
    ##   expe1 ~~ expe2      0.5353      0.0582    9.1925    0.0000 [ 0.4169; 0.6416 ] 
    ##   expe1 ~~ expe3      0.4694      0.0591    7.9483    0.0000 [ 0.3478; 0.5767 ] 
    ##   expe2 ~~ expe3      0.5467      0.0573    9.5364    0.0000 [ 0.4360; 0.6627 ] 
    ##   qual1 ~~ qual2      0.6053      0.0587   10.3071    0.0000 [ 0.4864; 0.7202 ] 
    ##   qual1 ~~ qual3      0.5406      0.0628    8.6043    0.0000 [ 0.4184; 0.6654 ] 
    ##   qual1 ~~ qual4      0.5662      0.0658    8.6047    0.0000 [ 0.4355; 0.6942 ] 
    ##   qual1 ~~ qual5      0.5180      0.0647    8.0017    0.0000 [ 0.3801; 0.6407 ] 
    ##   qual2 ~~ qual3      0.6187      0.0579   10.6904    0.0000 [ 0.4969; 0.7160 ] 
    ##   qual2 ~~ qual4      0.6517      0.0608   10.7251    0.0000 [ 0.5246; 0.7599 ] 
    ##   qual2 ~~ qual5      0.6291      0.0562   11.2026    0.0000 [ 0.5216; 0.7357 ] 
    ##   qual3 ~~ qual4      0.4752      0.0645    7.3715    0.0000 [ 0.3452; 0.6051 ] 
    ##   qual3 ~~ qual5      0.5074      0.0632    8.0296    0.0000 [ 0.3812; 0.6227 ] 
    ##   qual4 ~~ qual5      0.6402      0.0553   11.5702    0.0000 [ 0.5281; 0.7380 ] 
    ##   val1 ~~ val2        0.6344      0.0570   11.1234    0.0000 [ 0.5243; 0.7373 ] 
    ##   val1 ~~ val3        0.4602      0.0711    6.4733    0.0000 [ 0.3188; 0.5954 ] 
    ##   val2 ~~ val3        0.6288      0.0676    9.3088    0.0000 [ 0.4967; 0.7697 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0669    7.0466    0.0000 [ 0.3387; 0.6042 ] 
    ##   QUAL ~ IMAG       0.3933      0.0627    6.2769    0.0000 [ 0.2767; 0.5196 ] 
    ##   QUAL ~ EXPE       0.8344      0.0239   34.9215    0.0000 [ 0.7840; 0.8765 ] 
    ##   VAL ~ IMAG        0.2974      0.0612    4.8561    0.0000 [ 0.1829; 0.4284 ] 
    ##   VAL ~ EXPE        0.6309      0.0512   12.3308    0.0000 [ 0.5271; 0.7217 ] 
    ##   VAL ~ QUAL        0.7013      0.0836    8.3887    0.0000 [ 0.5288; 0.8514 ] 
    ##   SAT ~ IMAG        0.4807      0.0687    6.9951    0.0000 [ 0.3522; 0.6189 ] 
    ##   SAT ~ EXPE        0.5001      0.0586    8.5291    0.0000 [ 0.3854; 0.6142 ] 
    ##   SAT ~ QUAL        0.5911      0.0968    6.1052    0.0000 [ 0.3993; 0.7571 ] 
    ##   SAT ~ VAL         0.5270      0.0876    6.0170    0.0000 [ 0.3472; 0.6681 ] 
    ##   LOY ~ IMAG        0.4840      0.0678    7.1346    0.0000 [ 0.3486; 0.6163 ] 
    ##   LOY ~ EXPE        0.3142      0.0554    5.6760    0.0000 [ 0.2089; 0.4346 ] 
    ##   LOY ~ QUAL        0.3714      0.0832    4.4633    0.0000 [ 0.2242; 0.5467 ] 
    ##   LOY ~ VAL         0.3311      0.0748    4.4237    0.0000 [ 0.1908; 0.4784 ] 
    ##   LOY ~ SAT         0.6283      0.0837    7.5063    0.0000 [ 0.4694; 0.7843 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0627    6.2769    0.0000 [ 0.2767; 0.5196 ] 
    ##   VAL ~ IMAG           0.2974      0.0612    4.8561    0.0000 [ 0.1829; 0.4284 ] 
    ##   VAL ~ EXPE           0.5852      0.0723    8.0913    0.0000 [ 0.4347; 0.7261 ] 
    ##   SAT ~ IMAG           0.2357      0.0478    4.9303    0.0000 [ 0.1480; 0.3356 ] 
    ##   SAT ~ EXPE           0.5173      0.0679    7.6230    0.0000 [ 0.3827; 0.6485 ] 
    ##   SAT ~ QUAL           0.3696      0.0617    5.9910    0.0000 [ 0.2506; 0.4849 ] 
    ##   LOY ~ IMAG           0.3020      0.0606    4.9801    0.0000 [ 0.2024; 0.4383 ] 
    ##   LOY ~ EXPE           0.3142      0.0554    5.6760    0.0000 [ 0.2089; 0.4346 ] 
    ##   LOY ~ QUAL           0.3714      0.0832    4.4633    0.0000 [ 0.2242; 0.5467 ] 
    ##   LOY ~ VAL            0.3311      0.0748    4.4237    0.0000 [ 0.1908; 0.4784 ] 
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
