
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

## News (2026-08-24):

- Extend calculateHTMT() to allow for asymptotic inference for the HTMT.
  Thanks to Jason Berger for his contribution.

- Fix bug in calculating the moments used to estimate non-linear models.

- Fix bug in testMICOM(). Thanks to manzi0 for this contribution!

- Fix issue with setting a seed in non-interactive session. Thanks to
  Kjell S. Slupphaug and Jason Berger for their contribution.

- Fix bug in PLS-PM estimation of non-linear models involing
  second-order constructs. Thanks Thanks to Kjell S. Slupphaug for this
  contribution.

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
    ##  dG                      0.6493      0.3261  
    ##  SRMR                    0.0940      0.0516  
    ##  dL                      2.2340      0.6737  
    ##  dML                     2.9219      1.6373  
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
    ##  Out of 499 bootstrap replications 469 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -1813629629
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
    ## 
    ## Variance accounted for (VAF):
    ## =============================
    ##   Effects        Estimate  Std. error   t-stat.   p-value
    ##   QUAL ~ IMAG      1.0000          NA        NA        NA
    ##   VAL ~ IMAG       1.0000          NA        NA        NA
    ##   VAL ~ EXPE       0.9275          NA        NA        NA
    ##   SAT ~ IMAG       0.4904          NA        NA        NA
    ##   SAT ~ EXPE       1.0345          NA        NA        NA
    ##   SAT ~ QUAL       0.6252          NA        NA        NA
    ##   LOY ~ IMAG       0.6241          NA        NA        NA
    ##   LOY ~ EXPE       1.0000          NA        NA        NA
    ##   LOY ~ QUAL       1.0000          NA        NA        NA
    ##   LOY ~ VAL        1.0000          NA        NA        NA
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
    ##   expe1         1.4523         1.5925       1.9077         2.0990       0.0546
    ##   expe2         1.4087         1.5073       1.9303         2.0466       0.2011
    ##   expe3         1.6258         1.7357       2.1181         2.2180       0.1273
    ##   qual1         1.4764         1.5623       1.9327         2.0718       0.1114
    ##   qual2         1.5780         1.5502       2.0363         2.0674       0.2197
    ##   qual3         1.7341         1.7403       2.2209         2.2788       0.1212
    ##   qual4         1.2325         1.1898       1.5937         1.6225       0.2355
    ##   qual5         1.5023         1.5073       1.9273         1.9540       0.2001
    ##   val1          1.4453         1.3713       1.8681         1.7705       0.2505
    ##   val2          1.2242         1.2241       1.6441         1.7223       0.1752
    ##   val3          1.4815         1.3939       1.9724         1.9436       0.1479
    ##   sat1          1.2492         1.2313       1.6457         1.6209       0.3385
    ##   sat2          1.2354         1.2047       1.6367         1.6316       0.3104
    ##   sat3          1.3415         1.2959       1.6746         1.7265       0.2066
    ##   sat4          1.3213         1.2727       1.6746         1.6473       0.2708
    ##   loy1          1.6980         1.6627       2.2338         2.2242       0.2690
    ##   loy2          1.4861         1.4796       1.9121         1.9805       0.1316
    ##   loy3          1.7037         1.6757       2.2755         2.2685       0.2709
    ##   loy4          1.6922         1.6791       2.1765         2.2961       0.0921
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
    ##  Number of admissible results       = 481
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 912017757
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
    ##   EXPE ~ IMAG      0.4714      0.0648    7.2693    0.0000 [ 0.3464; 0.5971 ] 
    ##   QUAL ~ EXPE      0.8344      0.0237   35.1602    0.0000 [ 0.7818; 0.8739 ] 
    ##   VAL ~ EXPE       0.0457      0.0850    0.5378    0.5907 [-0.1243; 0.2264 ] 
    ##   VAL ~ QUAL       0.7013      0.0831    8.4353    0.0000 [ 0.5314; 0.8517 ] 
    ##   SAT ~ IMAG       0.2450      0.0562    4.3576    0.0000 [ 0.1369; 0.3548 ] 
    ##   SAT ~ EXPE      -0.0172      0.0732   -0.2354    0.8139 [-0.1671; 0.1154 ] 
    ##   SAT ~ QUAL       0.2215      0.1013    2.1860    0.0288 [ 0.0533; 0.4428 ] 
    ##   SAT ~ VAL        0.5270      0.0895    5.8905    0.0000 [ 0.3510; 0.6986 ] 
    ##   LOY ~ IMAG       0.1819      0.0748    2.4312    0.0151 [ 0.0169; 0.3188 ] 
    ##   LOY ~ SAT        0.6283      0.0776    8.0992    0.0000 [ 0.4882; 0.8051 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1027    6.1384    0.0000 [ 0.3981; 0.8092 ] 
    ##   IMAG =~ imag2      0.9246      0.0409   22.6120    0.0000 [ 0.8098; 0.9769 ] 
    ##   IMAG =~ imag3      0.9577      0.0297   32.2019    0.0000 [ 0.8818; 0.9931 ] 
    ##   EXPE =~ expe1      0.7525      0.0792    9.5063    0.0000 [ 0.5728; 0.8694 ] 
    ##   EXPE =~ expe2      0.9348      0.0273   34.3038    0.0000 [ 0.8674; 0.9698 ] 
    ##   EXPE =~ expe3      0.7295      0.0733    9.9483    0.0000 [ 0.5705; 0.8581 ] 
    ##   QUAL =~ qual1      0.7861      0.0711   11.0589    0.0000 [ 0.6264; 0.8866 ] 
    ##   QUAL =~ qual2      0.9244      0.0233   39.6958    0.0000 [ 0.8724; 0.9596 ] 
    ##   QUAL =~ qual3      0.7560      0.0620   12.1938    0.0000 [ 0.6046; 0.8590 ] 
    ##   QUAL =~ qual4      0.7632      0.0543   14.0584    0.0000 [ 0.6435; 0.8579 ] 
    ##   QUAL =~ qual5      0.7834      0.0492   15.9210    0.0000 [ 0.6720; 0.8566 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.5805    0.0000 [ 0.8964; 0.9860 ] 
    ##   VAL =~ val2        0.8056      0.0645   12.4920    0.0000 [ 0.6475; 0.9092 ] 
    ##   VAL =~ val3        0.6763      0.0749    9.0277    0.0000 [ 0.5209; 0.8022 ] 
    ##   SAT =~ sat1        0.9243      0.0220   41.9813    0.0000 [ 0.8732; 0.9632 ] 
    ##   SAT =~ sat2        0.8813      0.0291   30.2887    0.0000 [ 0.8150; 0.9261 ] 
    ##   SAT =~ sat3        0.7127      0.0501   14.2306    0.0000 [ 0.6093; 0.8115 ] 
    ##   SAT =~ sat4        0.7756      0.0503   15.4171    0.0000 [ 0.6714; 0.8667 ] 
    ##   LOY =~ loy1        0.9097      0.0492   18.5037    0.0000 [ 0.8021; 0.9863 ] 
    ##   LOY =~ loy2        0.5775      0.0820    7.0422    0.0000 [ 0.4194; 0.7288 ] 
    ##   LOY =~ loy3        0.9043      0.0430   21.0096    0.0000 [ 0.8081; 0.9774 ] 
    ##   LOY =~ loy4        0.4917      0.0956    5.1436    0.0000 [ 0.3080; 0.6708 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1188    0.1316    0.8953 [-0.2179; 0.2462 ] 
    ##   IMAG <~ imag2      0.4473      0.1570    2.8490    0.0044 [ 0.1242; 0.7339 ] 
    ##   IMAG <~ imag3      0.6020      0.1447    4.1620    0.0000 [ 0.2960; 0.8667 ] 
    ##   EXPE <~ expe1      0.2946      0.1157    2.5468    0.0109 [ 0.0544; 0.5153 ] 
    ##   EXPE <~ expe2      0.6473      0.0806    8.0327    0.0000 [ 0.4700; 0.7829 ] 
    ##   EXPE <~ expe3      0.2374      0.0945    2.5112    0.0120 [ 0.0616; 0.4222 ] 
    ##   QUAL <~ qual1      0.2370      0.0911    2.6018    0.0093 [ 0.0694; 0.4150 ] 
    ##   QUAL <~ qual2      0.4712      0.0797    5.9119    0.0000 [ 0.3185; 0.6099 ] 
    ##   QUAL <~ qual3      0.1831      0.0817    2.2396    0.0251 [ 0.0136; 0.3299 ] 
    ##   QUAL <~ qual4      0.1037      0.0621    1.6710    0.0947 [-0.0187; 0.2329 ] 
    ##   QUAL <~ qual5      0.2049      0.0631    3.2472    0.0012 [ 0.0704; 0.3163 ] 
    ##   VAL <~ val1        0.7163      0.0976    7.3394    0.0000 [ 0.5259; 0.8865 ] 
    ##   VAL <~ val2        0.2202      0.0935    2.3544    0.0186 [ 0.0422; 0.4166 ] 
    ##   VAL <~ val3        0.2082      0.0607    3.4266    0.0006 [ 0.0954; 0.3159 ] 
    ##   SAT <~ sat1        0.3209      0.0149   21.6049    0.0000 [ 0.2950; 0.3527 ] 
    ##   SAT <~ sat2        0.3059      0.0134   22.8650    0.0000 [ 0.2821; 0.3337 ] 
    ##   SAT <~ sat3        0.2474      0.0102   24.2334    0.0000 [ 0.2274; 0.2669 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.0126    0.0000 [ 0.2455; 0.2945 ] 
    ##   LOY <~ loy1        0.3834      0.0253   15.1296    0.0000 [ 0.3340; 0.4329 ] 
    ##   LOY <~ loy2        0.2434      0.0289    8.4171    0.0000 [ 0.1845; 0.2946 ] 
    ##   LOY <~ loy3        0.3812      0.0265   14.3721    0.0000 [ 0.3290; 0.4350 ] 
    ##   LOY <~ loy4        0.2073      0.0360    5.7504    0.0000 [ 0.1364; 0.2750 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0659    9.7707    0.0000 [ 0.5030; 0.7565 ] 
    ##   imag1 ~~ imag3      0.5433      0.0712    7.6323    0.0000 [ 0.3953; 0.6622 ] 
    ##   imag2 ~~ imag3      0.7761      0.0408   19.0282    0.0000 [ 0.6933; 0.8492 ] 
    ##   expe1 ~~ expe2      0.5353      0.0601    8.9049    0.0000 [ 0.4080; 0.6513 ] 
    ##   expe1 ~~ expe3      0.4694      0.0644    7.2843    0.0000 [ 0.3230; 0.5896 ] 
    ##   expe2 ~~ expe3      0.5467      0.0622    8.7937    0.0000 [ 0.4111; 0.6492 ] 
    ##   qual1 ~~ qual2      0.6053      0.0587   10.3038    0.0000 [ 0.4803; 0.7004 ] 
    ##   qual1 ~~ qual3      0.5406      0.0611    8.8426    0.0000 [ 0.4156; 0.6609 ] 
    ##   qual1 ~~ qual4      0.5662      0.0723    7.8276    0.0000 [ 0.4052; 0.6880 ] 
    ##   qual1 ~~ qual5      0.5180      0.0703    7.3706    0.0000 [ 0.3770; 0.6458 ] 
    ##   qual2 ~~ qual3      0.6187      0.0585   10.5809    0.0000 [ 0.4957; 0.7211 ] 
    ##   qual2 ~~ qual4      0.6517      0.0637   10.2281    0.0000 [ 0.5048; 0.7603 ] 
    ##   qual2 ~~ qual5      0.6291      0.0613   10.2539    0.0000 [ 0.5055; 0.7425 ] 
    ##   qual3 ~~ qual4      0.4752      0.0674    7.0533    0.0000 [ 0.3399; 0.6013 ] 
    ##   qual3 ~~ qual5      0.5074      0.0632    8.0295    0.0000 [ 0.3864; 0.6275 ] 
    ##   qual4 ~~ qual5      0.6402      0.0562   11.3896    0.0000 [ 0.5186; 0.7429 ] 
    ##   val1 ~~ val2        0.6344      0.0527   12.0416    0.0000 [ 0.5234; 0.7300 ] 
    ##   val1 ~~ val3        0.4602      0.0725    6.3502    0.0000 [ 0.3144; 0.5928 ] 
    ##   val2 ~~ val3        0.6288      0.0640    9.8247    0.0000 [ 0.5033; 0.7620 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0648    7.2693    0.0000 [ 0.3464; 0.5971 ] 
    ##   QUAL ~ IMAG       0.3933      0.0603    6.5186    0.0000 [ 0.2772; 0.5166 ] 
    ##   QUAL ~ EXPE       0.8344      0.0237   35.1602    0.0000 [ 0.7818; 0.8739 ] 
    ##   VAL ~ IMAG        0.2974      0.0598    4.9772    0.0000 [ 0.1848; 0.4222 ] 
    ##   VAL ~ EXPE        0.6309      0.0504   12.5065    0.0000 [ 0.5243; 0.7182 ] 
    ##   VAL ~ QUAL        0.7013      0.0831    8.4353    0.0000 [ 0.5314; 0.8517 ] 
    ##   SAT ~ IMAG        0.4807      0.0674    7.1337    0.0000 [ 0.3446; 0.6081 ] 
    ##   SAT ~ EXPE        0.5001      0.0581    8.6040    0.0000 [ 0.3800; 0.6067 ] 
    ##   SAT ~ QUAL        0.5911      0.0969    6.0995    0.0000 [ 0.4189; 0.7841 ] 
    ##   SAT ~ VAL         0.5270      0.0895    5.8905    0.0000 [ 0.3510; 0.6986 ] 
    ##   LOY ~ IMAG        0.4840      0.0676    7.1618    0.0000 [ 0.3517; 0.6113 ] 
    ##   LOY ~ EXPE        0.3142      0.0528    5.9505    0.0000 [ 0.2169; 0.4317 ] 
    ##   LOY ~ QUAL        0.3714      0.0869    4.2735    0.0000 [ 0.2319; 0.5624 ] 
    ##   LOY ~ VAL         0.3311      0.0751    4.4096    0.0000 [ 0.1954; 0.4895 ] 
    ##   LOY ~ SAT         0.6283      0.0776    8.0992    0.0000 [ 0.4882; 0.8051 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0603    6.5186    0.0000 [ 0.2772; 0.5166 ] 
    ##   VAL ~ IMAG           0.2974      0.0598    4.9772    0.0000 [ 0.1848; 0.4222 ] 
    ##   VAL ~ EXPE           0.5852      0.0721    8.1192    0.0000 [ 0.4395; 0.7184 ] 
    ##   SAT ~ IMAG           0.2357      0.0490    4.8123    0.0000 [ 0.1433; 0.3359 ] 
    ##   SAT ~ EXPE           0.5173      0.0685    7.5553    0.0000 [ 0.3986; 0.6612 ] 
    ##   SAT ~ QUAL           0.3696      0.0663    5.5726    0.0000 [ 0.2313; 0.5045 ] 
    ##   LOY ~ IMAG           0.3020      0.0545    5.5436    0.0000 [ 0.2064; 0.4196 ] 
    ##   LOY ~ EXPE           0.3142      0.0528    5.9505    0.0000 [ 0.2169; 0.4317 ] 
    ##   LOY ~ QUAL           0.3714      0.0869    4.2735    0.0000 [ 0.2319; 0.5624 ] 
    ##   LOY ~ VAL            0.3311      0.0751    4.4096    0.0000 [ 0.1954; 0.4895 ] 
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
