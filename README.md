
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

## News (2026-08-02):

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
    ##  dG                      0.6493      0.3103  
    ##  SRMR                    0.0940      0.0517  
    ##  dL                      2.2340      0.6756  
    ##  dML                     2.9219      1.5606  
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
    ##  Out of 499 bootstrap replications 466 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 1433336226
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
    ##   expe1         1.4551         1.5902       1.9036         2.0943       0.0552
    ##   expe2         1.4144         1.4995       1.9343         2.0287       0.2027
    ##   expe3         1.6290         1.7379       2.1227         2.2130       0.1279
    ##   qual1         1.4715         1.5536       1.9196         2.0531       0.1218
    ##   qual2         1.5773         1.5589       2.0445         2.0725       0.2151
    ##   qual3         1.7432         1.7588       2.2336         2.2907       0.1193
    ##   qual4         1.2382         1.1946       1.6018         1.6216       0.2326
    ##   qual5         1.5026         1.5013       1.9363         1.9452       0.1980
    ##   val1          1.4479         1.3696       1.8788         1.7644       0.2455
    ##   val2          1.2266         1.2148       1.6518         1.7063       0.1710
    ##   val3          1.4842         1.3938       1.9709         1.9381       0.1457
    ##   sat1          1.2469         1.2321       1.6465         1.6195       0.3403
    ##   sat2          1.2348         1.2028       1.6450         1.6241       0.3092
    ##   sat3          1.3379         1.2834       1.6724         1.7142       0.2114
    ##   sat4          1.3270         1.2579       1.6770         1.6307       0.2727
    ##   loy1          1.7000         1.6735       2.2468         2.2390       0.2635
    ##   loy2          1.4869         1.4918       1.9172         1.9882       0.1302
    ##   loy3          1.7006         1.6612       2.2862         2.2584       0.2686
    ##   loy4          1.6820         1.6804       2.1649         2.2936       0.0956
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
    ##  Number of admissible results       = 485
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 275567512
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
    ##   EXPE ~ IMAG      0.4714      0.0659    7.1538    0.0000 [ 0.3360; 0.5976 ] 
    ##   QUAL ~ EXPE      0.8344      0.0232   35.9310    0.0000 [ 0.7845; 0.8760 ] 
    ##   VAL ~ EXPE       0.0457      0.0846    0.5402    0.5890 [-0.0982; 0.2220 ] 
    ##   VAL ~ QUAL       0.7013      0.0823    8.5233    0.0000 [ 0.5153; 0.8526 ] 
    ##   SAT ~ IMAG       0.2450      0.0524    4.6729    0.0000 [ 0.1278; 0.3402 ] 
    ##   SAT ~ EXPE      -0.0172      0.0705   -0.2444    0.8069 [-0.1434; 0.1177 ] 
    ##   SAT ~ QUAL       0.2215      0.1048    2.1150    0.0344 [ 0.0260; 0.4359 ] 
    ##   SAT ~ VAL        0.5270      0.0870    6.0540    0.0000 [ 0.3525; 0.6978 ] 
    ##   LOY ~ IMAG       0.1819      0.0801    2.2708    0.0232 [ 0.0431; 0.3400 ] 
    ##   LOY ~ SAT        0.6283      0.0818    7.6814    0.0000 [ 0.4612; 0.7828 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0986    6.3937    0.0000 [ 0.4235; 0.7974 ] 
    ##   IMAG =~ imag2      0.9246      0.0421   21.9543    0.0000 [ 0.8285; 0.9781 ] 
    ##   IMAG =~ imag3      0.9577      0.0291   32.8991    0.0000 [ 0.8814; 0.9923 ] 
    ##   EXPE =~ expe1      0.7525      0.0739   10.1860    0.0000 [ 0.5836; 0.8718 ] 
    ##   EXPE =~ expe2      0.9348      0.0285   32.7953    0.0000 [ 0.8638; 0.9734 ] 
    ##   EXPE =~ expe3      0.7295      0.0725   10.0680    0.0000 [ 0.5523; 0.8522 ] 
    ##   QUAL =~ qual1      0.7861      0.0656   11.9912    0.0000 [ 0.6284; 0.8920 ] 
    ##   QUAL =~ qual2      0.9244      0.0246   37.5551    0.0000 [ 0.8669; 0.9588 ] 
    ##   QUAL =~ qual3      0.7560      0.0619   12.2083    0.0000 [ 0.6077; 0.8540 ] 
    ##   QUAL =~ qual4      0.7632      0.0540   14.1330    0.0000 [ 0.6390; 0.8562 ] 
    ##   QUAL =~ qual5      0.7834      0.0471   16.6457    0.0000 [ 0.6804; 0.8563 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.5971    0.0000 [ 0.8946; 0.9839 ] 
    ##   VAL =~ val2        0.8056      0.0685   11.7665    0.0000 [ 0.6522; 0.9096 ] 
    ##   VAL =~ val3        0.6763      0.0781    8.6544    0.0000 [ 0.5080; 0.8127 ] 
    ##   SAT =~ sat1        0.9243      0.0232   39.7748    0.0000 [ 0.8753; 0.9675 ] 
    ##   SAT =~ sat2        0.8813      0.0296   29.7484    0.0000 [ 0.8208; 0.9299 ] 
    ##   SAT =~ sat3        0.7127      0.0521   13.6815    0.0000 [ 0.6017; 0.8064 ] 
    ##   SAT =~ sat4        0.7756      0.0480   16.1485    0.0000 [ 0.6763; 0.8567 ] 
    ##   LOY =~ loy1        0.9097      0.0521   17.4436    0.0000 [ 0.7807; 0.9892 ] 
    ##   LOY =~ loy2        0.5775      0.0874    6.6095    0.0000 [ 0.4050; 0.7333 ] 
    ##   LOY =~ loy3        0.9043      0.0436   20.7318    0.0000 [ 0.8134; 0.9753 ] 
    ##   LOY =~ loy4        0.4917      0.0945    5.2062    0.0000 [ 0.3303; 0.6997 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1176    0.1330    0.8942 [-0.1988; 0.2499 ] 
    ##   IMAG <~ imag2      0.4473      0.1580    2.8317    0.0046 [ 0.1036; 0.7385 ] 
    ##   IMAG <~ imag3      0.6020      0.1449    4.1559    0.0000 [ 0.3206; 0.8558 ] 
    ##   EXPE <~ expe1      0.2946      0.1130    2.6069    0.0091 [ 0.0703; 0.5154 ] 
    ##   EXPE <~ expe2      0.6473      0.0855    7.5690    0.0000 [ 0.4715; 0.8015 ] 
    ##   EXPE <~ expe3      0.2374      0.0919    2.5837    0.0098 [ 0.0451; 0.4034 ] 
    ##   QUAL <~ qual1      0.2370      0.0890    2.6632    0.0077 [ 0.0862; 0.4429 ] 
    ##   QUAL <~ qual2      0.4712      0.0852    5.5330    0.0000 [ 0.3155; 0.6280 ] 
    ##   QUAL <~ qual3      0.1831      0.0820    2.2336    0.0255 [ 0.0049; 0.3360 ] 
    ##   QUAL <~ qual4      0.1037      0.0617    1.6816    0.0926 [-0.0197; 0.2124 ] 
    ##   QUAL <~ qual5      0.2049      0.0621    3.2969    0.0010 [ 0.0615; 0.3123 ] 
    ##   VAL <~ val1        0.7163      0.1000    7.1621    0.0000 [ 0.5202; 0.8867 ] 
    ##   VAL <~ val2        0.2202      0.0942    2.3367    0.0195 [ 0.0502; 0.4023 ] 
    ##   VAL <~ val3        0.2082      0.0608    3.4254    0.0006 [ 0.0929; 0.3207 ] 
    ##   SAT <~ sat1        0.3209      0.0149   21.4743    0.0000 [ 0.2945; 0.3533 ] 
    ##   SAT <~ sat2        0.3059      0.0138   22.1657    0.0000 [ 0.2824; 0.3366 ] 
    ##   SAT <~ sat3        0.2474      0.0111   22.3297    0.0000 [ 0.2237; 0.2644 ] 
    ##   SAT <~ sat4        0.2692      0.0121   22.2430    0.0000 [ 0.2452; 0.2929 ] 
    ##   LOY <~ loy1        0.3834      0.0267   14.3848    0.0000 [ 0.3324; 0.4343 ] 
    ##   LOY <~ loy2        0.2434      0.0296    8.2142    0.0000 [ 0.1767; 0.2921 ] 
    ##   LOY <~ loy3        0.3812      0.0281   13.5507    0.0000 [ 0.3275; 0.4341 ] 
    ##   LOY <~ loy4        0.2073      0.0353    5.8724    0.0000 [ 0.1433; 0.2815 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0610   10.5457    0.0000 [ 0.5069; 0.7492 ] 
    ##   imag1 ~~ imag3      0.5433      0.0658    8.2558    0.0000 [ 0.4182; 0.6803 ] 
    ##   imag2 ~~ imag3      0.7761      0.0403   19.2362    0.0000 [ 0.6913; 0.8522 ] 
    ##   expe1 ~~ expe2      0.5353      0.0572    9.3556    0.0000 [ 0.4183; 0.6484 ] 
    ##   expe1 ~~ expe3      0.4694      0.0611    7.6774    0.0000 [ 0.3543; 0.5901 ] 
    ##   expe2 ~~ expe3      0.5467      0.0602    9.0825    0.0000 [ 0.4214; 0.6535 ] 
    ##   qual1 ~~ qual2      0.6053      0.0570   10.6270    0.0000 [ 0.4913; 0.7041 ] 
    ##   qual1 ~~ qual3      0.5406      0.0560    9.6551    0.0000 [ 0.4283; 0.6471 ] 
    ##   qual1 ~~ qual4      0.5662      0.0666    8.4970    0.0000 [ 0.4406; 0.6880 ] 
    ##   qual1 ~~ qual5      0.5180      0.0667    7.7721    0.0000 [ 0.3876; 0.6431 ] 
    ##   qual2 ~~ qual3      0.6187      0.0556   11.1275    0.0000 [ 0.5010; 0.7170 ] 
    ##   qual2 ~~ qual4      0.6517      0.0625   10.4192    0.0000 [ 0.5232; 0.7626 ] 
    ##   qual2 ~~ qual5      0.6291      0.0602   10.4548    0.0000 [ 0.5022; 0.7274 ] 
    ##   qual3 ~~ qual4      0.4752      0.0664    7.1536    0.0000 [ 0.3423; 0.5960 ] 
    ##   qual3 ~~ qual5      0.5074      0.0644    7.8807    0.0000 [ 0.3821; 0.6224 ] 
    ##   qual4 ~~ qual5      0.6402      0.0586   10.9336    0.0000 [ 0.5264; 0.7379 ] 
    ##   val1 ~~ val2        0.6344      0.0589   10.7641    0.0000 [ 0.4968; 0.7377 ] 
    ##   val1 ~~ val3        0.4602      0.0753    6.1108    0.0000 [ 0.3008; 0.5907 ] 
    ##   val2 ~~ val3        0.6288      0.0635    9.8952    0.0000 [ 0.5059; 0.7453 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0659    7.1538    0.0000 [ 0.3360; 0.5976 ] 
    ##   QUAL ~ IMAG       0.3933      0.0614    6.4089    0.0000 [ 0.2691; 0.5144 ] 
    ##   QUAL ~ EXPE       0.8344      0.0232   35.9310    0.0000 [ 0.7845; 0.8760 ] 
    ##   VAL ~ IMAG        0.2974      0.0610    4.8720    0.0000 [ 0.1868; 0.4238 ] 
    ##   VAL ~ EXPE        0.6309      0.0521   12.1168    0.0000 [ 0.5374; 0.7267 ] 
    ##   VAL ~ QUAL        0.7013      0.0823    8.5233    0.0000 [ 0.5153; 0.8526 ] 
    ##   SAT ~ IMAG        0.4807      0.0651    7.3835    0.0000 [ 0.3561; 0.6045 ] 
    ##   SAT ~ EXPE        0.5001      0.0578    8.6524    0.0000 [ 0.3855; 0.6122 ] 
    ##   SAT ~ QUAL        0.5911      0.0959    6.1636    0.0000 [ 0.3930; 0.7691 ] 
    ##   SAT ~ VAL         0.5270      0.0870    6.0540    0.0000 [ 0.3525; 0.6978 ] 
    ##   LOY ~ IMAG        0.4840      0.0673    7.1935    0.0000 [ 0.3466; 0.6113 ] 
    ##   LOY ~ EXPE        0.3142      0.0556    5.6506    0.0000 [ 0.2078; 0.4213 ] 
    ##   LOY ~ QUAL        0.3714      0.0825    4.5019    0.0000 [ 0.2158; 0.5337 ] 
    ##   LOY ~ VAL         0.3311      0.0727    4.5544    0.0000 [ 0.2004; 0.4853 ] 
    ##   LOY ~ SAT         0.6283      0.0818    7.6814    0.0000 [ 0.4612; 0.7828 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0614    6.4089    0.0000 [ 0.2691; 0.5144 ] 
    ##   VAL ~ IMAG           0.2974      0.0610    4.8720    0.0000 [ 0.1868; 0.4238 ] 
    ##   VAL ~ EXPE           0.5852      0.0717    8.1671    0.0000 [ 0.4288; 0.7201 ] 
    ##   SAT ~ IMAG           0.2357      0.0484    4.8663    0.0000 [ 0.1454; 0.3351 ] 
    ##   SAT ~ EXPE           0.5173      0.0711    7.2784    0.0000 [ 0.3828; 0.6571 ] 
    ##   SAT ~ QUAL           0.3696      0.0662    5.5847    0.0000 [ 0.2445; 0.4981 ] 
    ##   LOY ~ IMAG           0.3020      0.0573    5.2722    0.0000 [ 0.1945; 0.4156 ] 
    ##   LOY ~ EXPE           0.3142      0.0556    5.6506    0.0000 [ 0.2078; 0.4213 ] 
    ##   LOY ~ QUAL           0.3714      0.0825    4.5019    0.0000 [ 0.2158; 0.5337 ] 
    ##   LOY ~ VAL            0.3311      0.0727    4.5544    0.0000 [ 0.2004; 0.4853 ] 
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
