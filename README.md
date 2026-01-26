
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

## News (2025-05-15):

- Implementation of doModelSearch to perform AGAS-PLS. Thanks to Gloria.

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
# install.packages("devtools")
devtools::install_github("FloSchuberth/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://floschuberth.github.io/cSEM/).

## Basic usage

The basic usage is illustrated below.

<img src="man/figures/api.png" width="80%" style="display: block; margin: auto;" />

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
    ##  dG                      0.6493      0.3275  
    ##  SRMR                    0.0940      0.0516  
    ##  dL                      2.2340      0.6745  
    ##  dML                     2.9219      1.6552  
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
    ##  Out of 499 bootstrap replications 475 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 30148742
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
    ##   expe1         1.4555         1.6089       1.9023         2.1231       0.0550
    ##   expe2         1.4141         1.5070       1.9311         2.0347       0.2028
    ##   expe3         1.6307         1.7444       2.1271         2.2279       0.1273
    ##   qual1         1.4774         1.5784       1.9289         2.0977       0.1140
    ##   qual2         1.5841         1.5479       2.0461         2.0611       0.2154
    ##   qual3         1.7307         1.7416       2.2191         2.2831       0.1227
    ##   qual4         1.2359         1.1979       1.6049         1.6430       0.2279
    ##   qual5         1.5006         1.5105       1.9288         1.9508       0.1992
    ##   val1          1.4525         1.3765       1.8780         1.7709       0.2455
    ##   val2          1.2242         1.2209       1.6472         1.7199       0.1710
    ##   val3          1.4812         1.3992       1.9674         1.9484       0.1474
    ##   sat1          1.2483         1.2334       1.6453         1.6251       0.3367
    ##   sat2          1.2282         1.1984       1.6372         1.6255       0.3094
    ##   sat3          1.3383         1.2892       1.6680         1.7209       0.2107
    ##   sat4          1.3095         1.2591       1.6604         1.6342       0.2793
    ##   loy1          1.6887         1.6566       2.2265         2.2189       0.2699
    ##   loy2          1.4892         1.4970       1.9216         1.9944       0.1281
    ##   loy3          1.7121         1.6799       2.2905         2.2769       0.2645
    ##   loy4          1.6818         1.6850       2.1732         2.3035       0.0901
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
    ##  Number of admissible results       = 484
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 1469621969
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
    ##   EXPE ~ IMAG      0.4714      0.0636    7.4090    0.0000 [ 0.3496; 0.6008 ] 
    ##   QUAL ~ EXPE      0.8344      0.0247   33.8101    0.0000 [ 0.7757; 0.8727 ] 
    ##   VAL ~ EXPE       0.0457      0.0831    0.5502    0.5822 [-0.1095; 0.2074 ] 
    ##   VAL ~ QUAL       0.7013      0.0798    8.7913    0.0000 [ 0.5468; 0.8490 ] 
    ##   SAT ~ IMAG       0.2450      0.0569    4.3072    0.0000 [ 0.1367; 0.3635 ] 
    ##   SAT ~ EXPE      -0.0172      0.0689   -0.2503    0.8023 [-0.1527; 0.1075 ] 
    ##   SAT ~ QUAL       0.2215      0.1018    2.1753    0.0296 [ 0.0372; 0.4233 ] 
    ##   SAT ~ VAL        0.5270      0.0913    5.7716    0.0000 [ 0.3470; 0.6991 ] 
    ##   LOY ~ IMAG       0.1819      0.0792    2.2962    0.0217 [ 0.0367; 0.3469 ] 
    ##   LOY ~ SAT        0.6283      0.0793    7.9211    0.0000 [ 0.4853; 0.7816 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0971    6.4950    0.0000 [ 0.4319; 0.7954 ] 
    ##   IMAG =~ imag2      0.9246      0.0419   22.0528    0.0000 [ 0.8271; 0.9783 ] 
    ##   IMAG =~ imag3      0.9577      0.0276   34.6501    0.0000 [ 0.8845; 0.9906 ] 
    ##   EXPE =~ expe1      0.7525      0.0798    9.4294    0.0000 [ 0.5643; 0.8773 ] 
    ##   EXPE =~ expe2      0.9348      0.0291   32.0725    0.0000 [ 0.8579; 0.9722 ] 
    ##   EXPE =~ expe3      0.7295      0.0746    9.7764    0.0000 [ 0.5516; 0.8521 ] 
    ##   QUAL =~ qual1      0.7861      0.0684   11.5004    0.0000 [ 0.6380; 0.8862 ] 
    ##   QUAL =~ qual2      0.9244      0.0227   40.8032    0.0000 [ 0.8713; 0.9561 ] 
    ##   QUAL =~ qual3      0.7560      0.0595   12.7111    0.0000 [ 0.6197; 0.8566 ] 
    ##   QUAL =~ qual4      0.7632      0.0538   14.1798    0.0000 [ 0.6425; 0.8537 ] 
    ##   QUAL =~ qual5      0.7834      0.0504   15.5363    0.0000 [ 0.6700; 0.8590 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.7082    0.0000 [ 0.8971; 0.9865 ] 
    ##   VAL =~ val2        0.8056      0.0648   12.4253    0.0000 [ 0.6577; 0.9055 ] 
    ##   VAL =~ val3        0.6763      0.0692    9.7765    0.0000 [ 0.5280; 0.8072 ] 
    ##   SAT =~ sat1        0.9243      0.0217   42.5420    0.0000 [ 0.8749; 0.9607 ] 
    ##   SAT =~ sat2        0.8813      0.0293   30.1082    0.0000 [ 0.8029; 0.9279 ] 
    ##   SAT =~ sat3        0.7127      0.0523   13.6335    0.0000 [ 0.6061; 0.8128 ] 
    ##   SAT =~ sat4        0.7756      0.0502   15.4621    0.0000 [ 0.6759; 0.8606 ] 
    ##   LOY =~ loy1        0.9097      0.0508   17.8991    0.0000 [ 0.7892; 0.9808 ] 
    ##   LOY =~ loy2        0.5775      0.0850    6.7955    0.0000 [ 0.3866; 0.7301 ] 
    ##   LOY =~ loy3        0.9043      0.0429   21.0872    0.0000 [ 0.8056; 0.9722 ] 
    ##   LOY =~ loy4        0.4917      0.0973    5.0559    0.0000 [ 0.3051; 0.6913 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1161    0.1347    0.8928 [-0.2177; 0.2395 ] 
    ##   IMAG <~ imag2      0.4473      0.1481    3.0210    0.0025 [ 0.1723; 0.7353 ] 
    ##   IMAG <~ imag3      0.6020      0.1377    4.3715    0.0000 [ 0.3089; 0.8329 ] 
    ##   EXPE <~ expe1      0.2946      0.1166    2.5260    0.0115 [ 0.0553; 0.5255 ] 
    ##   EXPE <~ expe2      0.6473      0.0893    7.2509    0.0000 [ 0.4360; 0.7928 ] 
    ##   EXPE <~ expe3      0.2374      0.0915    2.5956    0.0094 [ 0.0622; 0.4262 ] 
    ##   QUAL <~ qual1      0.2370      0.0876    2.7055    0.0068 [ 0.0733; 0.4245 ] 
    ##   QUAL <~ qual2      0.4712      0.0784    6.0138    0.0000 [ 0.3099; 0.6107 ] 
    ##   QUAL <~ qual3      0.1831      0.0764    2.3965    0.0166 [ 0.0369; 0.3347 ] 
    ##   QUAL <~ qual4      0.1037      0.0611    1.6973    0.0896 [-0.0041; 0.2236 ] 
    ##   QUAL <~ qual5      0.2049      0.0644    3.1811    0.0015 [ 0.0678; 0.3181 ] 
    ##   VAL <~ val1        0.7163      0.0970    7.3836    0.0000 [ 0.5246; 0.8857 ] 
    ##   VAL <~ val2        0.2202      0.0950    2.3175    0.0205 [ 0.0410; 0.4113 ] 
    ##   VAL <~ val3        0.2082      0.0623    3.3400    0.0008 [ 0.0873; 0.3264 ] 
    ##   SAT <~ sat1        0.3209      0.0152   21.1486    0.0000 [ 0.2953; 0.3502 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5035    0.0000 [ 0.2823; 0.3363 ] 
    ##   SAT <~ sat3        0.2474      0.0109   22.7922    0.0000 [ 0.2275; 0.2677 ] 
    ##   SAT <~ sat4        0.2692      0.0128   20.9781    0.0000 [ 0.2451; 0.2951 ] 
    ##   LOY <~ loy1        0.3834      0.0274   14.0083    0.0000 [ 0.3286; 0.4330 ] 
    ##   LOY <~ loy2        0.2434      0.0302    8.0512    0.0000 [ 0.1736; 0.2933 ] 
    ##   LOY <~ loy3        0.3812      0.0270   14.1044    0.0000 [ 0.3319; 0.4368 ] 
    ##   LOY <~ loy4        0.2073      0.0365    5.6732    0.0000 [ 0.1330; 0.2786 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0648    9.9284    0.0000 [ 0.5112; 0.7608 ] 
    ##   imag1 ~~ imag3      0.5433      0.0678    8.0170    0.0000 [ 0.4007; 0.6722 ] 
    ##   imag2 ~~ imag3      0.7761      0.0391   19.8351    0.0000 [ 0.6946; 0.8451 ] 
    ##   expe1 ~~ expe2      0.5353      0.0593    9.0241    0.0000 [ 0.4014; 0.6334 ] 
    ##   expe1 ~~ expe3      0.4694      0.0597    7.8626    0.0000 [ 0.3368; 0.5770 ] 
    ##   expe2 ~~ expe3      0.5467      0.0638    8.5705    0.0000 [ 0.4242; 0.6616 ] 
    ##   qual1 ~~ qual2      0.6053      0.0574   10.5391    0.0000 [ 0.4894; 0.7135 ] 
    ##   qual1 ~~ qual3      0.5406      0.0627    8.6260    0.0000 [ 0.4112; 0.6547 ] 
    ##   qual1 ~~ qual4      0.5662      0.0694    8.1526    0.0000 [ 0.4217; 0.7029 ] 
    ##   qual1 ~~ qual5      0.5180      0.0685    7.5579    0.0000 [ 0.3757; 0.6563 ] 
    ##   qual2 ~~ qual3      0.6187      0.0583   10.6170    0.0000 [ 0.4929; 0.7136 ] 
    ##   qual2 ~~ qual4      0.6517      0.0616   10.5806    0.0000 [ 0.5242; 0.7639 ] 
    ##   qual2 ~~ qual5      0.6291      0.0578   10.8873    0.0000 [ 0.5011; 0.7276 ] 
    ##   qual3 ~~ qual4      0.4752      0.0658    7.2168    0.0000 [ 0.3446; 0.5888 ] 
    ##   qual3 ~~ qual5      0.5074      0.0643    7.8940    0.0000 [ 0.3751; 0.6189 ] 
    ##   qual4 ~~ qual5      0.6402      0.0592   10.8102    0.0000 [ 0.5149; 0.7413 ] 
    ##   val1 ~~ val2        0.6344      0.0559   11.3544    0.0000 [ 0.5098; 0.7485 ] 
    ##   val1 ~~ val3        0.4602      0.0686    6.7117    0.0000 [ 0.3222; 0.5947 ] 
    ##   val2 ~~ val3        0.6288      0.0637    9.8707    0.0000 [ 0.5014; 0.7462 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0636    7.4090    0.0000 [ 0.3496; 0.6008 ] 
    ##   QUAL ~ IMAG       0.3933      0.0589    6.6772    0.0000 [ 0.2783; 0.5088 ] 
    ##   QUAL ~ EXPE       0.8344      0.0247   33.8101    0.0000 [ 0.7757; 0.8727 ] 
    ##   VAL ~ IMAG        0.2974      0.0572    5.1957    0.0000 [ 0.1942; 0.4166 ] 
    ##   VAL ~ EXPE        0.6309      0.0480   13.1354    0.0000 [ 0.5313; 0.7193 ] 
    ##   VAL ~ QUAL        0.7013      0.0798    8.7913    0.0000 [ 0.5468; 0.8490 ] 
    ##   SAT ~ IMAG        0.4807      0.0668    7.1949    0.0000 [ 0.3552; 0.6118 ] 
    ##   SAT ~ EXPE        0.5001      0.0560    8.9308    0.0000 [ 0.3870; 0.6011 ] 
    ##   SAT ~ QUAL        0.5911      0.0932    6.3424    0.0000 [ 0.4125; 0.7654 ] 
    ##   SAT ~ VAL         0.5270      0.0913    5.7716    0.0000 [ 0.3470; 0.6991 ] 
    ##   LOY ~ IMAG        0.4840      0.0654    7.3987    0.0000 [ 0.3533; 0.6118 ] 
    ##   LOY ~ EXPE        0.3142      0.0525    5.9824    0.0000 [ 0.2147; 0.4217 ] 
    ##   LOY ~ QUAL        0.3714      0.0819    4.5337    0.0000 [ 0.2299; 0.5349 ] 
    ##   LOY ~ VAL         0.3311      0.0779    4.2521    0.0000 [ 0.1902; 0.4973 ] 
    ##   LOY ~ SAT         0.6283      0.0793    7.9211    0.0000 [ 0.4853; 0.7816 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0589    6.6772    0.0000 [ 0.2783; 0.5088 ] 
    ##   VAL ~ IMAG           0.2974      0.0572    5.1957    0.0000 [ 0.1942; 0.4166 ] 
    ##   VAL ~ EXPE           0.5852      0.0709    8.2595    0.0000 [ 0.4489; 0.7260 ] 
    ##   SAT ~ IMAG           0.2357      0.0464    5.0787    0.0000 [ 0.1531; 0.3295 ] 
    ##   SAT ~ EXPE           0.5173      0.0659    7.8498    0.0000 [ 0.3959; 0.6466 ] 
    ##   SAT ~ QUAL           0.3696      0.0646    5.7229    0.0000 [ 0.2427; 0.4944 ] 
    ##   LOY ~ IMAG           0.3020      0.0555    5.4387    0.0000 [ 0.2010; 0.4165 ] 
    ##   LOY ~ EXPE           0.3142      0.0525    5.9824    0.0000 [ 0.2147; 0.4217 ] 
    ##   LOY ~ QUAL           0.3714      0.0819    4.5337    0.0000 [ 0.2299; 0.5349 ] 
    ##   LOY ~ VAL            0.3311      0.0779    4.2521    0.0000 [ 0.1902; 0.4973 ] 
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
