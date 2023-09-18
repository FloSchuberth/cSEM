
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

## News (2022-05-09):

- `predict()` function is now able to predict categorical indicators (a
  procedure known as OrdPLScPredict).

- Use singular value decomposition in GSCAm to deal with large datasets,
  which allows for large datasets (thanks to Heungsun Hwang)

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
- `summarize()` : summarize the results
- `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

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
    ##  SRMR                    0.0940      0.0530  
    ##  dL                      2.2340      0.7099  
    ##  dML                     2.9219      1.6308  
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
    ##  The seed used was: 980464158
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
    ##   expe1         1.4681         1.6013       1.9215         2.1197       0.0435
    ##   expe2         1.4233         1.5085       1.9455         2.0427       0.1907
    ##   expe3         1.6278         1.7411       2.1223         2.2187       0.1224
    ##   qual1         1.4835         1.5636       1.9342         2.0818       0.1115
    ##   qual2         1.5822         1.5459       2.0474         2.0710       0.2111
    ##   qual3         1.7312         1.7355       2.2205         2.2806       0.1176
    ##   qual4         1.2308         1.2007       1.5951         1.6432       0.2345
    ##   qual5         1.5032         1.5178       1.9349         1.9653       0.1958
    ##   val1          1.4477         1.3711       1.8707         1.7731       0.2481
    ##   val2          1.2276         1.2270       1.6515         1.7353       0.1713
    ##   val3          1.4746         1.3913       1.9619         1.9351       0.1520
    ##   sat1          1.2565         1.2362       1.6547         1.6278       0.3319
    ##   sat2          1.2372         1.2051       1.6465         1.6348       0.3035
    ##   sat3          1.3398         1.2889       1.6715         1.7256       0.2105
    ##   sat4          1.3171         1.2679       1.6718         1.6510       0.2734
    ##   loy1          1.7100         1.6798       2.2523         2.2520       0.2581
    ##   loy2          1.4907         1.4904       1.9127         1.9860       0.1354
    ##   loy3          1.7170         1.6789       2.2952         2.2773       0.2609
    ##   loy4          1.7045         1.7085       2.1899         2.3266       0.0850
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
    ##  Number of admissible results       = 487
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = -558304815
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
    ##   EXPE ~ IMAG      0.4714      0.0645    7.3038    0.0000 [ 0.3572; 0.5912 ] 
    ##   QUAL ~ EXPE      0.8344      0.0239   34.8885    0.0000 [ 0.7835; 0.8732 ] 
    ##   VAL ~ EXPE       0.0457      0.0830    0.5505    0.5820 [-0.1025; 0.2206 ] 
    ##   VAL ~ QUAL       0.7013      0.0796    8.8114    0.0000 [ 0.5327; 0.8562 ] 
    ##   SAT ~ IMAG       0.2450      0.0541    4.5315    0.0000 [ 0.1378; 0.3406 ] 
    ##   SAT ~ EXPE      -0.0172      0.0702   -0.2457    0.8059 [-0.1644; 0.1206 ] 
    ##   SAT ~ QUAL       0.2215      0.1034    2.1427    0.0321 [ 0.0491; 0.4451 ] 
    ##   SAT ~ VAL        0.5270      0.0878    5.9992    0.0000 [ 0.3340; 0.6845 ] 
    ##   LOY ~ IMAG       0.1819      0.0799    2.2764    0.0228 [ 0.0531; 0.3645 ] 
    ##   LOY ~ SAT        0.6283      0.0827    7.5987    0.0000 [ 0.4545; 0.7812 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0982    6.4240    0.0000 [ 0.4163; 0.8001 ] 
    ##   IMAG =~ imag2      0.9246      0.0414   22.3433    0.0000 [ 0.8208; 0.9761 ] 
    ##   IMAG =~ imag3      0.9577      0.0286   33.4955    0.0000 [ 0.8813; 0.9945 ] 
    ##   EXPE =~ expe1      0.7525      0.0759    9.9200    0.0000 [ 0.5820; 0.8621 ] 
    ##   EXPE =~ expe2      0.9348      0.0277   33.7516    0.0000 [ 0.8571; 0.9717 ] 
    ##   EXPE =~ expe3      0.7295      0.0700   10.4149    0.0000 [ 0.5649; 0.8461 ] 
    ##   QUAL =~ qual1      0.7861      0.0657   11.9741    0.0000 [ 0.6252; 0.8846 ] 
    ##   QUAL =~ qual2      0.9244      0.0218   42.4976    0.0000 [ 0.8691; 0.9566 ] 
    ##   QUAL =~ qual3      0.7560      0.0595   12.7096    0.0000 [ 0.6145; 0.8500 ] 
    ##   QUAL =~ qual4      0.7632      0.0531   14.3660    0.0000 [ 0.6453; 0.8510 ] 
    ##   QUAL =~ qual5      0.7834      0.0440   17.7936    0.0000 [ 0.6891; 0.8574 ] 
    ##   VAL =~ val1        0.9518      0.0246   38.6484    0.0000 [ 0.8941; 0.9854 ] 
    ##   VAL =~ val2        0.8056      0.0647   12.4583    0.0000 [ 0.6711; 0.9056 ] 
    ##   VAL =~ val3        0.6763      0.0733    9.2279    0.0000 [ 0.5271; 0.8104 ] 
    ##   SAT =~ sat1        0.9243      0.0240   38.5351    0.0000 [ 0.8710; 0.9662 ] 
    ##   SAT =~ sat2        0.8813      0.0308   28.6399    0.0000 [ 0.8099; 0.9273 ] 
    ##   SAT =~ sat3        0.7127      0.0550   12.9677    0.0000 [ 0.5923; 0.8009 ] 
    ##   SAT =~ sat4        0.7756      0.0500   15.5128    0.0000 [ 0.6697; 0.8675 ] 
    ##   LOY =~ loy1        0.9097      0.0501   18.1412    0.0000 [ 0.7959; 0.9838 ] 
    ##   LOY =~ loy2        0.5775      0.0833    6.9336    0.0000 [ 0.4093; 0.7281 ] 
    ##   LOY =~ loy3        0.9043      0.0428   21.1431    0.0000 [ 0.8068; 0.9688 ] 
    ##   LOY =~ loy4        0.4917      0.1016    4.8403    0.0000 [ 0.2938; 0.6902 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1192    0.1312    0.8956 [-0.2039; 0.2429 ] 
    ##   IMAG <~ imag2      0.4473      0.1540    2.9037    0.0037 [ 0.1171; 0.7296 ] 
    ##   IMAG <~ imag3      0.6020      0.1389    4.3349    0.0000 [ 0.3265; 0.8840 ] 
    ##   EXPE <~ expe1      0.2946      0.1121    2.6281    0.0086 [ 0.0595; 0.4870 ] 
    ##   EXPE <~ expe2      0.6473      0.0785    8.2469    0.0000 [ 0.4702; 0.7858 ] 
    ##   EXPE <~ expe3      0.2374      0.0902    2.6320    0.0085 [ 0.0560; 0.4173 ] 
    ##   QUAL <~ qual1      0.2370      0.0889    2.6667    0.0077 [ 0.0688; 0.4078 ] 
    ##   QUAL <~ qual2      0.4712      0.0739    6.3735    0.0000 [ 0.3344; 0.6058 ] 
    ##   QUAL <~ qual3      0.1831      0.0783    2.3383    0.0194 [ 0.0237; 0.3205 ] 
    ##   QUAL <~ qual4      0.1037      0.0585    1.7722    0.0764 [-0.0022; 0.2202 ] 
    ##   QUAL <~ qual5      0.2049      0.0603    3.3952    0.0007 [ 0.0849; 0.3109 ] 
    ##   VAL <~ val1        0.7163      0.0987    7.2574    0.0000 [ 0.5053; 0.8754 ] 
    ##   VAL <~ val2        0.2202      0.0939    2.3457    0.0190 [ 0.0601; 0.4145 ] 
    ##   VAL <~ val3        0.2082      0.0634    3.2811    0.0010 [ 0.0952; 0.3370 ] 
    ##   SAT <~ sat1        0.3209      0.0156   20.5742    0.0000 [ 0.2963; 0.3535 ] 
    ##   SAT <~ sat2        0.3059      0.0141   21.6609    0.0000 [ 0.2808; 0.3365 ] 
    ##   SAT <~ sat3        0.2474      0.0119   20.8686    0.0000 [ 0.2233; 0.2707 ] 
    ##   SAT <~ sat4        0.2692      0.0125   21.5568    0.0000 [ 0.2471; 0.2966 ] 
    ##   LOY <~ loy1        0.3834      0.0271   14.1384    0.0000 [ 0.3314; 0.4345 ] 
    ##   LOY <~ loy2        0.2434      0.0288    8.4639    0.0000 [ 0.1803; 0.2913 ] 
    ##   LOY <~ loy3        0.3812      0.0273   13.9809    0.0000 [ 0.3275; 0.4349 ] 
    ##   LOY <~ loy4        0.2073      0.0377    5.4925    0.0000 [ 0.1325; 0.2795 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0660    9.7459    0.0000 [ 0.5046; 0.7537 ] 
    ##   imag1 ~~ imag3      0.5433      0.0692    7.8542    0.0000 [ 0.3989; 0.6627 ] 
    ##   imag2 ~~ imag3      0.7761      0.0393   19.7338    0.0000 [ 0.6929; 0.8385 ] 
    ##   expe1 ~~ expe2      0.5353      0.0584    9.1611    0.0000 [ 0.4000; 0.6434 ] 
    ##   expe1 ~~ expe3      0.4694      0.0587    8.0001    0.0000 [ 0.3451; 0.5749 ] 
    ##   expe2 ~~ expe3      0.5467      0.0606    9.0141    0.0000 [ 0.4078; 0.6554 ] 
    ##   qual1 ~~ qual2      0.6053      0.0568   10.6494    0.0000 [ 0.4852; 0.7079 ] 
    ##   qual1 ~~ qual3      0.5406      0.0568    9.5226    0.0000 [ 0.4302; 0.6485 ] 
    ##   qual1 ~~ qual4      0.5662      0.0671    8.4422    0.0000 [ 0.4362; 0.6877 ] 
    ##   qual1 ~~ qual5      0.5180      0.0642    8.0657    0.0000 [ 0.3943; 0.6359 ] 
    ##   qual2 ~~ qual3      0.6187      0.0552   11.2007    0.0000 [ 0.5043; 0.7127 ] 
    ##   qual2 ~~ qual4      0.6517      0.0613   10.6365    0.0000 [ 0.5252; 0.7516 ] 
    ##   qual2 ~~ qual5      0.6291      0.0545   11.5384    0.0000 [ 0.5168; 0.7250 ] 
    ##   qual3 ~~ qual4      0.4752      0.0627    7.5821    0.0000 [ 0.3551; 0.5928 ] 
    ##   qual3 ~~ qual5      0.5074      0.0609    8.3332    0.0000 [ 0.3957; 0.6226 ] 
    ##   qual4 ~~ qual5      0.6402      0.0579   11.0481    0.0000 [ 0.5187; 0.7365 ] 
    ##   val1 ~~ val2        0.6344      0.0539   11.7782    0.0000 [ 0.5247; 0.7355 ] 
    ##   val1 ~~ val3        0.4602      0.0702    6.5561    0.0000 [ 0.3208; 0.5949 ] 
    ##   val2 ~~ val3        0.6288      0.0611   10.2953    0.0000 [ 0.5148; 0.7399 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0645    7.3038    0.0000 [ 0.3572; 0.5912 ] 
    ##   QUAL ~ IMAG       0.3933      0.0605    6.5021    0.0000 [ 0.2854; 0.5030 ] 
    ##   QUAL ~ EXPE       0.8344      0.0239   34.8885    0.0000 [ 0.7835; 0.8732 ] 
    ##   VAL ~ IMAG        0.2974      0.0600    4.9526    0.0000 [ 0.1930; 0.4247 ] 
    ##   VAL ~ EXPE        0.6309      0.0507   12.4388    0.0000 [ 0.5356; 0.7312 ] 
    ##   VAL ~ QUAL        0.7013      0.0796    8.8114    0.0000 [ 0.5327; 0.8562 ] 
    ##   SAT ~ IMAG        0.4807      0.0667    7.2079    0.0000 [ 0.3445; 0.6084 ] 
    ##   SAT ~ EXPE        0.5001      0.0569    8.7964    0.0000 [ 0.3856; 0.6051 ] 
    ##   SAT ~ QUAL        0.5911      0.0934    6.3297    0.0000 [ 0.4312; 0.7866 ] 
    ##   SAT ~ VAL         0.5270      0.0878    5.9992    0.0000 [ 0.3340; 0.6845 ] 
    ##   LOY ~ IMAG        0.4840      0.0671    7.2176    0.0000 [ 0.3701; 0.6248 ] 
    ##   LOY ~ EXPE        0.3142      0.0540    5.8235    0.0000 [ 0.2088; 0.4148 ] 
    ##   LOY ~ QUAL        0.3714      0.0800    4.6407    0.0000 [ 0.2269; 0.5264 ] 
    ##   LOY ~ VAL         0.3311      0.0769    4.3039    0.0000 [ 0.1797; 0.4755 ] 
    ##   LOY ~ SAT         0.6283      0.0827    7.5987    0.0000 [ 0.4545; 0.7812 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0605    6.5021    0.0000 [ 0.2854; 0.5030 ] 
    ##   VAL ~ IMAG           0.2974      0.0600    4.9526    0.0000 [ 0.1930; 0.4247 ] 
    ##   VAL ~ EXPE           0.5852      0.0680    8.6095    0.0000 [ 0.4469; 0.7188 ] 
    ##   SAT ~ IMAG           0.2357      0.0487    4.8406    0.0000 [ 0.1469; 0.3346 ] 
    ##   SAT ~ EXPE           0.5173      0.0663    7.7999    0.0000 [ 0.3961; 0.6611 ] 
    ##   SAT ~ QUAL           0.3696      0.0630    5.8700    0.0000 [ 0.2426; 0.4801 ] 
    ##   LOY ~ IMAG           0.3020      0.0554    5.4563    0.0000 [ 0.1934; 0.4106 ] 
    ##   LOY ~ EXPE           0.3142      0.0540    5.8235    0.0000 [ 0.2088; 0.4148 ] 
    ##   LOY ~ QUAL           0.3714      0.0800    4.6407    0.0000 [ 0.2269; 0.5264 ] 
    ##   LOY ~ VAL            0.3311      0.0769    4.3039    0.0000 [ 0.1797; 0.4755 ] 
    ## ________________________________________________________________________________

Several bootstrap-based confidence intervals are implemented, see
`?infer()`:

``` r
infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet
```

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as setting random seeds via the [future
framework](https://github.com/HenrikBengtsson/future/). For
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
