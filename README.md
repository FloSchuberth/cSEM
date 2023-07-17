
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
devtools::install_github("M-E-Rademaker/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://m-e-rademaker.github.io/cSEM/).

## Basic usage

The basic usage is illustrated below.

<div class="white">

<img src="man/figures/api.png" width="80%" style="display: block; margin: auto;" />

</div>

Usually, using `cSEM` is the same 3 step procedure:

> 1.  Pick a dataset and specify a model using [lavaan
>     syntax](https://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use `csem()`
> 3.  Apply one of the post-estimation functions listed below on the
>     resulting object.

## Post-estimation functions

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
    ##  dG                      0.6493      0.3256  
    ##  SRMR                    0.0940      0.0519  
    ##  dL                      2.2340      0.6808  
    ##  dML                     2.9219      1.6331  
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
    ##  The seed used was: -2112241140
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
    ##   expe1         1.4546         1.5962       1.9060         2.1068       0.0554
    ##   expe2         1.4073         1.4993       1.9294         2.0298       0.2041
    ##   expe3         1.6310         1.7420       2.1213         2.2193       0.1300
    ##   qual1         1.4759         1.5539       1.9312         2.0714       0.1137
    ##   qual2         1.5732         1.5345       2.0327         2.0531       0.2223
    ##   qual3         1.7306         1.7399       2.2199         2.2807       0.1251
    ##   qual4         1.2347         1.1973       1.5962         1.6310       0.2340
    ##   qual5         1.5061         1.5205       1.9325         1.9652       0.2044
    ##   val1          1.4409         1.3664       1.8632         1.7681       0.2540
    ##   val2          1.2297         1.2220       1.6503         1.7220       0.1712
    ##   val3          1.4853         1.3955       1.9709         1.9349       0.1453
    ##   sat1          1.2470         1.2262       1.6463         1.6163       0.3382
    ##   sat2          1.2312         1.1984       1.6397         1.6252       0.3083
    ##   sat3          1.3356         1.2841       1.6708         1.7178       0.2078
    ##   sat4          1.3192         1.2601       1.6723         1.6363       0.2751
    ##   loy1          1.6965         1.6688       2.2367         2.2395       0.2724
    ##   loy2          1.4925         1.4952       1.9201         1.9945       0.1267
    ##   loy3          1.7018         1.6768       2.2800         2.2863       0.2740
    ##   loy4          1.6985         1.6855       2.1985         2.3090       0.0754
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
    ##  Number of admissible results       = 489
    ##  Approach to handle inadmissibles   = "drop"
    ##  Sign change option                 = "none"
    ##  Random seed                        = 519363542
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
    ##   EXPE ~ IMAG      0.4714      0.0661    7.1282    0.0000 [ 0.3445; 0.5955 ] 
    ##   QUAL ~ EXPE      0.8344      0.0234   35.7261    0.0000 [ 0.7859; 0.8751 ] 
    ##   VAL ~ EXPE       0.0457      0.0849    0.5383    0.5904 [-0.1059; 0.2227 ] 
    ##   VAL ~ QUAL       0.7013      0.0817    8.5796    0.0000 [ 0.5332; 0.8494 ] 
    ##   SAT ~ IMAG       0.2450      0.0564    4.3455    0.0000 [ 0.1460; 0.3597 ] 
    ##   SAT ~ EXPE      -0.0172      0.0716   -0.2407    0.8098 [-0.1607; 0.1145 ] 
    ##   SAT ~ QUAL       0.2215      0.0927    2.3901    0.0168 [ 0.0438; 0.4028 ] 
    ##   SAT ~ VAL        0.5270      0.0849    6.2061    0.0000 [ 0.3648; 0.7003 ] 
    ##   LOY ~ IMAG       0.1819      0.0789    2.3057    0.0211 [ 0.0334; 0.3328 ] 
    ##   LOY ~ SAT        0.6283      0.0792    7.9350    0.0000 [ 0.4735; 0.7819 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1053    5.9864    0.0000 [ 0.3941; 0.8087 ] 
    ##   IMAG =~ imag2      0.9246      0.0393   23.5084    0.0000 [ 0.8301; 0.9755 ] 
    ##   IMAG =~ imag3      0.9577      0.0283   33.8577    0.0000 [ 0.8852; 0.9902 ] 
    ##   EXPE =~ expe1      0.7525      0.0838    8.9797    0.0000 [ 0.5571; 0.8797 ] 
    ##   EXPE =~ expe2      0.9348      0.0288   32.5033    0.0000 [ 0.8674; 0.9727 ] 
    ##   EXPE =~ expe3      0.7295      0.0694   10.5075    0.0000 [ 0.5781; 0.8401 ] 
    ##   QUAL =~ qual1      0.7861      0.0719   10.9288    0.0000 [ 0.6220; 0.8880 ] 
    ##   QUAL =~ qual2      0.9244      0.0246   37.5815    0.0000 [ 0.8623; 0.9613 ] 
    ##   QUAL =~ qual3      0.7560      0.0579   13.0671    0.0000 [ 0.6337; 0.8499 ] 
    ##   QUAL =~ qual4      0.7632      0.0522   14.6247    0.0000 [ 0.6518; 0.8607 ] 
    ##   QUAL =~ qual5      0.7834      0.0457   17.1414    0.0000 [ 0.6838; 0.8630 ] 
    ##   VAL =~ val1        0.9518      0.0229   41.6273    0.0000 [ 0.8997; 0.9847 ] 
    ##   VAL =~ val2        0.8056      0.0644   12.5061    0.0000 [ 0.6658; 0.9081 ] 
    ##   VAL =~ val3        0.6763      0.0746    9.0610    0.0000 [ 0.5204; 0.8152 ] 
    ##   SAT =~ sat1        0.9243      0.0218   42.3875    0.0000 [ 0.8799; 0.9635 ] 
    ##   SAT =~ sat2        0.8813      0.0293   30.0744    0.0000 [ 0.8225; 0.9299 ] 
    ##   SAT =~ sat3        0.7127      0.0511   13.9377    0.0000 [ 0.6036; 0.8015 ] 
    ##   SAT =~ sat4        0.7756      0.0490   15.8259    0.0000 [ 0.6723; 0.8587 ] 
    ##   LOY =~ loy1        0.9097      0.0488   18.6534    0.0000 [ 0.7973; 0.9839 ] 
    ##   LOY =~ loy2        0.5775      0.0843    6.8475    0.0000 [ 0.4074; 0.7303 ] 
    ##   LOY =~ loy3        0.9043      0.0423   21.3785    0.0000 [ 0.8034; 0.9769 ] 
    ##   LOY =~ loy4        0.4917      0.0965    5.0952    0.0000 [ 0.2982; 0.6658 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1157    0.1352    0.8925 [-0.2148; 0.2560 ] 
    ##   IMAG <~ imag2      0.4473      0.1410    3.1722    0.0015 [ 0.1823; 0.7096 ] 
    ##   IMAG <~ imag3      0.6020      0.1371    4.3907    0.0000 [ 0.3201; 0.8393 ] 
    ##   EXPE <~ expe1      0.2946      0.1203    2.4481    0.0144 [ 0.0584; 0.5040 ] 
    ##   EXPE <~ expe2      0.6473      0.0861    7.5205    0.0000 [ 0.4817; 0.7951 ] 
    ##   EXPE <~ expe3      0.2374      0.0905    2.6234    0.0087 [ 0.0548; 0.4085 ] 
    ##   QUAL <~ qual1      0.2370      0.0918    2.5820    0.0098 [ 0.0596; 0.4236 ] 
    ##   QUAL <~ qual2      0.4712      0.0793    5.9438    0.0000 [ 0.3077; 0.6258 ] 
    ##   QUAL <~ qual3      0.1831      0.0787    2.3268    0.0200 [ 0.0277; 0.3486 ] 
    ##   QUAL <~ qual4      0.1037      0.0634    1.6354    0.1020 [-0.0117; 0.2397 ] 
    ##   QUAL <~ qual5      0.2049      0.0643    3.1854    0.0014 [ 0.0767; 0.3241 ] 
    ##   VAL <~ val1        0.7163      0.0959    7.4732    0.0000 [ 0.5239; 0.8807 ] 
    ##   VAL <~ val2        0.2202      0.0937    2.3497    0.0188 [ 0.0495; 0.4041 ] 
    ##   VAL <~ val3        0.2082      0.0630    3.3017    0.0010 [ 0.0879; 0.3251 ] 
    ##   SAT <~ sat1        0.3209      0.0150   21.4186    0.0000 [ 0.2967; 0.3517 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5695    0.0000 [ 0.2831; 0.3349 ] 
    ##   SAT <~ sat3        0.2474      0.0107   23.1406    0.0000 [ 0.2263; 0.2681 ] 
    ##   SAT <~ sat4        0.2692      0.0119   22.5414    0.0000 [ 0.2448; 0.2915 ] 
    ##   LOY <~ loy1        0.3834      0.0255   15.0082    0.0000 [ 0.3343; 0.4385 ] 
    ##   LOY <~ loy2        0.2434      0.0284    8.5797    0.0000 [ 0.1803; 0.2913 ] 
    ##   LOY <~ loy3        0.3812      0.0269   14.1789    0.0000 [ 0.3297; 0.4370 ] 
    ##   LOY <~ loy4        0.2073      0.0365    5.6779    0.0000 [ 0.1339; 0.2735 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0676    9.5193    0.0000 [ 0.4970; 0.7531 ] 
    ##   imag1 ~~ imag3      0.5433      0.0740    7.3459    0.0000 [ 0.3743; 0.6705 ] 
    ##   imag2 ~~ imag3      0.7761      0.0408   19.0049    0.0000 [ 0.6891; 0.8417 ] 
    ##   expe1 ~~ expe2      0.5353      0.0648    8.2662    0.0000 [ 0.4035; 0.6500 ] 
    ##   expe1 ~~ expe3      0.4694      0.0627    7.4823    0.0000 [ 0.3441; 0.5904 ] 
    ##   expe2 ~~ expe3      0.5467      0.0589    9.2790    0.0000 [ 0.4240; 0.6423 ] 
    ##   qual1 ~~ qual2      0.6053      0.0590   10.2613    0.0000 [ 0.4836; 0.7075 ] 
    ##   qual1 ~~ qual3      0.5406      0.0595    9.0869    0.0000 [ 0.4187; 0.6513 ] 
    ##   qual1 ~~ qual4      0.5662      0.0653    8.6642    0.0000 [ 0.4328; 0.6835 ] 
    ##   qual1 ~~ qual5      0.5180      0.0673    7.6929    0.0000 [ 0.3690; 0.6375 ] 
    ##   qual2 ~~ qual3      0.6187      0.0543   11.3839    0.0000 [ 0.4982; 0.7051 ] 
    ##   qual2 ~~ qual4      0.6517      0.0584   11.1607    0.0000 [ 0.5375; 0.7555 ] 
    ##   qual2 ~~ qual5      0.6291      0.0561   11.2052    0.0000 [ 0.5184; 0.7305 ] 
    ##   qual3 ~~ qual4      0.4752      0.0613    7.7578    0.0000 [ 0.3582; 0.5944 ] 
    ##   qual3 ~~ qual5      0.5074      0.0588    8.6233    0.0000 [ 0.3991; 0.6225 ] 
    ##   qual4 ~~ qual5      0.6402      0.0548   11.6818    0.0000 [ 0.5268; 0.7397 ] 
    ##   val1 ~~ val2        0.6344      0.0547   11.5910    0.0000 [ 0.5285; 0.7358 ] 
    ##   val1 ~~ val3        0.4602      0.0728    6.3240    0.0000 [ 0.3229; 0.6042 ] 
    ##   val2 ~~ val3        0.6288      0.0610   10.3169    0.0000 [ 0.5067; 0.7353 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0661    7.1282    0.0000 [ 0.3445; 0.5955 ] 
    ##   QUAL ~ IMAG       0.3933      0.0615    6.3989    0.0000 [ 0.2810; 0.5081 ] 
    ##   QUAL ~ EXPE       0.8344      0.0234   35.7261    0.0000 [ 0.7859; 0.8751 ] 
    ##   VAL ~ IMAG        0.2974      0.0619    4.8051    0.0000 [ 0.1919; 0.4317 ] 
    ##   VAL ~ EXPE        0.6309      0.0520   12.1281    0.0000 [ 0.5261; 0.7348 ] 
    ##   VAL ~ QUAL        0.7013      0.0817    8.5796    0.0000 [ 0.5332; 0.8494 ] 
    ##   SAT ~ IMAG        0.4807      0.0666    7.2187    0.0000 [ 0.3491; 0.6111 ] 
    ##   SAT ~ EXPE        0.5001      0.0607    8.2431    0.0000 [ 0.3732; 0.6101 ] 
    ##   SAT ~ QUAL        0.5911      0.0901    6.5623    0.0000 [ 0.4151; 0.7572 ] 
    ##   SAT ~ VAL         0.5270      0.0849    6.2061    0.0000 [ 0.3648; 0.7003 ] 
    ##   LOY ~ IMAG        0.4840      0.0685    7.0650    0.0000 [ 0.3574; 0.6176 ] 
    ##   LOY ~ EXPE        0.3142      0.0562    5.5897    0.0000 [ 0.2116; 0.4382 ] 
    ##   LOY ~ QUAL        0.3714      0.0805    4.6127    0.0000 [ 0.2220; 0.5554 ] 
    ##   LOY ~ VAL         0.3311      0.0721    4.5913    0.0000 [ 0.2036; 0.4808 ] 
    ##   LOY ~ SAT         0.6283      0.0792    7.9350    0.0000 [ 0.4735; 0.7819 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0615    6.3989    0.0000 [ 0.2810; 0.5081 ] 
    ##   VAL ~ IMAG           0.2974      0.0619    4.8051    0.0000 [ 0.1919; 0.4317 ] 
    ##   VAL ~ EXPE           0.5852      0.0701    8.3434    0.0000 [ 0.4412; 0.7155 ] 
    ##   SAT ~ IMAG           0.2357      0.0504    4.6801    0.0000 [ 0.1468; 0.3386 ] 
    ##   SAT ~ EXPE           0.5173      0.0657    7.8700    0.0000 [ 0.3849; 0.6444 ] 
    ##   SAT ~ QUAL           0.3696      0.0620    5.9642    0.0000 [ 0.2507; 0.4975 ] 
    ##   LOY ~ IMAG           0.3020      0.0555    5.4383    0.0000 [ 0.2088; 0.4245 ] 
    ##   LOY ~ EXPE           0.3142      0.0562    5.5897    0.0000 [ 0.2116; 0.4382 ] 
    ##   LOY ~ QUAL           0.3714      0.0805    4.6127    0.0000 [ 0.2220; 0.5554 ] 
    ##   LOY ~ VAL            0.3311      0.0721    4.5913    0.0000 [ 0.2036; 0.4808 ] 
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
