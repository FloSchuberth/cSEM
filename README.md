
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

## News (2021-04-16):

-   New function `exportToExcel()`. The function conveniently exports
    the results from `assess()`, `predict()`, `summarize()` and
    `testOMF()` to an .xlsx file.

-   Several bug fixes, most notably: `calculateVifModeB()` did not
    calculate the VIFs for modeB constructs correctly because of a bug
    in the calculation of the R^2. PLEASE REVIEW YOUR CALCULATIONS in
    cSEM version &lt; 0.3.1:9000! (thanks to @Benjamin Liengaard for
    pointing it out)

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
    ##  dG                      0.6493      0.3136  
    ##  SRMR                    0.0940      0.0530  
    ##  dL                      2.2340      0.7109  
    ##  dML                     2.9219      1.5660  
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
    ##  The seed used was: 370797356
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
    ## ------------ Variance inflation factors (VIFs) for modeB constructs ------------
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
    ##   expe1         1.4521         1.5696       1.9042         2.1021       0.0559
    ##   expe2         1.4135         1.4755       1.9373         2.0310       0.1955
    ##   expe3         1.6332         1.7189       2.1300         2.2166       0.1205
    ##   qual1         1.4740         1.5559       1.9276         2.0808       0.1154
    ##   qual2         1.5766         1.5361       2.0393         2.0756       0.2145
    ##   qual3         1.7296         1.7219       2.2232         2.2799       0.1182
    ##   qual4         1.2322         1.1903       1.5938         1.6249       0.2328
    ##   qual5         1.5033         1.5069       1.9328         1.9660       0.1981
    ##   val1          1.4472         1.3627       1.8739         1.7686       0.2455
    ##   val2          1.2332         1.2106       1.6560         1.7266       0.1716
    ##   val3          1.4774         1.3708       1.9667         1.9255       0.1476
    ##   sat1          1.2491         1.2364       1.6512         1.6223       0.3383
    ##   sat2          1.2363         1.1923       1.6444         1.6209       0.3077
    ##   sat3          1.3404         1.2777       1.6696         1.7204       0.2141
    ##   sat4          1.3187         1.2598       1.6685         1.6347       0.2771
    ##   loy1          1.6891         1.6535       2.2237         2.2102       0.2736
    ##   loy2          1.4790         1.4682       1.9122         1.9845       0.1314
    ##   loy3          1.6946         1.6597       2.2762         2.2544       0.2706
    ##   loy4          1.6866         1.6656       2.1718         2.2857       0.0880
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential quantities are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways to compute resamples:

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

Now `summarize()` shows inferential quantities as well:

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
    ##  Number of admissible results     = 493
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = -1252687134
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
    ##   EXPE ~ IMAG      0.4714      0.0658    7.1624    0.0000 [ 0.3428; 0.5925 ] 
    ##   QUAL ~ EXPE      0.8344      0.0231   36.0770    0.0000 [ 0.7868; 0.8747 ] 
    ##   VAL ~ EXPE       0.0457      0.0837    0.5465    0.5847 [-0.1017; 0.2212 ] 
    ##   VAL ~ QUAL       0.7013      0.0832    8.4266    0.0000 [ 0.5298; 0.8549 ] 
    ##   SAT ~ IMAG       0.2450      0.0556    4.4020    0.0000 [ 0.1382; 0.3537 ] 
    ##   SAT ~ EXPE      -0.0172      0.0675   -0.2555    0.7983 [-0.1616; 0.1064 ] 
    ##   SAT ~ QUAL       0.2215      0.0980    2.2616    0.0237 [ 0.0240; 0.4135 ] 
    ##   SAT ~ VAL        0.5270      0.0826    6.3797    0.0000 [ 0.3691; 0.6911 ] 
    ##   LOY ~ IMAG       0.1819      0.0777    2.3418    0.0192 [ 0.0272; 0.3293 ] 
    ##   LOY ~ SAT        0.6283      0.0756    8.3131    0.0000 [ 0.4740; 0.7700 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0944    6.6801    0.0000 [ 0.4093; 0.7940 ] 
    ##   IMAG =~ imag2      0.9246      0.0402   22.9880    0.0000 [ 0.8290; 0.9753 ] 
    ##   IMAG =~ imag3      0.9577      0.0304   31.5136    0.0000 [ 0.8797; 0.9925 ] 
    ##   EXPE =~ expe1      0.7525      0.0755    9.9608    0.0000 [ 0.5707; 0.8741 ] 
    ##   EXPE =~ expe2      0.9348      0.0275   33.9457    0.0000 [ 0.8617; 0.9714 ] 
    ##   EXPE =~ expe3      0.7295      0.0791    9.2264    0.0000 [ 0.5354; 0.8433 ] 
    ##   QUAL =~ qual1      0.7861      0.0680   11.5553    0.0000 [ 0.6278; 0.8880 ] 
    ##   QUAL =~ qual2      0.9244      0.0229   40.2920    0.0000 [ 0.8722; 0.9568 ] 
    ##   QUAL =~ qual3      0.7560      0.0657   11.5020    0.0000 [ 0.5912; 0.8463 ] 
    ##   QUAL =~ qual4      0.7632      0.0530   14.3877    0.0000 [ 0.6374; 0.8454 ] 
    ##   QUAL =~ qual5      0.7834      0.0490   15.9865    0.0000 [ 0.6769; 0.8653 ] 
    ##   VAL =~ val1        0.9518      0.0238   39.9598    0.0000 [ 0.8910; 0.9846 ] 
    ##   VAL =~ val2        0.8056      0.0638   12.6226    0.0000 [ 0.6602; 0.9077 ] 
    ##   VAL =~ val3        0.6763      0.0719    9.4036    0.0000 [ 0.5209; 0.7956 ] 
    ##   SAT =~ sat1        0.9243      0.0239   38.7242    0.0000 [ 0.8684; 0.9623 ] 
    ##   SAT =~ sat2        0.8813      0.0299   29.4419    0.0000 [ 0.8179; 0.9288 ] 
    ##   SAT =~ sat3        0.7127      0.0531   13.4215    0.0000 [ 0.5972; 0.8014 ] 
    ##   SAT =~ sat4        0.7756      0.0477   16.2679    0.0000 [ 0.6740; 0.8633 ] 
    ##   LOY =~ loy1        0.9097      0.0506   17.9616    0.0000 [ 0.7979; 0.9829 ] 
    ##   LOY =~ loy2        0.5775      0.0850    6.7953    0.0000 [ 0.4106; 0.7251 ] 
    ##   LOY =~ loy3        0.9043      0.0436   20.7197    0.0000 [ 0.8037; 0.9723 ] 
    ##   LOY =~ loy4        0.4917      0.1021    4.8179    0.0000 [ 0.2956; 0.6853 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1094    0.1429    0.8863 [-0.1924; 0.2393 ] 
    ##   IMAG <~ imag2      0.4473      0.1524    2.9348    0.0033 [ 0.1538; 0.7666 ] 
    ##   IMAG <~ imag3      0.6020      0.1387    4.3405    0.0000 [ 0.3289; 0.8573 ] 
    ##   EXPE <~ expe1      0.2946      0.1179    2.4996    0.0124 [ 0.0678; 0.5140 ] 
    ##   EXPE <~ expe2      0.6473      0.0804    8.0539    0.0000 [ 0.4650; 0.7804 ] 
    ##   EXPE <~ expe3      0.2374      0.0959    2.4755    0.0133 [ 0.0331; 0.4197 ] 
    ##   QUAL <~ qual1      0.2370      0.0913    2.5957    0.0094 [ 0.0768; 0.4164 ] 
    ##   QUAL <~ qual2      0.4712      0.0764    6.1686    0.0000 [ 0.3168; 0.6089 ] 
    ##   QUAL <~ qual3      0.1831      0.0816    2.2443    0.0248 [ 0.0113; 0.3312 ] 
    ##   QUAL <~ qual4      0.1037      0.0595    1.7422    0.0815 [-0.0186; 0.2190 ] 
    ##   QUAL <~ qual5      0.2049      0.0636    3.2218    0.0013 [ 0.0796; 0.3229 ] 
    ##   VAL <~ val1        0.7163      0.0934    7.6718    0.0000 [ 0.5085; 0.8802 ] 
    ##   VAL <~ val2        0.2202      0.0925    2.3804    0.0173 [ 0.0584; 0.4182 ] 
    ##   VAL <~ val3        0.2082      0.0645    3.2272    0.0013 [ 0.0838; 0.3298 ] 
    ##   SAT <~ sat1        0.3209      0.0149   21.5705    0.0000 [ 0.2957; 0.3539 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5559    0.0000 [ 0.2827; 0.3346 ] 
    ##   SAT <~ sat3        0.2474      0.0112   22.0758    0.0000 [ 0.2244; 0.2679 ] 
    ##   SAT <~ sat4        0.2692      0.0124   21.7438    0.0000 [ 0.2477; 0.2970 ] 
    ##   LOY <~ loy1        0.3834      0.0282   13.5938    0.0000 [ 0.3317; 0.4415 ] 
    ##   LOY <~ loy2        0.2434      0.0303    8.0415    0.0000 [ 0.1781; 0.2948 ] 
    ##   LOY <~ loy3        0.3812      0.0264   14.4140    0.0000 [ 0.3334; 0.4341 ] 
    ##   LOY <~ loy4        0.2073      0.0381    5.4443    0.0000 [ 0.1335; 0.2794 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0628   10.2423    0.0000 [ 0.5180; 0.7518 ] 
    ##   imag1 ~~ imag3      0.5433      0.0688    7.8973    0.0000 [ 0.3958; 0.6709 ] 
    ##   imag2 ~~ imag3      0.7761      0.0406   19.1337    0.0000 [ 0.6817; 0.8411 ] 
    ##   expe1 ~~ expe2      0.5353      0.0584    9.1589    0.0000 [ 0.4156; 0.6431 ] 
    ##   expe1 ~~ expe3      0.4694      0.0635    7.3956    0.0000 [ 0.3395; 0.5835 ] 
    ##   expe2 ~~ expe3      0.5467      0.0640    8.5447    0.0000 [ 0.4116; 0.6566 ] 
    ##   qual1 ~~ qual2      0.6053      0.0575   10.5328    0.0000 [ 0.4794; 0.7099 ] 
    ##   qual1 ~~ qual3      0.5406      0.0614    8.8112    0.0000 [ 0.4109; 0.6527 ] 
    ##   qual1 ~~ qual4      0.5662      0.0683    8.2903    0.0000 [ 0.4110; 0.6799 ] 
    ##   qual1 ~~ qual5      0.5180      0.0698    7.4224    0.0000 [ 0.3690; 0.6395 ] 
    ##   qual2 ~~ qual3      0.6187      0.0574   10.7758    0.0000 [ 0.4942; 0.7118 ] 
    ##   qual2 ~~ qual4      0.6517      0.0611   10.6635    0.0000 [ 0.5226; 0.7655 ] 
    ##   qual2 ~~ qual5      0.6291      0.0608   10.3453    0.0000 [ 0.5039; 0.7267 ] 
    ##   qual3 ~~ qual4      0.4752      0.0685    6.9396    0.0000 [ 0.3282; 0.5893 ] 
    ##   qual3 ~~ qual5      0.5074      0.0662    7.6640    0.0000 [ 0.3563; 0.6224 ] 
    ##   qual4 ~~ qual5      0.6402      0.0568   11.2781    0.0000 [ 0.5162; 0.7376 ] 
    ##   val1 ~~ val2        0.6344      0.0555   11.4293    0.0000 [ 0.5042; 0.7325 ] 
    ##   val1 ~~ val3        0.4602      0.0694    6.6352    0.0000 [ 0.3227; 0.5861 ] 
    ##   val2 ~~ val3        0.6288      0.0654    9.6125    0.0000 [ 0.4908; 0.7404 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0658    7.1624    0.0000 [ 0.3428; 0.5925 ] 
    ##   QUAL ~ IMAG       0.3933      0.0614    6.4097    0.0000 [ 0.2748; 0.5100 ] 
    ##   QUAL ~ EXPE       0.8344      0.0231   36.0770    0.0000 [ 0.7868; 0.8747 ] 
    ##   VAL ~ IMAG        0.2974      0.0602    4.9432    0.0000 [ 0.1895; 0.4223 ] 
    ##   VAL ~ EXPE        0.6309      0.0501   12.6047    0.0000 [ 0.5285; 0.7275 ] 
    ##   VAL ~ QUAL        0.7013      0.0832    8.4266    0.0000 [ 0.5298; 0.8549 ] 
    ##   SAT ~ IMAG        0.4807      0.0656    7.3223    0.0000 [ 0.3348; 0.6019 ] 
    ##   SAT ~ EXPE        0.5001      0.0559    8.9470    0.0000 [ 0.3939; 0.6039 ] 
    ##   SAT ~ QUAL        0.5911      0.0909    6.5028    0.0000 [ 0.3986; 0.7472 ] 
    ##   SAT ~ VAL         0.5270      0.0826    6.3797    0.0000 [ 0.3691; 0.6911 ] 
    ##   LOY ~ IMAG        0.4840      0.0664    7.2835    0.0000 [ 0.3556; 0.6073 ] 
    ##   LOY ~ EXPE        0.3142      0.0523    6.0091    0.0000 [ 0.2125; 0.4170 ] 
    ##   LOY ~ QUAL        0.3714      0.0790    4.7029    0.0000 [ 0.2341; 0.5399 ] 
    ##   LOY ~ VAL         0.3311      0.0684    4.8402    0.0000 [ 0.1995; 0.4655 ] 
    ##   LOY ~ SAT         0.6283      0.0756    8.3131    0.0000 [ 0.4740; 0.7700 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0614    6.4097    0.0000 [ 0.2748; 0.5100 ] 
    ##   VAL ~ IMAG           0.2974      0.0602    4.9432    0.0000 [ 0.1895; 0.4223 ] 
    ##   VAL ~ EXPE           0.5852      0.0716    8.1743    0.0000 [ 0.4303; 0.7108 ] 
    ##   SAT ~ IMAG           0.2357      0.0469    5.0262    0.0000 [ 0.1549; 0.3281 ] 
    ##   SAT ~ EXPE           0.5173      0.0657    7.8778    0.0000 [ 0.3948; 0.6471 ] 
    ##   SAT ~ QUAL           0.3696      0.0601    6.1441    0.0000 [ 0.2447; 0.4972 ] 
    ##   LOY ~ IMAG           0.3020      0.0535    5.6455    0.0000 [ 0.2067; 0.4081 ] 
    ##   LOY ~ EXPE           0.3142      0.0523    6.0091    0.0000 [ 0.2125; 0.4170 ] 
    ##   LOY ~ QUAL           0.3714      0.0790    4.7029    0.0000 [ 0.2341; 0.5399 ] 
    ##   LOY ~ VAL            0.3311      0.0684    4.8402    0.0000 [ 0.1995; 0.4655 ] 
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
  .eval_plan       = "multiprocess")
```
