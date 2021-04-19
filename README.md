
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
<!-- WARNING: THIS IS WORK IN PROGRESS. BREAKING CHANGES TO THE API ARE VERY LIKELY.  -->
<!--          Use the package with caution and please report bugs to [the package developers](mailto:manuel.rademaker@uni-wuerzburg.de;f.schuberth@utwente.nl).  -->
<!--          The first stable relase will be version 0.0.1, most likely towards the end -->
<!--          of 2019. -->

## Purpose

Estimate, analyse, test, and study linear, nonlinear, hierachical and
multi-group structural equation models using composite-based approaches
and procedures, including estimation techniques such as partial least
squares path modeling (PLS-PM) and its derivatives (PLSc, OrdPLSc,
robustPLSc), generalized structured component analysis (GSCA),
generalized structured component analysis with uniqueness terms (GSCAm),
generalized canonical correlation analysis (GCCA), principal component
analysis (PCA), factor score regression (FSR) using sum score,
regression or bartlett scores (including bias correction using Croon’s
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
-   `infer()` : calculate common inferencial quantities (e.g., standard
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
    ##  dG                      0.6493      0.3250  
    ##  SRMR                    0.0940      0.0520  
    ##  dL                      2.2340      0.6852  
    ##  dML                     2.9219      1.6110  
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
    ##  The seed used was: -1196039766
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
    ##  Squared Euclidian distance  = 2.23402
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
    ## ------------------------------ Validity assessment -----------------------------
    ## 
    ##  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
    ## 
    ##           SAT LOY
    ## SAT 1.0000000   0
    ## LOY 0.7432489   1
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
    ##  Number of repetitions            = 10
    ##  Handle inadmissibles             = stop
    ##  Target                           = 'PLS-PM'
    ##  Benchmark                        = 'lm'
    ## 
    ## ------------------------------ Prediction metrics ------------------------------
    ## 
    ## 
    ##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
    ##   expe1         1.4563         1.5698       1.9071         2.0985       0.0538
    ##   expe2         1.4097         1.4771       1.9310         2.0265       0.2017
    ##   expe3         1.6318         1.7260       2.1258         2.2215       0.1238
    ##   qual1         1.4757         1.5489       1.9290         2.0677       0.1155
    ##   qual2         1.5754         1.5320       2.0348         2.0536       0.2193
    ##   qual3         1.7284         1.7240       2.2189         2.2776       0.1213
    ##   qual4         1.2328         1.1960       1.5939         1.6281       0.2346
    ##   qual5         1.5050         1.5029       1.9338         1.9566       0.1975
    ##   val1          1.4449         1.3622       1.8690         1.7661       0.2497
    ##   val2          1.2266         1.2060       1.6475         1.7141       0.1743
    ##   val3          1.4791         1.3808       1.9658         1.9352       0.1496
    ##   sat1          1.2451         1.2318       1.6431         1.6177       0.3406
    ##   sat2          1.2317         1.1954       1.6399         1.6261       0.3089
    ##   sat3          1.3410         1.2726       1.6732         1.7178       0.2088
    ##   sat4          1.3195         1.2623       1.6697         1.6349       0.2753
    ##   loy1          1.6899         1.6574       2.2304         2.2238       0.2699
    ##   loy2          1.4858         1.4702       1.9145         1.9783       0.1304
    ##   loy3          1.7040         1.6675       2.2810         2.2680       0.2693
    ##   loy4          1.6883         1.6693       2.1772         2.3036       0.0862
    ## ________________________________________________________________________________

#### Resampling and Inference

By default no inferential quantities are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways to compute resamples:

1.  Setting `.resample_method` in `csem()` to `"jackkinfe"` or
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

Now `summarize()` shows inferencial quantities as well:

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
    ##  Number of admissible results     = 490
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = -429071268
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
    ##   EXPE ~ IMAG      0.4714      0.0678    6.9506    0.0000 [ 0.3345; 0.6033 ] 
    ##   QUAL ~ EXPE      0.8344      0.0233   35.7851    0.0000 [ 0.7845; 0.8764 ] 
    ##   VAL ~ EXPE       0.0457      0.0903    0.5063    0.6127 [-0.1310; 0.2231 ] 
    ##   VAL ~ QUAL       0.7013      0.0872    8.0439    0.0000 [ 0.5318; 0.8690 ] 
    ##   SAT ~ IMAG       0.2450      0.0563    4.3527    0.0000 [ 0.1461; 0.3618 ] 
    ##   SAT ~ EXPE      -0.0172      0.0728   -0.2367    0.8129 [-0.1548; 0.1146 ] 
    ##   SAT ~ QUAL       0.2215      0.1040    2.1307    0.0331 [ 0.0511; 0.4375 ] 
    ##   SAT ~ VAL        0.5270      0.0915    5.7561    0.0000 [ 0.3393; 0.6936 ] 
    ##   LOY ~ IMAG       0.1819      0.0800    2.2748    0.0229 [ 0.0217; 0.3341 ] 
    ##   LOY ~ SAT        0.6283      0.0807    7.7822    0.0000 [ 0.4803; 0.7882 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1001    6.2995    0.0000 [ 0.4114; 0.7973 ] 
    ##   IMAG =~ imag2      0.9246      0.0399   23.1620    0.0000 [ 0.8273; 0.9782 ] 
    ##   IMAG =~ imag3      0.9577      0.0292   32.7890    0.0000 [ 0.8754; 0.9919 ] 
    ##   EXPE =~ expe1      0.7525      0.0827    9.0997    0.0000 [ 0.5700; 0.8720 ] 
    ##   EXPE =~ expe2      0.9348      0.0299   31.2152    0.0000 [ 0.8563; 0.9724 ] 
    ##   EXPE =~ expe3      0.7295      0.0781    9.3409    0.0000 [ 0.5531; 0.8471 ] 
    ##   QUAL =~ qual1      0.7861      0.0738   10.6524    0.0000 [ 0.6189; 0.8879 ] 
    ##   QUAL =~ qual2      0.9244      0.0221   41.8371    0.0000 [ 0.8696; 0.9551 ] 
    ##   QUAL =~ qual3      0.7560      0.0648   11.6719    0.0000 [ 0.6156; 0.8620 ] 
    ##   QUAL =~ qual4      0.7632      0.0530   14.4015    0.0000 [ 0.6560; 0.8533 ] 
    ##   QUAL =~ qual5      0.7834      0.0475   16.5036    0.0000 [ 0.6773; 0.8695 ] 
    ##   VAL =~ val1        0.9518      0.0241   39.5078    0.0000 [ 0.8907; 0.9870 ] 
    ##   VAL =~ val2        0.8056      0.0665   12.1111    0.0000 [ 0.6502; 0.9058 ] 
    ##   VAL =~ val3        0.6763      0.0749    9.0338    0.0000 [ 0.5236; 0.8035 ] 
    ##   SAT =~ sat1        0.9243      0.0238   38.8846    0.0000 [ 0.8684; 0.9661 ] 
    ##   SAT =~ sat2        0.8813      0.0292   30.1631    0.0000 [ 0.8151; 0.9304 ] 
    ##   SAT =~ sat3        0.7127      0.0524   13.5931    0.0000 [ 0.6077; 0.8108 ] 
    ##   SAT =~ sat4        0.7756      0.0497   15.6019    0.0000 [ 0.6705; 0.8646 ] 
    ##   LOY =~ loy1        0.9097      0.0502   18.1153    0.0000 [ 0.7940; 0.9890 ] 
    ##   LOY =~ loy2        0.5775      0.0869    6.6442    0.0000 [ 0.3938; 0.7307 ] 
    ##   LOY =~ loy3        0.9043      0.0416   21.7641    0.0000 [ 0.8148; 0.9733 ] 
    ##   LOY =~ loy4        0.4917      0.1040    4.7264    0.0000 [ 0.2705; 0.6854 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1140    0.1372    0.8909 [-0.2153; 0.2367 ] 
    ##   IMAG <~ imag2      0.4473      0.1514    2.9549    0.0031 [ 0.1383; 0.7293 ] 
    ##   IMAG <~ imag3      0.6020      0.1413    4.2614    0.0000 [ 0.3154; 0.8564 ] 
    ##   EXPE <~ expe1      0.2946      0.1260    2.3383    0.0194 [ 0.0416; 0.5217 ] 
    ##   EXPE <~ expe2      0.6473      0.0883    7.3300    0.0000 [ 0.4589; 0.8059 ] 
    ##   EXPE <~ expe3      0.2374      0.0993    2.3910    0.0168 [ 0.0543; 0.4311 ] 
    ##   QUAL <~ qual1      0.2370      0.0930    2.5492    0.0108 [ 0.0667; 0.4115 ] 
    ##   QUAL <~ qual2      0.4712      0.0765    6.1594    0.0000 [ 0.3053; 0.6072 ] 
    ##   QUAL <~ qual3      0.1831      0.0875    2.0912    0.0365 [ 0.0144; 0.3656 ] 
    ##   QUAL <~ qual4      0.1037      0.0597    1.7361    0.0825 [-0.0120; 0.2297 ] 
    ##   QUAL <~ qual5      0.2049      0.0656    3.1244    0.0018 [ 0.0795; 0.3371 ] 
    ##   VAL <~ val1        0.7163      0.0987    7.2545    0.0000 [ 0.5044; 0.8928 ] 
    ##   VAL <~ val2        0.2202      0.0941    2.3399    0.0193 [ 0.0390; 0.4080 ] 
    ##   VAL <~ val3        0.2082      0.0628    3.3134    0.0009 [ 0.0818; 0.3298 ] 
    ##   SAT <~ sat1        0.3209      0.0149   21.5022    0.0000 [ 0.2936; 0.3511 ] 
    ##   SAT <~ sat2        0.3059      0.0144   21.1972    0.0000 [ 0.2813; 0.3374 ] 
    ##   SAT <~ sat3        0.2474      0.0111   22.1904    0.0000 [ 0.2243; 0.2690 ] 
    ##   SAT <~ sat4        0.2692      0.0113   23.7985    0.0000 [ 0.2477; 0.2905 ] 
    ##   LOY <~ loy1        0.3834      0.0253   15.1745    0.0000 [ 0.3360; 0.4332 ] 
    ##   LOY <~ loy2        0.2434      0.0302    8.0587    0.0000 [ 0.1764; 0.2951 ] 
    ##   LOY <~ loy3        0.3812      0.0295   12.9071    0.0000 [ 0.3268; 0.4399 ] 
    ##   LOY <~ loy4        0.2073      0.0387    5.3495    0.0000 [ 0.1239; 0.2766 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0666    9.6712    0.0000 [ 0.5048; 0.7551 ] 
    ##   imag1 ~~ imag3      0.5433      0.0720    7.5432    0.0000 [ 0.3925; 0.6732 ] 
    ##   imag2 ~~ imag3      0.7761      0.0398   19.4873    0.0000 [ 0.6982; 0.8483 ] 
    ##   expe1 ~~ expe2      0.5353      0.0603    8.8807    0.0000 [ 0.4201; 0.6460 ] 
    ##   expe1 ~~ expe3      0.4694      0.0625    7.5043    0.0000 [ 0.3405; 0.5786 ] 
    ##   expe2 ~~ expe3      0.5467      0.0633    8.6320    0.0000 [ 0.4180; 0.6588 ] 
    ##   qual1 ~~ qual2      0.6053      0.0601   10.0803    0.0000 [ 0.4777; 0.7032 ] 
    ##   qual1 ~~ qual3      0.5406      0.0629    8.5920    0.0000 [ 0.4082; 0.6542 ] 
    ##   qual1 ~~ qual4      0.5662      0.0704    8.0467    0.0000 [ 0.4190; 0.6927 ] 
    ##   qual1 ~~ qual5      0.5180      0.0690    7.5103    0.0000 [ 0.3788; 0.6417 ] 
    ##   qual2 ~~ qual3      0.6187      0.0574   10.7827    0.0000 [ 0.4982; 0.7214 ] 
    ##   qual2 ~~ qual4      0.6517      0.0607   10.7448    0.0000 [ 0.5355; 0.7575 ] 
    ##   qual2 ~~ qual5      0.6291      0.0581   10.8324    0.0000 [ 0.5087; 0.7313 ] 
    ##   qual3 ~~ qual4      0.4752      0.0634    7.4917    0.0000 [ 0.3366; 0.5879 ] 
    ##   qual3 ~~ qual5      0.5074      0.0648    7.8278    0.0000 [ 0.3679; 0.6328 ] 
    ##   qual4 ~~ qual5      0.6402      0.0605   10.5890    0.0000 [ 0.4970; 0.7381 ] 
    ##   val1 ~~ val2        0.6344      0.0556   11.4198    0.0000 [ 0.5182; 0.7325 ] 
    ##   val1 ~~ val3        0.4602      0.0699    6.5859    0.0000 [ 0.3147; 0.5937 ] 
    ##   val2 ~~ val3        0.6288      0.0633    9.9383    0.0000 [ 0.4981; 0.7509 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0678    6.9506    0.0000 [ 0.3345; 0.6033 ] 
    ##   QUAL ~ IMAG       0.3933      0.0630    6.2459    0.0000 [ 0.2689; 0.5178 ] 
    ##   QUAL ~ EXPE       0.8344      0.0233   35.7851    0.0000 [ 0.7845; 0.8764 ] 
    ##   VAL ~ IMAG        0.2974      0.0630    4.7185    0.0000 [ 0.1842; 0.4282 ] 
    ##   VAL ~ EXPE        0.6309      0.0525   12.0070    0.0000 [ 0.5208; 0.7276 ] 
    ##   VAL ~ QUAL        0.7013      0.0872    8.0439    0.0000 [ 0.5318; 0.8690 ] 
    ##   SAT ~ IMAG        0.4807      0.0690    6.9620    0.0000 [ 0.3435; 0.6133 ] 
    ##   SAT ~ EXPE        0.5001      0.0579    8.6395    0.0000 [ 0.3795; 0.6094 ] 
    ##   SAT ~ QUAL        0.5911      0.0933    6.3385    0.0000 [ 0.4191; 0.7698 ] 
    ##   SAT ~ VAL         0.5270      0.0915    5.7561    0.0000 [ 0.3393; 0.6936 ] 
    ##   LOY ~ IMAG        0.4840      0.0718    6.7363    0.0000 [ 0.3345; 0.6147 ] 
    ##   LOY ~ EXPE        0.3142      0.0555    5.6663    0.0000 [ 0.2140; 0.4231 ] 
    ##   LOY ~ QUAL        0.3714      0.0826    4.4954    0.0000 [ 0.2363; 0.5558 ] 
    ##   LOY ~ VAL         0.3311      0.0775    4.2693    0.0000 [ 0.1855; 0.5003 ] 
    ##   LOY ~ SAT         0.6283      0.0807    7.7822    0.0000 [ 0.4803; 0.7882 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0630    6.2459    0.0000 [ 0.2689; 0.5178 ] 
    ##   VAL ~ IMAG           0.2974      0.0630    4.7185    0.0000 [ 0.1842; 0.4282 ] 
    ##   VAL ~ EXPE           0.5852      0.0741    7.8986    0.0000 [ 0.4401; 0.7363 ] 
    ##   SAT ~ IMAG           0.2357      0.0493    4.7836    0.0000 [ 0.1464; 0.3355 ] 
    ##   SAT ~ EXPE           0.5173      0.0663    7.8068    0.0000 [ 0.3898; 0.6384 ] 
    ##   SAT ~ QUAL           0.3696      0.0667    5.5435    0.0000 [ 0.2360; 0.5016 ] 
    ##   LOY ~ IMAG           0.3020      0.0556    5.4327    0.0000 [ 0.2058; 0.4225 ] 
    ##   LOY ~ EXPE           0.3142      0.0555    5.6663    0.0000 [ 0.2140; 0.4231 ] 
    ##   LOY ~ QUAL           0.3714      0.0826    4.4954    0.0000 [ 0.2363; 0.5558 ] 
    ##   LOY ~ VAL            0.3311      0.0775    4.2693    0.0000 [ 0.1855; 0.5003 ] 
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
Windows as many separate R instances are opened in the backround as
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
