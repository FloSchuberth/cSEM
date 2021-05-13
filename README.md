
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
    ##  dG                      0.6493      0.3307  
    ##  SRMR                    0.0940      0.0536  
    ##  dL                      2.2340      0.7269  
    ##  dML                     2.9219      1.6564  
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
    ##  Out of 499 bootstrap replications 481 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 808553412
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
    ##   expe1         1.4579         1.5723       1.9094         2.0994       0.0537
    ##   expe2         1.4108         1.4776       1.9328         2.0253       0.2016
    ##   expe3         1.6295         1.7223       2.1235         2.2166       0.1252
    ##   qual1         1.4760         1.5466       1.9283         2.0633       0.1157
    ##   qual2         1.5757         1.5329       2.0363         2.0547       0.2208
    ##   qual3         1.7305         1.7261       2.2220         2.2807       0.1198
    ##   qual4         1.2341         1.1991       1.5969         1.6316       0.2337
    ##   qual5         1.5051         1.4978       1.9360         1.9527       0.1974
    ##   val1          1.4470         1.3648       1.8700         1.7669       0.2507
    ##   val2          1.2262         1.2066       1.6481         1.7146       0.1735
    ##   val3          1.4816         1.3808       1.9681         1.9348       0.1478
    ##   sat1          1.2454         1.2325       1.6449         1.6194       0.3411
    ##   sat2          1.2335         1.1987       1.6409         1.6297       0.3096
    ##   sat3          1.3414         1.2761       1.6725         1.7212       0.2108
    ##   sat4          1.3194         1.2600       1.6700         1.6335       0.2769
    ##   loy1          1.6933         1.6581       2.2350         2.2238       0.2696
    ##   loy2          1.4821         1.4748       1.9104         1.9815       0.1326
    ##   loy3          1.7035         1.6671       2.2814         2.2686       0.2712
    ##   loy4          1.6915         1.6713       2.1807         2.3053       0.0861
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
    ##  Number of admissible results     = 478
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 2041383428
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
    ##   EXPE ~ IMAG      0.4714      0.0663    7.1049    0.0000 [ 0.3471; 0.6041 ] 
    ##   QUAL ~ EXPE      0.8344      0.0232   35.9447    0.0000 [ 0.7840; 0.8752 ] 
    ##   VAL ~ EXPE       0.0457      0.0835    0.5473    0.5841 [-0.1120; 0.1989 ] 
    ##   VAL ~ QUAL       0.7013      0.0798    8.7848    0.0000 [ 0.5359; 0.8563 ] 
    ##   SAT ~ IMAG       0.2450      0.0587    4.1759    0.0000 [ 0.1297; 0.3607 ] 
    ##   SAT ~ EXPE      -0.0172      0.0704   -0.2449    0.8066 [-0.1655; 0.1188 ] 
    ##   SAT ~ QUAL       0.2215      0.0914    2.4243    0.0153 [ 0.0613; 0.4155 ] 
    ##   SAT ~ VAL        0.5270      0.0828    6.3663    0.0000 [ 0.3712; 0.6824 ] 
    ##   LOY ~ IMAG       0.1819      0.0738    2.4665    0.0136 [ 0.0453; 0.3331 ] 
    ##   LOY ~ SAT        0.6283      0.0778    8.0721    0.0000 [ 0.4608; 0.7764 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1017    6.2025    0.0000 [ 0.4158; 0.8128 ] 
    ##   IMAG =~ imag2      0.9246      0.0404   22.8597    0.0000 [ 0.8253; 0.9768 ] 
    ##   IMAG =~ imag3      0.9577      0.0283   33.8654    0.0000 [ 0.8847; 0.9899 ] 
    ##   EXPE =~ expe1      0.7525      0.0801    9.3895    0.0000 [ 0.5831; 0.8800 ] 
    ##   EXPE =~ expe2      0.9348      0.0293   31.9506    0.0000 [ 0.8654; 0.9745 ] 
    ##   EXPE =~ expe3      0.7295      0.0704   10.3627    0.0000 [ 0.5625; 0.8441 ] 
    ##   QUAL =~ qual1      0.7861      0.0680   11.5541    0.0000 [ 0.6202; 0.8905 ] 
    ##   QUAL =~ qual2      0.9244      0.0216   42.8210    0.0000 [ 0.8713; 0.9554 ] 
    ##   QUAL =~ qual3      0.7560      0.0595   12.7113    0.0000 [ 0.6090; 0.8530 ] 
    ##   QUAL =~ qual4      0.7632      0.0516   14.7795    0.0000 [ 0.6561; 0.8502 ] 
    ##   QUAL =~ qual5      0.7834      0.0509   15.4057    0.0000 [ 0.6649; 0.8654 ] 
    ##   VAL =~ val1        0.9518      0.0220   43.2600    0.0000 [ 0.9035; 0.9848 ] 
    ##   VAL =~ val2        0.8056      0.0628   12.8227    0.0000 [ 0.6517; 0.9077 ] 
    ##   VAL =~ val3        0.6763      0.0728    9.2942    0.0000 [ 0.5366; 0.8019 ] 
    ##   SAT =~ sat1        0.9243      0.0231   40.0258    0.0000 [ 0.8726; 0.9649 ] 
    ##   SAT =~ sat2        0.8813      0.0308   28.6037    0.0000 [ 0.8108; 0.9296 ] 
    ##   SAT =~ sat3        0.7127      0.0501   14.2290    0.0000 [ 0.6049; 0.7993 ] 
    ##   SAT =~ sat4        0.7756      0.0489   15.8669    0.0000 [ 0.6736; 0.8640 ] 
    ##   LOY =~ loy1        0.9097      0.0478   19.0456    0.0000 [ 0.8021; 0.9883 ] 
    ##   LOY =~ loy2        0.5775      0.0822    7.0248    0.0000 [ 0.4190; 0.7249 ] 
    ##   LOY =~ loy3        0.9043      0.0413   21.8864    0.0000 [ 0.8140; 0.9788 ] 
    ##   LOY =~ loy4        0.4917      0.0952    5.1633    0.0000 [ 0.3090; 0.6806 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1205    0.1298    0.8967 [-0.2217; 0.2640 ] 
    ##   IMAG <~ imag2      0.4473      0.1449    3.0866    0.0020 [ 0.1560; 0.7148 ] 
    ##   IMAG <~ imag3      0.6020      0.1391    4.3284    0.0000 [ 0.3112; 0.8267 ] 
    ##   EXPE <~ expe1      0.2946      0.1222    2.4106    0.0159 [ 0.0529; 0.5267 ] 
    ##   EXPE <~ expe2      0.6473      0.0901    7.1843    0.0000 [ 0.4559; 0.8001 ] 
    ##   EXPE <~ expe3      0.2374      0.0900    2.6371    0.0084 [ 0.0503; 0.4127 ] 
    ##   QUAL <~ qual1      0.2370      0.0877    2.7028    0.0069 [ 0.0768; 0.4167 ] 
    ##   QUAL <~ qual2      0.4712      0.0749    6.2923    0.0000 [ 0.3183; 0.6169 ] 
    ##   QUAL <~ qual3      0.1831      0.0746    2.4546    0.0141 [ 0.0301; 0.3190 ] 
    ##   QUAL <~ qual4      0.1037      0.0623    1.6651    0.0959 [-0.0122; 0.2309 ] 
    ##   QUAL <~ qual5      0.2049      0.0673    3.0451    0.0023 [ 0.0491; 0.3234 ] 
    ##   VAL <~ val1        0.7163      0.0926    7.7379    0.0000 [ 0.5194; 0.8823 ] 
    ##   VAL <~ val2        0.2202      0.0907    2.4272    0.0152 [ 0.0552; 0.3951 ] 
    ##   VAL <~ val3        0.2082      0.0586    3.5542    0.0004 [ 0.0940; 0.3186 ] 
    ##   SAT <~ sat1        0.3209      0.0147   21.7783    0.0000 [ 0.2961; 0.3532 ] 
    ##   SAT <~ sat2        0.3059      0.0138   22.2332    0.0000 [ 0.2829; 0.3383 ] 
    ##   SAT <~ sat3        0.2474      0.0107   23.1441    0.0000 [ 0.2249; 0.2669 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.1596    0.0000 [ 0.2474; 0.2942 ] 
    ##   LOY <~ loy1        0.3834      0.0250   15.3151    0.0000 [ 0.3326; 0.4311 ] 
    ##   LOY <~ loy2        0.2434      0.0284    8.5798    0.0000 [ 0.1816; 0.2949 ] 
    ##   LOY <~ loy3        0.3812      0.0260   14.6496    0.0000 [ 0.3298; 0.4308 ] 
    ##   LOY <~ loy4        0.2073      0.0357    5.8071    0.0000 [ 0.1367; 0.2760 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0629   10.2281    0.0000 [ 0.5201; 0.7534 ] 
    ##   imag1 ~~ imag3      0.5433      0.0685    7.9276    0.0000 [ 0.4118; 0.6700 ] 
    ##   imag2 ~~ imag3      0.7761      0.0405   19.1428    0.0000 [ 0.6821; 0.8437 ] 
    ##   expe1 ~~ expe2      0.5353      0.0604    8.8565    0.0000 [ 0.4233; 0.6585 ] 
    ##   expe1 ~~ expe3      0.4694      0.0608    7.7262    0.0000 [ 0.3448; 0.5751 ] 
    ##   expe2 ~~ expe3      0.5467      0.0585    9.3474    0.0000 [ 0.4224; 0.6487 ] 
    ##   qual1 ~~ qual2      0.6053      0.0588   10.2917    0.0000 [ 0.4774; 0.7124 ] 
    ##   qual1 ~~ qual3      0.5406      0.0616    8.7786    0.0000 [ 0.4092; 0.6482 ] 
    ##   qual1 ~~ qual4      0.5662      0.0685    8.2599    0.0000 [ 0.4280; 0.6934 ] 
    ##   qual1 ~~ qual5      0.5180      0.0689    7.5216    0.0000 [ 0.3819; 0.6422 ] 
    ##   qual2 ~~ qual3      0.6187      0.0542   11.4122    0.0000 [ 0.5015; 0.7036 ] 
    ##   qual2 ~~ qual4      0.6517      0.0611   10.6738    0.0000 [ 0.5234; 0.7623 ] 
    ##   qual2 ~~ qual5      0.6291      0.0570   11.0361    0.0000 [ 0.5134; 0.7399 ] 
    ##   qual3 ~~ qual4      0.4752      0.0631    7.5355    0.0000 [ 0.3474; 0.6012 ] 
    ##   qual3 ~~ qual5      0.5074      0.0626    8.1067    0.0000 [ 0.3819; 0.6253 ] 
    ##   qual4 ~~ qual5      0.6402      0.0543   11.7793    0.0000 [ 0.5307; 0.7336 ] 
    ##   val1 ~~ val2        0.6344      0.0541   11.7250    0.0000 [ 0.5233; 0.7342 ] 
    ##   val1 ~~ val3        0.4602      0.0712    6.4643    0.0000 [ 0.3185; 0.5985 ] 
    ##   val2 ~~ val3        0.6288      0.0618   10.1826    0.0000 [ 0.5052; 0.7439 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0663    7.1049    0.0000 [ 0.3471; 0.6041 ] 
    ##   QUAL ~ IMAG       0.3933      0.0615    6.4000    0.0000 [ 0.2772; 0.5122 ] 
    ##   QUAL ~ EXPE       0.8344      0.0232   35.9447    0.0000 [ 0.7840; 0.8752 ] 
    ##   VAL ~ IMAG        0.2974      0.0604    4.9276    0.0000 [ 0.1879; 0.4265 ] 
    ##   VAL ~ EXPE        0.6309      0.0502   12.5678    0.0000 [ 0.5283; 0.7222 ] 
    ##   VAL ~ QUAL        0.7013      0.0798    8.7848    0.0000 [ 0.5359; 0.8563 ] 
    ##   SAT ~ IMAG        0.4807      0.0682    7.0465    0.0000 [ 0.3513; 0.6103 ] 
    ##   SAT ~ EXPE        0.5001      0.0559    8.9422    0.0000 [ 0.3903; 0.6019 ] 
    ##   SAT ~ QUAL        0.5911      0.0911    6.4860    0.0000 [ 0.4225; 0.7816 ] 
    ##   SAT ~ VAL         0.5270      0.0828    6.3663    0.0000 [ 0.3712; 0.6824 ] 
    ##   LOY ~ IMAG        0.4840      0.0634    7.6290    0.0000 [ 0.3653; 0.6114 ] 
    ##   LOY ~ EXPE        0.3142      0.0535    5.8750    0.0000 [ 0.2203; 0.4286 ] 
    ##   LOY ~ QUAL        0.3714      0.0802    4.6331    0.0000 [ 0.2365; 0.5438 ] 
    ##   LOY ~ VAL         0.3311      0.0710    4.6623    0.0000 [ 0.2072; 0.4803 ] 
    ##   LOY ~ SAT         0.6283      0.0778    8.0721    0.0000 [ 0.4608; 0.7764 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0615    6.4000    0.0000 [ 0.2772; 0.5122 ] 
    ##   VAL ~ IMAG           0.2974      0.0604    4.9276    0.0000 [ 0.1879; 0.4265 ] 
    ##   VAL ~ EXPE           0.5852      0.0694    8.4313    0.0000 [ 0.4416; 0.7173 ] 
    ##   SAT ~ IMAG           0.2357      0.0462    5.1049    0.0000 [ 0.1515; 0.3315 ] 
    ##   SAT ~ EXPE           0.5173      0.0690    7.5018    0.0000 [ 0.4025; 0.6618 ] 
    ##   SAT ~ QUAL           0.3696      0.0618    5.9824    0.0000 [ 0.2539; 0.4949 ] 
    ##   LOY ~ IMAG           0.3020      0.0576    5.2469    0.0000 [ 0.2006; 0.4268 ] 
    ##   LOY ~ EXPE           0.3142      0.0535    5.8750    0.0000 [ 0.2203; 0.4286 ] 
    ##   LOY ~ QUAL           0.3714      0.0802    4.6331    0.0000 [ 0.2365; 0.5438 ] 
    ##   LOY ~ VAL            0.3311      0.0710    4.6623    0.0000 [ 0.2072; 0.4803 ] 
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
