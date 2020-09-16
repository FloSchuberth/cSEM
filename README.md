
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![Build
Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)
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

<!-- ## Philosophy -->

<!-- - First and foremost: `cSEM` has a user-centered design!. "User-centered" mainly  -->

<!--   boils down to: `cSEM` is easy, i.e. intuitive to use by non-R experts!  -->

<!-- - Modern in a sense that the package integrates modern developments within  -->

<!--   the R community. This mainly includes ideas/recommendations/design choices that -->

<!--   fead into the packages of the [tidyverse](https://github.com/tidyverse/tidyverse). -->

<!-- - State of the art in a sense that we seek to quickly implement recent methodological -->

<!--   developments in composite-based SEM.  -->

## Basic usage

The basic usage is illustrated below.

<img src="man/figures/api.png" width="80%" style="display: block; margin: auto;" />

Usully, using `cSEM` is the same 3 step procedure:

> 1.  Pick a dataset and specify a model using [lavaan
>     syntax](http://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use `csem()`
> 3.  Apply one of the postestimation functions listed below on the
>     resulting object.

## Postestimation functions

There are five major postestimation verbs, three test family functions
and four do-family of function:

  - `assess()` : assess the model using common quality criteria
  - `infer()` : calculate common inferencial quantities (e.g., standard
    errors, confidence intervals)
  - `predict()` : predict endogenous indicator values
  - `summarize()` : summarize the results
  - `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

  - `testOMF()` : performs a test for overall model fit
  - `testMICOM()` : performs a test for composite measurement invariance
  - `testMGD()` : performs several tests to assess multi-group
    differences
  - `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently, three do functions are implemented:

  - `doIPMA()`: performs an importance-performance matrix analysis
  - `doNonlinearEffectsAnalysis()`: performs a nonlinear effects
    analysis such as floodlight and surface analysis
  - `doRedundancyAnalysis()`: performs a redundancy analysis

All functions require a `cSEMResults` object.

## Example

Models are defined using [lavaan
syntax](http://lavaan.ugent.be/tutorial/syntax1.html) with some slight
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
    ##                                       +------------------------------------------------------------------+
    ##                                       |                                                                  |
    ##                                       |   H0: The model-implied indicator covariance matrix equals the   |
    ##                                       |   population indicator covariance matrix.                        |
    ##                                       |                                                                  |
    ##                                       +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3143  
    ##  SRMR                    0.0940      0.0516  
    ##  dL                      2.2340      0.6740  
    ##  dML                     2.9219      1.6030  
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
    ##  Out of 499 bootstrap replications 472 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 1069960655
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
    ##  imag1      2.8359   
    ##  imag2      5.5535   
    ##  imag3      4.5088   
    ## 
    ##   Construct: 'EXPE'
    ## 
    ##  Weight    VIF value 
    ##  expe1      2.3551   
    ##  expe2      2.7116   
    ##  expe3      2.4116   
    ## 
    ##   Construct: 'QUAL'
    ## 
    ##  Weight    VIF value 
    ##  qual1      3.0835   
    ##  qual2      4.4376   
    ##  qual3      2.9575   
    ##  qual4      3.7341   
    ##  qual5      3.4566   
    ## 
    ##   Construct: 'VAL'
    ## 
    ##  Weight    VIF value 
    ##  val1       2.7725   
    ##  val2       3.8349   
    ##  val3       2.7307   
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
    ##   expe1         1.4555         1.5701       1.9097         2.0989       0.0537
    ##   expe2         1.4112         1.4782       1.9334         2.0290       0.2013
    ##   expe3         1.6315         1.7254       2.1267         2.2225       0.1238
    ##   qual1         1.4781         1.5534       1.9325         2.0736       0.1128
    ##   qual2         1.5780         1.5381       2.0381         2.0617       0.2185
    ##   qual3         1.7303         1.7256       2.2212         2.2800       0.1194
    ##   qual4         1.2325         1.2001       1.5945         1.6335       0.2342
    ##   qual5         1.5021         1.5023       1.9304         1.9534       0.1993
    ##   val1          1.4462         1.3674       1.8684         1.7684       0.2511
    ##   val2          1.2280         1.2075       1.6493         1.7166       0.1741
    ##   val3          1.4808         1.3794       1.9680         1.9327       0.1488
    ##   sat1          1.2460         1.2338       1.6440         1.6188       0.3405
    ##   sat2          1.2293         1.1940       1.6375         1.6257       0.3105
    ##   sat3          1.3398         1.2755       1.6711         1.7226       0.2121
    ##   sat4          1.3159         1.2605       1.6651         1.6332       0.2781
    ##   loy1          1.6887         1.6576       2.2295         2.2240       0.2700
    ##   loy2          1.4844         1.4738       1.9127         1.9802       0.1297
    ##   loy3          1.6989         1.6657       2.2764         2.2675       0.2715
    ##   loy4          1.6920         1.6737       2.1801         2.3048       0.0853
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

<!-- end list -->

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
    ##  Number of admissible results     = 486
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = -613951561
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
    ##   EXPE ~ IMAG      0.4714      0.0670    7.0345    0.0000 [ 0.3312; 0.5962 ] 
    ##   QUAL ~ EXPE      0.8344      0.0229   36.3622    0.0000 [ 0.7804; 0.8723 ] 
    ##   VAL ~ EXPE       0.0457      0.0880    0.5195    0.6034 [-0.1084; 0.2349 ] 
    ##   VAL ~ QUAL       0.7013      0.0837    8.3832    0.0000 [ 0.5190; 0.8596 ] 
    ##   SAT ~ IMAG       0.2450      0.0556    4.4037    0.0000 [ 0.1292; 0.3496 ] 
    ##   SAT ~ EXPE      -0.0172      0.0761   -0.2265    0.8208 [-0.1693; 0.1372 ] 
    ##   SAT ~ QUAL       0.2215      0.1052    2.1062    0.0352 [ 0.0352; 0.4398 ] 
    ##   SAT ~ VAL        0.5270      0.0830    6.3512    0.0000 [ 0.3677; 0.6932 ] 
    ##   LOY ~ IMAG       0.1819      0.0784    2.3195    0.0204 [ 0.0394; 0.3576 ] 
    ##   LOY ~ SAT        0.6283      0.0790    7.9553    0.0000 [ 0.4779; 0.7782 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0982    6.4210    0.0000 [ 0.4204; 0.7999 ] 
    ##   IMAG =~ imag2      0.9246      0.0403   22.9631    0.0000 [ 0.8252; 0.9723 ] 
    ##   IMAG =~ imag3      0.9577      0.0296   32.3592    0.0000 [ 0.8798; 0.9917 ] 
    ##   EXPE =~ expe1      0.7525      0.0769    9.7870    0.0000 [ 0.5695; 0.8665 ] 
    ##   EXPE =~ expe2      0.9348      0.0277   33.7041    0.0000 [ 0.8683; 0.9731 ] 
    ##   EXPE =~ expe3      0.7295      0.0722   10.1051    0.0000 [ 0.5624; 0.8439 ] 
    ##   QUAL =~ qual1      0.7861      0.0652   12.0548    0.0000 [ 0.6300; 0.8900 ] 
    ##   QUAL =~ qual2      0.9244      0.0228   40.5945    0.0000 [ 0.8716; 0.9575 ] 
    ##   QUAL =~ qual3      0.7560      0.0614   12.3061    0.0000 [ 0.6164; 0.8585 ] 
    ##   QUAL =~ qual4      0.7632      0.0535   14.2664    0.0000 [ 0.6444; 0.8505 ] 
    ##   QUAL =~ qual5      0.7834      0.0485   16.1379    0.0000 [ 0.6764; 0.8639 ] 
    ##   VAL =~ val1        0.9518      0.0235   40.5762    0.0000 [ 0.8975; 0.9852 ] 
    ##   VAL =~ val2        0.8056      0.0642   12.5501    0.0000 [ 0.6685; 0.9070 ] 
    ##   VAL =~ val3        0.6763      0.0707    9.5687    0.0000 [ 0.5363; 0.7985 ] 
    ##   SAT =~ sat1        0.9243      0.0216   42.7368    0.0000 [ 0.8758; 0.9605 ] 
    ##   SAT =~ sat2        0.8813      0.0280   31.5230    0.0000 [ 0.8219; 0.9285 ] 
    ##   SAT =~ sat3        0.7127      0.0519   13.7198    0.0000 [ 0.6111; 0.8045 ] 
    ##   SAT =~ sat4        0.7756      0.0461   16.8338    0.0000 [ 0.6685; 0.8579 ] 
    ##   LOY =~ loy1        0.9097      0.0488   18.6419    0.0000 [ 0.7942; 0.9802 ] 
    ##   LOY =~ loy2        0.5775      0.0859    6.7258    0.0000 [ 0.4033; 0.7306 ] 
    ##   LOY =~ loy3        0.9043      0.0440   20.5354    0.0000 [ 0.8061; 0.9794 ] 
    ##   LOY =~ loy4        0.4917      0.0979    5.0220    0.0000 [ 0.3136; 0.6900 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1203    0.1301    0.8965 [-0.2363; 0.2561 ] 
    ##   IMAG <~ imag2      0.4473      0.1453    3.0785    0.0021 [ 0.1587; 0.7138 ] 
    ##   IMAG <~ imag3      0.6020      0.1396    4.3134    0.0000 [ 0.3214; 0.8474 ] 
    ##   EXPE <~ expe1      0.2946      0.1141    2.5831    0.0098 [ 0.0623; 0.5137 ] 
    ##   EXPE <~ expe2      0.6473      0.0822    7.8772    0.0000 [ 0.4694; 0.7901 ] 
    ##   EXPE <~ expe3      0.2374      0.0909    2.6102    0.0090 [ 0.0598; 0.4076 ] 
    ##   QUAL <~ qual1      0.2370      0.0903    2.6259    0.0086 [ 0.0710; 0.4288 ] 
    ##   QUAL <~ qual2      0.4712      0.0787    5.9862    0.0000 [ 0.3188; 0.6038 ] 
    ##   QUAL <~ qual3      0.1831      0.0772    2.3706    0.0178 [ 0.0164; 0.3185 ] 
    ##   QUAL <~ qual4      0.1037      0.0604    1.7180    0.0858 [-0.0142; 0.2231 ] 
    ##   QUAL <~ qual5      0.2049      0.0614    3.3352    0.0009 [ 0.0708; 0.3165 ] 
    ##   VAL <~ val1        0.7163      0.0955    7.5023    0.0000 [ 0.5067; 0.8739 ] 
    ##   VAL <~ val2        0.2202      0.0927    2.3745    0.0176 [ 0.0493; 0.4144 ] 
    ##   VAL <~ val3        0.2082      0.0615    3.3839    0.0007 [ 0.0849; 0.3279 ] 
    ##   SAT <~ sat1        0.3209      0.0143   22.3644    0.0000 [ 0.2961; 0.3521 ] 
    ##   SAT <~ sat2        0.3059      0.0135   22.6694    0.0000 [ 0.2836; 0.3368 ] 
    ##   SAT <~ sat3        0.2474      0.0108   22.9732    0.0000 [ 0.2251; 0.2683 ] 
    ##   SAT <~ sat4        0.2692      0.0114   23.6275    0.0000 [ 0.2463; 0.2913 ] 
    ##   LOY <~ loy1        0.3834      0.0261   14.7146    0.0000 [ 0.3320; 0.4369 ] 
    ##   LOY <~ loy2        0.2434      0.0306    7.9619    0.0000 [ 0.1765; 0.2932 ] 
    ##   LOY <~ loy3        0.3812      0.0268   14.2389    0.0000 [ 0.3225; 0.4328 ] 
    ##   LOY <~ loy4        0.2073      0.0364    5.6927    0.0000 [ 0.1415; 0.2782 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0649    9.9226    0.0000 [ 0.4961; 0.7486 ] 
    ##   imag1 ~~ imag3      0.5433      0.0699    7.7761    0.0000 [ 0.3972; 0.6708 ] 
    ##   imag2 ~~ imag3      0.7761      0.0397   19.5605    0.0000 [ 0.6863; 0.8431 ] 
    ##   expe1 ~~ expe2      0.5353      0.0602    8.8915    0.0000 [ 0.4187; 0.6439 ] 
    ##   expe1 ~~ expe3      0.4694      0.0614    7.6442    0.0000 [ 0.3336; 0.5786 ] 
    ##   expe2 ~~ expe3      0.5467      0.0560    9.7545    0.0000 [ 0.4384; 0.6508 ] 
    ##   qual1 ~~ qual2      0.6053      0.0546   11.0796    0.0000 [ 0.5040; 0.7070 ] 
    ##   qual1 ~~ qual3      0.5406      0.0573    9.4349    0.0000 [ 0.4290; 0.6465 ] 
    ##   qual1 ~~ qual4      0.5662      0.0671    8.4416    0.0000 [ 0.4334; 0.6911 ] 
    ##   qual1 ~~ qual5      0.5180      0.0678    7.6462    0.0000 [ 0.3977; 0.6487 ] 
    ##   qual2 ~~ qual3      0.6187      0.0549   11.2770    0.0000 [ 0.5102; 0.7142 ] 
    ##   qual2 ~~ qual4      0.6517      0.0613   10.6301    0.0000 [ 0.5201; 0.7548 ] 
    ##   qual2 ~~ qual5      0.6291      0.0564   11.1487    0.0000 [ 0.5135; 0.7316 ] 
    ##   qual3 ~~ qual4      0.4752      0.0643    7.3911    0.0000 [ 0.3397; 0.5962 ] 
    ##   qual3 ~~ qual5      0.5074      0.0655    7.7515    0.0000 [ 0.3715; 0.6241 ] 
    ##   qual4 ~~ qual5      0.6402      0.0563   11.3706    0.0000 [ 0.5190; 0.7351 ] 
    ##   val1 ~~ val2        0.6344      0.0562   11.2903    0.0000 [ 0.5208; 0.7336 ] 
    ##   val1 ~~ val3        0.4602      0.0678    6.7869    0.0000 [ 0.3332; 0.5872 ] 
    ##   val2 ~~ val3        0.6288      0.0594   10.5777    0.0000 [ 0.5157; 0.7443 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0670    7.0345    0.0000 [ 0.3312; 0.5962 ] 
    ##   QUAL ~ IMAG       0.3933      0.0616    6.3838    0.0000 [ 0.2687; 0.5090 ] 
    ##   QUAL ~ EXPE       0.8344      0.0229   36.3622    0.0000 [ 0.7804; 0.8723 ] 
    ##   VAL ~ IMAG        0.2974      0.0607    4.8954    0.0000 [ 0.1868; 0.4235 ] 
    ##   VAL ~ EXPE        0.6309      0.0493   12.8027    0.0000 [ 0.5381; 0.7255 ] 
    ##   VAL ~ QUAL        0.7013      0.0837    8.3832    0.0000 [ 0.5190; 0.8596 ] 
    ##   SAT ~ IMAG        0.4807      0.0660    7.2823    0.0000 [ 0.3367; 0.6078 ] 
    ##   SAT ~ EXPE        0.5001      0.0557    8.9774    0.0000 [ 0.3978; 0.6112 ] 
    ##   SAT ~ QUAL        0.5911      0.1026    5.7592    0.0000 [ 0.3767; 0.7948 ] 
    ##   SAT ~ VAL         0.5270      0.0830    6.3512    0.0000 [ 0.3677; 0.6932 ] 
    ##   LOY ~ IMAG        0.4840      0.0663    7.3008    0.0000 [ 0.3582; 0.6171 ] 
    ##   LOY ~ EXPE        0.3142      0.0531    5.9178    0.0000 [ 0.2192; 0.4299 ] 
    ##   LOY ~ QUAL        0.3714      0.0849    4.3719    0.0000 [ 0.2086; 0.5412 ] 
    ##   LOY ~ VAL         0.3311      0.0709    4.6701    0.0000 [ 0.1967; 0.4677 ] 
    ##   LOY ~ SAT         0.6283      0.0790    7.9553    0.0000 [ 0.4779; 0.7782 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0616    6.3838    0.0000 [ 0.2687; 0.5090 ] 
    ##   VAL ~ IMAG           0.2974      0.0607    4.8954    0.0000 [ 0.1868; 0.4235 ] 
    ##   VAL ~ EXPE           0.5852      0.0710    8.2376    0.0000 [ 0.4333; 0.7159 ] 
    ##   SAT ~ IMAG           0.2357      0.0485    4.8640    0.0000 [ 0.1452; 0.3346 ] 
    ##   SAT ~ EXPE           0.5173      0.0724    7.1460    0.0000 [ 0.4008; 0.6699 ] 
    ##   SAT ~ QUAL           0.3696      0.0633    5.8411    0.0000 [ 0.2430; 0.4977 ] 
    ##   LOY ~ IMAG           0.3020      0.0560    5.3934    0.0000 [ 0.1966; 0.4206 ] 
    ##   LOY ~ EXPE           0.3142      0.0531    5.9178    0.0000 [ 0.2192; 0.4299 ] 
    ##   LOY ~ QUAL           0.3714      0.0849    4.3719    0.0000 [ 0.2086; 0.5412 ] 
    ##   LOY ~ VAL            0.3311      0.0709    4.6701    0.0000 [ 0.1967; 0.4677 ] 
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
