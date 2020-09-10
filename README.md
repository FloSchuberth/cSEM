
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
    ##                              +------------------------------------------------------------------+
    ##                              |                                                                  |
    ##                              |   H0: The model-implied indicator covariance matrix equals the   |
    ##                              |   population indicator covariance matrix.                        |
    ##                              |                                                                  |
    ##                              +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3220  
    ##  SRMR                    0.0940      0.0521  
    ##  dL                      2.2340      0.6872  
    ##  dML                     2.9219      1.5823  
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
    ##  The seed used was: 2132921959
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
    ##   expe1         1.4560         1.5678       1.9087         2.0968       0.0545
    ##   expe2         1.4138         1.4790       1.9358         2.0266       0.2014
    ##   expe3         1.6359         1.7289       2.1300         2.2245       0.1250
    ##   qual1         1.4754         1.5430       1.9302         2.0601       0.1144
    ##   qual2         1.5798         1.5340       2.0393         2.0549       0.2191
    ##   qual3         1.7317         1.7235       2.2225         2.2771       0.1205
    ##   qual4         1.2334         1.1934       1.5955         1.6223       0.2353
    ##   qual5         1.5057         1.5000       1.9337         1.9517       0.1994
    ##   val1          1.4471         1.3632       1.8700         1.7655       0.2511
    ##   val2          1.2267         1.2059       1.6485         1.7152       0.1746
    ##   val3          1.4793         1.3808       1.9696         1.9364       0.1484
    ##   sat1          1.2459         1.2320       1.6454         1.6165       0.3399
    ##   sat2          1.2340         1.1966       1.6426         1.6265       0.3088
    ##   sat3          1.3428         1.2751       1.6742         1.7217       0.2102
    ##   sat4          1.3173         1.2597       1.6669         1.6333       0.2772
    ##   loy1          1.6903         1.6598       2.2329         2.2269       0.2702
    ##   loy2          1.4824         1.4701       1.9096         1.9773       0.1325
    ##   loy3          1.7029         1.6682       2.2810         2.2707       0.2720
    ##   loy4          1.6868         1.6661       2.1741         2.2989       0.0884
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
    ##  Number of admissible results     = 487
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 794470647
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
    ##   EXPE ~ IMAG      0.4714      0.0682    6.9071    0.0000 [ 0.3341; 0.5995 ] 
    ##   QUAL ~ EXPE      0.8344      0.0243   34.3907    0.0000 [ 0.7812; 0.8776 ] 
    ##   VAL ~ EXPE       0.0457      0.0864    0.5290    0.5968 [-0.1188; 0.2115 ] 
    ##   VAL ~ QUAL       0.7013      0.0807    8.6873    0.0000 [ 0.5409; 0.8519 ] 
    ##   SAT ~ IMAG       0.2450      0.0569    4.3023    0.0000 [ 0.1372; 0.3686 ] 
    ##   SAT ~ EXPE      -0.0172      0.0703   -0.2453    0.8063 [-0.1618; 0.1121 ] 
    ##   SAT ~ QUAL       0.2215      0.1036    2.1390    0.0324 [ 0.0471; 0.4309 ] 
    ##   SAT ~ VAL        0.5270      0.0890    5.9184    0.0000 [ 0.3296; 0.6852 ] 
    ##   LOY ~ IMAG       0.1819      0.0777    2.3403    0.0193 [ 0.0473; 0.3450 ] 
    ##   LOY ~ SAT        0.6283      0.0810    7.7561    0.0000 [ 0.4523; 0.7801 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1024    6.1563    0.0000 [ 0.4147; 0.8174 ] 
    ##   IMAG =~ imag2      0.9246      0.0396   23.3210    0.0000 [ 0.8334; 0.9807 ] 
    ##   IMAG =~ imag3      0.9577      0.0298   32.1116    0.0000 [ 0.8763; 0.9917 ] 
    ##   EXPE =~ expe1      0.7525      0.0781    9.6412    0.0000 [ 0.5662; 0.8671 ] 
    ##   EXPE =~ expe2      0.9348      0.0286   32.7036    0.0000 [ 0.8568; 0.9724 ] 
    ##   EXPE =~ expe3      0.7295      0.0739    9.8667    0.0000 [ 0.5695; 0.8626 ] 
    ##   QUAL =~ qual1      0.7861      0.0675   11.6461    0.0000 [ 0.6314; 0.8867 ] 
    ##   QUAL =~ qual2      0.9244      0.0227   40.7635    0.0000 [ 0.8726; 0.9593 ] 
    ##   QUAL =~ qual3      0.7560      0.0598   12.6398    0.0000 [ 0.6248; 0.8628 ] 
    ##   QUAL =~ qual4      0.7632      0.0554   13.7648    0.0000 [ 0.6372; 0.8574 ] 
    ##   QUAL =~ qual5      0.7834      0.0456   17.1696    0.0000 [ 0.6719; 0.8525 ] 
    ##   VAL =~ val1        0.9518      0.0241   39.5580    0.0000 [ 0.8947; 0.9861 ] 
    ##   VAL =~ val2        0.8056      0.0650   12.3882    0.0000 [ 0.6550; 0.9092 ] 
    ##   VAL =~ val3        0.6763      0.0758    8.9259    0.0000 [ 0.5170; 0.8121 ] 
    ##   SAT =~ sat1        0.9243      0.0228   40.6176    0.0000 [ 0.8752; 0.9666 ] 
    ##   SAT =~ sat2        0.8813      0.0282   31.2774    0.0000 [ 0.8189; 0.9245 ] 
    ##   SAT =~ sat3        0.7127      0.0545   13.0792    0.0000 [ 0.6029; 0.8152 ] 
    ##   SAT =~ sat4        0.7756      0.0547   14.1854    0.0000 [ 0.6531; 0.8676 ] 
    ##   LOY =~ loy1        0.9097      0.0507   17.9292    0.0000 [ 0.7927; 0.9890 ] 
    ##   LOY =~ loy2        0.5775      0.0839    6.8831    0.0000 [ 0.4031; 0.7199 ] 
    ##   LOY =~ loy3        0.9043      0.0401   22.5686    0.0000 [ 0.8153; 0.9711 ] 
    ##   LOY =~ loy4        0.4917      0.1019    4.8267    0.0000 [ 0.2825; 0.6862 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1121    0.1396    0.8890 [-0.1922; 0.2633 ] 
    ##   IMAG <~ imag2      0.4473      0.1562    2.8644    0.0042 [ 0.1390; 0.7359 ] 
    ##   IMAG <~ imag3      0.6020      0.1422    4.2348    0.0000 [ 0.2894; 0.8289 ] 
    ##   EXPE <~ expe1      0.2946      0.1131    2.6053    0.0092 [ 0.0619; 0.5000 ] 
    ##   EXPE <~ expe2      0.6473      0.0869    7.4517    0.0000 [ 0.4624; 0.7964 ] 
    ##   EXPE <~ expe3      0.2374      0.0914    2.5960    0.0094 [ 0.0692; 0.4235 ] 
    ##   QUAL <~ qual1      0.2370      0.0850    2.7874    0.0053 [ 0.0890; 0.4137 ] 
    ##   QUAL <~ qual2      0.4712      0.0759    6.2110    0.0000 [ 0.3286; 0.6242 ] 
    ##   QUAL <~ qual3      0.1831      0.0792    2.3128    0.0207 [ 0.0214; 0.3458 ] 
    ##   QUAL <~ qual4      0.1037      0.0616    1.6844    0.0921 [-0.0043; 0.2297 ] 
    ##   QUAL <~ qual5      0.2049      0.0613    3.3396    0.0008 [ 0.0719; 0.3096 ] 
    ##   VAL <~ val1        0.7163      0.0965    7.4210    0.0000 [ 0.5036; 0.8865 ] 
    ##   VAL <~ val2        0.2202      0.0946    2.3284    0.0199 [ 0.0483; 0.4206 ] 
    ##   VAL <~ val3        0.2082      0.0609    3.4173    0.0006 [ 0.0956; 0.3189 ] 
    ##   SAT <~ sat1        0.3209      0.0165   19.4167    0.0000 [ 0.2941; 0.3587 ] 
    ##   SAT <~ sat2        0.3059      0.0145   21.0498    0.0000 [ 0.2807; 0.3399 ] 
    ##   SAT <~ sat3        0.2474      0.0113   21.8224    0.0000 [ 0.2245; 0.2677 ] 
    ##   SAT <~ sat4        0.2692      0.0124   21.7190    0.0000 [ 0.2461; 0.2949 ] 
    ##   LOY <~ loy1        0.3834      0.0267   14.3828    0.0000 [ 0.3353; 0.4378 ] 
    ##   LOY <~ loy2        0.2434      0.0291    8.3733    0.0000 [ 0.1783; 0.2964 ] 
    ##   LOY <~ loy3        0.3812      0.0280   13.6193    0.0000 [ 0.3303; 0.4349 ] 
    ##   LOY <~ loy4        0.2073      0.0376    5.5194    0.0000 [ 0.1265; 0.2794 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0661    9.7314    0.0000 [ 0.4994; 0.7637 ] 
    ##   imag1 ~~ imag3      0.5433      0.0717    7.5810    0.0000 [ 0.4056; 0.6766 ] 
    ##   imag2 ~~ imag3      0.7761      0.0402   19.3065    0.0000 [ 0.6912; 0.8501 ] 
    ##   expe1 ~~ expe2      0.5353      0.0597    8.9720    0.0000 [ 0.4092; 0.6429 ] 
    ##   expe1 ~~ expe3      0.4694      0.0626    7.5029    0.0000 [ 0.3454; 0.5762 ] 
    ##   expe2 ~~ expe3      0.5467      0.0647    8.4548    0.0000 [ 0.4028; 0.6609 ] 
    ##   qual1 ~~ qual2      0.6053      0.0594   10.1856    0.0000 [ 0.4816; 0.7161 ] 
    ##   qual1 ~~ qual3      0.5406      0.0576    9.3880    0.0000 [ 0.4243; 0.6483 ] 
    ##   qual1 ~~ qual4      0.5662      0.0731    7.7408    0.0000 [ 0.4072; 0.6888 ] 
    ##   qual1 ~~ qual5      0.5180      0.0690    7.5047    0.0000 [ 0.3735; 0.6348 ] 
    ##   qual2 ~~ qual3      0.6187      0.0567   10.9014    0.0000 [ 0.4993; 0.7207 ] 
    ##   qual2 ~~ qual4      0.6517      0.0630   10.3404    0.0000 [ 0.5170; 0.7622 ] 
    ##   qual2 ~~ qual5      0.6291      0.0605   10.4051    0.0000 [ 0.5019; 0.7343 ] 
    ##   qual3 ~~ qual4      0.4752      0.0660    7.2026    0.0000 [ 0.3459; 0.6010 ] 
    ##   qual3 ~~ qual5      0.5074      0.0629    8.0620    0.0000 [ 0.3784; 0.6252 ] 
    ##   qual4 ~~ qual5      0.6402      0.0558   11.4752    0.0000 [ 0.5224; 0.7427 ] 
    ##   val1 ~~ val2        0.6344      0.0555   11.4338    0.0000 [ 0.5225; 0.7396 ] 
    ##   val1 ~~ val3        0.4602      0.0741    6.2086    0.0000 [ 0.3188; 0.6071 ] 
    ##   val2 ~~ val3        0.6288      0.0646    9.7350    0.0000 [ 0.5000; 0.7443 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0682    6.9071    0.0000 [ 0.3341; 0.5995 ] 
    ##   QUAL ~ IMAG       0.3933      0.0632    6.2228    0.0000 [ 0.2704; 0.5139 ] 
    ##   QUAL ~ EXPE       0.8344      0.0243   34.3907    0.0000 [ 0.7812; 0.8776 ] 
    ##   VAL ~ IMAG        0.2974      0.0631    4.7100    0.0000 [ 0.1781; 0.4277 ] 
    ##   VAL ~ EXPE        0.6309      0.0542   11.6333    0.0000 [ 0.5169; 0.7280 ] 
    ##   VAL ~ QUAL        0.7013      0.0807    8.6873    0.0000 [ 0.5409; 0.8519 ] 
    ##   SAT ~ IMAG        0.4807      0.0700    6.8673    0.0000 [ 0.3393; 0.6270 ] 
    ##   SAT ~ EXPE        0.5001      0.0582    8.5923    0.0000 [ 0.3916; 0.6112 ] 
    ##   SAT ~ QUAL        0.5911      0.0930    6.3581    0.0000 [ 0.4101; 0.7734 ] 
    ##   SAT ~ VAL         0.5270      0.0890    5.9184    0.0000 [ 0.3296; 0.6852 ] 
    ##   LOY ~ IMAG        0.4840      0.0664    7.2910    0.0000 [ 0.3661; 0.6211 ] 
    ##   LOY ~ EXPE        0.3142      0.0563    5.5845    0.0000 [ 0.2131; 0.4317 ] 
    ##   LOY ~ QUAL        0.3714      0.0815    4.5579    0.0000 [ 0.2244; 0.5540 ] 
    ##   LOY ~ VAL         0.3311      0.0749    4.4229    0.0000 [ 0.1964; 0.4834 ] 
    ##   LOY ~ SAT         0.6283      0.0810    7.7561    0.0000 [ 0.4523; 0.7801 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0632    6.2228    0.0000 [ 0.2704; 0.5139 ] 
    ##   VAL ~ IMAG           0.2974      0.0631    4.7100    0.0000 [ 0.1781; 0.4277 ] 
    ##   VAL ~ EXPE           0.5852      0.0685    8.5424    0.0000 [ 0.4496; 0.7190 ] 
    ##   SAT ~ IMAG           0.2357      0.0485    4.8590    0.0000 [ 0.1403; 0.3323 ] 
    ##   SAT ~ EXPE           0.5173      0.0673    7.6823    0.0000 [ 0.3781; 0.6586 ] 
    ##   SAT ~ QUAL           0.3696      0.0645    5.7253    0.0000 [ 0.2429; 0.4979 ] 
    ##   LOY ~ IMAG           0.3020      0.0585    5.1632    0.0000 [ 0.1992; 0.4237 ] 
    ##   LOY ~ EXPE           0.3142      0.0563    5.5845    0.0000 [ 0.2131; 0.4317 ] 
    ##   LOY ~ QUAL           0.3714      0.0815    4.5579    0.0000 [ 0.2244; 0.5540 ] 
    ##   LOY ~ VAL            0.3311      0.0749    4.4229    0.0000 [ 0.1964; 0.4834 ] 
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
