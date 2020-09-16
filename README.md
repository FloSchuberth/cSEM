
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
    ##  dG                      0.6493      0.3242  
    ##  SRMR                    0.0940      0.0519  
    ##  dL                      2.2340      0.6809  
    ##  dML                     2.9219      1.6079  
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
    ##  Out of 499 bootstrap replications 482 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -477201325
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
    ##   expe1         1.4547         1.5685       1.9070         2.0943       0.0541
    ##   expe2         1.4144         1.4807       1.9368         2.0281       0.1980
    ##   expe3         1.6361         1.7319       2.1320         2.2286       0.1237
    ##   qual1         1.4762         1.5456       1.9282         2.0609       0.1153
    ##   qual2         1.5805         1.5388       2.0427         2.0619       0.2156
    ##   qual3         1.7355         1.7294       2.2285         2.2864       0.1169
    ##   qual4         1.2355         1.1982       1.5993         1.6302       0.2319
    ##   qual5         1.5072         1.4995       1.9369         1.9533       0.1971
    ##   val1          1.4485         1.3667       1.8730         1.7683       0.2483
    ##   val2          1.2257         1.2059       1.6496         1.7142       0.1714
    ##   val3          1.4812         1.3804       1.9698         1.9346       0.1472
    ##   sat1          1.2469         1.2361       1.6446         1.6210       0.3403
    ##   sat2          1.2320         1.1951       1.6405         1.6257       0.3086
    ##   sat3          1.3417         1.2781       1.6744         1.7216       0.2094
    ##   sat4          1.3212         1.2648       1.6723         1.6392       0.2749
    ##   loy1          1.6931         1.6584       2.2336         2.2236       0.2683
    ##   loy2          1.4855         1.4725       1.9128         1.9792       0.1304
    ##   loy3          1.7023         1.6670       2.2784         2.2679       0.2724
    ##   loy4          1.6858         1.6689       2.1741         2.3026       0.0880
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
    ##  Number of admissible results     = 488
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 1912439257
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
    ##   EXPE ~ IMAG      0.4714      0.0653    7.2221    0.0000 [ 0.3436; 0.5969 ] 
    ##   QUAL ~ EXPE      0.8344      0.0230   36.3345    0.0000 [ 0.7880; 0.8744 ] 
    ##   VAL ~ EXPE       0.0457      0.0835    0.5478    0.5838 [-0.0909; 0.2276 ] 
    ##   VAL ~ QUAL       0.7013      0.0817    8.5792    0.0000 [ 0.5237; 0.8513 ] 
    ##   SAT ~ IMAG       0.2450      0.0511    4.7974    0.0000 [ 0.1545; 0.3417 ] 
    ##   SAT ~ EXPE      -0.0172      0.0703   -0.2450    0.8064 [-0.1470; 0.1336 ] 
    ##   SAT ~ QUAL       0.2215      0.0924    2.3969    0.0165 [ 0.0674; 0.4125 ] 
    ##   SAT ~ VAL        0.5270      0.0822    6.4132    0.0000 [ 0.3538; 0.6724 ] 
    ##   LOY ~ IMAG       0.1819      0.0793    2.2951    0.0217 [ 0.0418; 0.3520 ] 
    ##   LOY ~ SAT        0.6283      0.0787    7.9811    0.0000 [ 0.4713; 0.7795 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1028    6.1330    0.0000 [ 0.4120; 0.8041 ] 
    ##   IMAG =~ imag2      0.9246      0.0428   21.5954    0.0000 [ 0.8116; 0.9782 ] 
    ##   IMAG =~ imag3      0.9577      0.0282   33.9190    0.0000 [ 0.8846; 0.9920 ] 
    ##   EXPE =~ expe1      0.7525      0.0767    9.8115    0.0000 [ 0.5662; 0.8657 ] 
    ##   EXPE =~ expe2      0.9348      0.0280   33.4198    0.0000 [ 0.8652; 0.9728 ] 
    ##   EXPE =~ expe3      0.7295      0.0712   10.2468    0.0000 [ 0.5642; 0.8493 ] 
    ##   QUAL =~ qual1      0.7861      0.0680   11.5640    0.0000 [ 0.6201; 0.8861 ] 
    ##   QUAL =~ qual2      0.9244      0.0210   43.9264    0.0000 [ 0.8751; 0.9576 ] 
    ##   QUAL =~ qual3      0.7560      0.0598   12.6411    0.0000 [ 0.6336; 0.8547 ] 
    ##   QUAL =~ qual4      0.7632      0.0519   14.7023    0.0000 [ 0.6481; 0.8545 ] 
    ##   QUAL =~ qual5      0.7834      0.0452   17.3208    0.0000 [ 0.6857; 0.8535 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.6264    0.0000 [ 0.8982; 0.9861 ] 
    ##   VAL =~ val2        0.8056      0.0664   12.1331    0.0000 [ 0.6510; 0.9039 ] 
    ##   VAL =~ val3        0.6763      0.0763    8.8587    0.0000 [ 0.5081; 0.8075 ] 
    ##   SAT =~ sat1        0.9243      0.0230   40.1495    0.0000 [ 0.8749; 0.9641 ] 
    ##   SAT =~ sat2        0.8813      0.0268   32.8578    0.0000 [ 0.8179; 0.9264 ] 
    ##   SAT =~ sat3        0.7127      0.0532   13.3964    0.0000 [ 0.6074; 0.8089 ] 
    ##   SAT =~ sat4        0.7756      0.0508   15.2745    0.0000 [ 0.6661; 0.8645 ] 
    ##   LOY =~ loy1        0.9097      0.0498   18.2515    0.0000 [ 0.7973; 0.9811 ] 
    ##   LOY =~ loy2        0.5775      0.0841    6.8664    0.0000 [ 0.3947; 0.7293 ] 
    ##   LOY =~ loy3        0.9043      0.0418   21.6282    0.0000 [ 0.8152; 0.9694 ] 
    ##   LOY =~ loy4        0.4917      0.0963    5.1076    0.0000 [ 0.3160; 0.6669 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1150    0.1360    0.8918 [-0.2075; 0.2435 ] 
    ##   IMAG <~ imag2      0.4473      0.1487    3.0076    0.0026 [ 0.1315; 0.7355 ] 
    ##   IMAG <~ imag3      0.6020      0.1418    4.2449    0.0000 [ 0.2873; 0.8664 ] 
    ##   EXPE <~ expe1      0.2946      0.1108    2.6595    0.0078 [ 0.0916; 0.5141 ] 
    ##   EXPE <~ expe2      0.6473      0.0849    7.6255    0.0000 [ 0.4660; 0.7878 ] 
    ##   EXPE <~ expe3      0.2374      0.0939    2.5281    0.0115 [ 0.0561; 0.4120 ] 
    ##   QUAL <~ qual1      0.2370      0.0857    2.7650    0.0057 [ 0.0765; 0.4262 ] 
    ##   QUAL <~ qual2      0.4712      0.0775    6.0783    0.0000 [ 0.3121; 0.6104 ] 
    ##   QUAL <~ qual3      0.1831      0.0806    2.2725    0.0231 [ 0.0238; 0.3304 ] 
    ##   QUAL <~ qual4      0.1037      0.0594    1.7453    0.0809 [-0.0080; 0.2179 ] 
    ##   QUAL <~ qual5      0.2049      0.0587    3.4925    0.0005 [ 0.0913; 0.3146 ] 
    ##   VAL <~ val1        0.7163      0.0996    7.1887    0.0000 [ 0.5130; 0.8864 ] 
    ##   VAL <~ val2        0.2202      0.0940    2.3417    0.0192 [ 0.0536; 0.4084 ] 
    ##   VAL <~ val3        0.2082      0.0598    3.4802    0.0005 [ 0.0982; 0.3307 ] 
    ##   SAT <~ sat1        0.3209      0.0153   20.9687    0.0000 [ 0.2963; 0.3537 ] 
    ##   SAT <~ sat2        0.3059      0.0139   22.0672    0.0000 [ 0.2826; 0.3357 ] 
    ##   SAT <~ sat3        0.2474      0.0112   21.9973    0.0000 [ 0.2253; 0.2704 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.1146    0.0000 [ 0.2455; 0.2934 ] 
    ##   LOY <~ loy1        0.3834      0.0262   14.6296    0.0000 [ 0.3300; 0.4310 ] 
    ##   LOY <~ loy2        0.2434      0.0291    8.3685    0.0000 [ 0.1786; 0.2915 ] 
    ##   LOY <~ loy3        0.3812      0.0281   13.5465    0.0000 [ 0.3292; 0.4394 ] 
    ##   LOY <~ loy4        0.2073      0.0354    5.8618    0.0000 [ 0.1401; 0.2736 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0670    9.6029    0.0000 [ 0.4988; 0.7645 ] 
    ##   imag1 ~~ imag3      0.5433      0.0730    7.4449    0.0000 [ 0.4052; 0.6748 ] 
    ##   imag2 ~~ imag3      0.7761      0.0404   19.1908    0.0000 [ 0.6901; 0.8458 ] 
    ##   expe1 ~~ expe2      0.5353      0.0599    8.9361    0.0000 [ 0.4140; 0.6431 ] 
    ##   expe1 ~~ expe3      0.4694      0.0599    7.8402    0.0000 [ 0.3479; 0.5805 ] 
    ##   expe2 ~~ expe3      0.5467      0.0631    8.6567    0.0000 [ 0.4182; 0.6588 ] 
    ##   qual1 ~~ qual2      0.6053      0.0578   10.4643    0.0000 [ 0.4780; 0.7067 ] 
    ##   qual1 ~~ qual3      0.5406      0.0617    8.7627    0.0000 [ 0.4084; 0.6493 ] 
    ##   qual1 ~~ qual4      0.5662      0.0681    8.3144    0.0000 [ 0.4234; 0.6928 ] 
    ##   qual1 ~~ qual5      0.5180      0.0669    7.7461    0.0000 [ 0.3787; 0.6457 ] 
    ##   qual2 ~~ qual3      0.6187      0.0555   11.1449    0.0000 [ 0.5036; 0.7147 ] 
    ##   qual2 ~~ qual4      0.6517      0.0607   10.7451    0.0000 [ 0.5161; 0.7585 ] 
    ##   qual2 ~~ qual5      0.6291      0.0553   11.3793    0.0000 [ 0.5043; 0.7264 ] 
    ##   qual3 ~~ qual4      0.4752      0.0638    7.4440    0.0000 [ 0.3556; 0.5909 ] 
    ##   qual3 ~~ qual5      0.5074      0.0632    8.0243    0.0000 [ 0.3775; 0.6158 ] 
    ##   qual4 ~~ qual5      0.6402      0.0576   11.1059    0.0000 [ 0.5129; 0.7397 ] 
    ##   val1 ~~ val2        0.6344      0.0526   12.0511    0.0000 [ 0.5213; 0.7273 ] 
    ##   val1 ~~ val3        0.4602      0.0728    6.3237    0.0000 [ 0.3116; 0.5953 ] 
    ##   val2 ~~ val3        0.6288      0.0626   10.0510    0.0000 [ 0.4946; 0.7339 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0653    7.2221    0.0000 [ 0.3436; 0.5969 ] 
    ##   QUAL ~ IMAG       0.3933      0.0602    6.5356    0.0000 [ 0.2851; 0.5136 ] 
    ##   QUAL ~ EXPE       0.8344      0.0230   36.3345    0.0000 [ 0.7880; 0.8744 ] 
    ##   VAL ~ IMAG        0.2974      0.0590    5.0366    0.0000 [ 0.1937; 0.4187 ] 
    ##   VAL ~ EXPE        0.6309      0.0489   12.8968    0.0000 [ 0.5367; 0.7270 ] 
    ##   VAL ~ QUAL        0.7013      0.0817    8.5792    0.0000 [ 0.5237; 0.8513 ] 
    ##   SAT ~ IMAG        0.4807      0.0627    7.6638    0.0000 [ 0.3667; 0.6014 ] 
    ##   SAT ~ EXPE        0.5001      0.0557    8.9708    0.0000 [ 0.3859; 0.6161 ] 
    ##   SAT ~ QUAL        0.5911      0.0872    6.7808    0.0000 [ 0.4205; 0.7687 ] 
    ##   SAT ~ VAL         0.5270      0.0822    6.4132    0.0000 [ 0.3538; 0.6724 ] 
    ##   LOY ~ IMAG        0.4840      0.0647    7.4840    0.0000 [ 0.3651; 0.6062 ] 
    ##   LOY ~ EXPE        0.3142      0.0531    5.9142    0.0000 [ 0.2171; 0.4287 ] 
    ##   LOY ~ QUAL        0.3714      0.0775    4.7931    0.0000 [ 0.2443; 0.5462 ] 
    ##   LOY ~ VAL         0.3311      0.0698    4.7429    0.0000 [ 0.1930; 0.4666 ] 
    ##   LOY ~ SAT         0.6283      0.0787    7.9811    0.0000 [ 0.4713; 0.7795 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0602    6.5356    0.0000 [ 0.2851; 0.5136 ] 
    ##   VAL ~ IMAG           0.2974      0.0590    5.0366    0.0000 [ 0.1937; 0.4187 ] 
    ##   VAL ~ EXPE           0.5852      0.0704    8.3070    0.0000 [ 0.4357; 0.7181 ] 
    ##   SAT ~ IMAG           0.2357      0.0477    4.9373    0.0000 [ 0.1529; 0.3360 ] 
    ##   SAT ~ EXPE           0.5173      0.0640    8.0858    0.0000 [ 0.3958; 0.6490 ] 
    ##   SAT ~ QUAL           0.3696      0.0608    6.0761    0.0000 [ 0.2413; 0.4901 ] 
    ##   LOY ~ IMAG           0.3020      0.0545    5.5431    0.0000 [ 0.2098; 0.4185 ] 
    ##   LOY ~ EXPE           0.3142      0.0531    5.9142    0.0000 [ 0.2171; 0.4287 ] 
    ##   LOY ~ QUAL           0.3714      0.0775    4.7931    0.0000 [ 0.2443; 0.5462 ] 
    ##   LOY ~ VAL            0.3311      0.0698    4.7429    0.0000 [ 0.1930; 0.4666 ] 
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
