
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
    ##  dG                      0.6493      0.3265  
    ##  SRMR                    0.0940      0.0512  
    ##  dL                      2.2340      0.6633  
    ##  dML                     2.9219      1.6345  
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
    ##  The seed used was: -1531845184
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
    ##  IMAG                      0.0000   
    ##  VAL                       0.0000   
    ##  SAT                       0.0000   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      3.2985   
    ##  QUAL                      4.4151   
    ##  IMAG                      1.7280   
    ##  VAL                       2.6726   
    ##  SAT                       0.0000   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct    VIF value 
    ##  EXPE                      0.0000   
    ##  QUAL                      0.0000   
    ##  IMAG                      1.9345   
    ##  VAL                       0.0000   
    ##  SAT                       1.9345   
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
    ##   expe1         1.4553         1.5678       1.9101         2.0981       0.0531
    ##   expe2         1.4085         1.4775       1.9312         2.0257       0.2039
    ##   expe3         1.6303         1.7280       2.1245         2.2239       0.1271
    ##   qual1         1.4763         1.5470       1.9309         2.0654       0.1147
    ##   qual2         1.5761         1.5339       2.0363         2.0553       0.2211
    ##   qual3         1.7325         1.7274       2.2215         2.2835       0.1210
    ##   qual4         1.2338         1.2001       1.5966         1.6301       0.2340
    ##   qual5         1.5041         1.5011       1.9328         1.9499       0.1988
    ##   val1          1.4434         1.3643       1.8665         1.7654       0.2531
    ##   val2          1.2258         1.2065       1.6483         1.7145       0.1747
    ##   val3          1.4786         1.3802       1.9704         1.9375       0.1470
    ##   sat1          1.2459         1.2297       1.6440         1.6138       0.3413
    ##   sat2          1.2336         1.1954       1.6407         1.6237       0.3105
    ##   sat3          1.3397         1.2743       1.6711         1.7176       0.2137
    ##   sat4          1.3175         1.2584       1.6667         1.6314       0.2795
    ##   loy1          1.6892         1.6594       2.2325         2.2249       0.2694
    ##   loy2          1.4860         1.4719       1.9133         1.9783       0.1319
    ##   loy3          1.7011         1.6660       2.2803         2.2681       0.2713
    ##   loy4          1.6876         1.6736       2.1772         2.3097       0.0867
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
    ##  Number of admissible results     = 490
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 1597651558
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
    ##   EXPE ~ IMAG      0.4714      0.0684    6.8895    0.0000 [ 0.3489; 0.6028 ] 
    ##   QUAL ~ EXPE      0.8344      0.0235   35.4788    0.0000 [ 0.7811; 0.8728 ] 
    ##   VAL ~ EXPE       0.0457      0.0842    0.5429    0.5872 [-0.1072; 0.2134 ] 
    ##   VAL ~ QUAL       0.7013      0.0840    8.3500    0.0000 [ 0.5247; 0.8523 ] 
    ##   SAT ~ IMAG       0.2450      0.0562    4.3551    0.0000 [ 0.1324; 0.3529 ] 
    ##   SAT ~ EXPE      -0.0172      0.0730   -0.2361    0.8133 [-0.1676; 0.1142 ] 
    ##   SAT ~ QUAL       0.2215      0.1020    2.1710    0.0299 [ 0.0381; 0.4395 ] 
    ##   SAT ~ VAL        0.5270      0.0878    6.0028    0.0000 [ 0.3642; 0.6958 ] 
    ##   LOY ~ IMAG       0.1819      0.0769    2.3656    0.0180 [ 0.0316; 0.3289 ] 
    ##   LOY ~ SAT        0.6283      0.0754    8.3339    0.0000 [ 0.4842; 0.7795 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1032    6.1102    0.0000 [ 0.4227; 0.8202 ] 
    ##   IMAG =~ imag2      0.9246      0.0395   23.3906    0.0000 [ 0.8301; 0.9761 ] 
    ##   IMAG =~ imag3      0.9577      0.0285   33.6603    0.0000 [ 0.8808; 0.9920 ] 
    ##   EXPE =~ expe1      0.7525      0.0804    9.3543    0.0000 [ 0.5878; 0.8812 ] 
    ##   EXPE =~ expe2      0.9348      0.0293   31.9527    0.0000 [ 0.8651; 0.9712 ] 
    ##   EXPE =~ expe3      0.7295      0.0707   10.3168    0.0000 [ 0.5760; 0.8389 ] 
    ##   QUAL =~ qual1      0.7861      0.0709   11.0876    0.0000 [ 0.6168; 0.8945 ] 
    ##   QUAL =~ qual2      0.9244      0.0237   39.0298    0.0000 [ 0.8650; 0.9589 ] 
    ##   QUAL =~ qual3      0.7560      0.0604   12.5118    0.0000 [ 0.6191; 0.8433 ] 
    ##   QUAL =~ qual4      0.7632      0.0533   14.3275    0.0000 [ 0.6455; 0.8516 ] 
    ##   QUAL =~ qual5      0.7834      0.0497   15.7609    0.0000 [ 0.6867; 0.8738 ] 
    ##   VAL =~ val1        0.9518      0.0227   41.9967    0.0000 [ 0.8995; 0.9854 ] 
    ##   VAL =~ val2        0.8056      0.0636   12.6619    0.0000 [ 0.6492; 0.9023 ] 
    ##   VAL =~ val3        0.6763      0.0744    9.0905    0.0000 [ 0.5195; 0.8010 ] 
    ##   SAT =~ sat1        0.9243      0.0224   41.2400    0.0000 [ 0.8724; 0.9588 ] 
    ##   SAT =~ sat2        0.8813      0.0307   28.7036    0.0000 [ 0.8128; 0.9303 ] 
    ##   SAT =~ sat3        0.7127      0.0541   13.1712    0.0000 [ 0.5901; 0.8010 ] 
    ##   SAT =~ sat4        0.7756      0.0467   16.5970    0.0000 [ 0.6852; 0.8639 ] 
    ##   LOY =~ loy1        0.9097      0.0500   18.2083    0.0000 [ 0.7877; 0.9884 ] 
    ##   LOY =~ loy2        0.5775      0.0814    7.0981    0.0000 [ 0.3935; 0.7232 ] 
    ##   LOY =~ loy3        0.9043      0.0429   21.0944    0.0000 [ 0.8082; 0.9688 ] 
    ##   LOY =~ loy4        0.4917      0.0968    5.0803    0.0000 [ 0.3001; 0.6763 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1175    0.1331    0.8941 [-0.2164; 0.2402 ] 
    ##   IMAG <~ imag2      0.4473      0.1519    2.9455    0.0032 [ 0.1692; 0.7306 ] 
    ##   IMAG <~ imag3      0.6020      0.1375    4.3774    0.0000 [ 0.3110; 0.8534 ] 
    ##   EXPE <~ expe1      0.2946      0.1212    2.4314    0.0150 [ 0.0612; 0.5510 ] 
    ##   EXPE <~ expe2      0.6473      0.0839    7.7184    0.0000 [ 0.4751; 0.7898 ] 
    ##   EXPE <~ expe3      0.2374      0.0932    2.5469    0.0109 [ 0.0411; 0.4223 ] 
    ##   QUAL <~ qual1      0.2370      0.0951    2.4923    0.0127 [ 0.0665; 0.4426 ] 
    ##   QUAL <~ qual2      0.4712      0.0784    6.0104    0.0000 [ 0.3047; 0.6140 ] 
    ##   QUAL <~ qual3      0.1831      0.0769    2.3819    0.0172 [ 0.0163; 0.3225 ] 
    ##   QUAL <~ qual4      0.1037      0.0600    1.7290    0.0838 [-0.0215; 0.2286 ] 
    ##   QUAL <~ qual5      0.2049      0.0690    2.9696    0.0030 [ 0.0558; 0.3388 ] 
    ##   VAL <~ val1        0.7163      0.0935    7.6605    0.0000 [ 0.5323; 0.8769 ] 
    ##   VAL <~ val2        0.2202      0.0901    2.4434    0.0145 [ 0.0484; 0.3991 ] 
    ##   VAL <~ val3        0.2082      0.0582    3.5772    0.0003 [ 0.0876; 0.3154 ] 
    ##   SAT <~ sat1        0.3209      0.0159   20.1887    0.0000 [ 0.2914; 0.3559 ] 
    ##   SAT <~ sat2        0.3059      0.0130   23.5684    0.0000 [ 0.2835; 0.3331 ] 
    ##   SAT <~ sat3        0.2474      0.0113   21.8170    0.0000 [ 0.2233; 0.2680 ] 
    ##   SAT <~ sat4        0.2692      0.0117   23.0163    0.0000 [ 0.2461; 0.2922 ] 
    ##   LOY <~ loy1        0.3834      0.0251   15.2532    0.0000 [ 0.3360; 0.4337 ] 
    ##   LOY <~ loy2        0.2434      0.0294    8.2663    0.0000 [ 0.1747; 0.2930 ] 
    ##   LOY <~ loy3        0.3812      0.0270   14.1066    0.0000 [ 0.3311; 0.4313 ] 
    ##   LOY <~ loy4        0.2073      0.0357    5.8117    0.0000 [ 0.1346; 0.2748 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0714    9.0126    0.0000 [ 0.4882; 0.7608 ] 
    ##   imag1 ~~ imag3      0.5433      0.0722    7.5234    0.0000 [ 0.3913; 0.6811 ] 
    ##   imag2 ~~ imag3      0.7761      0.0381   20.3858    0.0000 [ 0.6969; 0.8386 ] 
    ##   expe1 ~~ expe2      0.5353      0.0633    8.4594    0.0000 [ 0.3994; 0.6510 ] 
    ##   expe1 ~~ expe3      0.4694      0.0645    7.2827    0.0000 [ 0.3432; 0.5883 ] 
    ##   expe2 ~~ expe3      0.5467      0.0597    9.1605    0.0000 [ 0.4223; 0.6513 ] 
    ##   qual1 ~~ qual2      0.6053      0.0606    9.9863    0.0000 [ 0.4673; 0.7137 ] 
    ##   qual1 ~~ qual3      0.5406      0.0605    8.9295    0.0000 [ 0.4117; 0.6467 ] 
    ##   qual1 ~~ qual4      0.5662      0.0675    8.3871    0.0000 [ 0.4300; 0.6947 ] 
    ##   qual1 ~~ qual5      0.5180      0.0696    7.4405    0.0000 [ 0.3770; 0.6452 ] 
    ##   qual2 ~~ qual3      0.6187      0.0531   11.6407    0.0000 [ 0.5069; 0.7069 ] 
    ##   qual2 ~~ qual4      0.6517      0.0605   10.7741    0.0000 [ 0.5189; 0.7509 ] 
    ##   qual2 ~~ qual5      0.6291      0.0584   10.7646    0.0000 [ 0.5091; 0.7338 ] 
    ##   qual3 ~~ qual4      0.4752      0.0638    7.4477    0.0000 [ 0.3369; 0.5905 ] 
    ##   qual3 ~~ qual5      0.5074      0.0611    8.3003    0.0000 [ 0.3892; 0.6225 ] 
    ##   qual4 ~~ qual5      0.6402      0.0573   11.1735    0.0000 [ 0.5263; 0.7357 ] 
    ##   val1 ~~ val2        0.6344      0.0541   11.7332    0.0000 [ 0.5224; 0.7316 ] 
    ##   val1 ~~ val3        0.4602      0.0692    6.6516    0.0000 [ 0.3279; 0.5947 ] 
    ##   val2 ~~ val3        0.6288      0.0640    9.8254    0.0000 [ 0.4937; 0.7495 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0684    6.8895    0.0000 [ 0.3489; 0.6028 ] 
    ##   QUAL ~ IMAG       0.3933      0.0633    6.2134    0.0000 [ 0.2816; 0.5171 ] 
    ##   QUAL ~ EXPE       0.8344      0.0235   35.4788    0.0000 [ 0.7811; 0.8728 ] 
    ##   VAL ~ IMAG        0.2974      0.0618    4.8134    0.0000 [ 0.1898; 0.4232 ] 
    ##   VAL ~ EXPE        0.6309      0.0494   12.7694    0.0000 [ 0.5240; 0.7147 ] 
    ##   VAL ~ QUAL        0.7013      0.0840    8.3500    0.0000 [ 0.5247; 0.8523 ] 
    ##   SAT ~ IMAG        0.4807      0.0683    7.0386    0.0000 [ 0.3473; 0.6145 ] 
    ##   SAT ~ EXPE        0.5001      0.0589    8.4913    0.0000 [ 0.3782; 0.6150 ] 
    ##   SAT ~ QUAL        0.5911      0.0946    6.2504    0.0000 [ 0.4026; 0.7713 ] 
    ##   SAT ~ VAL         0.5270      0.0878    6.0028    0.0000 [ 0.3642; 0.6958 ] 
    ##   LOY ~ IMAG        0.4840      0.0679    7.1251    0.0000 [ 0.3418; 0.6216 ] 
    ##   LOY ~ EXPE        0.3142      0.0538    5.8371    0.0000 [ 0.2261; 0.4367 ] 
    ##   LOY ~ QUAL        0.3714      0.0777    4.7797    0.0000 [ 0.2352; 0.5343 ] 
    ##   LOY ~ VAL         0.3311      0.0753    4.3958    0.0000 [ 0.2048; 0.4893 ] 
    ##   LOY ~ SAT         0.6283      0.0754    8.3339    0.0000 [ 0.4842; 0.7795 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0633    6.2134    0.0000 [ 0.2816; 0.5171 ] 
    ##   VAL ~ IMAG           0.2974      0.0618    4.8134    0.0000 [ 0.1898; 0.4232 ] 
    ##   VAL ~ EXPE           0.5852      0.0721    8.1182    0.0000 [ 0.4268; 0.7255 ] 
    ##   SAT ~ IMAG           0.2357      0.0483    4.8822    0.0000 [ 0.1523; 0.3357 ] 
    ##   SAT ~ EXPE           0.5173      0.0689    7.5065    0.0000 [ 0.3927; 0.6584 ] 
    ##   SAT ~ QUAL           0.3696      0.0642    5.7590    0.0000 [ 0.2471; 0.4846 ] 
    ##   LOY ~ IMAG           0.3020      0.0569    5.3124    0.0000 [ 0.2062; 0.4222 ] 
    ##   LOY ~ EXPE           0.3142      0.0538    5.8371    0.0000 [ 0.2261; 0.4367 ] 
    ##   LOY ~ QUAL           0.3714      0.0777    4.7797    0.0000 [ 0.2352; 0.5343 ] 
    ##   LOY ~ VAL            0.3311      0.0753    4.3958    0.0000 [ 0.2048; 0.4893 ] 
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
