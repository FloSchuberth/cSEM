
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cSEM: Composite-based SEM <img src='man/figures/cSEMsticker.svg' align="right" height="200" /></a>

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![Build
Status](https://travis-ci.com/M-E-Rademaker/cSEM.svg?branch=master)](https://travis-ci.com/M-E-Rademaker/cSEM)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)

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

There are five major postestimation verbs, four test family functions
and three do-family of function:

  - `assess()` : assess the model using common quality criteria
  - `infer()` : calculate common inferencial quantities (e.g., standard
    errors, confidence intervals)
  - `predict()` : predict endogenous indicator values
  - `summarize()` : summarize the results
  - `verify()` : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently the
following tests are implemented:

  - `testOMF()` : performs a test for overall model fit
  - `testMICOM()` : performs a test for composite measurement invariance
  - `testMGD()` : performs several tests to assess multi-group
    differences
  - `testHausman()` : performs the regression-based Hausman test to test
    for endogeneity

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently two do functions are implemented:

  - `doFloodlightAnalysis()`: performs a floodlight analysis
  - `doSurfaceAnalysis()`: performs a surface analysis
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
    ##  Inner weighting scheme           = path
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
    ##  IMAG  Composite      First order   modeB
    ##  EXPE  Composite      First order   modeB
    ##  QUAL  Composite      First order   modeB
    ##  VAL   Composite      First order   modeB
    ##  SAT   Common factor  First order   modeA
    ##  LOY   Common factor  First order   modeA
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
    ##   Weights          Estimate  Std. error   t-stat.   p-value
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
testOMF(res, .verbose = FALSE)
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
    ##  dG                      0.6493      0.3160  
    ##  SRMR                    0.0940      0.0524  
    ##  dL                      2.2340      0.6945  
    ##  dML                     2.9219      1.5723  
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
    ## Additonal information:
    ## 
    ##  Out of 499 bootstrap replications 478 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 884410290
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
    ## 
    ## Variance accounted for (VAF):
    ## =============================
    ##   Effects        Estimate  Std. error   t-stat.   p-value
    ##   QUAL ~ IMAG      1.0000          NA        NA        NA
    ##   VAL ~ IMAG       1.0000          NA        NA        NA
    ##   VAL ~ EXPE       0.9275          NA        NA        NA
    ##   SAT ~ IMAG       0.4904          NA        NA        NA
    ##   SAT ~ EXPE       1.0345          NA        NA        NA
    ##   SAT ~ QUAL       0.6252          NA        NA        NA
    ##   LOY ~ IMAG       0.6241          NA        NA        NA
    ##   LOY ~ EXPE       1.0000          NA        NA        NA
    ##   LOY ~ QUAL       1.0000          NA        NA        NA
    ##   LOY ~ VAL        1.0000          NA        NA        NA
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
    ##   expe1         1.4546         1.5691       1.9075         2.0978       0.0538
    ##   expe2         1.4149         1.4833       1.9392         2.0335       0.1959
    ##   expe3         1.6326         1.7257       2.1280         2.2219       0.1229
    ##   qual1         1.4750         1.5484       1.9282         2.0677       0.1147
    ##   qual2         1.5801         1.5392       2.0426         2.0611       0.2143
    ##   qual3         1.7338         1.7282       2.2247         2.2831       0.1168
    ##   qual4         1.2347         1.1996       1.5973         1.6324       0.2324
    ##   qual5         1.5039         1.5009       1.9345         1.9542       0.1960
    ##   val1          1.4465         1.3679       1.8730         1.7728       0.2467
    ##   val2          1.2254         1.2071       1.6485         1.7151       0.1733
    ##   val3          1.4807         1.3804       1.9677         1.9324       0.1482
    ##   sat1          1.2473         1.2355       1.6474         1.6231       0.3370
    ##   sat2          1.2308         1.1957       1.6403         1.6279       0.3085
    ##   sat3          1.3405         1.2761       1.6730         1.7221       0.2095
    ##   sat4          1.3182         1.2601       1.6687         1.6356       0.2759
    ##   loy1          1.6917         1.6630       2.2318         2.2271       0.2684
    ##   loy2          1.4848         1.4773       1.9111         1.9845       0.1313
    ##   loy3          1.7027         1.6703       2.2789         2.2694       0.2706
    ##   loy4          1.6908         1.6741       2.1804         2.3046       0.0840
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
    ##  Inner weighting scheme           = path
    ##  Type of indicator correlation    = Pearson
    ##  Path model estimator             = OLS
    ##  Second-order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample method                  = bootstrap
    ##  Number of resamples              = 499
    ##  Number of admissible results     = 488
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -377489653
    ## 
    ##  Construct details:
    ##  ------------------
    ##  Name  Modeled as     Order         Mode 
    ## 
    ##  IMAG  Composite      First order   modeB
    ##  EXPE  Composite      First order   modeB
    ##  QUAL  Composite      First order   modeB
    ##  VAL   Composite      First order   modeB
    ##  SAT   Common factor  First order   modeA
    ##  LOY   Common factor  First order   modeA
    ## 
    ## ----------------------------------- Estimates ----------------------------------
    ## 
    ## Estimated path coefficients:
    ## ============================
    ##                                                              CI_percentile   
    ##   Path           Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG      0.4714      0.0680    6.9267    0.0000 [ 0.3375; 0.5995 ] 
    ##   QUAL ~ EXPE      0.8344      0.0248   33.7149    0.0000 [ 0.7789; 0.8750 ] 
    ##   VAL ~ EXPE       0.0457      0.0884    0.5169    0.6052 [-0.1109; 0.2306 ] 
    ##   VAL ~ QUAL       0.7013      0.0819    8.5620    0.0000 [ 0.5243; 0.8372 ] 
    ##   SAT ~ IMAG       0.2450      0.0569    4.3072    0.0000 [ 0.1269; 0.3487 ] 
    ##   SAT ~ EXPE      -0.0172      0.0690   -0.2497    0.8028 [-0.1562; 0.1190 ] 
    ##   SAT ~ QUAL       0.2215      0.1011    2.1909    0.0285 [ 0.0378; 0.4223 ] 
    ##   SAT ~ VAL        0.5270      0.0849    6.2049    0.0000 [ 0.3587; 0.7033 ] 
    ##   LOY ~ IMAG       0.1819      0.0761    2.3916    0.0168 [ 0.0395; 0.3302 ] 
    ##   LOY ~ SAT        0.6283      0.0810    7.7552    0.0000 [ 0.4723; 0.7825 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1001    6.2996    0.0000 [ 0.4242; 0.8089 ] 
    ##   IMAG =~ imag2      0.9246      0.0393   23.5103    0.0000 [ 0.8323; 0.9782 ] 
    ##   IMAG =~ imag3      0.9577      0.0291   32.8627    0.0000 [ 0.8792; 0.9911 ] 
    ##   EXPE =~ expe1      0.7525      0.0758    9.9290    0.0000 [ 0.5668; 0.8758 ] 
    ##   EXPE =~ expe2      0.9348      0.0285   32.8337    0.0000 [ 0.8686; 0.9724 ] 
    ##   EXPE =~ expe3      0.7295      0.0718   10.1617    0.0000 [ 0.5623; 0.8460 ] 
    ##   QUAL =~ qual1      0.7861      0.0684   11.4850    0.0000 [ 0.6197; 0.8920 ] 
    ##   QUAL =~ qual2      0.9244      0.0222   41.6871    0.0000 [ 0.8691; 0.9585 ] 
    ##   QUAL =~ qual3      0.7560      0.0616   12.2650    0.0000 [ 0.6168; 0.8547 ] 
    ##   QUAL =~ qual4      0.7632      0.0514   14.8365    0.0000 [ 0.6584; 0.8516 ] 
    ##   QUAL =~ qual5      0.7834      0.0449   17.4564    0.0000 [ 0.6851; 0.8611 ] 
    ##   VAL =~ val1        0.9518      0.0239   39.7799    0.0000 [ 0.8938; 0.9853 ] 
    ##   VAL =~ val2        0.8056      0.0641   12.5666    0.0000 [ 0.6536; 0.9102 ] 
    ##   VAL =~ val3        0.6763      0.0724    9.3471    0.0000 [ 0.5225; 0.8128 ] 
    ##   SAT =~ sat1        0.9243      0.0221   41.7831    0.0000 [ 0.8734; 0.9622 ] 
    ##   SAT =~ sat2        0.8813      0.0284   31.0611    0.0000 [ 0.8177; 0.9212 ] 
    ##   SAT =~ sat3        0.7127      0.0519   13.7267    0.0000 [ 0.6022; 0.7971 ] 
    ##   SAT =~ sat4        0.7756      0.0483   16.0429    0.0000 [ 0.6728; 0.8527 ] 
    ##   LOY =~ loy1        0.9097      0.0503   18.0716    0.0000 [ 0.7998; 0.9821 ] 
    ##   LOY =~ loy2        0.5775      0.0866    6.6663    0.0000 [ 0.3932; 0.7238 ] 
    ##   LOY =~ loy3        0.9043      0.0424   21.3410    0.0000 [ 0.8108; 0.9743 ] 
    ##   LOY =~ loy4        0.4917      0.0961    5.1155    0.0000 [ 0.2999; 0.6732 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1215    0.1288    0.8975 [-0.2214; 0.2561 ] 
    ##   IMAG <~ imag2      0.4473      0.1529    2.9260    0.0034 [ 0.1290; 0.7530 ] 
    ##   IMAG <~ imag3      0.6020      0.1413    4.2606    0.0000 [ 0.2955; 0.8330 ] 
    ##   EXPE <~ expe1      0.2946      0.1138    2.5899    0.0096 [ 0.0750; 0.5179 ] 
    ##   EXPE <~ expe2      0.6473      0.0832    7.7844    0.0000 [ 0.4785; 0.7856 ] 
    ##   EXPE <~ expe3      0.2374      0.0935    2.5389    0.0111 [ 0.0564; 0.4281 ] 
    ##   QUAL <~ qual1      0.2370      0.0894    2.6502    0.0080 [ 0.0722; 0.4213 ] 
    ##   QUAL <~ qual2      0.4712      0.0766    6.1541    0.0000 [ 0.3120; 0.5985 ] 
    ##   QUAL <~ qual3      0.1831      0.0819    2.2348    0.0254 [ 0.0218; 0.3338 ] 
    ##   QUAL <~ qual4      0.1037      0.0610    1.6994    0.0892 [-0.0065; 0.2361 ] 
    ##   QUAL <~ qual5      0.2049      0.0570    3.5946    0.0003 [ 0.0851; 0.2997 ] 
    ##   VAL <~ val1        0.7163      0.0964    7.4320    0.0000 [ 0.5154; 0.8755 ] 
    ##   VAL <~ val2        0.2202      0.0913    2.4119    0.0159 [ 0.0642; 0.4140 ] 
    ##   VAL <~ val3        0.2082      0.0615    3.3854    0.0007 [ 0.0758; 0.3235 ] 
    ##   SAT <~ sat1        0.3209      0.0151   21.2657    0.0000 [ 0.2975; 0.3552 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.5689    0.0000 [ 0.2826; 0.3337 ] 
    ##   SAT <~ sat3        0.2474      0.0108   22.9138    0.0000 [ 0.2242; 0.2674 ] 
    ##   SAT <~ sat4        0.2692      0.0117   23.0826    0.0000 [ 0.2466; 0.2946 ] 
    ##   LOY <~ loy1        0.3834      0.0264   14.5241    0.0000 [ 0.3337; 0.4352 ] 
    ##   LOY <~ loy2        0.2434      0.0307    7.9209    0.0000 [ 0.1735; 0.2966 ] 
    ##   LOY <~ loy3        0.3812      0.0277   13.7505    0.0000 [ 0.3310; 0.4369 ] 
    ##   LOY <~ loy4        0.2073      0.0357    5.8078    0.0000 [ 0.1320; 0.2769 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0620   10.3865    0.0000 [ 0.5124; 0.7552 ] 
    ##   imag1 ~~ imag3      0.5433      0.0667    8.1423    0.0000 [ 0.4092; 0.6712 ] 
    ##   imag2 ~~ imag3      0.7761      0.0373   20.8231    0.0000 [ 0.7002; 0.8438 ] 
    ##   expe1 ~~ expe2      0.5353      0.0577    9.2722    0.0000 [ 0.3978; 0.6362 ] 
    ##   expe1 ~~ expe3      0.4694      0.0578    8.1201    0.0000 [ 0.3471; 0.5729 ] 
    ##   expe2 ~~ expe3      0.5467      0.0595    9.1829    0.0000 [ 0.4102; 0.6501 ] 
    ##   qual1 ~~ qual2      0.6053      0.0574   10.5479    0.0000 [ 0.4787; 0.7003 ] 
    ##   qual1 ~~ qual3      0.5406      0.0593    9.1215    0.0000 [ 0.4278; 0.6502 ] 
    ##   qual1 ~~ qual4      0.5662      0.0662    8.5556    0.0000 [ 0.4358; 0.6886 ] 
    ##   qual1 ~~ qual5      0.5180      0.0684    7.5778    0.0000 [ 0.3777; 0.6448 ] 
    ##   qual2 ~~ qual3      0.6187      0.0536   11.5357    0.0000 [ 0.5034; 0.7055 ] 
    ##   qual2 ~~ qual4      0.6517      0.0595   10.9556    0.0000 [ 0.5300; 0.7571 ] 
    ##   qual2 ~~ qual5      0.6291      0.0576   10.9291    0.0000 [ 0.5070; 0.7314 ] 
    ##   qual3 ~~ qual4      0.4752      0.0607    7.8244    0.0000 [ 0.3543; 0.5825 ] 
    ##   qual3 ~~ qual5      0.5074      0.0614    8.2612    0.0000 [ 0.3883; 0.6228 ] 
    ##   qual4 ~~ qual5      0.6402      0.0573   11.1640    0.0000 [ 0.5162; 0.7371 ] 
    ##   val1 ~~ val2        0.6344      0.0556   11.4045    0.0000 [ 0.5193; 0.7405 ] 
    ##   val1 ~~ val3        0.4602      0.0712    6.4647    0.0000 [ 0.3315; 0.6051 ] 
    ##   val2 ~~ val3        0.6288      0.0611   10.2897    0.0000 [ 0.5062; 0.7455 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0680    6.9267    0.0000 [ 0.3375; 0.5995 ] 
    ##   QUAL ~ IMAG       0.3933      0.0640    6.1433    0.0000 [ 0.2703; 0.5149 ] 
    ##   QUAL ~ EXPE       0.8344      0.0248   33.7149    0.0000 [ 0.7789; 0.8750 ] 
    ##   VAL ~ IMAG        0.2974      0.0633    4.6962    0.0000 [ 0.1795; 0.4164 ] 
    ##   VAL ~ EXPE        0.6309      0.0540   11.6921    0.0000 [ 0.5210; 0.7301 ] 
    ##   VAL ~ QUAL        0.7013      0.0819    8.5620    0.0000 [ 0.5243; 0.8372 ] 
    ##   SAT ~ IMAG        0.4807      0.0694    6.9270    0.0000 [ 0.3505; 0.6139 ] 
    ##   SAT ~ EXPE        0.5001      0.0582    8.5892    0.0000 [ 0.3831; 0.6095 ] 
    ##   SAT ~ QUAL        0.5911      0.0967    6.1138    0.0000 [ 0.4087; 0.7756 ] 
    ##   SAT ~ VAL         0.5270      0.0849    6.2049    0.0000 [ 0.3587; 0.7033 ] 
    ##   LOY ~ IMAG        0.4840      0.0651    7.4357    0.0000 [ 0.3584; 0.6177 ] 
    ##   LOY ~ EXPE        0.3142      0.0564    5.5669    0.0000 [ 0.2170; 0.4351 ] 
    ##   LOY ~ QUAL        0.3714      0.0868    4.2775    0.0000 [ 0.2062; 0.5548 ] 
    ##   LOY ~ VAL         0.3311      0.0719    4.6018    0.0000 [ 0.2026; 0.4811 ] 
    ##   LOY ~ SAT         0.6283      0.0810    7.7552    0.0000 [ 0.4723; 0.7825 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0640    6.1433    0.0000 [ 0.2703; 0.5149 ] 
    ##   VAL ~ IMAG           0.2974      0.0633    4.6962    0.0000 [ 0.1795; 0.4164 ] 
    ##   VAL ~ EXPE           0.5852      0.0689    8.4928    0.0000 [ 0.4402; 0.7008 ] 
    ##   SAT ~ IMAG           0.2357      0.0508    4.6426    0.0000 [ 0.1419; 0.3353 ] 
    ##   SAT ~ EXPE           0.5173      0.0669    7.7305    0.0000 [ 0.3959; 0.6542 ] 
    ##   SAT ~ QUAL           0.3696      0.0619    5.9663    0.0000 [ 0.2505; 0.5021 ] 
    ##   LOY ~ IMAG           0.3020      0.0571    5.2856    0.0000 [ 0.2025; 0.4211 ] 
    ##   LOY ~ EXPE           0.3142      0.0564    5.5669    0.0000 [ 0.2170; 0.4351 ] 
    ##   LOY ~ QUAL           0.3714      0.0868    4.2775    0.0000 [ 0.2062; 0.5548 ] 
    ##   LOY ~ VAL            0.3311      0.0719    4.6018    0.0000 [ 0.2026; 0.4811 ] 
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
