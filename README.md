
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
    ##             +------------------------------------------------------------------+
    ##             |                                                                  |
    ##             |   H0: The model-implied indicator covariance matrix equals the   |
    ##             |   population indicator covariance matrix.                        |
    ##             |                                                                  |
    ##             +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3261  
    ##  SRMR                    0.0940      0.0534  
    ##  dL                      2.2340      0.7220  
    ##  dML                     2.9219      1.6189  
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
    ##  Out of 499 bootstrap replications 471 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -284148939
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
    ##   expe1         1.4557         1.5655       1.9070         2.0931       0.0556
    ##   expe2         1.4150         1.4792       1.9372         2.0290       0.1982
    ##   expe3         1.6286         1.7270       2.1223         2.2226       0.1257
    ##   qual1         1.4767         1.5456       1.9295         2.0577       0.1149
    ##   qual2         1.5803         1.5386       2.0422         2.0601       0.2154
    ##   qual3         1.7319         1.7293       2.2233         2.2862       0.1183
    ##   qual4         1.2370         1.1999       1.6002         1.6346       0.2308
    ##   qual5         1.5051         1.5026       1.9354         1.9564       0.1951
    ##   val1          1.4472         1.3666       1.8708         1.7704       0.2490
    ##   val2          1.2278         1.2053       1.6490         1.7132       0.1733
    ##   val3          1.4814         1.3806       1.9704         1.9373       0.1469
    ##   sat1          1.2479         1.2344       1.6465         1.6211       0.3393
    ##   sat2          1.2327         1.1969       1.6419         1.6300       0.3087
    ##   sat3          1.3413         1.2781       1.6735         1.7239       0.2106
    ##   sat4          1.3183         1.2640       1.6677         1.6370       0.2767
    ##   loy1          1.6897         1.6600       2.2338         2.2309       0.2684
    ##   loy2          1.4836         1.4757       1.9107         1.9848       0.1321
    ##   loy3          1.7044         1.6717       2.2836         2.2755       0.2689
    ##   loy4          1.6892         1.6677       2.1777         2.2998       0.0874
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
    ##  Resample methode                 = bootstrap
    ##  Number of resamples              = 499
    ##  Number of admissible results     = 489
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -2095294199
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
    ##   EXPE ~ IMAG      0.4714      0.0611    7.7112    0.0000 [ 0.3551; 0.5896 ] 
    ##   QUAL ~ EXPE      0.8344      0.0227   36.6984    0.0000 [ 0.7861; 0.8723 ] 
    ##   VAL ~ EXPE       0.0457      0.0812    0.5629    0.5735 [-0.1043; 0.1976 ] 
    ##   VAL ~ QUAL       0.7013      0.0797    8.7965    0.0000 [ 0.5479; 0.8471 ] 
    ##   SAT ~ IMAG       0.2450      0.0554    4.4257    0.0000 [ 0.1411; 0.3450 ] 
    ##   SAT ~ EXPE      -0.0172      0.0765   -0.2252    0.8218 [-0.1711; 0.1342 ] 
    ##   SAT ~ QUAL       0.2215      0.1046    2.1175    0.0342 [ 0.0356; 0.4333 ] 
    ##   SAT ~ VAL        0.5270      0.0889    5.9298    0.0000 [ 0.3531; 0.6965 ] 
    ##   LOY ~ IMAG       0.1819      0.0818    2.2243    0.0261 [ 0.0272; 0.3521 ] 
    ##   LOY ~ SAT        0.6283      0.0809    7.7711    0.0000 [ 0.4646; 0.7959 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0971    6.4953    0.0000 [ 0.4265; 0.7972 ] 
    ##   IMAG =~ imag2      0.9246      0.0367   25.1648    0.0000 [ 0.8272; 0.9756 ] 
    ##   IMAG =~ imag3      0.9577      0.0276   34.7497    0.0000 [ 0.8849; 0.9909 ] 
    ##   EXPE =~ expe1      0.7525      0.0742   10.1410    0.0000 [ 0.5831; 0.8717 ] 
    ##   EXPE =~ expe2      0.9348      0.0289   32.3491    0.0000 [ 0.8601; 0.9708 ] 
    ##   EXPE =~ expe3      0.7295      0.0654   11.1587    0.0000 [ 0.5843; 0.8417 ] 
    ##   QUAL =~ qual1      0.7861      0.0657   11.9720    0.0000 [ 0.6298; 0.8906 ] 
    ##   QUAL =~ qual2      0.9244      0.0231   40.0967    0.0000 [ 0.8686; 0.9550 ] 
    ##   QUAL =~ qual3      0.7560      0.0573   13.1992    0.0000 [ 0.6316; 0.8607 ] 
    ##   QUAL =~ qual4      0.7632      0.0520   14.6888    0.0000 [ 0.6460; 0.8472 ] 
    ##   QUAL =~ qual5      0.7834      0.0450   17.4097    0.0000 [ 0.6860; 0.8606 ] 
    ##   VAL =~ val1        0.9518      0.0223   42.7410    0.0000 [ 0.9004; 0.9850 ] 
    ##   VAL =~ val2        0.8056      0.0606   13.3014    0.0000 [ 0.6656; 0.8988 ] 
    ##   VAL =~ val3        0.6763      0.0715    9.4631    0.0000 [ 0.5256; 0.8087 ] 
    ##   SAT =~ sat1        0.9243      0.0223   41.5277    0.0000 [ 0.8700; 0.9614 ] 
    ##   SAT =~ sat2        0.8813      0.0283   31.1729    0.0000 [ 0.8198; 0.9256 ] 
    ##   SAT =~ sat3        0.7127      0.0513   13.8974    0.0000 [ 0.6096; 0.8030 ] 
    ##   SAT =~ sat4        0.7756      0.0510   15.1962    0.0000 [ 0.6693; 0.8647 ] 
    ##   LOY =~ loy1        0.9097      0.0475   19.1534    0.0000 [ 0.8001; 0.9861 ] 
    ##   LOY =~ loy2        0.5775      0.0832    6.9418    0.0000 [ 0.4178; 0.7318 ] 
    ##   LOY =~ loy3        0.9043      0.0420   21.5517    0.0000 [ 0.8193; 0.9741 ] 
    ##   LOY =~ loy4        0.4917      0.0953    5.1589    0.0000 [ 0.3142; 0.6754 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1166    0.1342    0.8933 [-0.2293; 0.2327 ] 
    ##   IMAG <~ imag2      0.4473      0.1464    3.0542    0.0023 [ 0.1634; 0.7113 ] 
    ##   IMAG <~ imag3      0.6020      0.1321    4.5588    0.0000 [ 0.3324; 0.8414 ] 
    ##   EXPE <~ expe1      0.2946      0.1113    2.6460    0.0081 [ 0.0756; 0.5169 ] 
    ##   EXPE <~ expe2      0.6473      0.0827    7.8304    0.0000 [ 0.4568; 0.7821 ] 
    ##   EXPE <~ expe3      0.2374      0.0850    2.7920    0.0052 [ 0.0754; 0.4034 ] 
    ##   QUAL <~ qual1      0.2370      0.0833    2.8446    0.0044 [ 0.0801; 0.4128 ] 
    ##   QUAL <~ qual2      0.4712      0.0777    6.0633    0.0000 [ 0.2988; 0.6087 ] 
    ##   QUAL <~ qual3      0.1831      0.0749    2.4430    0.0146 [ 0.0372; 0.3282 ] 
    ##   QUAL <~ qual4      0.1037      0.0577    1.7984    0.0721 [-0.0090; 0.2106 ] 
    ##   QUAL <~ qual5      0.2049      0.0600    3.4127    0.0006 [ 0.0859; 0.3089 ] 
    ##   VAL <~ val1        0.7163      0.0897    7.9814    0.0000 [ 0.5326; 0.8648 ] 
    ##   VAL <~ val2        0.2202      0.0871    2.5277    0.0115 [ 0.0669; 0.3952 ] 
    ##   VAL <~ val3        0.2082      0.0589    3.5347    0.0004 [ 0.1006; 0.3264 ] 
    ##   SAT <~ sat1        0.3209      0.0156   20.5824    0.0000 [ 0.2945; 0.3521 ] 
    ##   SAT <~ sat2        0.3059      0.0136   22.4447    0.0000 [ 0.2821; 0.3331 ] 
    ##   SAT <~ sat3        0.2474      0.0103   24.0424    0.0000 [ 0.2280; 0.2673 ] 
    ##   SAT <~ sat4        0.2692      0.0122   22.0832    0.0000 [ 0.2463; 0.2945 ] 
    ##   LOY <~ loy1        0.3834      0.0261   14.7097    0.0000 [ 0.3318; 0.4347 ] 
    ##   LOY <~ loy2        0.2434      0.0294    8.2714    0.0000 [ 0.1788; 0.2967 ] 
    ##   LOY <~ loy3        0.3812      0.0256   14.8944    0.0000 [ 0.3295; 0.4285 ] 
    ##   LOY <~ loy4        0.2073      0.0352    5.8910    0.0000 [ 0.1316; 0.2743 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0633   10.1735    0.0000 [ 0.5121; 0.7505 ] 
    ##   imag1 ~~ imag3      0.5433      0.0692    7.8511    0.0000 [ 0.4029; 0.6705 ] 
    ##   imag2 ~~ imag3      0.7761      0.0382   20.3096    0.0000 [ 0.6919; 0.8474 ] 
    ##   expe1 ~~ expe2      0.5353      0.0614    8.7244    0.0000 [ 0.4075; 0.6438 ] 
    ##   expe1 ~~ expe3      0.4694      0.0603    7.7892    0.0000 [ 0.3479; 0.5854 ] 
    ##   expe2 ~~ expe3      0.5467      0.0584    9.3551    0.0000 [ 0.4235; 0.6560 ] 
    ##   qual1 ~~ qual2      0.6053      0.0579   10.4586    0.0000 [ 0.4898; 0.7091 ] 
    ##   qual1 ~~ qual3      0.5406      0.0581    9.3074    0.0000 [ 0.4283; 0.6378 ] 
    ##   qual1 ~~ qual4      0.5662      0.0659    8.5855    0.0000 [ 0.4250; 0.6838 ] 
    ##   qual1 ~~ qual5      0.5180      0.0652    7.9460    0.0000 [ 0.3845; 0.6405 ] 
    ##   qual2 ~~ qual3      0.6187      0.0554   11.1761    0.0000 [ 0.4985; 0.7262 ] 
    ##   qual2 ~~ qual4      0.6517      0.0609   10.6941    0.0000 [ 0.5116; 0.7581 ] 
    ##   qual2 ~~ qual5      0.6291      0.0601   10.4707    0.0000 [ 0.5051; 0.7311 ] 
    ##   qual3 ~~ qual4      0.4752      0.0656    7.2473    0.0000 [ 0.3427; 0.5967 ] 
    ##   qual3 ~~ qual5      0.5074      0.0616    8.2413    0.0000 [ 0.3933; 0.6331 ] 
    ##   qual4 ~~ qual5      0.6402      0.0564   11.3541    0.0000 [ 0.5177; 0.7396 ] 
    ##   val1 ~~ val2        0.6344      0.0535   11.8597    0.0000 [ 0.5189; 0.7241 ] 
    ##   val1 ~~ val3        0.4602      0.0700    6.5744    0.0000 [ 0.3250; 0.5860 ] 
    ##   val2 ~~ val3        0.6288      0.0627   10.0297    0.0000 [ 0.5019; 0.7491 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0611    7.7112    0.0000 [ 0.3551; 0.5896 ] 
    ##   QUAL ~ IMAG       0.3933      0.0574    6.8553    0.0000 [ 0.2897; 0.5092 ] 
    ##   QUAL ~ EXPE       0.8344      0.0227   36.6984    0.0000 [ 0.7861; 0.8723 ] 
    ##   VAL ~ IMAG        0.2974      0.0575    5.1725    0.0000 [ 0.1990; 0.4176 ] 
    ##   VAL ~ EXPE        0.6309      0.0503   12.5332    0.0000 [ 0.5295; 0.7304 ] 
    ##   VAL ~ QUAL        0.7013      0.0797    8.7965    0.0000 [ 0.5479; 0.8471 ] 
    ##   SAT ~ IMAG        0.4807      0.0654    7.3475    0.0000 [ 0.3522; 0.6140 ] 
    ##   SAT ~ EXPE        0.5001      0.0600    8.3348    0.0000 [ 0.3717; 0.6051 ] 
    ##   SAT ~ QUAL        0.5911      0.0967    6.1140    0.0000 [ 0.4121; 0.7801 ] 
    ##   SAT ~ VAL         0.5270      0.0889    5.9298    0.0000 [ 0.3531; 0.6965 ] 
    ##   LOY ~ IMAG        0.4840      0.0704    6.8778    0.0000 [ 0.3505; 0.6187 ] 
    ##   LOY ~ EXPE        0.3142      0.0536    5.8573    0.0000 [ 0.2193; 0.4180 ] 
    ##   LOY ~ QUAL        0.3714      0.0847    4.3822    0.0000 [ 0.2247; 0.5601 ] 
    ##   LOY ~ VAL         0.3311      0.0790    4.1888    0.0000 [ 0.1844; 0.4891 ] 
    ##   LOY ~ SAT         0.6283      0.0809    7.7711    0.0000 [ 0.4646; 0.7959 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0574    6.8553    0.0000 [ 0.2897; 0.5092 ] 
    ##   VAL ~ IMAG           0.2974      0.0575    5.1725    0.0000 [ 0.1990; 0.4176 ] 
    ##   VAL ~ EXPE           0.5852      0.0692    8.4512    0.0000 [ 0.4451; 0.7197 ] 
    ##   SAT ~ IMAG           0.2357      0.0460    5.1210    0.0000 [ 0.1573; 0.3337 ] 
    ##   SAT ~ EXPE           0.5173      0.0717    7.2121    0.0000 [ 0.3831; 0.6669 ] 
    ##   SAT ~ QUAL           0.3696      0.0610    6.0617    0.0000 [ 0.2564; 0.4844 ] 
    ##   LOY ~ IMAG           0.3020      0.0497    6.0788    0.0000 [ 0.2021; 0.3919 ] 
    ##   LOY ~ EXPE           0.3142      0.0536    5.8573    0.0000 [ 0.2193; 0.4180 ] 
    ##   LOY ~ QUAL           0.3714      0.0847    4.3822    0.0000 [ 0.2247; 0.5601 ] 
    ##   LOY ~ VAL            0.3311      0.0790    4.1888    0.0000 [ 0.1844; 0.4891 ] 
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
