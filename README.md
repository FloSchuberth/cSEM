
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
and two do-family of function:

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
    ##                                            +------------------------------------------------------------------+
    ##                                            |                                                                  |
    ##                                            |   H0: The model-implied indicator covariance matrix equals the   |
    ##                                            |   population indicator covariance matrix.                        |
    ##                                            |                                                                  |
    ##                                            +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3309  
    ##  SRMR                    0.0940      0.0531  
    ##  dL                      2.2340      0.7137  
    ##  dML                     2.9219      1.6844  
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
    ##  Out of 499 bootstrap replications 472 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 644455280
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
    ##   expe1         1.4554         1.5684       1.9068         2.0938       0.0545
    ##   expe2         1.4155         1.4817       1.9363         2.0313       0.1986
    ##   expe3         1.6301         1.7233       2.1253         2.2193       0.1256
    ##   qual1         1.4759         1.5458       1.9284         2.0576       0.1148
    ##   qual2         1.5776         1.5366       2.0389         2.0605       0.2174
    ##   qual3         1.7320         1.7278       2.2234         2.2830       0.1204
    ##   qual4         1.2324         1.1947       1.5948         1.6264       0.2336
    ##   qual5         1.5060         1.5004       1.9355         1.9529       0.1965
    ##   val1          1.4487         1.3692       1.8736         1.7742       0.2478
    ##   val2          1.2288         1.2080       1.6521         1.7180       0.1702
    ##   val3          1.4822         1.3811       1.9705         1.9355       0.1470
    ##   sat1          1.2487         1.2336       1.6465         1.6195       0.3374
    ##   sat2          1.2308         1.1948       1.6392         1.6244       0.3092
    ##   sat3          1.3397         1.2763       1.6697         1.7188       0.2128
    ##   sat4          1.3184         1.2595       1.6677         1.6330       0.2777
    ##   loy1          1.6928         1.6591       2.2343         2.2240       0.2672
    ##   loy2          1.4829         1.4724       1.9110         1.9808       0.1322
    ##   loy3          1.7022         1.6679       2.2785         2.2666       0.2716
    ##   loy4          1.6896         1.6726       2.1787         2.3072       0.0850
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
    ##  Number of admissible results     = 487
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 664237232
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
    ##   EXPE ~ IMAG      0.4714      0.0656    7.1867    0.0000 [ 0.3360; 0.6003 ] 
    ##   QUAL ~ EXPE      0.8344      0.0220   38.0123    0.0000 [ 0.7879; 0.8758 ] 
    ##   VAL ~ EXPE       0.0457      0.0894    0.5116    0.6089 [-0.1091; 0.2325 ] 
    ##   VAL ~ QUAL       0.7013      0.0845    8.3002    0.0000 [ 0.5216; 0.8630 ] 
    ##   SAT ~ IMAG       0.2450      0.0539    4.5423    0.0000 [ 0.1516; 0.3496 ] 
    ##   SAT ~ EXPE      -0.0172      0.0708   -0.2433    0.8078 [-0.1573; 0.1133 ] 
    ##   SAT ~ QUAL       0.2215      0.0967    2.2912    0.0220 [ 0.0553; 0.4130 ] 
    ##   SAT ~ VAL        0.5270      0.0847    6.2203    0.0000 [ 0.3531; 0.6738 ] 
    ##   LOY ~ IMAG       0.1819      0.0813    2.2374    0.0253 [ 0.0285; 0.3418 ] 
    ##   LOY ~ SAT        0.6283      0.0808    7.7771    0.0000 [ 0.4756; 0.7909 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1035    6.0914    0.0000 [ 0.4127; 0.8086 ] 
    ##   IMAG =~ imag2      0.9246      0.0403   22.9566    0.0000 [ 0.8221; 0.9767 ] 
    ##   IMAG =~ imag3      0.9577      0.0285   33.5784    0.0000 [ 0.8809; 0.9893 ] 
    ##   EXPE =~ expe1      0.7525      0.0811    9.2810    0.0000 [ 0.5625; 0.8805 ] 
    ##   EXPE =~ expe2      0.9348      0.0303   30.8551    0.0000 [ 0.8607; 0.9717 ] 
    ##   EXPE =~ expe3      0.7295      0.0764    9.5433    0.0000 [ 0.5546; 0.8528 ] 
    ##   QUAL =~ qual1      0.7861      0.0696   11.2924    0.0000 [ 0.6162; 0.8888 ] 
    ##   QUAL =~ qual2      0.9244      0.0238   38.7672    0.0000 [ 0.8664; 0.9573 ] 
    ##   QUAL =~ qual3      0.7560      0.0633   11.9409    0.0000 [ 0.6108; 0.8576 ] 
    ##   QUAL =~ qual4      0.7632      0.0513   14.8852    0.0000 [ 0.6458; 0.8420 ] 
    ##   QUAL =~ qual5      0.7834      0.0467   16.7644    0.0000 [ 0.6795; 0.8558 ] 
    ##   VAL =~ val1        0.9518      0.0225   42.2838    0.0000 [ 0.8954; 0.9841 ] 
    ##   VAL =~ val2        0.8056      0.0628   12.8292    0.0000 [ 0.6619; 0.9051 ] 
    ##   VAL =~ val3        0.6763      0.0718    9.4153    0.0000 [ 0.5255; 0.7993 ] 
    ##   SAT =~ sat1        0.9243      0.0229   40.4119    0.0000 [ 0.8724; 0.9615 ] 
    ##   SAT =~ sat2        0.8813      0.0284   30.9871    0.0000 [ 0.8161; 0.9287 ] 
    ##   SAT =~ sat3        0.7127      0.0520   13.6999    0.0000 [ 0.5979; 0.8039 ] 
    ##   SAT =~ sat4        0.7756      0.0500   15.5207    0.0000 [ 0.6735; 0.8589 ] 
    ##   LOY =~ loy1        0.9097      0.0503   18.0999    0.0000 [ 0.7954; 0.9914 ] 
    ##   LOY =~ loy2        0.5775      0.0800    7.2199    0.0000 [ 0.4058; 0.7199 ] 
    ##   LOY =~ loy3        0.9043      0.0412   21.9748    0.0000 [ 0.8172; 0.9723 ] 
    ##   LOY =~ loy4        0.4917      0.0988    4.9787    0.0000 [ 0.3044; 0.6912 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1212    0.1291    0.8973 [-0.2147; 0.2524 ] 
    ##   IMAG <~ imag2      0.4473      0.1429    3.1302    0.0017 [ 0.1854; 0.7127 ] 
    ##   IMAG <~ imag3      0.6020      0.1395    4.3172    0.0000 [ 0.3113; 0.8444 ] 
    ##   EXPE <~ expe1      0.2946      0.1176    2.5055    0.0122 [ 0.0563; 0.5161 ] 
    ##   EXPE <~ expe2      0.6473      0.0910    7.1150    0.0000 [ 0.4433; 0.7999 ] 
    ##   EXPE <~ expe3      0.2374      0.0937    2.5324    0.0113 [ 0.0634; 0.4217 ] 
    ##   QUAL <~ qual1      0.2370      0.0881    2.6905    0.0071 [ 0.0702; 0.4071 ] 
    ##   QUAL <~ qual2      0.4712      0.0811    5.8120    0.0000 [ 0.3041; 0.6193 ] 
    ##   QUAL <~ qual3      0.1831      0.0804    2.2772    0.0228 [ 0.0239; 0.3416 ] 
    ##   QUAL <~ qual4      0.1037      0.0589    1.7604    0.0783 [-0.0137; 0.2178 ] 
    ##   QUAL <~ qual5      0.2049      0.0629    3.2577    0.0011 [ 0.0716; 0.3313 ] 
    ##   VAL <~ val1        0.7163      0.0946    7.5751    0.0000 [ 0.5177; 0.8880 ] 
    ##   VAL <~ val2        0.2202      0.0914    2.4086    0.0160 [ 0.0546; 0.4105 ] 
    ##   VAL <~ val3        0.2082      0.0565    3.6859    0.0002 [ 0.0860; 0.3112 ] 
    ##   SAT <~ sat1        0.3209      0.0154   20.8703    0.0000 [ 0.2959; 0.3567 ] 
    ##   SAT <~ sat2        0.3059      0.0138   22.2443    0.0000 [ 0.2826; 0.3351 ] 
    ##   SAT <~ sat3        0.2474      0.0111   22.3187    0.0000 [ 0.2243; 0.2667 ] 
    ##   SAT <~ sat4        0.2692      0.0117   23.0286    0.0000 [ 0.2472; 0.2922 ] 
    ##   LOY <~ loy1        0.3834      0.0262   14.6298    0.0000 [ 0.3311; 0.4341 ] 
    ##   LOY <~ loy2        0.2434      0.0286    8.4970    0.0000 [ 0.1769; 0.2949 ] 
    ##   LOY <~ loy3        0.3812      0.0263   14.5195    0.0000 [ 0.3294; 0.4323 ] 
    ##   LOY <~ loy4        0.2073      0.0359    5.7813    0.0000 [ 0.1360; 0.2818 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0648    9.9405    0.0000 [ 0.5069; 0.7532 ] 
    ##   imag1 ~~ imag3      0.5433      0.0710    7.6506    0.0000 [ 0.3894; 0.6726 ] 
    ##   imag2 ~~ imag3      0.7761      0.0394   19.7078    0.0000 [ 0.6886; 0.8381 ] 
    ##   expe1 ~~ expe2      0.5353      0.0606    8.8396    0.0000 [ 0.4091; 0.6424 ] 
    ##   expe1 ~~ expe3      0.4694      0.0598    7.8513    0.0000 [ 0.3577; 0.5810 ] 
    ##   expe2 ~~ expe3      0.5467      0.0646    8.4683    0.0000 [ 0.4048; 0.6583 ] 
    ##   qual1 ~~ qual2      0.6053      0.0619    9.7871    0.0000 [ 0.4728; 0.7060 ] 
    ##   qual1 ~~ qual3      0.5406      0.0625    8.6559    0.0000 [ 0.4075; 0.6477 ] 
    ##   qual1 ~~ qual4      0.5662      0.0659    8.5977    0.0000 [ 0.4302; 0.6805 ] 
    ##   qual1 ~~ qual5      0.5180      0.0678    7.6444    0.0000 [ 0.3780; 0.6448 ] 
    ##   qual2 ~~ qual3      0.6187      0.0581   10.6505    0.0000 [ 0.4874; 0.7126 ] 
    ##   qual2 ~~ qual4      0.6517      0.0589   11.0673    0.0000 [ 0.5341; 0.7528 ] 
    ##   qual2 ~~ qual5      0.6291      0.0575   10.9335    0.0000 [ 0.5071; 0.7238 ] 
    ##   qual3 ~~ qual4      0.4752      0.0676    7.0333    0.0000 [ 0.3242; 0.5930 ] 
    ##   qual3 ~~ qual5      0.5074      0.0639    7.9347    0.0000 [ 0.3689; 0.6252 ] 
    ##   qual4 ~~ qual5      0.6402      0.0534   11.9808    0.0000 [ 0.5251; 0.7283 ] 
    ##   val1 ~~ val2        0.6344      0.0540   11.7487    0.0000 [ 0.5323; 0.7370 ] 
    ##   val1 ~~ val3        0.4602      0.0678    6.7871    0.0000 [ 0.3233; 0.5929 ] 
    ##   val2 ~~ val3        0.6288      0.0615   10.2208    0.0000 [ 0.5102; 0.7433 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0656    7.1867    0.0000 [ 0.3360; 0.6003 ] 
    ##   QUAL ~ IMAG       0.3933      0.0599    6.5627    0.0000 [ 0.2777; 0.5130 ] 
    ##   QUAL ~ EXPE       0.8344      0.0220   38.0123    0.0000 [ 0.7879; 0.8758 ] 
    ##   VAL ~ IMAG        0.2974      0.0594    5.0062    0.0000 [ 0.1891; 0.4295 ] 
    ##   VAL ~ EXPE        0.6309      0.0480   13.1401    0.0000 [ 0.5294; 0.7194 ] 
    ##   VAL ~ QUAL        0.7013      0.0845    8.3002    0.0000 [ 0.5216; 0.8630 ] 
    ##   SAT ~ IMAG        0.4807      0.0644    7.4691    0.0000 [ 0.3571; 0.6151 ] 
    ##   SAT ~ EXPE        0.5001      0.0593    8.4402    0.0000 [ 0.3940; 0.6091 ] 
    ##   SAT ~ QUAL        0.5911      0.0918    6.4392    0.0000 [ 0.4173; 0.7662 ] 
    ##   SAT ~ VAL         0.5270      0.0847    6.2203    0.0000 [ 0.3531; 0.6738 ] 
    ##   LOY ~ IMAG        0.4840      0.0674    7.1765    0.0000 [ 0.3561; 0.6179 ] 
    ##   LOY ~ EXPE        0.3142      0.0557    5.6428    0.0000 [ 0.2054; 0.4183 ] 
    ##   LOY ~ QUAL        0.3714      0.0825    4.4993    0.0000 [ 0.2324; 0.5474 ] 
    ##   LOY ~ VAL         0.3311      0.0719    4.6073    0.0000 [ 0.1953; 0.4821 ] 
    ##   LOY ~ SAT         0.6283      0.0808    7.7771    0.0000 [ 0.4756; 0.7909 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0599    6.5627    0.0000 [ 0.2777; 0.5130 ] 
    ##   VAL ~ IMAG           0.2974      0.0594    5.0062    0.0000 [ 0.1891; 0.4295 ] 
    ##   VAL ~ EXPE           0.5852      0.0736    7.9526    0.0000 [ 0.4388; 0.7260 ] 
    ##   SAT ~ IMAG           0.2357      0.0478    4.9308    0.0000 [ 0.1465; 0.3315 ] 
    ##   SAT ~ EXPE           0.5173      0.0663    7.8080    0.0000 [ 0.3973; 0.6535 ] 
    ##   SAT ~ QUAL           0.3696      0.0620    5.9577    0.0000 [ 0.2422; 0.4789 ] 
    ##   LOY ~ IMAG           0.3020      0.0572    5.2781    0.0000 [ 0.1979; 0.4234 ] 
    ##   LOY ~ EXPE           0.3142      0.0557    5.6428    0.0000 [ 0.2054; 0.4183 ] 
    ##   LOY ~ QUAL           0.3714      0.0825    4.4993    0.0000 [ 0.2324; 0.5474 ] 
    ##   LOY ~ VAL            0.3311      0.0719    4.6073    0.0000 [ 0.1953; 0.4821 ] 
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
