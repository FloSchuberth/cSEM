
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

There are five major postestimation verbs, four test family functions
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
functions. Currently, four do functions are implemented:

  - `doIPMA()`: performs an importance-performance matrix analysis
  - `doNonLinearEffectsnalysis()`: performs a nonlinear effects analysis
    such as floodlight and surface analysis
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
    ##              +------------------------------------------------------------------+
    ##              |                                                                  |
    ##              |   H0: The model-implied indicator covariance matrix equals the   |
    ##              |   population indicator covariance matrix.                        |
    ##              |                                                                  |
    ##              +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3165  
    ##  SRMR                    0.0940      0.0525  
    ##  dL                      2.2340      0.6964  
    ##  dML                     2.9219      1.6237  
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
    ##  Out of 499 bootstrap replications 482 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: 118873338
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
    ##   expe1         1.4550         1.5652       1.9075         2.0931       0.0530
    ##   expe2         1.4106         1.4806       1.9323         2.0294       0.2006
    ##   expe3         1.6285         1.7236       2.1248         2.2203       0.1241
    ##   qual1         1.4746         1.5452       1.9280         2.0615       0.1149
    ##   qual2         1.5761         1.5350       2.0363         2.0552       0.2196
    ##   qual3         1.7305         1.7284       2.2220         2.2838       0.1187
    ##   qual4         1.2318         1.1973       1.5923         1.6257       0.2370
    ##   qual5         1.5029         1.4968       1.9307         1.9481       0.1991
    ##   val1          1.4445         1.3616       1.8686         1.7621       0.2515
    ##   val2          1.2273         1.2035       1.6479         1.7085       0.1736
    ##   val3          1.4795         1.3806       1.9673         1.9353       0.1494
    ##   sat1          1.2454         1.2315       1.6431         1.6167       0.3422
    ##   sat2          1.2301         1.1943       1.6381         1.6234       0.3115
    ##   sat3          1.3425         1.2730       1.6733         1.7159       0.2104
    ##   sat4          1.3211         1.2619       1.6711         1.6332       0.2759
    ##   loy1          1.6912         1.6588       2.2311         2.2252       0.2701
    ##   loy2          1.4865         1.4746       1.9146         1.9838       0.1310
    ##   loy3          1.7020         1.6667       2.2795         2.2671       0.2717
    ##   loy4          1.6901         1.6738       2.1808         2.3070       0.0853
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
    ##  Number of admissible results     = 481
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = -1541810141
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
    ##   EXPE ~ IMAG      0.4714      0.0689    6.8443    0.0000 [ 0.3449; 0.6041 ] 
    ##   QUAL ~ EXPE      0.8344      0.0227   36.7309    0.0000 [ 0.7885; 0.8769 ] 
    ##   VAL ~ EXPE       0.0457      0.0848    0.5388    0.5900 [-0.1094; 0.2404 ] 
    ##   VAL ~ QUAL       0.7013      0.0831    8.4368    0.0000 [ 0.5279; 0.8484 ] 
    ##   SAT ~ IMAG       0.2450      0.0556    4.4072    0.0000 [ 0.1320; 0.3538 ] 
    ##   SAT ~ EXPE      -0.0172      0.0708   -0.2434    0.8077 [-0.1526; 0.1222 ] 
    ##   SAT ~ QUAL       0.2215      0.1043    2.1239    0.0337 [ 0.0403; 0.4445 ] 
    ##   SAT ~ VAL        0.5270      0.0925    5.6984    0.0000 [ 0.3447; 0.6968 ] 
    ##   LOY ~ IMAG       0.1819      0.0789    2.3069    0.0211 [ 0.0382; 0.3471 ] 
    ##   LOY ~ SAT        0.6283      0.0847    7.4168    0.0000 [ 0.4692; 0.7805 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0986    6.3965    0.0000 [ 0.3983; 0.7893 ] 
    ##   IMAG =~ imag2      0.9246      0.0377   24.4994    0.0000 [ 0.8357; 0.9737 ] 
    ##   IMAG =~ imag3      0.9577      0.0279   34.2897    0.0000 [ 0.8794; 0.9918 ] 
    ##   EXPE =~ expe1      0.7525      0.0773    9.7320    0.0000 [ 0.5855; 0.8766 ] 
    ##   EXPE =~ expe2      0.9348      0.0303   30.8987    0.0000 [ 0.8585; 0.9742 ] 
    ##   EXPE =~ expe3      0.7295      0.0729   10.0038    0.0000 [ 0.5421; 0.8440 ] 
    ##   QUAL =~ qual1      0.7861      0.0671   11.7232    0.0000 [ 0.6312; 0.8979 ] 
    ##   QUAL =~ qual2      0.9244      0.0241   38.3748    0.0000 [ 0.8700; 0.9621 ] 
    ##   QUAL =~ qual3      0.7560      0.0621   12.1779    0.0000 [ 0.6125; 0.8482 ] 
    ##   QUAL =~ qual4      0.7632      0.0518   14.7413    0.0000 [ 0.6519; 0.8577 ] 
    ##   QUAL =~ qual5      0.7834      0.0466   16.8158    0.0000 [ 0.6950; 0.8611 ] 
    ##   VAL =~ val1        0.9518      0.0236   40.3836    0.0000 [ 0.8946; 0.9847 ] 
    ##   VAL =~ val2        0.8056      0.0653   12.3454    0.0000 [ 0.6545; 0.9059 ] 
    ##   VAL =~ val3        0.6763      0.0740    9.1435    0.0000 [ 0.5128; 0.7977 ] 
    ##   SAT =~ sat1        0.9243      0.0222   41.5426    0.0000 [ 0.8733; 0.9608 ] 
    ##   SAT =~ sat2        0.8813      0.0292   30.2078    0.0000 [ 0.8206; 0.9284 ] 
    ##   SAT =~ sat3        0.7127      0.0526   13.5483    0.0000 [ 0.6072; 0.8060 ] 
    ##   SAT =~ sat4        0.7756      0.0514   15.0898    0.0000 [ 0.6510; 0.8641 ] 
    ##   LOY =~ loy1        0.9097      0.0485   18.7619    0.0000 [ 0.8003; 0.9870 ] 
    ##   LOY =~ loy2        0.5775      0.0829    6.9708    0.0000 [ 0.4036; 0.7361 ] 
    ##   LOY =~ loy3        0.9043      0.0438   20.6387    0.0000 [ 0.8101; 0.9750 ] 
    ##   LOY =~ loy4        0.4917      0.0978    5.0297    0.0000 [ 0.3077; 0.6791 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1210    0.1293    0.8971 [-0.2241; 0.2666 ] 
    ##   IMAG <~ imag2      0.4473      0.1442    3.1011    0.0019 [ 0.1521; 0.7165 ] 
    ##   IMAG <~ imag3      0.6020      0.1318    4.5695    0.0000 [ 0.3288; 0.8443 ] 
    ##   EXPE <~ expe1      0.2946      0.1200    2.4550    0.0141 [ 0.0660; 0.5372 ] 
    ##   EXPE <~ expe2      0.6473      0.0875    7.3979    0.0000 [ 0.4597; 0.7997 ] 
    ##   EXPE <~ expe3      0.2374      0.0910    2.6087    0.0091 [ 0.0453; 0.4088 ] 
    ##   QUAL <~ qual1      0.2370      0.0901    2.6312    0.0085 [ 0.0759; 0.4446 ] 
    ##   QUAL <~ qual2      0.4712      0.0817    5.7645    0.0000 [ 0.3028; 0.6208 ] 
    ##   QUAL <~ qual3      0.1831      0.0797    2.2969    0.0216 [ 0.0093; 0.3192 ] 
    ##   QUAL <~ qual4      0.1037      0.0553    1.8747    0.0608 [ 0.0046; 0.2184 ] 
    ##   QUAL <~ qual5      0.2049      0.0634    3.2300    0.0012 [ 0.0855; 0.3188 ] 
    ##   VAL <~ val1        0.7163      0.0964    7.4306    0.0000 [ 0.5140; 0.8833 ] 
    ##   VAL <~ val2        0.2202      0.0904    2.4351    0.0149 [ 0.0645; 0.3926 ] 
    ##   VAL <~ val3        0.2082      0.0578    3.5997    0.0003 [ 0.0921; 0.3172 ] 
    ##   SAT <~ sat1        0.3209      0.0156   20.6254    0.0000 [ 0.2941; 0.3556 ] 
    ##   SAT <~ sat2        0.3059      0.0139   22.0213    0.0000 [ 0.2829; 0.3359 ] 
    ##   SAT <~ sat3        0.2474      0.0112   22.1548    0.0000 [ 0.2246; 0.2677 ] 
    ##   SAT <~ sat4        0.2692      0.0120   22.3454    0.0000 [ 0.2458; 0.2907 ] 
    ##   LOY <~ loy1        0.3834      0.0257   14.9266    0.0000 [ 0.3329; 0.4331 ] 
    ##   LOY <~ loy2        0.2434      0.0284    8.5766    0.0000 [ 0.1835; 0.2932 ] 
    ##   LOY <~ loy3        0.3812      0.0275   13.8622    0.0000 [ 0.3314; 0.4356 ] 
    ##   LOY <~ loy4        0.2073      0.0362    5.7235    0.0000 [ 0.1380; 0.2747 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0652    9.8792    0.0000 [ 0.5086; 0.7492 ] 
    ##   imag1 ~~ imag3      0.5433      0.0709    7.6643    0.0000 [ 0.3999; 0.6694 ] 
    ##   imag2 ~~ imag3      0.7761      0.0393   19.7647    0.0000 [ 0.6944; 0.8477 ] 
    ##   expe1 ~~ expe2      0.5353      0.0610    8.7698    0.0000 [ 0.4183; 0.6445 ] 
    ##   expe1 ~~ expe3      0.4694      0.0621    7.5580    0.0000 [ 0.3483; 0.5890 ] 
    ##   expe2 ~~ expe3      0.5467      0.0630    8.6791    0.0000 [ 0.4106; 0.6579 ] 
    ##   qual1 ~~ qual2      0.6053      0.0566   10.6981    0.0000 [ 0.4935; 0.7081 ] 
    ##   qual1 ~~ qual3      0.5406      0.0611    8.8494    0.0000 [ 0.4160; 0.6527 ] 
    ##   qual1 ~~ qual4      0.5662      0.0676    8.3744    0.0000 [ 0.4304; 0.6927 ] 
    ##   qual1 ~~ qual5      0.5180      0.0685    7.5633    0.0000 [ 0.3868; 0.6567 ] 
    ##   qual2 ~~ qual3      0.6187      0.0567   10.9153    0.0000 [ 0.4950; 0.7125 ] 
    ##   qual2 ~~ qual4      0.6517      0.0638   10.2161    0.0000 [ 0.5145; 0.7646 ] 
    ##   qual2 ~~ qual5      0.6291      0.0577   10.8932    0.0000 [ 0.5066; 0.7336 ] 
    ##   qual3 ~~ qual4      0.4752      0.0637    7.4605    0.0000 [ 0.3387; 0.5945 ] 
    ##   qual3 ~~ qual5      0.5074      0.0606    8.3745    0.0000 [ 0.3809; 0.6163 ] 
    ##   qual4 ~~ qual5      0.6402      0.0599   10.6798    0.0000 [ 0.5061; 0.7407 ] 
    ##   val1 ~~ val2        0.6344      0.0543   11.6767    0.0000 [ 0.5179; 0.7265 ] 
    ##   val1 ~~ val3        0.4602      0.0700    6.5777    0.0000 [ 0.3173; 0.5976 ] 
    ##   val2 ~~ val3        0.6288      0.0624   10.0741    0.0000 [ 0.5027; 0.7493 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0689    6.8443    0.0000 [ 0.3449; 0.6041 ] 
    ##   QUAL ~ IMAG       0.3933      0.0637    6.1775    0.0000 [ 0.2759; 0.5178 ] 
    ##   QUAL ~ EXPE       0.8344      0.0227   36.7309    0.0000 [ 0.7885; 0.8769 ] 
    ##   VAL ~ IMAG        0.2974      0.0632    4.7082    0.0000 [ 0.1929; 0.4290 ] 
    ##   VAL ~ EXPE        0.6309      0.0513   12.2969    0.0000 [ 0.5336; 0.7286 ] 
    ##   VAL ~ QUAL        0.7013      0.0831    8.4368    0.0000 [ 0.5279; 0.8484 ] 
    ##   SAT ~ IMAG        0.4807      0.0699    6.8741    0.0000 [ 0.3353; 0.6077 ] 
    ##   SAT ~ EXPE        0.5001      0.0584    8.5648    0.0000 [ 0.3827; 0.6192 ] 
    ##   SAT ~ QUAL        0.5911      0.0930    6.3550    0.0000 [ 0.4125; 0.7734 ] 
    ##   SAT ~ VAL         0.5270      0.0925    5.6984    0.0000 [ 0.3447; 0.6968 ] 
    ##   LOY ~ IMAG        0.4840      0.0649    7.4532    0.0000 [ 0.3597; 0.6269 ] 
    ##   LOY ~ EXPE        0.3142      0.0577    5.4423    0.0000 [ 0.2053; 0.4326 ] 
    ##   LOY ~ QUAL        0.3714      0.0828    4.4854    0.0000 [ 0.2265; 0.5414 ] 
    ##   LOY ~ VAL         0.3311      0.0786    4.2121    0.0000 [ 0.1942; 0.5002 ] 
    ##   LOY ~ SAT         0.6283      0.0847    7.4168    0.0000 [ 0.4692; 0.7805 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0637    6.1775    0.0000 [ 0.2759; 0.5178 ] 
    ##   VAL ~ IMAG           0.2974      0.0632    4.7082    0.0000 [ 0.1929; 0.4290 ] 
    ##   VAL ~ EXPE           0.5852      0.0719    8.1338    0.0000 [ 0.4449; 0.7136 ] 
    ##   SAT ~ IMAG           0.2357      0.0497    4.7400    0.0000 [ 0.1489; 0.3488 ] 
    ##   SAT ~ EXPE           0.5173      0.0664    7.7877    0.0000 [ 0.3891; 0.6535 ] 
    ##   SAT ~ QUAL           0.3696      0.0661    5.5912    0.0000 [ 0.2472; 0.4909 ] 
    ##   LOY ~ IMAG           0.3020      0.0594    5.0860    0.0000 [ 0.1962; 0.4243 ] 
    ##   LOY ~ EXPE           0.3142      0.0577    5.4423    0.0000 [ 0.2053; 0.4326 ] 
    ##   LOY ~ QUAL           0.3714      0.0828    4.4854    0.0000 [ 0.2265; 0.5414 ] 
    ##   LOY ~ VAL            0.3311      0.0786    4.2121    0.0000 [ 0.1942; 0.5002 ] 
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
