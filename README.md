
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
    ##                       +------------------------------------------------------------------+
    ##                       |                                                                  |
    ##                       |   H0: The model-implied indicator covariance matrix equals the   |
    ##                       |   population indicator covariance matrix.                        |
    ##                       |                                                                  |
    ##                       +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3180  
    ##  SRMR                    0.0940      0.0532  
    ##  dL                      2.2340      0.7166  
    ##  dML                     2.9219      1.6041  
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
    ##  Out of 499 bootstrap replications 483 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -538581867
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
    ##   expe1         1.4569         1.5683       1.9102         2.0964       0.0513
    ##   expe2         1.4150         1.4821       1.9365         2.0302       0.1972
    ##   expe3         1.6318         1.7227       2.1273         2.2192       0.1224
    ##   qual1         1.4792         1.5497       1.9325         2.0658       0.1111
    ##   qual2         1.5805         1.5376       2.0409         2.0587       0.2151
    ##   qual3         1.7333         1.7217       2.2241         2.2735       0.1176
    ##   qual4         1.2330         1.1968       1.5969         1.6266       0.2322
    ##   qual5         1.5054         1.4996       1.9353         1.9505       0.1953
    ##   val1          1.4468         1.3660       1.8701         1.7654       0.2481
    ##   val2          1.2262         1.2025       1.6481         1.7097       0.1730
    ##   val3          1.4804         1.3790       1.9684         1.9333       0.1481
    ##   sat1          1.2458         1.2313       1.6444         1.6174       0.3402
    ##   sat2          1.2330         1.1957       1.6412         1.6269       0.3082
    ##   sat3          1.3408         1.2760       1.6738         1.7197       0.2088
    ##   sat4          1.3190         1.2613       1.6692         1.6358       0.2752
    ##   loy1          1.6946         1.6643       2.2365         2.2309       0.2664
    ##   loy2          1.4873         1.4755       1.9143         1.9828       0.1298
    ##   loy3          1.7049         1.6701       2.2841         2.2727       0.2682
    ##   loy4          1.6896         1.6706       2.1788         2.3002       0.0859
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
    ##  Number of admissible results     = 483
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 2117431115
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
    ##   EXPE ~ IMAG      0.4714      0.0667    7.0668    0.0000 [ 0.3445; 0.6060 ] 
    ##   QUAL ~ EXPE      0.8344      0.0221   37.7484    0.0000 [ 0.7895; 0.8741 ] 
    ##   VAL ~ EXPE       0.0457      0.0857    0.5334    0.5938 [-0.1076; 0.2214 ] 
    ##   VAL ~ QUAL       0.7013      0.0812    8.6400    0.0000 [ 0.5256; 0.8566 ] 
    ##   SAT ~ IMAG       0.2450      0.0569    4.3031    0.0000 [ 0.1287; 0.3620 ] 
    ##   SAT ~ EXPE      -0.0172      0.0756   -0.2280    0.8197 [-0.1642; 0.1309 ] 
    ##   SAT ~ QUAL       0.2215      0.1048    2.1132    0.0346 [ 0.0183; 0.4252 ] 
    ##   SAT ~ VAL        0.5270      0.0883    5.9711    0.0000 [ 0.3516; 0.7049 ] 
    ##   LOY ~ IMAG       0.1819      0.0765    2.3771    0.0174 [ 0.0359; 0.3470 ] 
    ##   LOY ~ SAT        0.6283      0.0793    7.9243    0.0000 [ 0.4715; 0.7791 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0982    6.4233    0.0000 [ 0.4375; 0.8012 ] 
    ##   IMAG =~ imag2      0.9246      0.0398   23.2131    0.0000 [ 0.8225; 0.9773 ] 
    ##   IMAG =~ imag3      0.9577      0.0292   32.7473    0.0000 [ 0.8769; 0.9915 ] 
    ##   EXPE =~ expe1      0.7525      0.0769    9.7796    0.0000 [ 0.5554; 0.8717 ] 
    ##   EXPE =~ expe2      0.9348      0.0265   35.2804    0.0000 [ 0.8748; 0.9702 ] 
    ##   EXPE =~ expe3      0.7295      0.0726   10.0536    0.0000 [ 0.5608; 0.8453 ] 
    ##   QUAL =~ qual1      0.7861      0.0677   11.6162    0.0000 [ 0.6265; 0.8894 ] 
    ##   QUAL =~ qual2      0.9244      0.0218   42.3598    0.0000 [ 0.8681; 0.9579 ] 
    ##   QUAL =~ qual3      0.7560      0.0609   12.4104    0.0000 [ 0.6176; 0.8437 ] 
    ##   QUAL =~ qual4      0.7632      0.0517   14.7668    0.0000 [ 0.6505; 0.8454 ] 
    ##   QUAL =~ qual5      0.7834      0.0496   15.7861    0.0000 [ 0.6717; 0.8590 ] 
    ##   VAL =~ val1        0.9518      0.0236   40.3716    0.0000 [ 0.8977; 0.9851 ] 
    ##   VAL =~ val2        0.8056      0.0632   12.7495    0.0000 [ 0.6666; 0.9036 ] 
    ##   VAL =~ val3        0.6763      0.0735    9.1982    0.0000 [ 0.5149; 0.8047 ] 
    ##   SAT =~ sat1        0.9243      0.0229   40.4263    0.0000 [ 0.8704; 0.9604 ] 
    ##   SAT =~ sat2        0.8813      0.0286   30.8132    0.0000 [ 0.8091; 0.9229 ] 
    ##   SAT =~ sat3        0.7127      0.0534   13.3417    0.0000 [ 0.6011; 0.8184 ] 
    ##   SAT =~ sat4        0.7756      0.0531   14.6148    0.0000 [ 0.6627; 0.8720 ] 
    ##   LOY =~ loy1        0.9097      0.0484   18.8102    0.0000 [ 0.7973; 0.9812 ] 
    ##   LOY =~ loy2        0.5775      0.0870    6.6378    0.0000 [ 0.3861; 0.7207 ] 
    ##   LOY =~ loy3        0.9043      0.0413   21.9092    0.0000 [ 0.8089; 0.9723 ] 
    ##   LOY =~ loy4        0.4917      0.1004    4.8990    0.0000 [ 0.2879; 0.6791 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1161    0.1347    0.8928 [-0.1948; 0.2335 ] 
    ##   IMAG <~ imag2      0.4473      0.1493    2.9968    0.0027 [ 0.1615; 0.7375 ] 
    ##   IMAG <~ imag3      0.6020      0.1406    4.2822    0.0000 [ 0.3014; 0.8617 ] 
    ##   EXPE <~ expe1      0.2946      0.1173    2.5121    0.0120 [ 0.0657; 0.5251 ] 
    ##   EXPE <~ expe2      0.6473      0.0801    8.0836    0.0000 [ 0.4955; 0.7849 ] 
    ##   EXPE <~ expe3      0.2374      0.0924    2.5679    0.0102 [ 0.0443; 0.4082 ] 
    ##   QUAL <~ qual1      0.2370      0.0895    2.6472    0.0081 [ 0.0658; 0.4191 ] 
    ##   QUAL <~ qual2      0.4712      0.0760    6.1963    0.0000 [ 0.3118; 0.6167 ] 
    ##   QUAL <~ qual3      0.1831      0.0793    2.3090    0.0209 [ 0.0223; 0.3313 ] 
    ##   QUAL <~ qual4      0.1037      0.0599    1.7322    0.0832 [-0.0091; 0.2148 ] 
    ##   QUAL <~ qual5      0.2049      0.0617    3.3227    0.0009 [ 0.0672; 0.3122 ] 
    ##   VAL <~ val1        0.7163      0.0958    7.4774    0.0000 [ 0.5312; 0.8744 ] 
    ##   VAL <~ val2        0.2202      0.0940    2.3425    0.0192 [ 0.0529; 0.4030 ] 
    ##   VAL <~ val3        0.2082      0.0588    3.5377    0.0004 [ 0.0882; 0.3189 ] 
    ##   SAT <~ sat1        0.3209      0.0155   20.7562    0.0000 [ 0.2931; 0.3549 ] 
    ##   SAT <~ sat2        0.3059      0.0148   20.7158    0.0000 [ 0.2808; 0.3358 ] 
    ##   SAT <~ sat3        0.2474      0.0112   22.0214    0.0000 [ 0.2260; 0.2709 ] 
    ##   SAT <~ sat4        0.2692      0.0120   22.5268    0.0000 [ 0.2452; 0.2913 ] 
    ##   LOY <~ loy1        0.3834      0.0266   14.4067    0.0000 [ 0.3314; 0.4316 ] 
    ##   LOY <~ loy2        0.2434      0.0305    7.9715    0.0000 [ 0.1737; 0.2953 ] 
    ##   LOY <~ loy3        0.3812      0.0276   13.7990    0.0000 [ 0.3321; 0.4423 ] 
    ##   LOY <~ loy4        0.2073      0.0372    5.5702    0.0000 [ 0.1290; 0.2756 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0635   10.1386    0.0000 [ 0.5062; 0.7547 ] 
    ##   imag1 ~~ imag3      0.5433      0.0659    8.2433    0.0000 [ 0.4144; 0.6720 ] 
    ##   imag2 ~~ imag3      0.7761      0.0388   19.9935    0.0000 [ 0.6917; 0.8489 ] 
    ##   expe1 ~~ expe2      0.5353      0.0574    9.3318    0.0000 [ 0.4220; 0.6499 ] 
    ##   expe1 ~~ expe3      0.4694      0.0613    7.6537    0.0000 [ 0.3525; 0.5891 ] 
    ##   expe2 ~~ expe3      0.5467      0.0594    9.1966    0.0000 [ 0.4226; 0.6567 ] 
    ##   qual1 ~~ qual2      0.6053      0.0571   10.5971    0.0000 [ 0.4881; 0.7095 ] 
    ##   qual1 ~~ qual3      0.5406      0.0599    9.0328    0.0000 [ 0.4223; 0.6533 ] 
    ##   qual1 ~~ qual4      0.5662      0.0659    8.5982    0.0000 [ 0.4344; 0.6919 ] 
    ##   qual1 ~~ qual5      0.5180      0.0682    7.6014    0.0000 [ 0.3740; 0.6474 ] 
    ##   qual2 ~~ qual3      0.6187      0.0539   11.4767    0.0000 [ 0.5045; 0.7172 ] 
    ##   qual2 ~~ qual4      0.6517      0.0613   10.6268    0.0000 [ 0.5221; 0.7566 ] 
    ##   qual2 ~~ qual5      0.6291      0.0614   10.2380    0.0000 [ 0.4961; 0.7342 ] 
    ##   qual3 ~~ qual4      0.4752      0.0653    7.2719    0.0000 [ 0.3383; 0.5843 ] 
    ##   qual3 ~~ qual5      0.5074      0.0638    7.9523    0.0000 [ 0.3741; 0.6184 ] 
    ##   qual4 ~~ qual5      0.6402      0.0566   11.3196    0.0000 [ 0.5192; 0.7368 ] 
    ##   val1 ~~ val2        0.6344      0.0518   12.2471    0.0000 [ 0.5311; 0.7272 ] 
    ##   val1 ~~ val3        0.4602      0.0678    6.7915    0.0000 [ 0.3336; 0.5926 ] 
    ##   val2 ~~ val3        0.6288      0.0617   10.1961    0.0000 [ 0.5137; 0.7446 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0667    7.0668    0.0000 [ 0.3445; 0.6060 ] 
    ##   QUAL ~ IMAG       0.3933      0.0614    6.4025    0.0000 [ 0.2757; 0.5148 ] 
    ##   QUAL ~ EXPE       0.8344      0.0221   37.7484    0.0000 [ 0.7895; 0.8741 ] 
    ##   VAL ~ IMAG        0.2974      0.0611    4.8644    0.0000 [ 0.1886; 0.4305 ] 
    ##   VAL ~ EXPE        0.6309      0.0494   12.7681    0.0000 [ 0.5375; 0.7316 ] 
    ##   VAL ~ QUAL        0.7013      0.0812    8.6400    0.0000 [ 0.5256; 0.8566 ] 
    ##   SAT ~ IMAG        0.4807      0.0680    7.0674    0.0000 [ 0.3686; 0.6198 ] 
    ##   SAT ~ EXPE        0.5001      0.0571    8.7557    0.0000 [ 0.3943; 0.6127 ] 
    ##   SAT ~ QUAL        0.5911      0.0994    5.9490    0.0000 [ 0.4031; 0.7785 ] 
    ##   SAT ~ VAL         0.5270      0.0883    5.9711    0.0000 [ 0.3516; 0.7049 ] 
    ##   LOY ~ IMAG        0.4840      0.0623    7.7666    0.0000 [ 0.3682; 0.6150 ] 
    ##   LOY ~ EXPE        0.3142      0.0538    5.8350    0.0000 [ 0.2185; 0.4176 ] 
    ##   LOY ~ QUAL        0.3714      0.0827    4.4884    0.0000 [ 0.2240; 0.5350 ] 
    ##   LOY ~ VAL         0.3311      0.0734    4.5129    0.0000 [ 0.1984; 0.4718 ] 
    ##   LOY ~ SAT         0.6283      0.0793    7.9243    0.0000 [ 0.4715; 0.7791 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0614    6.4025    0.0000 [ 0.2757; 0.5148 ] 
    ##   VAL ~ IMAG           0.2974      0.0611    4.8644    0.0000 [ 0.1886; 0.4305 ] 
    ##   VAL ~ EXPE           0.5852      0.0684    8.5574    0.0000 [ 0.4441; 0.7215 ] 
    ##   SAT ~ IMAG           0.2357      0.0490    4.8126    0.0000 [ 0.1544; 0.3473 ] 
    ##   SAT ~ EXPE           0.5173      0.0688    7.5240    0.0000 [ 0.3834; 0.6480 ] 
    ##   SAT ~ QUAL           0.3696      0.0633    5.8387    0.0000 [ 0.2391; 0.4939 ] 
    ##   LOY ~ IMAG           0.3020      0.0577    5.2308    0.0000 [ 0.2097; 0.4268 ] 
    ##   LOY ~ EXPE           0.3142      0.0538    5.8350    0.0000 [ 0.2185; 0.4176 ] 
    ##   LOY ~ QUAL           0.3714      0.0827    4.4884    0.0000 [ 0.2240; 0.5350 ] 
    ##   LOY ~ VAL            0.3311      0.0734    4.5129    0.0000 [ 0.1984; 0.4718 ] 
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
