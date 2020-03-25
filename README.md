
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
    ##                                     +------------------------------------------------------------------+
    ##                                     |                                                                  |
    ##                                     |   H0: The model-implied indicator covariance matrix equals the   |
    ##                                     |   population indicator covariance matrix.                        |
    ##                                     |                                                                  |
    ##                                     +------------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3180  
    ##  SRMR                    0.0940      0.0530  
    ##  dL                      2.2340      0.7114  
    ##  dML                     2.9219      1.6095  
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
    ##  Out of 499 bootstrap replications 479 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -419531147
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
    ##  RMS_theta_mi   = 0.05069299
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
    ##   expe1         1.4571         1.5708       1.9083         2.1011       0.0531
    ##   expe2         1.4103         1.4841       1.9320         2.0346       0.2016
    ##   expe3         1.6306         1.7271       2.1266         2.2231       0.1256
    ##   qual1         1.4748         1.5493       1.9276         2.0720       0.1148
    ##   qual2         1.5757         1.5366       2.0365         2.0612       0.2198
    ##   qual3         1.7309         1.7297       2.2226         2.2842       0.1200
    ##   qual4         1.2330         1.2026       1.5952         1.6365       0.2338
    ##   qual5         1.5074         1.5057       1.9371         1.9595       0.1962
    ##   val1          1.4440         1.3654       1.8670         1.7697       0.2516
    ##   val2          1.2265         1.2075       1.6472         1.7176       0.1736
    ##   val3          1.4791         1.3823       1.9671         1.9378       0.1482
    ##   sat1          1.2446         1.2327       1.6406         1.6184       0.3423
    ##   sat2          1.2296         1.1946       1.6353         1.6254       0.3117
    ##   sat3          1.3393         1.2755       1.6692         1.7199       0.2123
    ##   sat4          1.3158         1.2600       1.6647         1.6338       0.2790
    ##   loy1          1.6903         1.6587       2.2331         2.2261       0.2681
    ##   loy2          1.4853         1.4710       1.9122         1.9787       0.1308
    ##   loy3          1.7022         1.6682       2.2809         2.2708       0.2693
    ##   loy4          1.6907         1.6682       2.1782         2.3004       0.0867
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
    ##  Number of admissible results     = 486
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 2073546134
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
    ##   EXPE ~ IMAG      0.4714      0.0649    7.2647    0.0000 [ 0.3485; 0.5893 ] 
    ##   QUAL ~ EXPE      0.8344      0.0227   36.8114    0.0000 [ 0.7858; 0.8781 ] 
    ##   VAL ~ EXPE       0.0457      0.0868    0.5265    0.5986 [-0.1084; 0.2356 ] 
    ##   VAL ~ QUAL       0.7013      0.0842    8.3304    0.0000 [ 0.5309; 0.8536 ] 
    ##   SAT ~ IMAG       0.2450      0.0561    4.3695    0.0000 [ 0.1408; 0.3531 ] 
    ##   SAT ~ EXPE      -0.0172      0.0691   -0.2496    0.8029 [-0.1529; 0.1084 ] 
    ##   SAT ~ QUAL       0.2215      0.0981    2.2589    0.0239 [ 0.0426; 0.4349 ] 
    ##   SAT ~ VAL        0.5270      0.0847    6.2222    0.0000 [ 0.3589; 0.6775 ] 
    ##   LOY ~ IMAG       0.1819      0.0813    2.2368    0.0253 [ 0.0286; 0.3370 ] 
    ##   LOY ~ SAT        0.6283      0.0781    8.0415    0.0000 [ 0.4629; 0.7707 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1018    6.1973    0.0000 [ 0.3973; 0.7988 ] 
    ##   IMAG =~ imag2      0.9246      0.0390   23.7099    0.0000 [ 0.8281; 0.9739 ] 
    ##   IMAG =~ imag3      0.9577      0.0259   37.0269    0.0000 [ 0.8910; 0.9914 ] 
    ##   EXPE =~ expe1      0.7525      0.0793    9.4844    0.0000 [ 0.5586; 0.8796 ] 
    ##   EXPE =~ expe2      0.9348      0.0281   33.2967    0.0000 [ 0.8650; 0.9720 ] 
    ##   EXPE =~ expe3      0.7295      0.0731    9.9775    0.0000 [ 0.5602; 0.8477 ] 
    ##   QUAL =~ qual1      0.7861      0.0699   11.2511    0.0000 [ 0.6301; 0.8924 ] 
    ##   QUAL =~ qual2      0.9244      0.0228   40.5854    0.0000 [ 0.8687; 0.9580 ] 
    ##   QUAL =~ qual3      0.7560      0.0613   12.3387    0.0000 [ 0.6214; 0.8571 ] 
    ##   QUAL =~ qual4      0.7632      0.0532   14.3364    0.0000 [ 0.6379; 0.8499 ] 
    ##   QUAL =~ qual5      0.7834      0.0468   16.7380    0.0000 [ 0.6678; 0.8580 ] 
    ##   VAL =~ val1        0.9518      0.0239   39.7474    0.0000 [ 0.8962; 0.9860 ] 
    ##   VAL =~ val2        0.8056      0.0656   12.2803    0.0000 [ 0.6545; 0.9038 ] 
    ##   VAL =~ val3        0.6763      0.0773    8.7510    0.0000 [ 0.4935; 0.7997 ] 
    ##   SAT =~ sat1        0.9243      0.0235   39.3024    0.0000 [ 0.8729; 0.9646 ] 
    ##   SAT =~ sat2        0.8813      0.0283   31.0933    0.0000 [ 0.8207; 0.9321 ] 
    ##   SAT =~ sat3        0.7127      0.0499   14.2696    0.0000 [ 0.6078; 0.8053 ] 
    ##   SAT =~ sat4        0.7756      0.0477   16.2600    0.0000 [ 0.6853; 0.8737 ] 
    ##   LOY =~ loy1        0.9097      0.0506   17.9822    0.0000 [ 0.7928; 0.9870 ] 
    ##   LOY =~ loy2        0.5775      0.0813    7.1021    0.0000 [ 0.4015; 0.7240 ] 
    ##   LOY =~ loy3        0.9043      0.0404   22.4083    0.0000 [ 0.8101; 0.9778 ] 
    ##   LOY =~ loy4        0.4917      0.1017    4.8360    0.0000 [ 0.2868; 0.6783 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1180    0.1325    0.8946 [-0.2156; 0.2330 ] 
    ##   IMAG <~ imag2      0.4473      0.1449    3.0874    0.0020 [ 0.1335; 0.6918 ] 
    ##   IMAG <~ imag3      0.6020      0.1349    4.4629    0.0000 [ 0.3508; 0.8661 ] 
    ##   EXPE <~ expe1      0.2946      0.1194    2.4680    0.0136 [ 0.0634; 0.5150 ] 
    ##   EXPE <~ expe2      0.6473      0.0853    7.5906    0.0000 [ 0.4773; 0.7989 ] 
    ##   EXPE <~ expe3      0.2374      0.0916    2.5915    0.0096 [ 0.0446; 0.4113 ] 
    ##   QUAL <~ qual1      0.2370      0.0886    2.6749    0.0075 [ 0.0842; 0.4148 ] 
    ##   QUAL <~ qual2      0.4712      0.0812    5.8019    0.0000 [ 0.2937; 0.6311 ] 
    ##   QUAL <~ qual3      0.1831      0.0775    2.3631    0.0181 [ 0.0214; 0.3283 ] 
    ##   QUAL <~ qual4      0.1037      0.0614    1.6895    0.0911 [-0.0189; 0.2186 ] 
    ##   QUAL <~ qual5      0.2049      0.0630    3.2530    0.0011 [ 0.0672; 0.3134 ] 
    ##   VAL <~ val1        0.7163      0.0977    7.3337    0.0000 [ 0.5194; 0.8851 ] 
    ##   VAL <~ val2        0.2202      0.0944    2.3327    0.0197 [ 0.0514; 0.4054 ] 
    ##   VAL <~ val3        0.2082      0.0602    3.4564    0.0005 [ 0.0843; 0.3211 ] 
    ##   SAT <~ sat1        0.3209      0.0145   22.1969    0.0000 [ 0.2918; 0.3473 ] 
    ##   SAT <~ sat2        0.3059      0.0133   22.9528    0.0000 [ 0.2823; 0.3333 ] 
    ##   SAT <~ sat3        0.2474      0.0107   23.0840    0.0000 [ 0.2250; 0.2662 ] 
    ##   SAT <~ sat4        0.2692      0.0113   23.7277    0.0000 [ 0.2495; 0.2929 ] 
    ##   LOY <~ loy1        0.3834      0.0262   14.6543    0.0000 [ 0.3300; 0.4311 ] 
    ##   LOY <~ loy2        0.2434      0.0282    8.6388    0.0000 [ 0.1806; 0.2946 ] 
    ##   LOY <~ loy3        0.3812      0.0275   13.8740    0.0000 [ 0.3296; 0.4351 ] 
    ##   LOY <~ loy4        0.2073      0.0379    5.4678    0.0000 [ 0.1279; 0.2721 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0674    9.5547    0.0000 [ 0.4919; 0.7520 ] 
    ##   imag1 ~~ imag3      0.5433      0.0717    7.5742    0.0000 [ 0.3885; 0.6799 ] 
    ##   imag2 ~~ imag3      0.7761      0.0381   20.3452    0.0000 [ 0.6948; 0.8452 ] 
    ##   expe1 ~~ expe2      0.5353      0.0608    8.8102    0.0000 [ 0.4124; 0.6476 ] 
    ##   expe1 ~~ expe3      0.4694      0.0594    7.8980    0.0000 [ 0.3526; 0.5781 ] 
    ##   expe2 ~~ expe3      0.5467      0.0636    8.5926    0.0000 [ 0.4159; 0.6704 ] 
    ##   qual1 ~~ qual2      0.6053      0.0612    9.8914    0.0000 [ 0.4766; 0.7117 ] 
    ##   qual1 ~~ qual3      0.5406      0.0609    8.8750    0.0000 [ 0.4145; 0.6535 ] 
    ##   qual1 ~~ qual4      0.5662      0.0703    8.0526    0.0000 [ 0.4185; 0.6831 ] 
    ##   qual1 ~~ qual5      0.5180      0.0691    7.5005    0.0000 [ 0.3843; 0.6441 ] 
    ##   qual2 ~~ qual3      0.6187      0.0582   10.6388    0.0000 [ 0.5030; 0.7332 ] 
    ##   qual2 ~~ qual4      0.6517      0.0627   10.3873    0.0000 [ 0.5168; 0.7672 ] 
    ##   qual2 ~~ qual5      0.6291      0.0604   10.4188    0.0000 [ 0.4980; 0.7315 ] 
    ##   qual3 ~~ qual4      0.4752      0.0656    7.2398    0.0000 [ 0.3420; 0.5929 ] 
    ##   qual3 ~~ qual5      0.5074      0.0626    8.1109    0.0000 [ 0.3838; 0.6293 ] 
    ##   qual4 ~~ qual5      0.6402      0.0591   10.8307    0.0000 [ 0.5071; 0.7414 ] 
    ##   val1 ~~ val2        0.6344      0.0549   11.5579    0.0000 [ 0.5077; 0.7268 ] 
    ##   val1 ~~ val3        0.4602      0.0696    6.6117    0.0000 [ 0.3138; 0.5782 ] 
    ##   val2 ~~ val3        0.6288      0.0620   10.1502    0.0000 [ 0.5010; 0.7415 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0649    7.2647    0.0000 [ 0.3485; 0.5893 ] 
    ##   QUAL ~ IMAG       0.3933      0.0607    6.4848    0.0000 [ 0.2821; 0.5102 ] 
    ##   QUAL ~ EXPE       0.8344      0.0227   36.8114    0.0000 [ 0.7858; 0.8781 ] 
    ##   VAL ~ IMAG        0.2974      0.0607    4.9029    0.0000 [ 0.1865; 0.4248 ] 
    ##   VAL ~ EXPE        0.6309      0.0513   12.2983    0.0000 [ 0.5318; 0.7268 ] 
    ##   VAL ~ QUAL        0.7013      0.0842    8.3304    0.0000 [ 0.5309; 0.8536 ] 
    ##   SAT ~ IMAG        0.4807      0.0670    7.1713    0.0000 [ 0.3519; 0.6134 ] 
    ##   SAT ~ EXPE        0.5001      0.0583    8.5826    0.0000 [ 0.3985; 0.6120 ] 
    ##   SAT ~ QUAL        0.5911      0.0917    6.4480    0.0000 [ 0.4075; 0.7624 ] 
    ##   SAT ~ VAL         0.5270      0.0847    6.2222    0.0000 [ 0.3589; 0.6775 ] 
    ##   LOY ~ IMAG        0.4840      0.0705    6.8684    0.0000 [ 0.3530; 0.6273 ] 
    ##   LOY ~ EXPE        0.3142      0.0546    5.7544    0.0000 [ 0.2095; 0.4160 ] 
    ##   LOY ~ QUAL        0.3714      0.0793    4.6838    0.0000 [ 0.2208; 0.5245 ] 
    ##   LOY ~ VAL         0.3311      0.0727    4.5548    0.0000 [ 0.1921; 0.4843 ] 
    ##   LOY ~ SAT         0.6283      0.0781    8.0415    0.0000 [ 0.4629; 0.7707 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0607    6.4848    0.0000 [ 0.2821; 0.5102 ] 
    ##   VAL ~ IMAG           0.2974      0.0607    4.9029    0.0000 [ 0.1865; 0.4248 ] 
    ##   VAL ~ EXPE           0.5852      0.0719    8.1345    0.0000 [ 0.4455; 0.7245 ] 
    ##   SAT ~ IMAG           0.2357      0.0494    4.7727    0.0000 [ 0.1512; 0.3462 ] 
    ##   SAT ~ EXPE           0.5173      0.0677    7.6458    0.0000 [ 0.3897; 0.6537 ] 
    ##   SAT ~ QUAL           0.3696      0.0637    5.8033    0.0000 [ 0.2413; 0.4871 ] 
    ##   LOY ~ IMAG           0.3020      0.0541    5.5813    0.0000 [ 0.2043; 0.4192 ] 
    ##   LOY ~ EXPE           0.3142      0.0546    5.7544    0.0000 [ 0.2095; 0.4160 ] 
    ##   LOY ~ QUAL           0.3714      0.0793    4.6838    0.0000 [ 0.2208; 0.5245 ] 
    ##   LOY ~ VAL            0.3311      0.0727    4.5548    0.0000 [ 0.1921; 0.4843 ] 
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
