
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

  - `doFloodlightAnalysis()`: performs a floodlight analysis
  - `doSurfaceAnalysis()`: performs a surface analysis
  - `doRedundancyAnalysis()`: performs a redundancy analysis
  - `doIPMA()`: performs an importance-performance matrix analysis

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
    ##  dG                      0.6493      0.3258  
    ##  SRMR                    0.0940      0.0523  
    ##  dL                      2.2340      0.6920  
    ##  dML                     2.9219      1.6431  
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
    ##  The seed used was: -1940569770
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
    ##   expe1         1.4556         1.5663       1.9073         2.0911       0.0554
    ##   expe2         1.4125         1.4783       1.9350         2.0271       0.2009
    ##   expe3         1.6311         1.7242       2.1242         2.2196       0.1253
    ##   qual1         1.4779         1.5457       1.9303         2.0596       0.1146
    ##   qual2         1.5808         1.5376       2.0420         2.0605       0.2162
    ##   qual3         1.7331         1.7280       2.2235         2.2829       0.1191
    ##   qual4         1.2338         1.1980       1.5970         1.6310       0.2322
    ##   qual5         1.5049         1.5016       1.9364         1.9563       0.1954
    ##   val1          1.4473         1.3681       1.8731         1.7718       0.2477
    ##   val2          1.2281         1.2070       1.6487         1.7115       0.1737
    ##   val3          1.4806         1.3801       1.9693         1.9352       0.1471
    ##   sat1          1.2497         1.2345       1.6487         1.6200       0.3380
    ##   sat2          1.2333         1.1965       1.6409         1.6251       0.3095
    ##   sat3          1.3424         1.2758       1.6747         1.7192       0.2100
    ##   sat4          1.3181         1.2604       1.6690         1.6338       0.2766
    ##   loy1          1.6931         1.6609       2.2333         2.2235       0.2694
    ##   loy2          1.4837         1.4721       1.9095         1.9800       0.1336
    ##   loy3          1.7001         1.6641       2.2761         2.2617       0.2727
    ##   loy4          1.6919         1.6739       2.1813         2.3066       0.0832
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
    ##  Number of admissible results     = 489
    ##  Approach to handle inadmissibles = "drop"
    ##  Sign change option               = "none"
    ##  Random seed                      = 837157285
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
    ##   EXPE ~ IMAG      0.4714      0.0659    7.1553    0.0000 [ 0.3461; 0.5946 ] 
    ##   QUAL ~ EXPE      0.8344      0.0242   34.4325    0.0000 [ 0.7812; 0.8741 ] 
    ##   VAL ~ EXPE       0.0457      0.0853    0.5361    0.5919 [-0.1101; 0.2204 ] 
    ##   VAL ~ QUAL       0.7013      0.0833    8.4195    0.0000 [ 0.5290; 0.8539 ] 
    ##   SAT ~ IMAG       0.2450      0.0569    4.3038    0.0000 [ 0.1381; 0.3532 ] 
    ##   SAT ~ EXPE      -0.0172      0.0725   -0.2376    0.8122 [-0.1591; 0.1209 ] 
    ##   SAT ~ QUAL       0.2215      0.1012    2.1886    0.0286 [ 0.0360; 0.4252 ] 
    ##   SAT ~ VAL        0.5270      0.0866    6.0822    0.0000 [ 0.3534; 0.6846 ] 
    ##   LOY ~ IMAG       0.1819      0.0805    2.2601    0.0238 [ 0.0325; 0.3558 ] 
    ##   LOY ~ SAT        0.6283      0.0827    7.5978    0.0000 [ 0.4530; 0.7917 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0950    6.6379    0.0000 [ 0.4246; 0.8064 ] 
    ##   IMAG =~ imag2      0.9246      0.0417   22.1561    0.0000 [ 0.8271; 0.9783 ] 
    ##   IMAG =~ imag3      0.9577      0.0288   33.2636    0.0000 [ 0.8838; 0.9899 ] 
    ##   EXPE =~ expe1      0.7525      0.0778    9.6770    0.0000 [ 0.5715; 0.8667 ] 
    ##   EXPE =~ expe2      0.9348      0.0278   33.6346    0.0000 [ 0.8691; 0.9727 ] 
    ##   EXPE =~ expe3      0.7295      0.0696   10.4767    0.0000 [ 0.5662; 0.8423 ] 
    ##   QUAL =~ qual1      0.7861      0.0691   11.3804    0.0000 [ 0.6169; 0.8831 ] 
    ##   QUAL =~ qual2      0.9244      0.0237   39.0137    0.0000 [ 0.8710; 0.9619 ] 
    ##   QUAL =~ qual3      0.7560      0.0573   13.2014    0.0000 [ 0.6266; 0.8471 ] 
    ##   QUAL =~ qual4      0.7632      0.0564   13.5385    0.0000 [ 0.6400; 0.8536 ] 
    ##   QUAL =~ qual5      0.7834      0.0442   17.7072    0.0000 [ 0.6790; 0.8589 ] 
    ##   VAL =~ val1        0.9518      0.0237   40.1938    0.0000 [ 0.8967; 0.9855 ] 
    ##   VAL =~ val2        0.8056      0.0636   12.6744    0.0000 [ 0.6626; 0.9078 ] 
    ##   VAL =~ val3        0.6763      0.0750    9.0194    0.0000 [ 0.5358; 0.8154 ] 
    ##   SAT =~ sat1        0.9243      0.0228   40.5148    0.0000 [ 0.8774; 0.9665 ] 
    ##   SAT =~ sat2        0.8813      0.0296   29.7528    0.0000 [ 0.8090; 0.9290 ] 
    ##   SAT =~ sat3        0.7127      0.0520   13.7117    0.0000 [ 0.6073; 0.8036 ] 
    ##   SAT =~ sat4        0.7756      0.0472   16.4183    0.0000 [ 0.6826; 0.8624 ] 
    ##   LOY =~ loy1        0.9097      0.0496   18.3336    0.0000 [ 0.7808; 0.9819 ] 
    ##   LOY =~ loy2        0.5775      0.0821    7.0341    0.0000 [ 0.4067; 0.7364 ] 
    ##   LOY =~ loy3        0.9043      0.0432   20.9283    0.0000 [ 0.8158; 0.9773 ] 
    ##   LOY =~ loy4        0.4917      0.1011    4.8635    0.0000 [ 0.2809; 0.6790 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1115    0.1403    0.8884 [-0.2015; 0.2365 ] 
    ##   IMAG <~ imag2      0.4473      0.1535    2.9147    0.0036 [ 0.1550; 0.7323 ] 
    ##   IMAG <~ imag3      0.6020      0.1460    4.1228    0.0000 [ 0.2821; 0.8310 ] 
    ##   EXPE <~ expe1      0.2946      0.1184    2.4879    0.0128 [ 0.0683; 0.5200 ] 
    ##   EXPE <~ expe2      0.6473      0.0846    7.6475    0.0000 [ 0.4781; 0.8004 ] 
    ##   EXPE <~ expe3      0.2374      0.0914    2.5964    0.0094 [ 0.0536; 0.4100 ] 
    ##   QUAL <~ qual1      0.2370      0.0901    2.6310    0.0085 [ 0.0568; 0.4035 ] 
    ##   QUAL <~ qual2      0.4712      0.0828    5.6893    0.0000 [ 0.3273; 0.6367 ] 
    ##   QUAL <~ qual3      0.1831      0.0752    2.4350    0.0149 [ 0.0342; 0.3180 ] 
    ##   QUAL <~ qual4      0.1037      0.0646    1.6052    0.1084 [-0.0146; 0.2468 ] 
    ##   QUAL <~ qual5      0.2049      0.0638    3.2123    0.0013 [ 0.0671; 0.3212 ] 
    ##   VAL <~ val1        0.7163      0.0976    7.3379    0.0000 [ 0.4989; 0.8804 ] 
    ##   VAL <~ val2        0.2202      0.0933    2.3605    0.0183 [ 0.0492; 0.4209 ] 
    ##   VAL <~ val3        0.2082      0.0601    3.4614    0.0005 [ 0.0882; 0.3289 ] 
    ##   SAT <~ sat1        0.3209      0.0150   21.3646    0.0000 [ 0.2952; 0.3524 ] 
    ##   SAT <~ sat2        0.3059      0.0134   22.7948    0.0000 [ 0.2814; 0.3339 ] 
    ##   SAT <~ sat3        0.2474      0.0107   23.0978    0.0000 [ 0.2262; 0.2672 ] 
    ##   SAT <~ sat4        0.2692      0.0119   22.6482    0.0000 [ 0.2459; 0.2940 ] 
    ##   LOY <~ loy1        0.3834      0.0255   15.0455    0.0000 [ 0.3347; 0.4314 ] 
    ##   LOY <~ loy2        0.2434      0.0287    8.4797    0.0000 [ 0.1842; 0.2980 ] 
    ##   LOY <~ loy3        0.3812      0.0274   13.9190    0.0000 [ 0.3298; 0.4376 ] 
    ##   LOY <~ loy4        0.2073      0.0382    5.4293    0.0000 [ 0.1265; 0.2766 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0623   10.3309    0.0000 [ 0.5183; 0.7611 ] 
    ##   imag1 ~~ imag3      0.5433      0.0688    7.8981    0.0000 [ 0.4089; 0.6774 ] 
    ##   imag2 ~~ imag3      0.7761      0.0394   19.7008    0.0000 [ 0.6867; 0.8483 ] 
    ##   expe1 ~~ expe2      0.5353      0.0568    9.4304    0.0000 [ 0.4069; 0.6391 ] 
    ##   expe1 ~~ expe3      0.4694      0.0596    7.8802    0.0000 [ 0.3546; 0.5793 ] 
    ##   expe2 ~~ expe3      0.5467      0.0571    9.5684    0.0000 [ 0.4327; 0.6541 ] 
    ##   qual1 ~~ qual2      0.6053      0.0566   10.6936    0.0000 [ 0.4879; 0.7032 ] 
    ##   qual1 ~~ qual3      0.5406      0.0586    9.2282    0.0000 [ 0.4245; 0.6487 ] 
    ##   qual1 ~~ qual4      0.5662      0.0671    8.4412    0.0000 [ 0.4346; 0.6977 ] 
    ##   qual1 ~~ qual5      0.5180      0.0666    7.7730    0.0000 [ 0.3857; 0.6474 ] 
    ##   qual2 ~~ qual3      0.6187      0.0524   11.8084    0.0000 [ 0.5064; 0.7169 ] 
    ##   qual2 ~~ qual4      0.6517      0.0646   10.0826    0.0000 [ 0.5204; 0.7619 ] 
    ##   qual2 ~~ qual5      0.6291      0.0572   10.9948    0.0000 [ 0.5138; 0.7335 ] 
    ##   qual3 ~~ qual4      0.4752      0.0627    7.5755    0.0000 [ 0.3584; 0.5884 ] 
    ##   qual3 ~~ qual5      0.5074      0.0617    8.2186    0.0000 [ 0.3960; 0.6269 ] 
    ##   qual4 ~~ qual5      0.6402      0.0536   11.9410    0.0000 [ 0.5227; 0.7358 ] 
    ##   val1 ~~ val2        0.6344      0.0527   12.0356    0.0000 [ 0.5292; 0.7346 ] 
    ##   val1 ~~ val3        0.4602      0.0717    6.4174    0.0000 [ 0.3256; 0.6100 ] 
    ##   val2 ~~ val3        0.6288      0.0622   10.1138    0.0000 [ 0.5084; 0.7414 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0659    7.1553    0.0000 [ 0.3461; 0.5946 ] 
    ##   QUAL ~ IMAG       0.3933      0.0610    6.4435    0.0000 [ 0.2824; 0.5051 ] 
    ##   QUAL ~ EXPE       0.8344      0.0242   34.4325    0.0000 [ 0.7812; 0.8741 ] 
    ##   VAL ~ IMAG        0.2974      0.0602    4.9413    0.0000 [ 0.1953; 0.4187 ] 
    ##   VAL ~ EXPE        0.6309      0.0490   12.8739    0.0000 [ 0.5366; 0.7235 ] 
    ##   VAL ~ QUAL        0.7013      0.0833    8.4195    0.0000 [ 0.5290; 0.8539 ] 
    ##   SAT ~ IMAG        0.4807      0.0675    7.1245    0.0000 [ 0.3521; 0.6026 ] 
    ##   SAT ~ EXPE        0.5001      0.0571    8.7571    0.0000 [ 0.3876; 0.6036 ] 
    ##   SAT ~ QUAL        0.5911      0.0965    6.1226    0.0000 [ 0.3954; 0.7850 ] 
    ##   SAT ~ VAL         0.5270      0.0866    6.0822    0.0000 [ 0.3534; 0.6846 ] 
    ##   LOY ~ IMAG        0.4840      0.0683    7.0853    0.0000 [ 0.3515; 0.6150 ] 
    ##   LOY ~ EXPE        0.3142      0.0543    5.7823    0.0000 [ 0.2111; 0.4163 ] 
    ##   LOY ~ QUAL        0.3714      0.0860    4.3172    0.0000 [ 0.2168; 0.5638 ] 
    ##   LOY ~ VAL         0.3311      0.0749    4.4207    0.0000 [ 0.1882; 0.4777 ] 
    ##   LOY ~ SAT         0.6283      0.0827    7.5978    0.0000 [ 0.4530; 0.7917 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0610    6.4435    0.0000 [ 0.2824; 0.5051 ] 
    ##   VAL ~ IMAG           0.2974      0.0602    4.9413    0.0000 [ 0.1953; 0.4187 ] 
    ##   VAL ~ EXPE           0.5852      0.0732    7.9990    0.0000 [ 0.4422; 0.7268 ] 
    ##   SAT ~ IMAG           0.2357      0.0484    4.8730    0.0000 [ 0.1529; 0.3384 ] 
    ##   SAT ~ EXPE           0.5173      0.0671    7.7148    0.0000 [ 0.3888; 0.6620 ] 
    ##   SAT ~ QUAL           0.3696      0.0640    5.7732    0.0000 [ 0.2497; 0.4890 ] 
    ##   LOY ~ IMAG           0.3020      0.0547    5.5260    0.0000 [ 0.2025; 0.4138 ] 
    ##   LOY ~ EXPE           0.3142      0.0543    5.7823    0.0000 [ 0.2111; 0.4163 ] 
    ##   LOY ~ QUAL           0.3714      0.0860    4.3172    0.0000 [ 0.2168; 0.5638 ] 
    ##   LOY ~ VAL            0.3311      0.0749    4.4207    0.0000 [ 0.1882; 0.4777 ] 
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
