
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
    ##  Second order approach            = NA
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
    ##                                   +------------------------------------------------------------+
    ##                                   |                                                            |
    ##                                   |   H0: Population indicator covariance matrix is equal to   |
    ##                                   |   model-implied indicator covariance matrix.               |
    ##                                   |                                                            |
    ##                                   +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3269  
    ##  SRMR                    0.0940      0.0522  
    ##  dL                      2.2340      0.6896  
    ##  
    ## 
    ## Decision: 
    ## 
    ##                          Significance level
    ##  Distance measure          95%   
    ##  dG                      reject  
    ##  SRMR                    reject  
    ##  dL                      reject  
    ##  
    ## Additonal information:
    ## 
    ##  Out of 499 bootstrap replications 482 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -702134351
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
    ##  Chi_square_df  = 4.181386
    ##  CFI            = 0.8573048
    ##  GFI            = 0.9642375
    ##  IFI            = 0.8593711
    ##  NFI            = 0.8229918
    ##  NNFI           = 0.8105598
    ##  RMSEA          = 0.1130338
    ##  RMS_theta      = 0.05069299
    ##  RMS_theta_mi   = 0.05069299
    ##  SRMR           = 0.09396871
    ## 
    ##  Degrees of freedom    = 174
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
    ##  Independent construct   Effect size
    ##  IMAG                      0.2856   
    ## 
    ##   Dependent construct: 'QUAL'
    ## 
    ##  Independent construct   Effect size
    ##  EXPE                      2.2928   
    ## 
    ##   Dependent construct: 'VAL'
    ## 
    ##  Independent construct   Effect size
    ##  EXPE                      0.0014   
    ##  QUAL                      0.3301   
    ## 
    ##   Dependent construct: 'SAT'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      0.1462   
    ##  EXPE                      0.0004   
    ##  QUAL                      0.0468   
    ##  VAL                       0.4373   
    ## 
    ##   Dependent construct: 'LOY'
    ## 
    ##  Independent construct   Effect size
    ##  IMAG                      0.0414   
    ##  SAT                       0.4938   
    ## 
    ## ------------------------------ Validity assessment -----------------------------
    ## 
    ##  Heterotrait-montrait ratio of correlation matrix (HTMT matrix)
    ## 
    ##           SAT LOY
    ## SAT 0.0000000   0
    ## LOY 0.7432489   0
    ## 
    ## 
    ##  Fornell-Larcker matrix
    ## 
    ##           SAT       LOY
    ## SAT 0.6851491 0.5696460
    ## LOY 0.5696460 0.5551718
    ## 
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
    ##   expe1         1.4566         1.5683       1.9071         2.0925       0.0536
    ##   expe2         1.4113         1.4788       1.9317         2.0267       0.2027
    ##   expe3         1.6341         1.7275       2.1285         2.2238       0.1252
    ##   qual1         1.4767         1.5473       1.9281         2.0602       0.1154
    ##   qual2         1.5769         1.5350       2.0369         2.0563       0.2202
    ##   qual3         1.7315         1.7282       2.2228         2.2811       0.1204
    ##   qual4         1.2346         1.1980       1.5972         1.6317       0.2336
    ##   qual5         1.5059         1.5020       1.9366         1.9568       0.1977
    ##   val1          1.4481         1.3653       1.8706         1.7669       0.2514
    ##   val2          1.2270         1.2039       1.6482         1.7105       0.1749
    ##   val3          1.4810         1.3803       1.9691         1.9346       0.1487
    ##   sat1          1.2449         1.2329       1.6451         1.6200       0.3409
    ##   sat2          1.2334         1.1967       1.6415         1.6276       0.3100
    ##   sat3          1.3413         1.2774       1.6735         1.7246       0.2112
    ##   sat4          1.3172         1.2621       1.6681         1.6363       0.2784
    ##   loy1          1.6915         1.6617       2.2341         2.2315       0.2684
    ##   loy2          1.4845         1.4765       1.9121         1.9826       0.1324
    ##   loy3          1.7018         1.6692       2.2818         2.2742       0.2703
    ##   loy4          1.6914         1.6687       2.1810         2.2964       0.0843
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
    ##  Second order approach            = NA
    ##  Type of path model               = Linear
    ##  Disattenuated                    = Yes (PLSc)
    ## 
    ##  Resample information:
    ##  ---------------------
    ##  Resample methode                 = bootstrap
    ##  Number of resamples              = 499
    ##  Number of admissible results     = 481
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = -1766939280
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
    ##   EXPE ~ IMAG      0.4714      0.0684    6.8917    0.0000 [ 0.3426; 0.6091 ] 
    ##   QUAL ~ EXPE      0.8344      0.0242   34.4109    0.0000 [ 0.7811; 0.8787 ] 
    ##   VAL ~ EXPE       0.0457      0.0850    0.5378    0.5907 [-0.1060; 0.2240 ] 
    ##   VAL ~ QUAL       0.7013      0.0795    8.8224    0.0000 [ 0.5300; 0.8448 ] 
    ##   SAT ~ IMAG       0.2450      0.0600    4.0828    0.0000 [ 0.1372; 0.3677 ] 
    ##   SAT ~ EXPE      -0.0172      0.0740   -0.2328    0.8159 [-0.1589; 0.1335 ] 
    ##   SAT ~ QUAL       0.2215      0.1020    2.1720    0.0299 [ 0.0451; 0.4224 ] 
    ##   SAT ~ VAL        0.5270      0.0870    6.0600    0.0000 [ 0.3577; 0.6931 ] 
    ##   LOY ~ IMAG       0.1819      0.0779    2.3360    0.0195 [ 0.0307; 0.3405 ] 
    ##   LOY ~ SAT        0.6283      0.0778    8.0712    0.0000 [ 0.4671; 0.7704 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.1015    6.2137    0.0000 [ 0.4161; 0.8085 ] 
    ##   IMAG =~ imag2      0.9246      0.0414   22.3411    0.0000 [ 0.8145; 0.9785 ] 
    ##   IMAG =~ imag3      0.9577      0.0300   31.9304    0.0000 [ 0.8776; 0.9939 ] 
    ##   EXPE =~ expe1      0.7525      0.0808    9.3145    0.0000 [ 0.5641; 0.8831 ] 
    ##   EXPE =~ expe2      0.9348      0.0290   32.2878    0.0000 [ 0.8588; 0.9744 ] 
    ##   EXPE =~ expe3      0.7295      0.0716   10.1948    0.0000 [ 0.5540; 0.8409 ] 
    ##   QUAL =~ qual1      0.7861      0.0709   11.0947    0.0000 [ 0.6178; 0.8903 ] 
    ##   QUAL =~ qual2      0.9244      0.0231   39.9508    0.0000 [ 0.8688; 0.9579 ] 
    ##   QUAL =~ qual3      0.7560      0.0598   12.6433    0.0000 [ 0.6207; 0.8543 ] 
    ##   QUAL =~ qual4      0.7632      0.0559   13.6628    0.0000 [ 0.6340; 0.8550 ] 
    ##   QUAL =~ qual5      0.7834      0.0476   16.4718    0.0000 [ 0.6756; 0.8667 ] 
    ##   VAL =~ val1        0.9518      0.0235   40.4248    0.0000 [ 0.8963; 0.9850 ] 
    ##   VAL =~ val2        0.8056      0.0639   12.6096    0.0000 [ 0.6723; 0.9096 ] 
    ##   VAL =~ val3        0.6763      0.0709    9.5446    0.0000 [ 0.5153; 0.8048 ] 
    ##   SAT =~ sat1        0.9243      0.0226   40.8414    0.0000 [ 0.8736; 0.9606 ] 
    ##   SAT =~ sat2        0.8813      0.0282   31.2392    0.0000 [ 0.8186; 0.9293 ] 
    ##   SAT =~ sat3        0.7127      0.0518   13.7549    0.0000 [ 0.6014; 0.8092 ] 
    ##   SAT =~ sat4        0.7756      0.0492   15.7656    0.0000 [ 0.6669; 0.8637 ] 
    ##   LOY =~ loy1        0.9097      0.0511   17.8088    0.0000 [ 0.7824; 0.9869 ] 
    ##   LOY =~ loy2        0.5775      0.0903    6.3952    0.0000 [ 0.3852; 0.7339 ] 
    ##   LOY =~ loy3        0.9043      0.0406   22.2853    0.0000 [ 0.8235; 0.9735 ] 
    ##   LOY =~ loy4        0.4917      0.0962    5.1100    0.0000 [ 0.2956; 0.6679 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1159    0.1350    0.8926 [-0.1971; 0.2500 ] 
    ##   IMAG <~ imag2      0.4473      0.1569    2.8500    0.0044 [ 0.1395; 0.7308 ] 
    ##   IMAG <~ imag3      0.6020      0.1453    4.1449    0.0000 [ 0.3026; 0.8703 ] 
    ##   EXPE <~ expe1      0.2946      0.1197    2.4619    0.0138 [ 0.0643; 0.5372 ] 
    ##   EXPE <~ expe2      0.6473      0.0868    7.4548    0.0000 [ 0.4508; 0.7990 ] 
    ##   EXPE <~ expe3      0.2374      0.0919    2.5840    0.0098 [ 0.0558; 0.4120 ] 
    ##   QUAL <~ qual1      0.2370      0.0889    2.6654    0.0077 [ 0.0737; 0.4154 ] 
    ##   QUAL <~ qual2      0.4712      0.0820    5.7487    0.0000 [ 0.2958; 0.6155 ] 
    ##   QUAL <~ qual3      0.1831      0.0793    2.3097    0.0209 [ 0.0199; 0.3394 ] 
    ##   QUAL <~ qual4      0.1037      0.0619    1.6769    0.0936 [-0.0155; 0.2267 ] 
    ##   QUAL <~ qual5      0.2049      0.0643    3.1871    0.0014 [ 0.0792; 0.3172 ] 
    ##   VAL <~ val1        0.7163      0.0959    7.4671    0.0000 [ 0.5211; 0.8888 ] 
    ##   VAL <~ val2        0.2202      0.0965    2.2827    0.0224 [ 0.0548; 0.4097 ] 
    ##   VAL <~ val3        0.2082      0.0587    3.5455    0.0004 [ 0.0860; 0.3341 ] 
    ##   SAT <~ sat1        0.3209      0.0145   22.0645    0.0000 [ 0.2969; 0.3509 ] 
    ##   SAT <~ sat2        0.3059      0.0138   22.1549    0.0000 [ 0.2814; 0.3368 ] 
    ##   SAT <~ sat3        0.2474      0.0111   22.2131    0.0000 [ 0.2228; 0.2683 ] 
    ##   SAT <~ sat4        0.2692      0.0123   21.9533    0.0000 [ 0.2455; 0.2937 ] 
    ##   LOY <~ loy1        0.3834      0.0266   14.4056    0.0000 [ 0.3318; 0.4381 ] 
    ##   LOY <~ loy2        0.2434      0.0319    7.6221    0.0000 [ 0.1713; 0.2945 ] 
    ##   LOY <~ loy3        0.3812      0.0261   14.5770    0.0000 [ 0.3306; 0.4350 ] 
    ##   LOY <~ loy4        0.2073      0.0369    5.6186    0.0000 [ 0.1308; 0.2781 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0651    9.8907    0.0000 [ 0.5132; 0.7608 ] 
    ##   imag1 ~~ imag3      0.5433      0.0718    7.5662    0.0000 [ 0.3904; 0.6742 ] 
    ##   imag2 ~~ imag3      0.7761      0.0392   19.8190    0.0000 [ 0.7007; 0.8504 ] 
    ##   expe1 ~~ expe2      0.5353      0.0637    8.4064    0.0000 [ 0.4077; 0.6462 ] 
    ##   expe1 ~~ expe3      0.4694      0.0635    7.3872    0.0000 [ 0.3269; 0.5901 ] 
    ##   expe2 ~~ expe3      0.5467      0.0605    9.0345    0.0000 [ 0.4094; 0.6488 ] 
    ##   qual1 ~~ qual2      0.6053      0.0611    9.9139    0.0000 [ 0.4838; 0.7076 ] 
    ##   qual1 ~~ qual3      0.5406      0.0629    8.5955    0.0000 [ 0.4096; 0.6573 ] 
    ##   qual1 ~~ qual4      0.5662      0.0734    7.7170    0.0000 [ 0.4078; 0.6851 ] 
    ##   qual1 ~~ qual5      0.5180      0.0688    7.5291    0.0000 [ 0.3785; 0.6382 ] 
    ##   qual2 ~~ qual3      0.6187      0.0545   11.3411    0.0000 [ 0.5052; 0.7121 ] 
    ##   qual2 ~~ qual4      0.6517      0.0630   10.3382    0.0000 [ 0.5029; 0.7597 ] 
    ##   qual2 ~~ qual5      0.6291      0.0580   10.8438    0.0000 [ 0.5065; 0.7333 ] 
    ##   qual3 ~~ qual4      0.4752      0.0666    7.1401    0.0000 [ 0.3364; 0.5848 ] 
    ##   qual3 ~~ qual5      0.5074      0.0626    8.1102    0.0000 [ 0.3815; 0.6260 ] 
    ##   qual4 ~~ qual5      0.6402      0.0574   11.1437    0.0000 [ 0.5038; 0.7405 ] 
    ##   val1 ~~ val2        0.6344      0.0555   11.4415    0.0000 [ 0.5252; 0.7414 ] 
    ##   val1 ~~ val3        0.4602      0.0717    6.4214    0.0000 [ 0.3065; 0.5945 ] 
    ##   val2 ~~ val3        0.6288      0.0637    9.8677    0.0000 [ 0.4982; 0.7455 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0684    6.8917    0.0000 [ 0.3426; 0.6091 ] 
    ##   QUAL ~ IMAG       0.3933      0.0637    6.1736    0.0000 [ 0.2763; 0.5232 ] 
    ##   QUAL ~ EXPE       0.8344      0.0242   34.4109    0.0000 [ 0.7811; 0.8787 ] 
    ##   VAL ~ IMAG        0.2974      0.0625    4.7615    0.0000 [ 0.1884; 0.4308 ] 
    ##   VAL ~ EXPE        0.6309      0.0508   12.4279    0.0000 [ 0.5341; 0.7270 ] 
    ##   VAL ~ QUAL        0.7013      0.0795    8.8224    0.0000 [ 0.5300; 0.8448 ] 
    ##   SAT ~ IMAG        0.4807      0.0707    6.7974    0.0000 [ 0.3368; 0.6253 ] 
    ##   SAT ~ EXPE        0.5001      0.0566    8.8317    0.0000 [ 0.3905; 0.6122 ] 
    ##   SAT ~ QUAL        0.5911      0.0951    6.2150    0.0000 [ 0.4015; 0.7764 ] 
    ##   SAT ~ VAL         0.5270      0.0870    6.0600    0.0000 [ 0.3577; 0.6931 ] 
    ##   LOY ~ IMAG        0.4840      0.0692    6.9894    0.0000 [ 0.3441; 0.6184 ] 
    ##   LOY ~ EXPE        0.3142      0.0539    5.8333    0.0000 [ 0.2190; 0.4359 ] 
    ##   LOY ~ QUAL        0.3714      0.0819    4.5352    0.0000 [ 0.2206; 0.5464 ] 
    ##   LOY ~ VAL         0.3311      0.0764    4.3327    0.0000 [ 0.1991; 0.4947 ] 
    ##   LOY ~ SAT         0.6283      0.0778    8.0712    0.0000 [ 0.4671; 0.7704 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0637    6.1736    0.0000 [ 0.2763; 0.5232 ] 
    ##   VAL ~ IMAG           0.2974      0.0625    4.7615    0.0000 [ 0.1884; 0.4308 ] 
    ##   VAL ~ EXPE           0.5852      0.0678    8.6369    0.0000 [ 0.4400; 0.7123 ] 
    ##   SAT ~ IMAG           0.2357      0.0491    4.8006    0.0000 [ 0.1485; 0.3410 ] 
    ##   SAT ~ EXPE           0.5173      0.0697    7.4271    0.0000 [ 0.3786; 0.6539 ] 
    ##   SAT ~ QUAL           0.3696      0.0622    5.9403    0.0000 [ 0.2503; 0.4919 ] 
    ##   LOY ~ IMAG           0.3020      0.0556    5.4325    0.0000 [ 0.2122; 0.4346 ] 
    ##   LOY ~ EXPE           0.3142      0.0539    5.8333    0.0000 [ 0.2190; 0.4359 ] 
    ##   LOY ~ QUAL           0.3714      0.0819    4.5352    0.0000 [ 0.2206; 0.5464 ] 
    ##   LOY ~ VAL            0.3311      0.0764    4.3327    0.0000 [ 0.1991; 0.4947 ] 
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
