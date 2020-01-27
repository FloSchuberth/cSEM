
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
    ##                                                 +------------------------------------------------------------+
    ##                                                 |                                                            |
    ##                                                 |   H0: Population indicator covariance matrix is equal to   |
    ##                                                 |   model-implied indicator covariance matrix.               |
    ##                                                 |                                                            |
    ##                                                 +------------------------------------------------------------+
    ## 
    ## Test statistic and critical value: 
    ## 
    ##                                      Critical value
    ##  Distance measure    Test statistic    95%   
    ##  dG                      0.6493      0.3183  
    ##  SRMR                    0.0940      0.0510  
    ##  dL                      2.2340      0.6580  
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
    ##  Out of 499 bootstrap replications 466 are admissible.
    ##  See ?verify() for what constitutes an inadmissible result.
    ## 
    ##  The seed used was: -896903287
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
    ##   expe1         1.4537         1.5644       1.9062         2.0935       0.0556
    ##   expe2         1.4136         1.4797       1.9359         2.0300       0.2002
    ##   expe3         1.6323         1.7263       2.1266         2.2216       0.1246
    ##   qual1         1.4744         1.5444       1.9288         2.0642       0.1153
    ##   qual2         1.5805         1.5375       2.0414         2.0606       0.2178
    ##   qual3         1.7326         1.7254       2.2234         2.2791       0.1192
    ##   qual4         1.2336         1.1957       1.5969         1.6277       0.2334
    ##   qual5         1.5075         1.5027       1.9356         1.9555       0.1981
    ##   val1          1.4485         1.3628       1.8727         1.7653       0.2488
    ##   val2          1.2290         1.2064       1.6516         1.7138       0.1719
    ##   val3          1.4793         1.3782       1.9681         1.9315       0.1486
    ##   sat1          1.2466         1.2319       1.6474         1.6191       0.3396
    ##   sat2          1.2337         1.1939       1.6430         1.6249       0.3078
    ##   sat3          1.3417         1.2751       1.6757         1.7204       0.2085
    ##   sat4          1.3186         1.2607       1.6706         1.6342       0.2756
    ##   loy1          1.6949         1.6582       2.2347         2.2231       0.2685
    ##   loy2          1.4851         1.4725       1.9102         1.9772       0.1335
    ##   loy3          1.7073         1.6696       2.2809         2.2691       0.2718
    ##   loy4          1.6890         1.6653       2.1766         2.2923       0.0879
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
    ##  Number of admissible results     = 483
    ##  Approach to handle inadmissibles = drop
    ##  Sign change option               = none
    ##  Random seed                      = 551525737
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
    ##   EXPE ~ IMAG      0.4714      0.0647    7.2864    0.0000 [ 0.3499; 0.5959 ] 
    ##   QUAL ~ EXPE      0.8344      0.0236   35.4025    0.0000 [ 0.7820; 0.8749 ] 
    ##   VAL ~ EXPE       0.0457      0.0882    0.5184    0.6042 [-0.1215; 0.2244 ] 
    ##   VAL ~ QUAL       0.7013      0.0866    8.1010    0.0000 [ 0.5238; 0.8788 ] 
    ##   SAT ~ IMAG       0.2450      0.0571    4.2909    0.0000 [ 0.1358; 0.3590 ] 
    ##   SAT ~ EXPE      -0.0172      0.0708   -0.2436    0.8076 [-0.1720; 0.1218 ] 
    ##   SAT ~ QUAL       0.2215      0.1084    2.0446    0.0409 [ 0.0207; 0.4494 ] 
    ##   SAT ~ VAL        0.5270      0.0904    5.8307    0.0000 [ 0.3523; 0.6989 ] 
    ##   LOY ~ IMAG       0.1819      0.0837    2.1742    0.0297 [ 0.0412; 0.3519 ] 
    ##   LOY ~ SAT        0.6283      0.0865    7.2650    0.0000 [ 0.4417; 0.7929 ] 
    ## 
    ## Estimated loadings:
    ## ===================
    ##                                                                CI_percentile   
    ##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG =~ imag1      0.6306      0.0984    6.4074    0.0000 [ 0.4341; 0.8148 ] 
    ##   IMAG =~ imag2      0.9246      0.0407   22.7049    0.0000 [ 0.8252; 0.9804 ] 
    ##   IMAG =~ imag3      0.9577      0.0298   32.1122    0.0000 [ 0.8818; 0.9937 ] 
    ##   EXPE =~ expe1      0.7525      0.0768    9.7932    0.0000 [ 0.5793; 0.8712 ] 
    ##   EXPE =~ expe2      0.9348      0.0280   33.4437    0.0000 [ 0.8635; 0.9727 ] 
    ##   EXPE =~ expe3      0.7295      0.0744    9.8088    0.0000 [ 0.5719; 0.8437 ] 
    ##   QUAL =~ qual1      0.7861      0.0664   11.8318    0.0000 [ 0.6358; 0.8901 ] 
    ##   QUAL =~ qual2      0.9244      0.0220   41.9460    0.0000 [ 0.8730; 0.9538 ] 
    ##   QUAL =~ qual3      0.7560      0.0637   11.8630    0.0000 [ 0.6043; 0.8527 ] 
    ##   QUAL =~ qual4      0.7632      0.0534   14.3019    0.0000 [ 0.6385; 0.8563 ] 
    ##   QUAL =~ qual5      0.7834      0.0451   17.3826    0.0000 [ 0.6817; 0.8603 ] 
    ##   VAL =~ val1        0.9518      0.0240   39.7358    0.0000 [ 0.8962; 0.9871 ] 
    ##   VAL =~ val2        0.8056      0.0703   11.4560    0.0000 [ 0.6485; 0.9164 ] 
    ##   VAL =~ val3        0.6763      0.0768    8.8047    0.0000 [ 0.5045; 0.8034 ] 
    ##   SAT =~ sat1        0.9243      0.0233   39.6095    0.0000 [ 0.8735; 0.9626 ] 
    ##   SAT =~ sat2        0.8813      0.0309   28.5575    0.0000 [ 0.8122; 0.9301 ] 
    ##   SAT =~ sat3        0.7127      0.0528   13.4981    0.0000 [ 0.6001; 0.8091 ] 
    ##   SAT =~ sat4        0.7756      0.0522   14.8627    0.0000 [ 0.6575; 0.8658 ] 
    ##   LOY =~ loy1        0.9097      0.0506   17.9929    0.0000 [ 0.7899; 0.9883 ] 
    ##   LOY =~ loy2        0.5775      0.0882    6.5443    0.0000 [ 0.3708; 0.7216 ] 
    ##   LOY =~ loy3        0.9043      0.0428   21.1455    0.0000 [ 0.8114; 0.9747 ] 
    ##   LOY =~ loy4        0.4917      0.0989    4.9746    0.0000 [ 0.2996; 0.6770 ] 
    ## 
    ## Estimated weights:
    ## ==================
    ##                                                                CI_percentile   
    ##   Weights          Estimate  Std. error   t-stat.   p-value         95%        
    ##   IMAG <~ imag1      0.0156      0.1115    0.1402    0.8885 [-0.1890; 0.2620 ] 
    ##   IMAG <~ imag2      0.4473      0.1544    2.8978    0.0038 [ 0.1559; 0.7227 ] 
    ##   IMAG <~ imag3      0.6020      0.1487    4.0476    0.0001 [ 0.2834; 0.8700 ] 
    ##   EXPE <~ expe1      0.2946      0.1159    2.5413    0.0110 [ 0.0735; 0.5269 ] 
    ##   EXPE <~ expe2      0.6473      0.0846    7.6514    0.0000 [ 0.4634; 0.7972 ] 
    ##   EXPE <~ expe3      0.2374      0.0937    2.5330    0.0113 [ 0.0491; 0.4012 ] 
    ##   QUAL <~ qual1      0.2370      0.0906    2.6174    0.0089 [ 0.0848; 0.4188 ] 
    ##   QUAL <~ qual2      0.4712      0.0741    6.3578    0.0000 [ 0.3333; 0.6128 ] 
    ##   QUAL <~ qual3      0.1831      0.0815    2.2472    0.0246 [ 0.0351; 0.3391 ] 
    ##   QUAL <~ qual4      0.1037      0.0562    1.8461    0.0649 [ 0.0053; 0.2260 ] 
    ##   QUAL <~ qual5      0.2049      0.0631    3.2470    0.0012 [ 0.0557; 0.3150 ] 
    ##   VAL <~ val1        0.7163      0.1012    7.0807    0.0000 [ 0.5054; 0.8953 ] 
    ##   VAL <~ val2        0.2202      0.0990    2.2246    0.0261 [ 0.0322; 0.4186 ] 
    ##   VAL <~ val3        0.2082      0.0591    3.5238    0.0004 [ 0.0801; 0.3167 ] 
    ##   SAT <~ sat1        0.3209      0.0161   19.9886    0.0000 [ 0.2931; 0.3561 ] 
    ##   SAT <~ sat2        0.3059      0.0139   21.9914    0.0000 [ 0.2824; 0.3349 ] 
    ##   SAT <~ sat3        0.2474      0.0113   21.9528    0.0000 [ 0.2263; 0.2697 ] 
    ##   SAT <~ sat4        0.2692      0.0124   21.6886    0.0000 [ 0.2460; 0.2931 ] 
    ##   LOY <~ loy1        0.3834      0.0269   14.2687    0.0000 [ 0.3333; 0.4377 ] 
    ##   LOY <~ loy2        0.2434      0.0310    7.8540    0.0000 [ 0.1705; 0.2913 ] 
    ##   LOY <~ loy3        0.3812      0.0282   13.5023    0.0000 [ 0.3275; 0.4379 ] 
    ##   LOY <~ loy4        0.2073      0.0368    5.6291    0.0000 [ 0.1385; 0.2745 ] 
    ## 
    ## Estimated indicator correlations:
    ## =================================
    ##                                                                 CI_percentile   
    ##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
    ##   imag1 ~~ imag2      0.6437      0.0680    9.4661    0.0000 [ 0.4899; 0.7569 ] 
    ##   imag1 ~~ imag3      0.5433      0.0659    8.2418    0.0000 [ 0.4183; 0.6765 ] 
    ##   imag2 ~~ imag3      0.7761      0.0374   20.7265    0.0000 [ 0.6987; 0.8408 ] 
    ##   expe1 ~~ expe2      0.5353      0.0586    9.1393    0.0000 [ 0.4116; 0.6378 ] 
    ##   expe1 ~~ expe3      0.4694      0.0623    7.5308    0.0000 [ 0.3330; 0.5806 ] 
    ##   expe2 ~~ expe3      0.5467      0.0609    8.9753    0.0000 [ 0.4087; 0.6550 ] 
    ##   qual1 ~~ qual2      0.6053      0.0561   10.7853    0.0000 [ 0.4891; 0.7017 ] 
    ##   qual1 ~~ qual3      0.5406      0.0578    9.3543    0.0000 [ 0.4074; 0.6462 ] 
    ##   qual1 ~~ qual4      0.5662      0.0664    8.5304    0.0000 [ 0.4237; 0.6753 ] 
    ##   qual1 ~~ qual5      0.5180      0.0660    7.8460    0.0000 [ 0.3673; 0.6343 ] 
    ##   qual2 ~~ qual3      0.6187      0.0566   10.9245    0.0000 [ 0.4981; 0.7163 ] 
    ##   qual2 ~~ qual4      0.6517      0.0610   10.6820    0.0000 [ 0.5110; 0.7603 ] 
    ##   qual2 ~~ qual5      0.6291      0.0570   11.0386    0.0000 [ 0.5093; 0.7353 ] 
    ##   qual3 ~~ qual4      0.4752      0.0632    7.5134    0.0000 [ 0.3431; 0.5883 ] 
    ##   qual3 ~~ qual5      0.5074      0.0603    8.4167    0.0000 [ 0.3932; 0.6125 ] 
    ##   qual4 ~~ qual5      0.6402      0.0551   11.6259    0.0000 [ 0.5196; 0.7297 ] 
    ##   val1 ~~ val2        0.6344      0.0569   11.1479    0.0000 [ 0.5165; 0.7448 ] 
    ##   val1 ~~ val3        0.4602      0.0722    6.3728    0.0000 [ 0.3255; 0.5935 ] 
    ##   val2 ~~ val3        0.6288      0.0637    9.8769    0.0000 [ 0.4890; 0.7460 ] 
    ## 
    ## ------------------------------------ Effects -----------------------------------
    ## 
    ## Estimated total effects:
    ## ========================
    ##                                                               CI_percentile   
    ##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   EXPE ~ IMAG       0.4714      0.0647    7.2864    0.0000 [ 0.3499; 0.5959 ] 
    ##   QUAL ~ IMAG       0.3933      0.0605    6.4973    0.0000 [ 0.2752; 0.5127 ] 
    ##   QUAL ~ EXPE       0.8344      0.0236   35.4025    0.0000 [ 0.7820; 0.8749 ] 
    ##   VAL ~ IMAG        0.2974      0.0597    4.9808    0.0000 [ 0.1885; 0.4171 ] 
    ##   VAL ~ EXPE        0.6309      0.0510   12.3804    0.0000 [ 0.5293; 0.7153 ] 
    ##   VAL ~ QUAL        0.7013      0.0866    8.1010    0.0000 [ 0.5238; 0.8788 ] 
    ##   SAT ~ IMAG        0.4807      0.0689    6.9762    0.0000 [ 0.3504; 0.6090 ] 
    ##   SAT ~ EXPE        0.5001      0.0582    8.5937    0.0000 [ 0.3814; 0.6061 ] 
    ##   SAT ~ QUAL        0.5911      0.0999    5.9191    0.0000 [ 0.4076; 0.7906 ] 
    ##   SAT ~ VAL         0.5270      0.0904    5.8307    0.0000 [ 0.3523; 0.6989 ] 
    ##   LOY ~ IMAG        0.4840      0.0675    7.1719    0.0000 [ 0.3512; 0.6224 ] 
    ##   LOY ~ EXPE        0.3142      0.0546    5.7592    0.0000 [ 0.2055; 0.4222 ] 
    ##   LOY ~ QUAL        0.3714      0.0876    4.2383    0.0000 [ 0.2181; 0.5695 ] 
    ##   LOY ~ VAL         0.3311      0.0778    4.2581    0.0000 [ 0.1874; 0.4755 ] 
    ##   LOY ~ SAT         0.6283      0.0865    7.2650    0.0000 [ 0.4417; 0.7929 ] 
    ## 
    ## Estimated indirect effects:
    ## ===========================
    ##                                                                  CI_percentile   
    ##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
    ##   QUAL ~ IMAG          0.3933      0.0605    6.4973    0.0000 [ 0.2752; 0.5127 ] 
    ##   VAL ~ IMAG           0.2974      0.0597    4.9808    0.0000 [ 0.1885; 0.4171 ] 
    ##   VAL ~ EXPE           0.5852      0.0748    7.8236    0.0000 [ 0.4308; 0.7232 ] 
    ##   SAT ~ IMAG           0.2357      0.0469    5.0243    0.0000 [ 0.1472; 0.3263 ] 
    ##   SAT ~ EXPE           0.5173      0.0702    7.3687    0.0000 [ 0.3815; 0.6458 ] 
    ##   SAT ~ QUAL           0.3696      0.0614    6.0187    0.0000 [ 0.2475; 0.4840 ] 
    ##   LOY ~ IMAG           0.3020      0.0595    5.0765    0.0000 [ 0.1987; 0.4232 ] 
    ##   LOY ~ EXPE           0.3142      0.0546    5.7592    0.0000 [ 0.2055; 0.4222 ] 
    ##   LOY ~ QUAL           0.3714      0.0876    4.2383    0.0000 [ 0.2181; 0.5695 ] 
    ##   LOY ~ VAL            0.3311      0.0778    4.2581    0.0000 [ 0.1874; 0.4755 ] 
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
