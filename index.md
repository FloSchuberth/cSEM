# cSEM: Composite-based SEM ![](reference/figures/cSEMsticker.svg)

[![CRAN
status](https://www.r-pkg.org/badges/version/cSEM)](https://cran.r-project.org/package=cSEM)
[![R-CMD-check](https://github.com/FloSchuberth/cSEM/workflows/R-CMD-check/badge.svg)](https://github.com/FloSchuberth/cSEM/actions)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/M-E-Rademaker/cSEM?branch=master&svg=true)](https://ci.appveyor.com/project/M-E-Rademaker/csem)
![Lifecycle
Status](https://img.shields.io/badge/lifecycle-maturing-blue.svg)[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/cSEM)](https://cran.r-project.org/package=cSEM)

## Purpose

Estimate, analyse, test, and study linear, nonlinear, hierarchical and
multi-group structural equation models using composite-based approaches
and procedures, including estimation techniques such as partial least
squares path modeling (PLS-PM) and its derivatives (PLSc, OrdPLSc,
robustPLSc), generalized structured component analysis (GSCA),
generalized structured component analysis with uniqueness terms (GSCAm),
generalized canonical correlation analysis (GCCA), principal component
analysis (PCA), factor score regression (FSR) using sum score,
regression or Bartlett scores (including bias correction using Croon’s
approach), as well as several tests and typical post-estimation
procedures (e.g., verify admissibility of the estimates, assess the
model fit, test the model fit, compute confidence intervals, compare
groups, etc.).

## News (2025-05-15):

- Implementation of doModelSearch to perform AGAS-PLS. Thanks to Gloria.

- Release of cSEM version 0.6.1

- Release of cSEM Version 0.6.0

- Implementation of a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) function to
  visualize cSEM models. Thanks to Nguyen.

- Enhancement of the
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  function

## Installation

The package is available on [CRAN](https://cran.r-project.org/):

``` r
install.packages("cSEM")
```

To install the development version, which is recommended, use:

``` r
# install.packages("pak")
pak::pak("FloSchuberth/cSEM")
```

## Getting started

The best place to get started is the
[cSEM-website](https://floschuberth.github.io/cSEM/).

## Basic usage

The basic usage is illustrated below.

![](reference/figures/api.png)

Usually, using `cSEM` is the same 3 step procedure:

> 1.  Pick a dataset and specify a model using [lavaan
>     syntax](https://lavaan.ugent.be/tutorial/syntax1.html)
> 2.  Use
>     [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
> 3.  Apply one of the post-estimation functions listed below on the
>     resulting object.

## Post-Estimation Functions

There are five major post-estimation verbs, three test family functions
and three do-family of function:

- [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  : assess the model using common quality criteria
- [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md) :
  calculate common inferential quantities (e.g., standard errors,
  confidence intervals)
- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  : predict endogenous indicator values
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) : Plot the
  cSEM model
- [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  : summarize the results
- [`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md)
  : verify admissibility of the estimates

Tests are performed by using the test family of functions. Currently,
the following tests are implemented:

- [`testCVPAT()`](https://floschuberth.github.io/cSEM/reference/testCVPAT.md)
  performs a cross-validated predictive ability test
- [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
  : performs a test for overall model fit
- [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md)
  : performs a test for composite measurement invariance
- [`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)
  : performs several tests to assess multi-group differences
- [`testHausman()`](https://floschuberth.github.io/cSEM/reference/testHausman.md)
  : performs the regression-based Hausman test to test for endogeneity

Other miscellaneous post-estimation functions belong do the do-family of
functions. Currently, three do functions are implemented:

- [`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md):
  performs an importance-performance matrix analysis
- [`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md):
  performs a nonlinear effects analysis such as floodlight and surface
  analysis
- [`doRedundancyAnalysis()`](https://floschuberth.github.io/cSEM/reference/doRedundancyAnalysis.md):
  performs a redundancy analysis

All functions require a `cSEMResults` object.

## Example

Models are defined using [lavaan
syntax](https://lavaan.ugent.be/tutorial/syntax1.html) with some slight
modifications (see the [Specifying a
model](https://floschuberth.github.io/cSEM/articles/cSEM.html#using-csem)
section on the [cSEM-website](https://floschuberth.github.io/cSEM/)).
For illustration we use the build-in and well-known `satisfaction`
dataset.

``` r
require(cSEM)
    
## Note: The operator "<~" tells cSEM that the construct to its left is modeled
##       as a composite.
##       The operator "=~" tells cSEM that the construct to its left is modeled
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

The estimation is conducted using the
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
function.

``` r
# Estimate using defaults
res <- csem(.data = satisfaction, .model = model)
res
```

``` R
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
```

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
package](https://github.com/timelyportfolio/listviewer/). If you are new
to `cSEM` this might be a good way to familiarize yourself with the
structure of a `cSEMResults` object.

``` r
listviewer::jsonedit(res, mode = "view") # requires the listviewer package.
```

Apply post-estimation functions:

``` r
## Get a summary
summarize(res) 
```

``` R
## ________________________________________________________________________________
## ----------------------------------- Overview -----------------------------------
## 
##  General information:
##  ------------------------
##  Estimation status                  = Ok
##  Number of observations             = 250
##  Weight estimator                   = PLS-PM
##  Inner weighting scheme             = "path"
##  Type of indicator correlation      = Pearson
##  Path model estimator               = OLS
##  Second-order approach              = NA
##  Type of path model                 = Linear
##  Disattenuated                      = Yes (PLSc)
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
```

``` r
## Verify admissibility of the results
verify(res) 
```

``` R
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
```

``` r
## Test overall model fit
testOMF(res)
```

``` R
## ________________________________________________________________________________
## --------- Test for overall model fit based on Beran & Srivastava (1985) --------
## 
## Null hypothesis:
## 
##        ┌──────────────────────────────────────────────────────────────────┐
##        │                                                                  │
##        │   H0: The model-implied indicator covariance matrix equals the   │
##        │   population indicator covariance matrix.                        │
##        │                                                                  │
##        └──────────────────────────────────────────────────────────────────┘
## 
## Test statistic and critical value: 
## 
##                                      Critical value
##  Distance measure    Test statistic    95%   
##  dG                      0.6493      0.3190  
##  SRMR                    0.0940      0.0517  
##  dL                      2.2340      0.6753  
##  dML                     2.9219      1.5478  
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
##  Out of 499 bootstrap replications 482 are admissible.
##  See ?verify() for what constitutes an inadmissible result.
## 
##  The seed used was: -1769142168
## ________________________________________________________________________________
```

``` r
## Assess the model
assess(res)
```

``` R
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
##  Geodesic distance             = 0.6493432
##  Squared Euclidean distance    = 2.23402
##  ML distance                   = 2.921932
## 
##  Chi_square       = 727.5611
##  Chi_square_df    = 3.954137
##  CFI              = 0.8598825
##  CN               = 75.14588
##  GFI              = 0.7280612
##  IFI              = 0.8615598
##  NFI              = 0.8229918
##  NNFI             = 0.8240917
##  RMSEA            = 0.108922
##  RMS_theta        = 0.05069299
##  SRMR             = 0.09396871
## 
##  Degrees of freedom       = 184
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
## 
##   Dependent construct: 'SAT'
## 
##  Independent construct    VIF value 
##  EXPE                      3.2985   
##  QUAL                      4.4151   
##  IMAG                      1.7280   
##  VAL                       2.6726   
## 
##   Dependent construct: 'LOY'
## 
##  Independent construct    VIF value 
##  IMAG                      1.9345   
##  SAT                       1.9345   
## 
## -------------- Variance inflation factors (VIFs) for modeB weights -------------
## 
##   Construct: 'IMAG'
## 
##  Weight    VIF value 
##  imag1      1.7215   
##  imag2      3.0515   
##  imag3      2.5356   
## 
##   Construct: 'EXPE'
## 
##  Weight    VIF value 
##  expe1      1.4949   
##  expe2      1.6623   
##  expe3      1.5212   
## 
##   Construct: 'QUAL'
## 
##  Weight    VIF value 
##  qual1      1.8401   
##  qual2      2.5005   
##  qual3      1.7796   
##  qual4      2.1557   
##  qual5      2.0206   
## 
##   Construct: 'VAL'
## 
##  Weight    VIF value 
##  val1       1.6912   
##  val2       2.2049   
##  val3       1.6714   
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
## ----------------------- Discriminant validity assessment -----------------------
## 
##  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
## 
##           SAT LOY
## SAT 1.0000000   0
## LOY 0.7432489   1
## 
## 
##  Advanced heterotrait-monotrait ratio of correlations matrix (HTMT2 matrix)
## 
##           SAT LOY
## SAT 1.0000000   0
## LOY 0.7140046   1
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
```

``` r
## Predict indicator scores of endogenous constructs
predict(res)
```

``` R
## ________________________________________________________________________________
## ----------------------------------- Overview -----------------------------------
## 
##  Number of obs. training            = 225
##  Number of obs. test                = 25
##  Number of cv folds                 = 10
##  Number of repetitions              = 1
##  Handle inadmissibles               = stop
##  Estimator target                   = 'PLS-PM'
##  Estimator benchmark                = 'lm'
##  Disattenuation target              = 'TRUE'
##  Disattenuation benchmark           = 'FALSE'
##  Approach to predict                = 'earliest'
## 
## ------------------------------ Prediction metrics ------------------------------
## 
## 
##   Name      MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
##   expe1         1.4485         1.5897       1.9030         2.0890       0.0564
##   expe2         1.4208         1.5024       1.9411         2.0314       0.1968
##   expe3         1.6302         1.7376       2.1286         2.2232       0.1222
##   qual1         1.4783         1.5520       1.9303         2.0396       0.1190
##   qual2         1.5790         1.5584       2.0449         2.0755       0.2138
##   qual3         1.7358         1.7511       2.2289         2.3008       0.1151
##   qual4         1.2363         1.1892       1.5985         1.6291       0.2282
##   qual5         1.5123         1.5188       1.9450         1.9635       0.1884
##   val1          1.4498         1.3687       1.8750         1.7690       0.2433
##   val2          1.2274         1.2204       1.6520         1.7076       0.1715
##   val3          1.4868         1.3927       1.9758         1.9371       0.1419
##   sat1          1.2551         1.2406       1.6542         1.6306       0.3365
##   sat2          1.2398         1.2131       1.6481         1.6373       0.3047
##   sat3          1.3384         1.2954       1.6739         1.7285       0.2096
##   sat4          1.3220         1.2675       1.6717         1.6393       0.2752
##   loy1          1.6955         1.6836       2.2394         2.2465       0.2703
##   loy2          1.4821         1.4905       1.9099         1.9848       0.1307
##   loy3          1.6993         1.6661       2.2771         2.2617       0.2725
##   loy4          1.6925         1.6744       2.1834         2.2952       0.0822
## ________________________________________________________________________________
```

#### Resampling and Inference

By default no inferential statistics are calculated since most
composite-based estimators have no closed-form expressions for standard
errors. Resampling is used instead. `cSEM` mostly relies on the
`bootstrap` procedure (although `jackknife` is implemented as well) to
estimate standard errors, test statistics, and critical quantiles.

`cSEM` offers two ways for resampling:

1.  Setting `.resample_method` in
    [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) to
    `"jackknife"` or `"bootstrap"` and subsequently using
    post-estimation functions
    [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
    or
    [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).
2.  The same result is achieved by passing a `cSEMResults` object to
    [`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
    and subsequently using post-estimation functions
    [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
    or
    [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).

``` r
# Setting `.resample_method`
b1 <- csem(.data = satisfaction, .model = model, .resample_method = "bootstrap")
# Using resamplecSEMResults()
b2 <- resamplecSEMResults(res)
```

The
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
function reports the inferential statistics:

``` r
summarize(b1)
```

``` R
## ________________________________________________________________________________
## ----------------------------------- Overview -----------------------------------
## 
##  General information:
##  ------------------------
##  Estimation status                  = Ok
##  Number of observations             = 250
##  Weight estimator                   = PLS-PM
##  Inner weighting scheme             = "path"
##  Type of indicator correlation      = Pearson
##  Path model estimator               = OLS
##  Second-order approach              = NA
##  Type of path model                 = Linear
##  Disattenuated                      = Yes (PLSc)
## 
##  Resample information:
##  ---------------------
##  Resample method                    = "bootstrap"
##  Number of resamples                = 499
##  Number of admissible results       = 487
##  Approach to handle inadmissibles   = "drop"
##  Sign change option                 = "none"
##  Random seed                        = -941137440
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
##   EXPE ~ IMAG      0.4714      0.0655    7.1980    0.0000 [ 0.3526; 0.6091 ] 
##   QUAL ~ EXPE      0.8344      0.0230   36.2822    0.0000 [ 0.7855; 0.8742 ] 
##   VAL ~ EXPE       0.0457      0.0893    0.5122    0.6085 [-0.1250; 0.2236 ] 
##   VAL ~ QUAL       0.7013      0.0851    8.2438    0.0000 [ 0.5419; 0.8596 ] 
##   SAT ~ IMAG       0.2450      0.0525    4.6645    0.0000 [ 0.1392; 0.3465 ] 
##   SAT ~ EXPE      -0.0172      0.0697   -0.2472    0.8048 [-0.1504; 0.1109 ] 
##   SAT ~ QUAL       0.2215      0.1079    2.0528    0.0401 [ 0.0420; 0.4618 ] 
##   SAT ~ VAL        0.5270      0.0862    6.1102    0.0000 [ 0.3471; 0.6923 ] 
##   LOY ~ IMAG       0.1819      0.0781    2.3307    0.0198 [ 0.0298; 0.3334 ] 
##   LOY ~ SAT        0.6283      0.0784    8.0177    0.0000 [ 0.4729; 0.7861 ] 
## 
## Estimated loadings:
## ===================
##                                                                CI_percentile   
##   Loading          Estimate  Std. error   t-stat.   p-value         95%        
##   IMAG =~ imag1      0.6306      0.0991    6.3651    0.0000 [ 0.4253; 0.8012 ] 
##   IMAG =~ imag2      0.9246      0.0412   22.4219    0.0000 [ 0.8155; 0.9789 ] 
##   IMAG =~ imag3      0.9577      0.0299   32.0520    0.0000 [ 0.8794; 0.9927 ] 
##   EXPE =~ expe1      0.7525      0.0788    9.5467    0.0000 [ 0.5586; 0.8705 ] 
##   EXPE =~ expe2      0.9348      0.0259   36.1542    0.0000 [ 0.8669; 0.9687 ] 
##   EXPE =~ expe3      0.7295      0.0743    9.8126    0.0000 [ 0.5384; 0.8351 ] 
##   QUAL =~ qual1      0.7861      0.0712   11.0371    0.0000 [ 0.6220; 0.8848 ] 
##   QUAL =~ qual2      0.9244      0.0219   42.2277    0.0000 [ 0.8703; 0.9559 ] 
##   QUAL =~ qual3      0.7560      0.0617   12.2518    0.0000 [ 0.6049; 0.8551 ] 
##   QUAL =~ qual4      0.7632      0.0541   14.0989    0.0000 [ 0.6501; 0.8570 ] 
##   QUAL =~ qual5      0.7834      0.0452   17.3154    0.0000 [ 0.6789; 0.8610 ] 
##   VAL =~ val1        0.9518      0.0236   40.3971    0.0000 [ 0.8925; 0.9863 ] 
##   VAL =~ val2        0.8056      0.0666   12.0888    0.0000 [ 0.6632; 0.9027 ] 
##   VAL =~ val3        0.6763      0.0747    9.0483    0.0000 [ 0.5199; 0.8020 ] 
##   SAT =~ sat1        0.9243      0.0223   41.3799    0.0000 [ 0.8761; 0.9654 ] 
##   SAT =~ sat2        0.8813      0.0282   31.2051    0.0000 [ 0.8206; 0.9252 ] 
##   SAT =~ sat3        0.7127      0.0539   13.2346    0.0000 [ 0.5914; 0.7962 ] 
##   SAT =~ sat4        0.7756      0.0492   15.7694    0.0000 [ 0.6664; 0.8589 ] 
##   LOY =~ loy1        0.9097      0.0512   17.7690    0.0000 [ 0.7983; 0.9908 ] 
##   LOY =~ loy2        0.5775      0.0852    6.7764    0.0000 [ 0.3946; 0.7282 ] 
##   LOY =~ loy3        0.9043      0.0417   21.6654    0.0000 [ 0.8097; 0.9720 ] 
##   LOY =~ loy4        0.4917      0.0974    5.0466    0.0000 [ 0.2889; 0.6607 ] 
## 
## Estimated weights:
## ==================
##                                                                CI_percentile   
##   Weight           Estimate  Std. error   t-stat.   p-value         95%        
##   IMAG <~ imag1      0.0156      0.1188    0.1316    0.8953 [-0.2130; 0.2534 ] 
##   IMAG <~ imag2      0.4473      0.1566    2.8562    0.0043 [ 0.1393; 0.7297 ] 
##   IMAG <~ imag3      0.6020      0.1463    4.1165    0.0000 [ 0.3011; 0.8501 ] 
##   EXPE <~ expe1      0.2946      0.1128    2.6122    0.0090 [ 0.0706; 0.5200 ] 
##   EXPE <~ expe2      0.6473      0.0773    8.3761    0.0000 [ 0.4805; 0.7806 ] 
##   EXPE <~ expe3      0.2374      0.0940    2.5260    0.0115 [ 0.0629; 0.4260 ] 
##   QUAL <~ qual1      0.2370      0.0906    2.6151    0.0089 [ 0.0720; 0.4268 ] 
##   QUAL <~ qual2      0.4712      0.0777    6.0633    0.0000 [ 0.3096; 0.6206 ] 
##   QUAL <~ qual3      0.1831      0.0838    2.1850    0.0289 [ 0.0222; 0.3407 ] 
##   QUAL <~ qual4      0.1037      0.0596    1.7411    0.0817 [-0.0020; 0.2367 ] 
##   QUAL <~ qual5      0.2049      0.0616    3.3255    0.0009 [ 0.0790; 0.3252 ] 
##   VAL <~ val1        0.7163      0.0976    7.3384    0.0000 [ 0.5129; 0.8851 ] 
##   VAL <~ val2        0.2202      0.0953    2.3115    0.0208 [ 0.0199; 0.4028 ] 
##   VAL <~ val3        0.2082      0.0593    3.5073    0.0005 [ 0.0870; 0.3302 ] 
##   SAT <~ sat1        0.3209      0.0156   20.6335    0.0000 [ 0.2955; 0.3591 ] 
##   SAT <~ sat2        0.3059      0.0134   22.8684    0.0000 [ 0.2822; 0.3372 ] 
##   SAT <~ sat3        0.2474      0.0112   22.0782    0.0000 [ 0.2213; 0.2647 ] 
##   SAT <~ sat4        0.2692      0.0123   21.9171    0.0000 [ 0.2469; 0.2941 ] 
##   LOY <~ loy1        0.3834      0.0273   14.0531    0.0000 [ 0.3336; 0.4426 ] 
##   LOY <~ loy2        0.2434      0.0300    8.1061    0.0000 [ 0.1776; 0.2978 ] 
##   LOY <~ loy3        0.3812      0.0280   13.6219    0.0000 [ 0.3305; 0.4420 ] 
##   LOY <~ loy4        0.2073      0.0358    5.7862    0.0000 [ 0.1330; 0.2689 ] 
## 
## Estimated indicator correlations:
## =================================
##                                                                 CI_percentile   
##   Correlation       Estimate  Std. error   t-stat.   p-value         95%        
##   imag1 ~~ imag2      0.6437      0.0666    9.6703    0.0000 [ 0.4992; 0.7478 ] 
##   imag1 ~~ imag3      0.5433      0.0698    7.7857    0.0000 [ 0.3990; 0.6717 ] 
##   imag2 ~~ imag3      0.7761      0.0394   19.7008    0.0000 [ 0.6922; 0.8508 ] 
##   expe1 ~~ expe2      0.5353      0.0619    8.6440    0.0000 [ 0.4096; 0.6395 ] 
##   expe1 ~~ expe3      0.4694      0.0657    7.1475    0.0000 [ 0.3439; 0.5874 ] 
##   expe2 ~~ expe3      0.5467      0.0613    8.9152    0.0000 [ 0.4150; 0.6457 ] 
##   qual1 ~~ qual2      0.6053      0.0600   10.0958    0.0000 [ 0.4693; 0.7032 ] 
##   qual1 ~~ qual3      0.5406      0.0608    8.8987    0.0000 [ 0.4098; 0.6533 ] 
##   qual1 ~~ qual4      0.5662      0.0705    8.0281    0.0000 [ 0.4080; 0.6785 ] 
##   qual1 ~~ qual5      0.5180      0.0699    7.4159    0.0000 [ 0.3659; 0.6363 ] 
##   qual2 ~~ qual3      0.6187      0.0549   11.2661    0.0000 [ 0.5031; 0.7168 ] 
##   qual2 ~~ qual4      0.6517      0.0613   10.6265    0.0000 [ 0.5220; 0.7632 ] 
##   qual2 ~~ qual5      0.6291      0.0566   11.1133    0.0000 [ 0.5113; 0.7320 ] 
##   qual3 ~~ qual4      0.4752      0.0649    7.3252    0.0000 [ 0.3400; 0.5879 ] 
##   qual3 ~~ qual5      0.5074      0.0603    8.4207    0.0000 [ 0.3926; 0.6214 ] 
##   qual4 ~~ qual5      0.6402      0.0566   11.3069    0.0000 [ 0.5245; 0.7411 ] 
##   val1 ~~ val2        0.6344      0.0543   11.6784    0.0000 [ 0.5181; 0.7354 ] 
##   val1 ~~ val3        0.4602      0.0704    6.5413    0.0000 [ 0.3209; 0.5976 ] 
##   val2 ~~ val3        0.6288      0.0656    9.5867    0.0000 [ 0.4918; 0.7545 ] 
## 
## ------------------------------------ Effects -----------------------------------
## 
## Estimated total effects:
## ========================
##                                                               CI_percentile   
##   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
##   EXPE ~ IMAG       0.4714      0.0655    7.1980    0.0000 [ 0.3526; 0.6091 ] 
##   QUAL ~ IMAG       0.3933      0.0603    6.5213    0.0000 [ 0.2782; 0.5175 ] 
##   QUAL ~ EXPE       0.8344      0.0230   36.2822    0.0000 [ 0.7855; 0.8742 ] 
##   VAL ~ IMAG        0.2974      0.0611    4.8690    0.0000 [ 0.1849; 0.4221 ] 
##   VAL ~ EXPE        0.6309      0.0513   12.3092    0.0000 [ 0.5301; 0.7258 ] 
##   VAL ~ QUAL        0.7013      0.0851    8.2438    0.0000 [ 0.5419; 0.8596 ] 
##   SAT ~ IMAG        0.4807      0.0659    7.2963    0.0000 [ 0.3446; 0.6086 ] 
##   SAT ~ EXPE        0.5001      0.0561    8.9180    0.0000 [ 0.3930; 0.6072 ] 
##   SAT ~ QUAL        0.5911      0.0969    6.0994    0.0000 [ 0.4162; 0.7705 ] 
##   SAT ~ VAL         0.5270      0.0862    6.1102    0.0000 [ 0.3471; 0.6923 ] 
##   LOY ~ IMAG        0.4840      0.0674    7.1819    0.0000 [ 0.3530; 0.6203 ] 
##   LOY ~ EXPE        0.3142      0.0532    5.9090    0.0000 [ 0.2201; 0.4234 ] 
##   LOY ~ QUAL        0.3714      0.0843    4.4075    0.0000 [ 0.2301; 0.5551 ] 
##   LOY ~ VAL         0.3311      0.0704    4.7038    0.0000 [ 0.2012; 0.4752 ] 
##   LOY ~ SAT         0.6283      0.0784    8.0177    0.0000 [ 0.4729; 0.7861 ] 
## 
## Estimated indirect effects:
## ===========================
##                                                                  CI_percentile   
##   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
##   QUAL ~ IMAG          0.3933      0.0603    6.5213    0.0000 [ 0.2782; 0.5175 ] 
##   VAL ~ IMAG           0.2974      0.0611    4.8690    0.0000 [ 0.1849; 0.4221 ] 
##   VAL ~ EXPE           0.5852      0.0722    8.1027    0.0000 [ 0.4479; 0.7153 ] 
##   SAT ~ IMAG           0.2357      0.0487    4.8393    0.0000 [ 0.1465; 0.3330 ] 
##   SAT ~ EXPE           0.5173      0.0690    7.4969    0.0000 [ 0.4031; 0.6489 ] 
##   SAT ~ QUAL           0.3696      0.0637    5.8009    0.0000 [ 0.2435; 0.4912 ] 
##   LOY ~ IMAG           0.3020      0.0558    5.4141    0.0000 [ 0.2021; 0.4226 ] 
##   LOY ~ EXPE           0.3142      0.0532    5.9090    0.0000 [ 0.2201; 0.4234 ] 
##   LOY ~ QUAL           0.3714      0.0843    4.4075    0.0000 [ 0.2301; 0.5551 ] 
##   LOY ~ VAL            0.3311      0.0704    4.7038    0.0000 [ 0.2012; 0.4752 ] 
## ________________________________________________________________________________
```

Several bootstrap-based confidence intervals are implemented, see
`?infer()`:

``` r
infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet
```

Both bootstrap and jackknife resampling support platform-independent
multiprocessing as well as setting random seeds via the [future
framework](https://github.com/futureverse/future/). For multiprocessing
simply set `.eval_plan = "multisession"` in which case the maximum
number of available cores is used if not on Windows. On Windows as many
separate R instances are opened in the background as there are cores
available instead. Note that this naturally has some overhead so for a
small number of resamples multiprocessing will not always be faster
compared to sequential (single core) processing (the default). Seeds are
set via the `.seed` argument.

``` r
b <- csem(
  .data            = satisfaction,
  .model           = model, 
  .resample_method = "bootstrap",
  .R               = 999,
  .seed            = 98234,
  .eval_plan       = "multisession")
```
