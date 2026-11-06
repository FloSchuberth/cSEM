# Summarize model

**\[stable\]**

## Usage

``` r
summarize(
 .object = NULL, 
 .alpha  = 0.05,
 .ci     = NULL,
 ...
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .ci:

  A vector of character strings naming the confidence interval to
  compute. For possible choices see
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).

- ...:

  Further arguments to `summarize()`. Currently ignored.

## Value

An object of class `cSEMSummarize`. A `cSEMSummarize` object has the
same structure as the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object with a couple differences:

1.  Elements `$Path_estimates`, `$Loadings_estimates`,
    `$Weight_estimates`, `$Weight_estimates`, and
    `$Residual_correlation` are standardized data frames instead of
    matrices.

2.  Data frames `$Effect_estimates`, `$Indicator_correlation`, and
    `$Exo_construct_correlation` are added to `$Estimates`.

The data frame format is usually much more convenient if users intend to
present the results in e.g., a paper or a presentation.

## Details

The summary is mainly focused on estimated parameters. For quality
criteria such as the average variance extracted (AVE), reliability
estimates, effect size estimates etc., use
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).

If `.object` contains resamples, standard errors, t-values and p-values
(assuming estimates are standard normally distributed) are printed as
well. By default the percentile confidence interval is given as well.
For other confidence intervals use the `.ci` argument. See
[`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md) for
possible choices and a description.

## See also

[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`exportToExcel()`](https://floschuberth.github.io/cSEM/reference/exportToExcel.md)

## Examples

``` r
## Take a look at the dataset
#?threecommonfactors

## Specify the (correct) model
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Estimate
res <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 40)

## Postestimation
res_summarize <- summarize(res)
res_summarize
#> ________________________________________________________________________________
#> ----------------------------------- Overview -----------------------------------
#> 
#>  General information:
#>  ------------------------
#>  Estimation status                  = Ok
#>  Number of observations             = 500
#>  Weight estimator                   = PLS-PM
#>  Inner weighting scheme             = "path"
#>  Type of indicator correlation      = Pearson
#>  Path model estimator               = OLS
#>  Second-order approach              = NA
#>  Type of path model                 = Linear
#>  Disattenuated                      = Yes (PLSc)
#> 
#>  Resample information:
#>  ---------------------
#>  Resample method                    = "bootstrap"
#>  Number of resamples                = 40
#>  Number of admissible results       = 40
#>  Approach to handle inadmissibles   = "drop"
#>  Sign change option                 = "none"
#>  Random seed                        = -411350815
#> 
#>  Construct details:
#>  ------------------
#>  Name  Modeled as     Order         Mode      
#> 
#>  eta1  Common factor  First order   "modeA"   
#>  eta2  Common factor  First order   "modeA"   
#>  eta3  Common factor  First order   "modeA"   
#> 
#> ----------------------------------- Estimates ----------------------------------
#> 
#> Estimated path coefficients:
#> ============================
#>                                                              CI_percentile   
#>   Path           Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1      0.6713      0.0363   18.4965    0.0000 [ 0.6021; 0.7474 ] 
#>   eta3 ~ eta1      0.4585      0.0733    6.2574    0.0000 [ 0.3141; 0.5965 ] 
#>   eta3 ~ eta2      0.3052      0.0696    4.3836    0.0000 [ 0.1800; 0.4273 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0424   15.6311    0.0000 [ 0.5777; 0.7275 ] 
#>   eta1 =~ y12      0.6493      0.0451   14.3911    0.0000 [ 0.5585; 0.7075 ] 
#>   eta1 =~ y13      0.7613      0.0307   24.7825    0.0000 [ 0.6991; 0.8090 ] 
#>   eta2 =~ y21      0.5165      0.0546    9.4630    0.0000 [ 0.3827; 0.6344 ] 
#>   eta2 =~ y22      0.7554      0.0390   19.3926    0.0000 [ 0.6886; 0.8138 ] 
#>   eta2 =~ y23      0.7997      0.0348   22.9489    0.0000 [ 0.7357; 0.8628 ] 
#>   eta3 =~ y31      0.8223      0.0379   21.7015    0.0000 [ 0.7495; 0.8761 ] 
#>   eta3 =~ y32      0.6581      0.0352   18.7004    0.0000 [ 0.5920; 0.7129 ] 
#>   eta3 =~ y33      0.7474      0.0310   24.0868    0.0000 [ 0.7019; 0.8179 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0248   15.9777    0.0000 [ 0.3509; 0.4300 ] 
#>   eta1 <~ y12      0.3873      0.0206   18.8277    0.0000 [ 0.3398; 0.4229 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.8230    0.0000 [ 0.4221; 0.4931 ] 
#>   eta2 <~ y21      0.3058      0.0286   10.7005    0.0000 [ 0.2352; 0.3539 ] 
#>   eta2 <~ y22      0.4473      0.0208   21.5187    0.0000 [ 0.4198; 0.4843 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1702    0.0000 [ 0.4323; 0.5121 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5962    0.0000 [ 0.4080; 0.4658 ] 
#>   eta3 <~ y32      0.3521      0.0174   20.1905    0.0000 [ 0.3180; 0.3810 ] 
#>   eta3 <~ y33      0.3999      0.0169   23.7238    0.0000 [ 0.3750; 0.4370 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0363   18.4965    0.0000 [ 0.6021; 0.7474 ] 
#>   eta3 ~ eta1       0.6634      0.0373   17.7984    0.0000 [ 0.5903; 0.7185 ] 
#>   eta3 ~ eta2       0.3052      0.0696    4.3836    0.0000 [ 0.1800; 0.4273 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0471    4.3489    0.0000 [ 0.1218; 0.2867 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04242002 15.631059  4.473037e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04511678 14.391054  5.889180e-47
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03072113 24.782477 1.385215e-135
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05457595  9.463048  2.990963e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03895243 19.392568  8.917147e-84
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03484538 22.948918 1.510574e-116
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03789030 21.701526 1.984643e-104
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03519014 18.700379  4.915368e-78
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03103048 24.086774 3.439662e-128
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5777019          0.7274622
#> 2          0.5585091          0.7075393
#> 3          0.6990908          0.8090061
#> 4          0.3826657          0.6343591
#> 5          0.6886303          0.8138418
#> 6          0.7356744          0.8628401
#> 7          0.7495169          0.8761028
#> 8          0.5919687          0.7128537
#> 9          0.7019498          0.8179341

## By default only the 95% percentile confidence interval is printed. User
## can have several confidence interval computed, however, only the first
## will be printed.

res_summarize <- summarize(res, .ci = c("CI_standard_t", "CI_percentile"), 
                           .alpha = c(0.05, 0.01))
res_summarize
#> ________________________________________________________________________________
#> ----------------------------------- Overview -----------------------------------
#> 
#>  General information:
#>  ------------------------
#>  Estimation status                  = Ok
#>  Number of observations             = 500
#>  Weight estimator                   = PLS-PM
#>  Inner weighting scheme             = "path"
#>  Type of indicator correlation      = Pearson
#>  Path model estimator               = OLS
#>  Second-order approach              = NA
#>  Type of path model                 = Linear
#>  Disattenuated                      = Yes (PLSc)
#> 
#>  Resample information:
#>  ---------------------
#>  Resample method                    = "bootstrap"
#>  Number of resamples                = 40
#>  Number of admissible results       = 40
#>  Approach to handle inadmissibles   = "drop"
#>  Sign change option                 = "none"
#>  Random seed                        = -411350815
#> 
#>  Construct details:
#>  ------------------
#>  Name  Modeled as     Order         Mode      
#> 
#>  eta1  Common factor  First order   "modeA"   
#>  eta2  Common factor  First order   "modeA"   
#>  eta3  Common factor  First order   "modeA"   
#> 
#> ----------------------------------- Estimates ----------------------------------By default, only one confidence interval supplied to `.ci` is printed.
#> Use `xxx` to print all confidence intervals (not yet implemented).
#> 
#> 
#> 
#> Estimated path coefficients:
#> ============================
#>                                                              CI_standard_t   
#>   Path           Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1      0.6713      0.0363   18.4965    0.0000 [ 0.5758; 0.7635 ] 
#>   eta3 ~ eta1      0.4585      0.0733    6.2574    0.0000 [ 0.2672; 0.6462 ] 
#>   eta3 ~ eta2      0.3052      0.0696    4.3836    0.0000 [ 0.1349; 0.4949 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0424   15.6311    0.0000 [ 0.5629; 0.7822 ] 
#>   eta1 =~ y12      0.6493      0.0451   14.3911    0.0000 [ 0.5339; 0.7673 ] 
#>   eta1 =~ y13      0.7613      0.0307   24.7825    0.0000 [ 0.6856; 0.8445 ] 
#>   eta2 =~ y21      0.5165      0.0546    9.4630    0.0000 [ 0.3819; 0.6642 ] 
#>   eta2 =~ y22      0.7554      0.0390   19.3926    0.0000 [ 0.6574; 0.8588 ] 
#>   eta2 =~ y23      0.7997      0.0348   22.9489    0.0000 [ 0.7141; 0.8943 ] 
#>   eta3 =~ y31      0.8223      0.0379   21.7015    0.0000 [ 0.7332; 0.9291 ] 
#>   eta3 =~ y32      0.6581      0.0352   18.7004    0.0000 [ 0.5706; 0.7526 ] 
#>   eta3 =~ y33      0.7474      0.0310   24.0868    0.0000 [ 0.6543; 0.8148 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0248   15.9777    0.0000 [ 0.3333; 0.4614 ] 
#>   eta1 <~ y12      0.3873      0.0206   18.8277    0.0000 [ 0.3315; 0.4378 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.8230    0.0000 [ 0.3979; 0.5056 ] 
#>   eta2 <~ y21      0.3058      0.0286   10.7005    0.0000 [ 0.2336; 0.3814 ] 
#>   eta2 <~ y22      0.4473      0.0208   21.5187    0.0000 [ 0.3914; 0.4989 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1702    0.0000 [ 0.4113; 0.5327 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5962    0.0000 [ 0.4007; 0.4895 ] 
#>   eta3 <~ y32      0.3521      0.0174   20.1905    0.0000 [ 0.3092; 0.3994 ] 
#>   eta3 <~ y33      0.3999      0.0169   23.7238    0.0000 [ 0.3496; 0.4368 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0363   18.4965    0.0000 [ 0.5758; 0.7635 ] 
#>   eta3 ~ eta1       0.6634      0.0373   17.7984    0.0000 [ 0.5714; 0.7641 ] 
#>   eta3 ~ eta2       0.3052      0.0696    4.3836    0.0000 [ 0.1349; 0.4949 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0471    4.3489    0.0000 [ 0.0892; 0.3328 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03629525 18.496454 2.205119e-76
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07327485  6.257355 3.915607e-10
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.06961131  4.383643 1.167110e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5757833          0.7634819          0.5983223          0.7409429
#> 2          0.2672230          0.6461591          0.3127258          0.6006563
#> 3          0.1348619          0.4948523          0.1780897          0.4516245
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5989118          0.7493150          0.6021363          0.7473643
#> 2          0.2921109          0.5968467          0.3140516          0.5964515
#> 3          0.1764926          0.4419382          0.1800462          0.4273031
```
