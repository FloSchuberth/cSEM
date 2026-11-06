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
#>  Random seed                        = 992949468
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
#>   eta2 ~ eta1      0.6713      0.0380   17.6619    0.0000 [ 0.5913; 0.7373 ] 
#>   eta3 ~ eta1      0.4585      0.0671    6.8376    0.0000 [ 0.3037; 0.5875 ] 
#>   eta3 ~ eta2      0.3052      0.0675    4.5223    0.0000 [ 0.1994; 0.4694 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0397   16.7121    0.0000 [ 0.5773; 0.7321 ] 
#>   eta1 =~ y12      0.6493      0.0385   16.8765    0.0000 [ 0.5656; 0.7149 ] 
#>   eta1 =~ y13      0.7613      0.0356   21.4040    0.0000 [ 0.6784; 0.8229 ] 
#>   eta2 =~ y21      0.5165      0.0514   10.0384    0.0000 [ 0.4390; 0.6244 ] 
#>   eta2 =~ y22      0.7554      0.0420   17.9657    0.0000 [ 0.6817; 0.8232 ] 
#>   eta2 =~ y23      0.7997      0.0290   27.5366    0.0000 [ 0.7453; 0.8504 ] 
#>   eta3 =~ y31      0.8223      0.0355   23.1853    0.0000 [ 0.7493; 0.8931 ] 
#>   eta3 =~ y32      0.6581      0.0370   17.8030    0.0000 [ 0.6083; 0.7402 ] 
#>   eta3 =~ y33      0.7474      0.0391   19.1268    0.0000 [ 0.6556; 0.8218 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0179   22.1333    0.0000 [ 0.3609; 0.4317 ] 
#>   eta1 <~ y12      0.3873      0.0185   20.8871    0.0000 [ 0.3586; 0.4294 ] 
#>   eta1 <~ y13      0.4542      0.0223   20.3507    0.0000 [ 0.4203; 0.5105 ] 
#>   eta2 <~ y21      0.3058      0.0255   11.9987    0.0000 [ 0.2676; 0.3567 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6075    0.0000 [ 0.4024; 0.4824 ] 
#>   eta2 <~ y23      0.4735      0.0176   26.9459    0.0000 [ 0.4380; 0.5079 ] 
#>   eta3 <~ y31      0.4400      0.0143   30.8205    0.0000 [ 0.4121; 0.4625 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9261    0.0000 [ 0.3274; 0.4004 ] 
#>   eta3 <~ y33      0.3999      0.0218   18.3819    0.0000 [ 0.3611; 0.4357 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0380   17.6619    0.0000 [ 0.5913; 0.7373 ] 
#>   eta3 ~ eta1       0.6634      0.0352   18.8697    0.0000 [ 0.6091; 0.7332 ] 
#>   eta3 ~ eta2       0.3052      0.0675    4.5223    0.0000 [ 0.1994; 0.4694 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0446    4.5909    0.0000 [ 0.1449; 0.3220 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03967613 16.71206  1.070699e-62
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03847229 16.87651  6.699151e-64
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03557030 21.40398 1.226832e-101
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05144816 10.03835  1.033850e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04204612 17.96569  3.617480e-72
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02904000 27.53663 6.398766e-167
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03546550 23.18527 6.411219e-119
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03696399 17.80297  6.701414e-71
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03907733 19.12680  1.510806e-81
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5772694          0.7320997
#> 2          0.5656036          0.7148632
#> 3          0.6784435          0.8228693
#> 4          0.4389976          0.6243546
#> 5          0.6816818          0.8232248
#> 6          0.7452517          0.8503611
#> 7          0.7493006          0.8931082
#> 8          0.6083350          0.7402329
#> 9          0.6556429          0.8218336

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
#>  Random seed                        = 992949468
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
#>   eta2 ~ eta1      0.6713      0.0380   17.6619    0.0000 [ 0.5770; 0.7735 ] 
#>   eta3 ~ eta1      0.4585      0.0671    6.8376    0.0000 [ 0.2906; 0.6374 ] 
#>   eta3 ~ eta2      0.3052      0.0675    4.5223    0.0000 [ 0.1248; 0.4737 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0397   16.7121    0.0000 [ 0.5607; 0.7659 ] 
#>   eta1 =~ y12      0.6493      0.0385   16.8765    0.0000 [ 0.5481; 0.7470 ] 
#>   eta1 =~ y13      0.7613      0.0356   21.4040    0.0000 [ 0.6611; 0.8451 ] 
#>   eta2 =~ y21      0.5165      0.0514   10.0384    0.0000 [ 0.3777; 0.6438 ] 
#>   eta2 =~ y22      0.7554      0.0420   17.9657    0.0000 [ 0.6463; 0.8637 ] 
#>   eta2 =~ y23      0.7997      0.0290   27.5366    0.0000 [ 0.7253; 0.8755 ] 
#>   eta3 =~ y31      0.8223      0.0355   23.1853    0.0000 [ 0.7386; 0.9220 ] 
#>   eta3 =~ y32      0.6581      0.0370   17.8030    0.0000 [ 0.5595; 0.7507 ] 
#>   eta3 =~ y33      0.7474      0.0391   19.1268    0.0000 [ 0.6490; 0.8511 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0179   22.1333    0.0000 [ 0.3526; 0.4450 ] 
#>   eta1 <~ y12      0.3873      0.0185   20.8871    0.0000 [ 0.3413; 0.4372 ] 
#>   eta1 <~ y13      0.4542      0.0223   20.3507    0.0000 [ 0.3948; 0.5102 ] 
#>   eta2 <~ y21      0.3058      0.0255   11.9987    0.0000 [ 0.2381; 0.3700 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6075    0.0000 [ 0.3898; 0.5078 ] 
#>   eta2 <~ y23      0.4735      0.0176   26.9459    0.0000 [ 0.4302; 0.5211 ] 
#>   eta3 <~ y31      0.4400      0.0143   30.8205    0.0000 [ 0.4053; 0.4791 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9261    0.0000 [ 0.2979; 0.3995 ] 
#>   eta3 <~ y33      0.3999      0.0218   18.3819    0.0000 [ 0.3430; 0.4555 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0380   17.6619    0.0000 [ 0.5770; 0.7735 ] 
#>   eta3 ~ eta1       0.6634      0.0352   18.8697    0.0000 [ 0.5759; 0.7578 ] 
#>   eta3 ~ eta2       0.3052      0.0675    4.5223    0.0000 [ 0.1248; 0.4737 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0446    4.5909    0.0000 [ 0.0874; 0.3182 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03801017 17.661943 8.235463e-70
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06705714  6.837553 8.055705e-12
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.06747748  4.522266 6.118104e-06
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5769652          0.7735324          0.6005691          0.7499285
#> 2          0.2906476          0.6374292          0.3322893          0.5957875
#> 3          0.1247589          0.4737143          0.1666616          0.4318115
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5912811          0.7595306          0.5913190          0.7372885
#> 2          0.2913931          0.6101472          0.3037439          0.5874601
#> 3          0.1625976          0.4983694          0.1994227          0.4694368
```
