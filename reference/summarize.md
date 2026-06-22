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
#>  Random seed                        = -527017318
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
#>   eta2 ~ eta1      0.6713      0.0425   15.7968    0.0000 [ 0.5912; 0.7537 ] 
#>   eta3 ~ eta1      0.4585      0.0811    5.6550    0.0000 [ 0.2994; 0.5837 ] 
#>   eta3 ~ eta2      0.3052      0.0847    3.6029    0.0003 [ 0.1829; 0.4695 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0343   19.3542    0.0000 [ 0.5883; 0.7084 ] 
#>   eta1 =~ y12      0.6493      0.0418   15.5213    0.0000 [ 0.5832; 0.7185 ] 
#>   eta1 =~ y13      0.7613      0.0353   21.5500    0.0000 [ 0.6952; 0.8324 ] 
#>   eta2 =~ y21      0.5165      0.0410   12.6040    0.0000 [ 0.4410; 0.5919 ] 
#>   eta2 =~ y22      0.7554      0.0326   23.1800    0.0000 [ 0.6872; 0.8164 ] 
#>   eta2 =~ y23      0.7997      0.0281   28.5016    0.0000 [ 0.7469; 0.8363 ] 
#>   eta3 =~ y31      0.8223      0.0339   24.2385    0.0000 [ 0.7665; 0.8741 ] 
#>   eta3 =~ y32      0.6581      0.0454   14.4919    0.0000 [ 0.5663; 0.7551 ] 
#>   eta3 =~ y33      0.7474      0.0420   17.8169    0.0000 [ 0.6719; 0.8169 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0204   19.3629    0.0000 [ 0.3540; 0.4249 ] 
#>   eta1 <~ y12      0.3873      0.0191   20.3031    0.0000 [ 0.3555; 0.4159 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9425    0.0000 [ 0.4112; 0.4958 ] 
#>   eta2 <~ y21      0.3058      0.0211   14.5225    0.0000 [ 0.2668; 0.3386 ] 
#>   eta2 <~ y22      0.4473      0.0180   24.8493    0.0000 [ 0.4157; 0.4855 ] 
#>   eta2 <~ y23      0.4735      0.0151   31.4112    0.0000 [ 0.4427; 0.5017 ] 
#>   eta3 <~ y31      0.4400      0.0200   22.0213    0.0000 [ 0.3944; 0.4736 ] 
#>   eta3 <~ y32      0.3521      0.0222   15.8615    0.0000 [ 0.3019; 0.3832 ] 
#>   eta3 <~ y33      0.3999      0.0180   22.2350    0.0000 [ 0.3704; 0.4320 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0425   15.7968    0.0000 [ 0.5912; 0.7537 ] 
#>   eta3 ~ eta1       0.6634      0.0348   19.0521    0.0000 [ 0.5892; 0.7185 ] 
#>   eta3 ~ eta2       0.3052      0.0847    3.6029    0.0003 [ 0.1829; 0.4695 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0596    3.4366    0.0006 [ 0.1337; 0.3176 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03425968 19.35424  1.877571e-83
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04183142 15.52130  2.489631e-54
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03532923 21.55003 5.291574e-103
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04097556 12.60397  2.007757e-36
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03258784 23.18005 7.238729e-119
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02805681 28.50159 1.119399e-178
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03392440 24.23852 8.738139e-130
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04540931 14.49194  1.362457e-47
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04195022 17.81693  5.222511e-71
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5882928          0.7084067
#> 2          0.5832213          0.7185418
#> 3          0.6951806          0.8323649
#> 4          0.4410489          0.5919163
#> 5          0.6872241          0.8164143
#> 6          0.7469057          0.8362926
#> 7          0.7665246          0.8741133
#> 8          0.5662910          0.7550555
#> 9          0.6719395          0.8168708

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
#>  Random seed                        = -527017318
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
#>   eta2 ~ eta1      0.6713      0.0425   15.7968    0.0000 [ 0.5556; 0.7753 ] 
#>   eta3 ~ eta1      0.4585      0.0811    5.6550    0.0000 [ 0.2898; 0.7091 ] 
#>   eta3 ~ eta2      0.3052      0.0847    3.6029    0.0003 [ 0.0423; 0.4803 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0343   19.3542    0.0000 [ 0.5764; 0.7535 ] 
#>   eta1 =~ y12      0.6493      0.0418   15.5213    0.0000 [ 0.5401; 0.7564 ] 
#>   eta1 =~ y13      0.7613      0.0353   21.5500    0.0000 [ 0.6719; 0.8546 ] 
#>   eta2 =~ y21      0.5165      0.0410   12.6040    0.0000 [ 0.4081; 0.6200 ] 
#>   eta2 =~ y22      0.7554      0.0326   23.1800    0.0000 [ 0.6701; 0.8386 ] 
#>   eta2 =~ y23      0.7997      0.0281   28.5016    0.0000 [ 0.7320; 0.8771 ] 
#>   eta3 =~ y31      0.8223      0.0339   24.2385    0.0000 [ 0.7370; 0.9125 ] 
#>   eta3 =~ y32      0.6581      0.0454   14.4919    0.0000 [ 0.5403; 0.7752 ] 
#>   eta3 =~ y33      0.7474      0.0420   17.8169    0.0000 [ 0.6362; 0.8532 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0204   19.3629    0.0000 [ 0.3432; 0.4488 ] 
#>   eta1 <~ y12      0.3873      0.0191   20.3031    0.0000 [ 0.3370; 0.4357 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9425    0.0000 [ 0.3984; 0.5106 ] 
#>   eta2 <~ y21      0.3058      0.0211   14.5225    0.0000 [ 0.2498; 0.3587 ] 
#>   eta2 <~ y22      0.4473      0.0180   24.8493    0.0000 [ 0.3997; 0.4928 ] 
#>   eta2 <~ y23      0.4735      0.0151   31.4112    0.0000 [ 0.4369; 0.5149 ] 
#>   eta3 <~ y31      0.4400      0.0200   22.0213    0.0000 [ 0.3899; 0.4932 ] 
#>   eta3 <~ y32      0.3521      0.0222   15.8615    0.0000 [ 0.2949; 0.4097 ] 
#>   eta3 <~ y33      0.3999      0.0180   22.2350    0.0000 [ 0.3525; 0.4455 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0425   15.7968    0.0000 [ 0.5556; 0.7753 ] 
#>   eta3 ~ eta1       0.6634      0.0348   19.0521    0.0000 [ 0.5829; 0.7630 ] 
#>   eta3 ~ eta2       0.3052      0.0847    3.6029    0.0003 [ 0.0423; 0.4803 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0596    3.4366    0.0006 [ 0.0193; 0.3276 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04249795 15.796844 3.270849e-56
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08107950  5.655027 1.558218e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08469564  3.602914 3.146694e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55556106          0.7753364         0.58195179          0.7489457
#> 2         0.28982308          0.7091205         0.34017250          0.6587711
#> 3         0.04227126          0.4802693         0.09486626          0.4276743
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5893819          0.7707784          0.5911853          0.7537144
#> 2          0.2992726          0.6400037          0.2994087          0.5836923
#> 3          0.1547092          0.4821009          0.1829247          0.4695489
```
