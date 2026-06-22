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
#>  Random seed                        = 947384364
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
#>   eta2 ~ eta1      0.6713      0.0425   15.7859    0.0000 [ 0.5884; 0.7378 ] 
#>   eta3 ~ eta1      0.4585      0.0798    5.7482    0.0000 [ 0.3459; 0.6145 ] 
#>   eta3 ~ eta2      0.3052      0.0964    3.1645    0.0016 [ 0.1005; 0.4522 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0322   20.5997    0.0000 [ 0.6035; 0.7087 ] 
#>   eta1 =~ y12      0.6493      0.0370   17.5678    0.0000 [ 0.5655; 0.7056 ] 
#>   eta1 =~ y13      0.7613      0.0311   24.4555    0.0000 [ 0.6991; 0.8092 ] 
#>   eta2 =~ y21      0.5165      0.0519    9.9538    0.0000 [ 0.4004; 0.5967 ] 
#>   eta2 =~ y22      0.7554      0.0313   24.1420    0.0000 [ 0.6999; 0.8064 ] 
#>   eta2 =~ y23      0.7997      0.0413   19.3516    0.0000 [ 0.7208; 0.8681 ] 
#>   eta3 =~ y31      0.8223      0.0355   23.1483    0.0000 [ 0.7491; 0.8817 ] 
#>   eta3 =~ y32      0.6581      0.0436   15.1103    0.0000 [ 0.5691; 0.7345 ] 
#>   eta3 =~ y33      0.7474      0.0391   19.0989    0.0000 [ 0.6798; 0.8226 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0182   21.7707    0.0000 [ 0.3737; 0.4180 ] 
#>   eta1 <~ y12      0.3873      0.0149   26.0755    0.0000 [ 0.3637; 0.4170 ] 
#>   eta1 <~ y13      0.4542      0.0194   23.3656    0.0000 [ 0.4258; 0.5033 ] 
#>   eta2 <~ y21      0.3058      0.0265   11.5429    0.0000 [ 0.2584; 0.3494 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.5998    0.0000 [ 0.4097; 0.4940 ] 
#>   eta2 <~ y23      0.4735      0.0211   22.4248    0.0000 [ 0.4444; 0.5088 ] 
#>   eta3 <~ y31      0.4400      0.0219   20.1212    0.0000 [ 0.4017; 0.4702 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.6025    0.0000 [ 0.3167; 0.3820 ] 
#>   eta3 <~ y33      0.3999      0.0192   20.7965    0.0000 [ 0.3605; 0.4270 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0425   15.7859    0.0000 [ 0.5884; 0.7378 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5529    0.0000 [ 0.6096; 0.7388 ] 
#>   eta3 ~ eta2       0.3052      0.0964    3.1645    0.0016 [ 0.1005; 0.4522 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0607    3.3776    0.0007 [ 0.0643; 0.2782 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03218828 20.599730  2.759751e-94
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03695851 17.567752  4.350703e-69
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03113192 24.455468 4.401537e-132
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05188523  9.953792  2.427529e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03128931 24.142036 9.052682e-129
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04132280 19.351635  1.974786e-83
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03552215 23.148301 1.512171e-118
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04355098 15.110310  1.384857e-51
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03913441 19.098897  2.578795e-81
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6035220          0.7086877
#> 2          0.5654981          0.7055909
#> 3          0.6991230          0.8092421
#> 4          0.4003596          0.5966745
#> 5          0.6998661          0.8063870
#> 6          0.7207815          0.8681480
#> 7          0.7491438          0.8816671
#> 8          0.5690642          0.7344508
#> 9          0.6798359          0.8226157

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
#>  Random seed                        = 947384364
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
#>   eta2 ~ eta1      0.6713      0.0425   15.7859    0.0000 [ 0.5608; 0.7808 ] 
#>   eta3 ~ eta1      0.4585      0.0798    5.7482    0.0000 [ 0.2225; 0.6350 ] 
#>   eta3 ~ eta2      0.3052      0.0964    3.1645    0.0016 [ 0.0850; 0.5837 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0322   20.5997    0.0000 [ 0.5872; 0.7537 ] 
#>   eta1 =~ y12      0.6493      0.0370   17.5678    0.0000 [ 0.5490; 0.7401 ] 
#>   eta1 =~ y13      0.7613      0.0311   24.4555    0.0000 [ 0.6911; 0.8521 ] 
#>   eta2 =~ y21      0.5165      0.0519    9.9538    0.0000 [ 0.3997; 0.6680 ] 
#>   eta2 =~ y22      0.7554      0.0313   24.1420    0.0000 [ 0.6752; 0.8370 ] 
#>   eta2 =~ y23      0.7997      0.0413   19.3516    0.0000 [ 0.6911; 0.9048 ] 
#>   eta3 =~ y31      0.8223      0.0355   23.1483    0.0000 [ 0.7305; 0.9142 ] 
#>   eta3 =~ y32      0.6581      0.0436   15.1103    0.0000 [ 0.5493; 0.7745 ] 
#>   eta3 =~ y33      0.7474      0.0391   19.0989    0.0000 [ 0.6490; 0.8514 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0182   21.7707    0.0000 [ 0.3493; 0.4433 ] 
#>   eta1 <~ y12      0.3873      0.0149   26.0755    0.0000 [ 0.3427; 0.4196 ] 
#>   eta1 <~ y13      0.4542      0.0194   23.3656    0.0000 [ 0.4058; 0.5063 ] 
#>   eta2 <~ y21      0.3058      0.0265   11.5429    0.0000 [ 0.2453; 0.3823 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.5998    0.0000 [ 0.3812; 0.5055 ] 
#>   eta2 <~ y23      0.4735      0.0211   22.4248    0.0000 [ 0.4137; 0.5229 ] 
#>   eta3 <~ y31      0.4400      0.0219   20.1212    0.0000 [ 0.3815; 0.4946 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.6025    0.0000 [ 0.3040; 0.4019 ] 
#>   eta3 <~ y33      0.3999      0.0192   20.7965    0.0000 [ 0.3501; 0.4496 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0425   15.7859    0.0000 [ 0.5608; 0.7808 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5529    0.0000 [ 0.5621; 0.7470 ] 
#>   eta3 ~ eta2       0.3052      0.0964    3.1645    0.0016 [ 0.0850; 0.5837 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0607    3.3776    0.0007 [ 0.0690; 0.3826 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04252750 15.78587 3.892705e-56
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07976486  5.74823 9.018231e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09643042  3.16447 1.553659e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5608326          0.7807609          0.5872417          0.7543518
#> 2          0.2225329          0.6350317          0.2720659          0.5854986
#> 3          0.0850419          0.5837256          0.1449241          0.5238434
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58393669          0.7472615          0.5884389          0.7377857
#> 2         0.33508175          0.6415317          0.3458841          0.6144842
#> 3         0.06550328          0.4881467          0.1004945          0.4521768
```
