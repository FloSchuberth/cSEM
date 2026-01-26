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
#>  Random seed                        = 180679782
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
#>   eta2 ~ eta1      0.6713      0.0342   19.6559    0.0000 [ 0.6298; 0.7324 ] 
#>   eta3 ~ eta1      0.4585      0.0845    5.4286    0.0000 [ 0.2820; 0.5993 ] 
#>   eta3 ~ eta2      0.3052      0.0841    3.6289    0.0003 [ 0.1583; 0.4504 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0408   16.2403    0.0000 [ 0.5932; 0.7456 ] 
#>   eta1 =~ y12      0.6493      0.0390   16.6317    0.0000 [ 0.5612; 0.6982 ] 
#>   eta1 =~ y13      0.7613      0.0338   22.5489    0.0000 [ 0.6980; 0.8092 ] 
#>   eta2 =~ y21      0.5165      0.0372   13.8778    0.0000 [ 0.4312; 0.5665 ] 
#>   eta2 =~ y22      0.7554      0.0374   20.1894    0.0000 [ 0.6871; 0.8041 ] 
#>   eta2 =~ y23      0.7997      0.0363   22.0342    0.0000 [ 0.7357; 0.8659 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.7202    0.0000 [ 0.7526; 0.8735 ] 
#>   eta3 =~ y32      0.6581      0.0387   17.0205    0.0000 [ 0.6006; 0.7457 ] 
#>   eta3 =~ y33      0.7474      0.0333   22.4461    0.0000 [ 0.6837; 0.7930 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0223   17.7030    0.0000 [ 0.3641; 0.4401 ] 
#>   eta1 <~ y12      0.3873      0.0173   22.4062    0.0000 [ 0.3498; 0.4214 ] 
#>   eta1 <~ y13      0.4542      0.0190   23.9183    0.0000 [ 0.4264; 0.4852 ] 
#>   eta2 <~ y21      0.3058      0.0197   15.5389    0.0000 [ 0.2633; 0.3318 ] 
#>   eta2 <~ y22      0.4473      0.0164   27.3241    0.0000 [ 0.4191; 0.4863 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.2336    0.0000 [ 0.4373; 0.5174 ] 
#>   eta3 <~ y31      0.4400      0.0174   25.2263    0.0000 [ 0.4083; 0.4631 ] 
#>   eta3 <~ y32      0.3521      0.0158   22.2679    0.0000 [ 0.3273; 0.3874 ] 
#>   eta3 <~ y33      0.3999      0.0195   20.4791    0.0000 [ 0.3648; 0.4361 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0342   19.6559    0.0000 [ 0.6298; 0.7324 ] 
#>   eta3 ~ eta1       0.6634      0.0413   16.0767    0.0000 [ 0.5735; 0.7236 ] 
#>   eta3 ~ eta2       0.3052      0.0841    3.6289    0.0003 [ 0.1583; 0.4504 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0568    3.6043    0.0003 [ 0.1054; 0.2966 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04082858 16.24034  2.614817e-59
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03903854 16.63172  4.106616e-62
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03376423 22.54889 1.376965e-112
#> 4 eta2 =~ y21  Common factor 0.5164548 0.03721450 13.87778  8.637437e-44
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03741513 20.18936  1.214190e-90
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03629196 22.03418 1.354644e-107
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03326342 24.72016 6.492651e-135
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03866327 17.02052  5.785759e-65
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03329861 22.44610 1.396915e-111
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5932325          0.7455876
#> 2          0.5612394          0.6982418
#> 3          0.6979840          0.8092481
#> 4          0.4311776          0.5665361
#> 5          0.6870898          0.8041198
#> 6          0.7356527          0.8658530
#> 7          0.7526452          0.8734555
#> 8          0.6006428          0.7457300
#> 9          0.6836845          0.7930326

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
#>  Random seed                        = 180679782
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
#>   eta2 ~ eta1      0.6713      0.0342   19.6559    0.0000 [ 0.5751; 0.7517 ] 
#>   eta3 ~ eta1      0.4585      0.0845    5.4286    0.0000 [ 0.2400; 0.6768 ] 
#>   eta3 ~ eta2      0.3052      0.0841    3.6289    0.0003 [ 0.0952; 0.5300 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0408   16.2403    0.0000 [ 0.5547; 0.7658 ] 
#>   eta1 =~ y12      0.6493      0.0390   16.6317    0.0000 [ 0.5514; 0.7533 ] 
#>   eta1 =~ y13      0.7613      0.0338   22.5489    0.0000 [ 0.6748; 0.8495 ] 
#>   eta2 =~ y21      0.5165      0.0372   13.8778    0.0000 [ 0.4340; 0.6265 ] 
#>   eta2 =~ y22      0.7554      0.0374   20.1894    0.0000 [ 0.6682; 0.8617 ] 
#>   eta2 =~ y23      0.7997      0.0363   22.0342    0.0000 [ 0.7093; 0.8970 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.7202    0.0000 [ 0.7453; 0.9174 ] 
#>   eta3 =~ y32      0.6581      0.0387   17.0205    0.0000 [ 0.5464; 0.7464 ] 
#>   eta3 =~ y33      0.7474      0.0333   22.4461    0.0000 [ 0.6670; 0.8392 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0223   17.7030    0.0000 [ 0.3359; 0.4515 ] 
#>   eta1 <~ y12      0.3873      0.0173   22.4062    0.0000 [ 0.3445; 0.4339 ] 
#>   eta1 <~ y13      0.4542      0.0190   23.9183    0.0000 [ 0.4053; 0.5035 ] 
#>   eta2 <~ y21      0.3058      0.0197   15.5389    0.0000 [ 0.2578; 0.3596 ] 
#>   eta2 <~ y22      0.4473      0.0164   27.3241    0.0000 [ 0.4027; 0.4874 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.2336    0.0000 [ 0.4142; 0.5196 ] 
#>   eta3 <~ y31      0.4400      0.0174   25.2263    0.0000 [ 0.3987; 0.4889 ] 
#>   eta3 <~ y32      0.3521      0.0158   22.2679    0.0000 [ 0.3044; 0.3862 ] 
#>   eta3 <~ y33      0.3999      0.0195   20.4791    0.0000 [ 0.3515; 0.4524 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0342   19.6559    0.0000 [ 0.5751; 0.7517 ] 
#>   eta3 ~ eta1       0.6634      0.0413   16.0767    0.0000 [ 0.5594; 0.7728 ] 
#>   eta3 ~ eta2       0.3052      0.0841    3.6289    0.0003 [ 0.0952; 0.5300 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0568    3.6043    0.0003 [ 0.0608; 0.3547 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03415428 19.655908 5.145995e-86
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08446157  5.428584 5.680285e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08408887  3.628912 2.846179e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5750851          0.7517117          0.5962945          0.7305023
#> 2          0.2399912          0.6767788          0.2924409          0.6243292
#> 3          0.0951752          0.5300353          0.1473934          0.4778171
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5879146          0.7467598          0.6297849          0.7323918
#> 2          0.2625687          0.6023628          0.2820494          0.5993464
#> 3          0.1148813          0.4517997          0.1582745          0.4503662
```
