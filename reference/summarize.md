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
#>  Random seed                        = -2019400107
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
#>   eta2 ~ eta1      0.6713      0.0497   13.5186    0.0000 [ 0.5873; 0.7667 ] 
#>   eta3 ~ eta1      0.4585      0.0613    7.4851    0.0000 [ 0.3455; 0.5562 ] 
#>   eta3 ~ eta2      0.3052      0.0703    4.3426    0.0000 [ 0.1726; 0.4265 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0422   15.7216    0.0000 [ 0.5905; 0.7221 ] 
#>   eta1 =~ y12      0.6493      0.0410   15.8206    0.0000 [ 0.5773; 0.7250 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.1829    0.0000 [ 0.6968; 0.8281 ] 
#>   eta2 =~ y21      0.5165      0.0562    9.1838    0.0000 [ 0.4219; 0.5970 ] 
#>   eta2 =~ y22      0.7554      0.0389   19.4393    0.0000 [ 0.6886; 0.8242 ] 
#>   eta2 =~ y23      0.7997      0.0388   20.6314    0.0000 [ 0.7481; 0.8683 ] 
#>   eta3 =~ y31      0.8223      0.0424   19.3791    0.0000 [ 0.7411; 0.8903 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1275    0.0000 [ 0.5795; 0.7313 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2281    0.0000 [ 0.6637; 0.8106 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0226   17.5023    0.0000 [ 0.3554; 0.4249 ] 
#>   eta1 <~ y12      0.3873      0.0198   19.5681    0.0000 [ 0.3608; 0.4232 ] 
#>   eta1 <~ y13      0.4542      0.0229   19.7941    0.0000 [ 0.4188; 0.5025 ] 
#>   eta2 <~ y21      0.3058      0.0300   10.1948    0.0000 [ 0.2571; 0.3565 ] 
#>   eta2 <~ y22      0.4473      0.0223   20.0174    0.0000 [ 0.4147; 0.4824 ] 
#>   eta2 <~ y23      0.4735      0.0226   20.9251    0.0000 [ 0.4258; 0.5135 ] 
#>   eta3 <~ y31      0.4400      0.0227   19.3868    0.0000 [ 0.4009; 0.4819 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.4923    0.0000 [ 0.3147; 0.3851 ] 
#>   eta3 <~ y33      0.3999      0.0207   19.3355    0.0000 [ 0.3670; 0.4448 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0497   13.5186    0.0000 [ 0.5873; 0.7667 ] 
#>   eta3 ~ eta1       0.6634      0.0356   18.6211    0.0000 [ 0.5912; 0.7382 ] 
#>   eta3 ~ eta2       0.3052      0.0703    4.3426    0.0000 [ 0.1726; 0.4265 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0501    4.0886    0.0000 [ 0.1178; 0.2999 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04217569 15.721612 1.075513e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04103997 15.820626 2.242529e-56
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03594158 21.182872 1.373949e-99
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05623556  9.183776 4.162344e-20
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03885882 19.439284 3.591415e-84
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03875947 20.631439 1.433218e-94
#> 7 eta3 =~ y31  Common factor 0.8222773 0.04243113 19.379106 1.158417e-83
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04080404 16.127542 1.633948e-58
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04100392 18.228113 3.087958e-74
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5904701          0.7221343
#> 2          0.5772653          0.7250272
#> 3          0.6968106          0.8281119
#> 4          0.4219402          0.5969587
#> 5          0.6886398          0.8242473
#> 6          0.7481103          0.8682984
#> 7          0.7411189          0.8902970
#> 8          0.5795189          0.7313403
#> 9          0.6637124          0.8105933

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
#>  Random seed                        = -2019400107
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
#>   eta2 ~ eta1      0.6713      0.0497   13.5186    0.0000 [ 0.5380; 0.7949 ] 
#>   eta3 ~ eta1      0.4585      0.0613    7.4851    0.0000 [ 0.2988; 0.6155 ] 
#>   eta3 ~ eta2      0.3052      0.0703    4.3426    0.0000 [ 0.1211; 0.4844 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0422   15.7216    0.0000 [ 0.5586; 0.7767 ] 
#>   eta1 =~ y12      0.6493      0.0410   15.8206    0.0000 [ 0.5326; 0.7448 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.1829    0.0000 [ 0.6732; 0.8590 ] 
#>   eta2 =~ y21      0.5165      0.0562    9.1838    0.0000 [ 0.3677; 0.6586 ] 
#>   eta2 =~ y22      0.7554      0.0389   19.4393    0.0000 [ 0.6618; 0.8628 ] 
#>   eta2 =~ y23      0.7997      0.0388   20.6314    0.0000 [ 0.6955; 0.8959 ] 
#>   eta3 =~ y31      0.8223      0.0424   19.3791    0.0000 [ 0.7136; 0.9330 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1275    0.0000 [ 0.5536; 0.7646 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2281    0.0000 [ 0.6483; 0.8603 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0226   17.5023    0.0000 [ 0.3404; 0.4573 ] 
#>   eta1 <~ y12      0.3873      0.0198   19.5681    0.0000 [ 0.3305; 0.4329 ] 
#>   eta1 <~ y13      0.4542      0.0229   19.7941    0.0000 [ 0.3981; 0.5167 ] 
#>   eta2 <~ y21      0.3058      0.0300   10.1948    0.0000 [ 0.2269; 0.3820 ] 
#>   eta2 <~ y22      0.4473      0.0223   20.0174    0.0000 [ 0.3940; 0.5096 ] 
#>   eta2 <~ y23      0.4735      0.0226   20.9251    0.0000 [ 0.4130; 0.5300 ] 
#>   eta3 <~ y31      0.4400      0.0227   19.3868    0.0000 [ 0.3793; 0.4966 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.4923    0.0000 [ 0.3016; 0.4000 ] 
#>   eta3 <~ y33      0.3999      0.0207   19.3355    0.0000 [ 0.3479; 0.4549 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0497   13.5186    0.0000 [ 0.5380; 0.7949 ] 
#>   eta3 ~ eta1       0.6634      0.0356   18.6211    0.0000 [ 0.5669; 0.7511 ] 
#>   eta3 ~ eta2       0.3052      0.0703    4.3426    0.0000 [ 0.1211; 0.4844 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0501    4.0886    0.0000 [ 0.0723; 0.3314 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04965990 13.518621 1.214309e-41
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06125622  7.485065 7.151195e-14
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07026865  4.342636 1.407836e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5380429          0.7948559          0.5688811          0.7640177
#> 2          0.2987571          0.6155397          0.3367965          0.5775003
#> 3          0.1210523          0.4844421          0.1646883          0.4408060
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5470602          0.7897117          0.5872631          0.7666756
#> 2          0.3448629          0.5691468          0.3454970          0.5561597
#> 3          0.1674138          0.4357728          0.1725656          0.4264921
```
