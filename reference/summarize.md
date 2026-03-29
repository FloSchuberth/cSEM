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
#>  Random seed                        = 1203308964
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
#>   eta2 ~ eta1      0.6713      0.0464   14.4693    0.0000 [ 0.5902; 0.7401 ] 
#>   eta3 ~ eta1      0.4585      0.0857    5.3478    0.0000 [ 0.3069; 0.6222 ] 
#>   eta3 ~ eta2      0.3052      0.0995    3.0654    0.0022 [ 0.0581; 0.4639 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0394   16.8312    0.0000 [ 0.5883; 0.7235 ] 
#>   eta1 =~ y12      0.6493      0.0406   16.0029    0.0000 [ 0.5563; 0.6949 ] 
#>   eta1 =~ y13      0.7613      0.0369   20.6292    0.0000 [ 0.6969; 0.8348 ] 
#>   eta2 =~ y21      0.5165      0.0553    9.3344    0.0000 [ 0.4180; 0.6065 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1421    0.0000 [ 0.6707; 0.7911 ] 
#>   eta2 =~ y23      0.7997      0.0460   17.3911    0.0000 [ 0.6943; 0.8661 ] 
#>   eta3 =~ y31      0.8223      0.0339   24.2912    0.0000 [ 0.7720; 0.8844 ] 
#>   eta3 =~ y32      0.6581      0.0472   13.9296    0.0000 [ 0.5816; 0.7432 ] 
#>   eta3 =~ y33      0.7474      0.0413   18.0886    0.0000 [ 0.6813; 0.8347 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0236   16.7538    0.0000 [ 0.3604; 0.4426 ] 
#>   eta1 <~ y12      0.3873      0.0173   22.3406    0.0000 [ 0.3592; 0.4104 ] 
#>   eta1 <~ y13      0.4542      0.0227   20.0295    0.0000 [ 0.4253; 0.5097 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5637    0.0000 [ 0.2536; 0.3788 ] 
#>   eta2 <~ y22      0.4473      0.0226   19.8085    0.0000 [ 0.4049; 0.4815 ] 
#>   eta2 <~ y23      0.4735      0.0225   21.0342    0.0000 [ 0.4299; 0.5105 ] 
#>   eta3 <~ y31      0.4400      0.0214   20.5803    0.0000 [ 0.4125; 0.4846 ] 
#>   eta3 <~ y32      0.3521      0.0219   16.0676    0.0000 [ 0.3116; 0.3837 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.8833    0.0000 [ 0.3578; 0.4381 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0464   14.4693    0.0000 [ 0.5902; 0.7401 ] 
#>   eta3 ~ eta1       0.6634      0.0318   20.8628    0.0000 [ 0.6090; 0.7272 ] 
#>   eta3 ~ eta2       0.3052      0.0995    3.0654    0.0022 [ 0.0581; 0.4639 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0701    2.9238    0.0035 [ 0.0426; 0.3376 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03939527 16.831206  1.441355e-63
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04057248 16.002916  1.219272e-57
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03690630 20.629157  1.502470e-94
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05532795  9.334428  1.015379e-20
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03572912 21.142075  3.264075e-99
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04598116 17.391118  9.633667e-68
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03385087 24.291173 2.430110e-130
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04724264 13.929553  4.189685e-44
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04132016 18.088607  3.918674e-73
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5882564          0.7235266
#> 2          0.5563273          0.6949384
#> 3          0.6969439          0.8348064
#> 4          0.4179581          0.6065165
#> 5          0.6707083          0.7911144
#> 6          0.6943186          0.8661270
#> 7          0.7719722          0.8843756
#> 8          0.5815783          0.7432356
#> 9          0.6812673          0.8347086

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
#>  Random seed                        = 1203308964
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
#>   eta2 ~ eta1      0.6713      0.0464   14.4693    0.0000 [ 0.5491; 0.7890 ] 
#>   eta3 ~ eta1      0.4585      0.0857    5.3478    0.0000 [ 0.2465; 0.6899 ] 
#>   eta3 ~ eta2      0.3052      0.0995    3.0654    0.0022 [ 0.0397; 0.5545 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0394   16.8312    0.0000 [ 0.5670; 0.7707 ] 
#>   eta1 =~ y12      0.6493      0.0406   16.0029    0.0000 [ 0.5560; 0.7658 ] 
#>   eta1 =~ y13      0.7613      0.0369   20.6292    0.0000 [ 0.6606; 0.8515 ] 
#>   eta2 =~ y21      0.5165      0.0553    9.3344    0.0000 [ 0.3711; 0.6573 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1421    0.0000 [ 0.6831; 0.8679 ] 
#>   eta2 =~ y23      0.7997      0.0460   17.3911    0.0000 [ 0.6906; 0.9284 ] 
#>   eta3 =~ y31      0.8223      0.0339   24.2912    0.0000 [ 0.7318; 0.9069 ] 
#>   eta3 =~ y32      0.6581      0.0472   13.9296    0.0000 [ 0.5362; 0.7805 ] 
#>   eta3 =~ y33      0.7474      0.0413   18.0886    0.0000 [ 0.6424; 0.8561 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0236   16.7538    0.0000 [ 0.3348; 0.4569 ] 
#>   eta1 <~ y12      0.3873      0.0173   22.3406    0.0000 [ 0.3468; 0.4364 ] 
#>   eta1 <~ y13      0.4542      0.0227   20.0295    0.0000 [ 0.3886; 0.5059 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5637    0.0000 [ 0.2157; 0.3811 ] 
#>   eta2 <~ y22      0.4473      0.0226   19.8085    0.0000 [ 0.3919; 0.5087 ] 
#>   eta2 <~ y23      0.4735      0.0225   21.0342    0.0000 [ 0.4119; 0.5283 ] 
#>   eta3 <~ y31      0.4400      0.0214   20.5803    0.0000 [ 0.3835; 0.4940 ] 
#>   eta3 <~ y32      0.3521      0.0219   16.0676    0.0000 [ 0.2962; 0.4096 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.8833    0.0000 [ 0.3494; 0.4534 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0464   14.4693    0.0000 [ 0.5491; 0.7890 ] 
#>   eta3 ~ eta1       0.6634      0.0318   20.8628    0.0000 [ 0.5849; 0.7493 ] 
#>   eta3 ~ eta2       0.3052      0.0995    3.0654    0.0022 [ 0.0397; 0.5545 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0701    2.9238    0.0035 [ 0.0178; 0.3801 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04639703 14.469318 1.893391e-47
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08573794  5.347770 8.904469e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09954728  3.065389 2.173872e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.54909342          0.7890327          0.5779054          0.7602207
#> 2         0.24646743          0.6898557          0.2997097          0.6366134
#> 3         0.03969181          0.5544941          0.1015095          0.4926764
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58316597          0.7600852         0.59018943          0.7401136
#> 2         0.27362751          0.6787858         0.30687818          0.6222216
#> 3         0.02888848          0.4897524         0.05806914          0.4638641
```
