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
#>  Random seed                        = -2062153737
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
#>   eta2 ~ eta1      0.6713      0.0492   13.6331    0.0000 [ 0.5853; 0.7583 ] 
#>   eta3 ~ eta1      0.4585      0.0825    5.5590    0.0000 [ 0.3238; 0.6262 ] 
#>   eta3 ~ eta2      0.3052      0.0910    3.3528    0.0008 [ 0.1079; 0.4575 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0397   16.7221    0.0000 [ 0.5877; 0.7036 ] 
#>   eta1 =~ y12      0.6493      0.0364   17.8404    0.0000 [ 0.5683; 0.7219 ] 
#>   eta1 =~ y13      0.7613      0.0340   22.3643    0.0000 [ 0.7124; 0.8289 ] 
#>   eta2 =~ y21      0.5165      0.0575    8.9784    0.0000 [ 0.4137; 0.5927 ] 
#>   eta2 =~ y22      0.7554      0.0392   19.2627    0.0000 [ 0.6894; 0.8221 ] 
#>   eta2 =~ y23      0.7997      0.0430   18.5987    0.0000 [ 0.7414; 0.8704 ] 
#>   eta3 =~ y31      0.8223      0.0337   24.4294    0.0000 [ 0.7685; 0.8799 ] 
#>   eta3 =~ y32      0.6581      0.0380   17.3011    0.0000 [ 0.5797; 0.7292 ] 
#>   eta3 =~ y33      0.7474      0.0401   18.6485    0.0000 [ 0.6492; 0.7964 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0220   17.9843    0.0000 [ 0.3401; 0.4294 ] 
#>   eta1 <~ y12      0.3873      0.0197   19.6831    0.0000 [ 0.3557; 0.4222 ] 
#>   eta1 <~ y13      0.4542      0.0192   23.6579    0.0000 [ 0.4287; 0.4900 ] 
#>   eta2 <~ y21      0.3058      0.0302   10.1148    0.0000 [ 0.2481; 0.3640 ] 
#>   eta2 <~ y22      0.4473      0.0231   19.3295    0.0000 [ 0.4171; 0.5036 ] 
#>   eta2 <~ y23      0.4735      0.0247   19.1675    0.0000 [ 0.4363; 0.5207 ] 
#>   eta3 <~ y31      0.4400      0.0181   24.3009    0.0000 [ 0.4094; 0.4740 ] 
#>   eta3 <~ y32      0.3521      0.0151   23.2752    0.0000 [ 0.3321; 0.3897 ] 
#>   eta3 <~ y33      0.3999      0.0219   18.2682    0.0000 [ 0.3499; 0.4289 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0492   13.6331    0.0000 [ 0.5853; 0.7583 ] 
#>   eta3 ~ eta1       0.6634      0.0414   16.0159    0.0000 [ 0.5763; 0.7240 ] 
#>   eta3 ~ eta2       0.3052      0.0910    3.3528    0.0008 [ 0.1079; 0.4575 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4292    0.0006 [ 0.0787; 0.3236 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03965238 16.722070  9.051926e-63
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03639370 17.840395  3.432594e-71
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03404291 22.364302 8.764479e-111
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05752187  8.978408  2.747129e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03921502 19.262714  1.104354e-82
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04299576 18.598663  3.294242e-77
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03365931 24.429416 8.329442e-132
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03803632 17.301066  4.617634e-67
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04007948 18.648548  1.297526e-77
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5877450          0.7035584
#> 2          0.5683123          0.7219321
#> 3          0.7123951          0.8288994
#> 4          0.4137029          0.5927347
#> 5          0.6893793          0.8220857
#> 6          0.7414355          0.8704154
#> 7          0.7685124          0.8799267
#> 8          0.5796707          0.7291575
#> 9          0.6491625          0.7964118

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
#>  Random seed                        = -2062153737
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
#>   eta2 ~ eta1      0.6713      0.0492   13.6331    0.0000 [ 0.5390; 0.7937 ] 
#>   eta3 ~ eta1      0.4585      0.0825    5.5590    0.0000 [ 0.2475; 0.6740 ] 
#>   eta3 ~ eta2      0.3052      0.0910    3.3528    0.0008 [ 0.0688; 0.5395 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0397   16.7221    0.0000 [ 0.5786; 0.7837 ] 
#>   eta1 =~ y12      0.6493      0.0364   17.8404    0.0000 [ 0.5555; 0.7437 ] 
#>   eta1 =~ y13      0.7613      0.0340   22.3643    0.0000 [ 0.6661; 0.8422 ] 
#>   eta2 =~ y21      0.5165      0.0575    8.9784    0.0000 [ 0.3744; 0.6719 ] 
#>   eta2 =~ y22      0.7554      0.0392   19.2627    0.0000 [ 0.6516; 0.8544 ] 
#>   eta2 =~ y23      0.7997      0.0430   18.5987    0.0000 [ 0.6865; 0.9089 ] 
#>   eta3 =~ y31      0.8223      0.0337   24.4294    0.0000 [ 0.7395; 0.9136 ] 
#>   eta3 =~ y32      0.6581      0.0380   17.3011    0.0000 [ 0.5552; 0.7519 ] 
#>   eta3 =~ y33      0.7474      0.0401   18.6485    0.0000 [ 0.6475; 0.8548 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0220   17.9843    0.0000 [ 0.3467; 0.4605 ] 
#>   eta1 <~ y12      0.3873      0.0197   19.6831    0.0000 [ 0.3339; 0.4357 ] 
#>   eta1 <~ y13      0.4542      0.0192   23.6579    0.0000 [ 0.3969; 0.4962 ] 
#>   eta2 <~ y21      0.3058      0.0302   10.1148    0.0000 [ 0.2320; 0.3884 ] 
#>   eta2 <~ y22      0.4473      0.0231   19.3295    0.0000 [ 0.3861; 0.5057 ] 
#>   eta2 <~ y23      0.4735      0.0247   19.1675    0.0000 [ 0.4085; 0.5362 ] 
#>   eta3 <~ y31      0.4400      0.0181   24.3009    0.0000 [ 0.3944; 0.4880 ] 
#>   eta3 <~ y32      0.3521      0.0151   23.2752    0.0000 [ 0.3100; 0.3882 ] 
#>   eta3 <~ y33      0.3999      0.0219   18.2682    0.0000 [ 0.3443; 0.4576 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0492   13.6331    0.0000 [ 0.5390; 0.7937 ] 
#>   eta3 ~ eta1       0.6634      0.0414   16.0159    0.0000 [ 0.5576; 0.7718 ] 
#>   eta3 ~ eta2       0.3052      0.0910    3.3528    0.0008 [ 0.0688; 0.5395 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4292    0.0006 [ 0.0495; 0.3584 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04924291 13.633099 2.545371e-42
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08248068  5.558959 2.713878e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09101504  3.352755 8.001149e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5390086          0.7936651          0.5695878          0.7630858
#> 2          0.2474780          0.6740215          0.2986975          0.6228019
#> 3          0.0687902          0.5394686          0.1253095          0.4829493
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58205074          0.7646893          0.5853074          0.7583212
#> 2         0.26072483          0.6536039          0.3238235          0.6261988
#> 3         0.07177521          0.5078281          0.1079038          0.4574874
```
