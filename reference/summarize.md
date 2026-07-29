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
#>  Random seed                        = -2072701448
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
#>   eta2 ~ eta1      0.6713      0.0453   14.8317    0.0000 [ 0.5937; 0.7488 ] 
#>   eta3 ~ eta1      0.4585      0.0774    5.9226    0.0000 [ 0.2806; 0.5911 ] 
#>   eta3 ~ eta2      0.3052      0.0842    3.6258    0.0003 [ 0.1814; 0.4638 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0421   15.7337    0.0000 [ 0.5610; 0.7142 ] 
#>   eta1 =~ y12      0.6493      0.0446   14.5711    0.0000 [ 0.5478; 0.7122 ] 
#>   eta1 =~ y13      0.7613      0.0372   20.4575    0.0000 [ 0.6919; 0.8181 ] 
#>   eta2 =~ y21      0.5165      0.0473   10.9156    0.0000 [ 0.4442; 0.5812 ] 
#>   eta2 =~ y22      0.7554      0.0332   22.7660    0.0000 [ 0.7171; 0.8338 ] 
#>   eta2 =~ y23      0.7997      0.0386   20.7141    0.0000 [ 0.7100; 0.8562 ] 
#>   eta3 =~ y31      0.8223      0.0334   24.6061    0.0000 [ 0.7595; 0.8764 ] 
#>   eta3 =~ y32      0.6581      0.0435   15.1211    0.0000 [ 0.5766; 0.7340 ] 
#>   eta3 =~ y33      0.7474      0.0365   20.4840    0.0000 [ 0.6590; 0.7986 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0215   18.4085    0.0000 [ 0.3398; 0.4197 ] 
#>   eta1 <~ y12      0.3873      0.0237   16.3270    0.0000 [ 0.3488; 0.4401 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8755    0.0000 [ 0.4133; 0.5009 ] 
#>   eta2 <~ y21      0.3058      0.0241   12.7104    0.0000 [ 0.2680; 0.3499 ] 
#>   eta2 <~ y22      0.4473      0.0190   23.5467    0.0000 [ 0.4210; 0.4860 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1666    0.0000 [ 0.4350; 0.5065 ] 
#>   eta3 <~ y31      0.4400      0.0187   23.5452    0.0000 [ 0.4010; 0.4662 ] 
#>   eta3 <~ y32      0.3521      0.0217   16.2461    0.0000 [ 0.3184; 0.3900 ] 
#>   eta3 <~ y33      0.3999      0.0170   23.5389    0.0000 [ 0.3670; 0.4263 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0453   14.8317    0.0000 [ 0.5937; 0.7488 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5042    0.0000 [ 0.5902; 0.7110 ] 
#>   eta3 ~ eta2       0.3052      0.0842    3.6258    0.0003 [ 0.1814; 0.4638 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0574    3.5701    0.0004 [ 0.1191; 0.3263 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04214318 15.73374  8.880411e-56
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04455941 14.57106  4.291715e-48
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03721599 20.45749  5.152411e-93
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04731329 10.91564  9.704326e-28
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03318056 22.76597 9.970534e-115
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03860474 20.71413  2.583147e-95
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03341766 24.60607 1.087643e-133
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04351999 15.12107  1.176135e-51
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03648820 20.48399  2.990959e-93
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5610224          0.7142424
#> 2          0.5477878          0.7121708
#> 3          0.6919330          0.8180943
#> 4          0.4442158          0.5812043
#> 5          0.7171327          0.8338434
#> 6          0.7099743          0.8562025
#> 7          0.7595001          0.8764086
#> 8          0.5766284          0.7340201
#> 9          0.6589659          0.7985669

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
#>  Random seed                        = -2072701448
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
#>   eta2 ~ eta1      0.6713      0.0453   14.8317    0.0000 [ 0.5533; 0.7874 ] 
#>   eta3 ~ eta1      0.4585      0.0774    5.9226    0.0000 [ 0.2758; 0.6761 ] 
#>   eta3 ~ eta2      0.3052      0.0842    3.6258    0.0003 [ 0.0714; 0.5066 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0421   15.7337    0.0000 [ 0.5670; 0.7849 ] 
#>   eta1 =~ y12      0.6493      0.0446   14.5711    0.0000 [ 0.5298; 0.7602 ] 
#>   eta1 =~ y13      0.7613      0.0372   20.4575    0.0000 [ 0.6655; 0.8579 ] 
#>   eta2 =~ y21      0.5165      0.0473   10.9156    0.0000 [ 0.3898; 0.6345 ] 
#>   eta2 =~ y22      0.7554      0.0332   22.7660    0.0000 [ 0.6626; 0.8341 ] 
#>   eta2 =~ y23      0.7997      0.0386   20.7141    0.0000 [ 0.7077; 0.9073 ] 
#>   eta3 =~ y31      0.8223      0.0334   24.6061    0.0000 [ 0.7365; 0.9093 ] 
#>   eta3 =~ y32      0.6581      0.0435   15.1211    0.0000 [ 0.5383; 0.7634 ] 
#>   eta3 =~ y33      0.7474      0.0365   20.4840    0.0000 [ 0.6577; 0.8464 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0215   18.4085    0.0000 [ 0.3455; 0.4567 ] 
#>   eta1 <~ y12      0.3873      0.0237   16.3270    0.0000 [ 0.3213; 0.4439 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8755    0.0000 [ 0.3954; 0.5079 ] 
#>   eta2 <~ y21      0.3058      0.0241   12.7104    0.0000 [ 0.2421; 0.3666 ] 
#>   eta2 <~ y22      0.4473      0.0190   23.5467    0.0000 [ 0.3951; 0.4934 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1666    0.0000 [ 0.4265; 0.5322 ] 
#>   eta3 <~ y31      0.4400      0.0187   23.5452    0.0000 [ 0.3926; 0.4892 ] 
#>   eta3 <~ y32      0.3521      0.0217   16.2461    0.0000 [ 0.2928; 0.4049 ] 
#>   eta3 <~ y33      0.3999      0.0170   23.5389    0.0000 [ 0.3592; 0.4470 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0453   14.8317    0.0000 [ 0.5533; 0.7874 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5042    0.0000 [ 0.5774; 0.7627 ] 
#>   eta3 ~ eta2       0.3052      0.0842    3.6258    0.0003 [ 0.0714; 0.5066 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0574    3.5701    0.0004 [ 0.0457; 0.3425 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04526349 14.831676 9.142099e-50
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07741616  5.922623 3.168464e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08416128  3.625790 2.880795e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55333295          0.7874102          0.5814411          0.7593021
#> 2         0.27578537          0.6761381          0.3238599          0.6280635
#> 3         0.07137717          0.5066118          0.1236403          0.4543486
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5810917          0.7552252          0.5936763          0.7487916
#> 2          0.2745248          0.6238684          0.2806473          0.5910573
#> 3          0.1044827          0.4823992          0.1813978          0.4638045
```
