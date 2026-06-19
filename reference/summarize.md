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
#>  Random seed                        = -807378051
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
#>   eta2 ~ eta1      0.6713      0.0400   16.7794    0.0000 [ 0.5787; 0.7236 ] 
#>   eta3 ~ eta1      0.4585      0.0938    4.8862    0.0000 [ 0.2258; 0.5818 ] 
#>   eta3 ~ eta2      0.3052      0.0922    3.3114    0.0009 [ 0.1831; 0.4971 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0432   15.3490    0.0000 [ 0.5602; 0.7339 ] 
#>   eta1 =~ y12      0.6493      0.0374   17.3660    0.0000 [ 0.5679; 0.7031 ] 
#>   eta1 =~ y13      0.7613      0.0351   21.7076    0.0000 [ 0.6928; 0.8213 ] 
#>   eta2 =~ y21      0.5165      0.0570    9.0618    0.0000 [ 0.3953; 0.5999 ] 
#>   eta2 =~ y22      0.7554      0.0352   21.4547    0.0000 [ 0.6994; 0.8213 ] 
#>   eta2 =~ y23      0.7997      0.0392   20.4081    0.0000 [ 0.7358; 0.8614 ] 
#>   eta3 =~ y31      0.8223      0.0361   22.7532    0.0000 [ 0.7642; 0.8935 ] 
#>   eta3 =~ y32      0.6581      0.0431   15.2623    0.0000 [ 0.5746; 0.7240 ] 
#>   eta3 =~ y33      0.7474      0.0426   17.5462    0.0000 [ 0.6651; 0.8168 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0190   20.7976    0.0000 [ 0.3631; 0.4269 ] 
#>   eta1 <~ y12      0.3873      0.0211   18.3145    0.0000 [ 0.3525; 0.4193 ] 
#>   eta1 <~ y13      0.4542      0.0221   20.5472    0.0000 [ 0.4224; 0.5012 ] 
#>   eta2 <~ y21      0.3058      0.0305   10.0198    0.0000 [ 0.2331; 0.3482 ] 
#>   eta2 <~ y22      0.4473      0.0216   20.7494    0.0000 [ 0.4193; 0.4981 ] 
#>   eta2 <~ y23      0.4735      0.0230   20.5628    0.0000 [ 0.4455; 0.5157 ] 
#>   eta3 <~ y31      0.4400      0.0211   20.8538    0.0000 [ 0.4183; 0.4851 ] 
#>   eta3 <~ y32      0.3521      0.0173   20.3609    0.0000 [ 0.3271; 0.3794 ] 
#>   eta3 <~ y33      0.3999      0.0233   17.1321    0.0000 [ 0.3605; 0.4499 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0400   16.7794    0.0000 [ 0.5787; 0.7236 ] 
#>   eta3 ~ eta1       0.6634      0.0469   14.1570    0.0000 [ 0.5737; 0.7098 ] 
#>   eta3 ~ eta2       0.3052      0.0922    3.3114    0.0009 [ 0.1831; 0.4971 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0626    3.2706    0.0011 [ 0.1147; 0.3501 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04319956 15.348996  3.597079e-53
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03738783 17.366024  1.492131e-67
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03507284 21.707558 1.740638e-104
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05699279  9.061757  1.283654e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03520843 21.454739 4.124224e-102
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03918370 20.408070  1.417744e-92
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03613896 22.753212 1.333683e-114
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04311722 15.262323  1.363070e-52
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04259738 17.546246  6.354330e-69
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5601910          0.7338518
#> 2          0.5679143          0.7030727
#> 3          0.6927774          0.8212662
#> 4          0.3952897          0.5998747
#> 5          0.6994431          0.8212748
#> 6          0.7358037          0.8614442
#> 7          0.7642171          0.8935082
#> 8          0.5746105          0.7240338
#> 9          0.6650835          0.8167918

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
#>  Random seed                        = -807378051
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
#>   eta2 ~ eta1      0.6713      0.0400   16.7794    0.0000 [ 0.5805; 0.7875 ] 
#>   eta3 ~ eta1      0.4585      0.0938    4.8862    0.0000 [ 0.2367; 0.7220 ] 
#>   eta3 ~ eta2      0.3052      0.0922    3.3114    0.0009 [ 0.0536; 0.5302 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0432   15.3490    0.0000 [ 0.5488; 0.7722 ] 
#>   eta1 =~ y12      0.6493      0.0374   17.3660    0.0000 [ 0.5625; 0.7558 ] 
#>   eta1 =~ y13      0.7613      0.0351   21.7076    0.0000 [ 0.6764; 0.8577 ] 
#>   eta2 =~ y21      0.5165      0.0570    9.0618    0.0000 [ 0.3875; 0.6823 ] 
#>   eta2 =~ y22      0.7554      0.0352   21.4547    0.0000 [ 0.6589; 0.8409 ] 
#>   eta2 =~ y23      0.7997      0.0392   20.4081    0.0000 [ 0.6991; 0.9017 ] 
#>   eta3 =~ y31      0.8223      0.0361   22.7532    0.0000 [ 0.7250; 0.9119 ] 
#>   eta3 =~ y32      0.6581      0.0431   15.2623    0.0000 [ 0.5512; 0.7742 ] 
#>   eta3 =~ y33      0.7474      0.0426   17.5462    0.0000 [ 0.6431; 0.8634 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0190   20.7976    0.0000 [ 0.3416; 0.4400 ] 
#>   eta1 <~ y12      0.3873      0.0211   18.3145    0.0000 [ 0.3351; 0.4445 ] 
#>   eta1 <~ y13      0.4542      0.0221   20.5472    0.0000 [ 0.3962; 0.5105 ] 
#>   eta2 <~ y21      0.3058      0.0305   10.0198    0.0000 [ 0.2362; 0.3940 ] 
#>   eta2 <~ y22      0.4473      0.0216   20.7494    0.0000 [ 0.3852; 0.4967 ] 
#>   eta2 <~ y23      0.4735      0.0230   20.5628    0.0000 [ 0.4112; 0.5303 ] 
#>   eta3 <~ y31      0.4400      0.0211   20.8538    0.0000 [ 0.3815; 0.4906 ] 
#>   eta3 <~ y32      0.3521      0.0173   20.3609    0.0000 [ 0.3088; 0.3983 ] 
#>   eta3 <~ y33      0.3999      0.0233   17.1321    0.0000 [ 0.3411; 0.4618 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0400   16.7794    0.0000 [ 0.5805; 0.7875 ] 
#>   eta3 ~ eta1       0.6634      0.0469   14.1570    0.0000 [ 0.5581; 0.8004 ] 
#>   eta3 ~ eta2       0.3052      0.0922    3.3114    0.0009 [ 0.0536; 0.5302 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0626    3.2706    0.0011 [ 0.0380; 0.3619 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04000941 16.779389 3.453700e-63
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09383762  4.886172 1.028152e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09215175  3.311398 9.283109e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.58054759          0.7874537          0.6053930          0.7626083
#> 2         0.23671904          0.7219942          0.2949911          0.6637222
#> 3         0.05359972          0.5301565          0.1108249          0.4729314
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5703566          0.7236612          0.5787391          0.7235854
#> 2          0.2246324          0.6050709          0.2258325          0.5817683
#> 3          0.1248904          0.5077666          0.1830599          0.4971393
```
