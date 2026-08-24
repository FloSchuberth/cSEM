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
#>  Random seed                        = -1190939065
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
#>   eta2 ~ eta1      0.6713      0.0402   16.6996    0.0000 [ 0.6097; 0.7319 ] 
#>   eta3 ~ eta1      0.4585      0.0808    5.6718    0.0000 [ 0.3241; 0.5862 ] 
#>   eta3 ~ eta2      0.3052      0.0728    4.1897    0.0000 [ 0.1720; 0.4360 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0363   18.2471    0.0000 [ 0.6024; 0.7184 ] 
#>   eta1 =~ y12      0.6493      0.0420   15.4645    0.0000 [ 0.5657; 0.7212 ] 
#>   eta1 =~ y13      0.7613      0.0336   22.6324    0.0000 [ 0.6942; 0.8236 ] 
#>   eta2 =~ y21      0.5165      0.0592    8.7295    0.0000 [ 0.3960; 0.6092 ] 
#>   eta2 =~ y22      0.7554      0.0332   22.7664    0.0000 [ 0.6922; 0.8060 ] 
#>   eta2 =~ y23      0.7997      0.0377   21.1998    0.0000 [ 0.7459; 0.8742 ] 
#>   eta3 =~ y31      0.8223      0.0340   24.1513    0.0000 [ 0.7559; 0.8858 ] 
#>   eta3 =~ y32      0.6581      0.0417   15.7958    0.0000 [ 0.5675; 0.7283 ] 
#>   eta3 =~ y33      0.7474      0.0300   24.9294    0.0000 [ 0.7054; 0.8130 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0178   22.1805    0.0000 [ 0.3646; 0.4284 ] 
#>   eta1 <~ y12      0.3873      0.0212   18.2991    0.0000 [ 0.3574; 0.4354 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0353    0.0000 [ 0.4056; 0.4886 ] 
#>   eta2 <~ y21      0.3058      0.0270   11.3103    0.0000 [ 0.2339; 0.3461 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.5822    0.0000 [ 0.4090; 0.4912 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3489    0.0000 [ 0.4481; 0.5110 ] 
#>   eta3 <~ y31      0.4400      0.0170   25.8602    0.0000 [ 0.4105; 0.4713 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5618    0.0000 [ 0.3162; 0.3855 ] 
#>   eta3 <~ y33      0.3999      0.0156   25.6755    0.0000 [ 0.3752; 0.4326 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0402   16.6996    0.0000 [ 0.6097; 0.7319 ] 
#>   eta3 ~ eta1       0.6634      0.0476   13.9416    0.0000 [ 0.5654; 0.7668 ] 
#>   eta3 ~ eta2       0.3052      0.0728    4.1897    0.0000 [ 0.1720; 0.4360 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0470    4.3555    0.0000 [ 0.1138; 0.2869 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03633828 18.247146  2.180056e-74
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04198493 15.464546  6.019586e-54
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03363963 22.632405 2.079482e-113
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05916204  8.729496  2.558100e-18
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03317988 22.766437 9.864668e-115
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03772034 21.199799 9.590575e-100
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03404698 24.151254 7.243465e-129
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04166091 15.795834  3.323666e-56
#> 9 eta3 =~ y33  Common factor 0.7474241 0.02998163 24.929401 3.572262e-137
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6024078          0.7184051
#> 2          0.5656868          0.7212030
#> 3          0.6942345          0.8235513
#> 4          0.3959860          0.6091544
#> 5          0.6922067          0.8060195
#> 6          0.7459269          0.8741754
#> 7          0.7558985          0.8857706
#> 8          0.5675293          0.7283006
#> 9          0.7054420          0.8130414

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
#>  Random seed                        = -1190939065
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
#>   eta2 ~ eta1      0.6713      0.0402   16.6996    0.0000 [ 0.5650; 0.7729 ] 
#>   eta3 ~ eta1      0.4585      0.0808    5.6718    0.0000 [ 0.2344; 0.6525 ] 
#>   eta3 ~ eta2      0.3052      0.0728    4.1897    0.0000 [ 0.1257; 0.5023 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0363   18.2471    0.0000 [ 0.5735; 0.7614 ] 
#>   eta1 =~ y12      0.6493      0.0420   15.4645    0.0000 [ 0.5474; 0.7646 ] 
#>   eta1 =~ y13      0.7613      0.0336   22.6324    0.0000 [ 0.6775; 0.8515 ] 
#>   eta2 =~ y21      0.5165      0.0592    8.7295    0.0000 [ 0.3761; 0.6820 ] 
#>   eta2 =~ y22      0.7554      0.0332   22.7664    0.0000 [ 0.6703; 0.8419 ] 
#>   eta2 =~ y23      0.7997      0.0377   21.1998    0.0000 [ 0.6962; 0.8912 ] 
#>   eta3 =~ y31      0.8223      0.0340   24.1513    0.0000 [ 0.7324; 0.9084 ] 
#>   eta3 =~ y32      0.6581      0.0417   15.7958    0.0000 [ 0.5532; 0.7686 ] 
#>   eta3 =~ y33      0.7474      0.0300   24.9294    0.0000 [ 0.6675; 0.8225 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0178   22.1805    0.0000 [ 0.3483; 0.4405 ] 
#>   eta1 <~ y12      0.3873      0.0212   18.2991    0.0000 [ 0.3330; 0.4424 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0353    0.0000 [ 0.3982; 0.5048 ] 
#>   eta2 <~ y21      0.3058      0.0270   11.3103    0.0000 [ 0.2432; 0.3830 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.5822    0.0000 [ 0.3871; 0.5052 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3489    0.0000 [ 0.4137; 0.5233 ] 
#>   eta3 <~ y31      0.4400      0.0170   25.8602    0.0000 [ 0.3956; 0.4836 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5618    0.0000 [ 0.3052; 0.4033 ] 
#>   eta3 <~ y33      0.3999      0.0156   25.6755    0.0000 [ 0.3589; 0.4394 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0402   16.6996    0.0000 [ 0.5650; 0.7729 ] 
#>   eta3 ~ eta1       0.6634      0.0476   13.9416    0.0000 [ 0.5313; 0.7773 ] 
#>   eta3 ~ eta2       0.3052      0.0728    4.1897    0.0000 [ 0.1257; 0.5023 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0470    4.3555    0.0000 [ 0.0892; 0.3324 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04020054 16.699612 1.319205e-62
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08084025  5.671763 1.413352e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07283341  4.189713 2.793071e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5650307          0.7729252          0.5899948          0.7479612
#> 2          0.2344451          0.6525052          0.2846459          0.6023043
#> 3          0.1256804          0.5023337          0.1709091          0.4571050
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.6057970          0.7760986          0.6097018          0.7319109
#> 2          0.2788930          0.6881041          0.3241342          0.5861825
#> 3          0.1075492          0.4419076          0.1720476          0.4359973
```
