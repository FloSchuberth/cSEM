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
#>  Random seed                        = 1509378279
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
#>   eta2 ~ eta1      0.6713      0.0496   13.5255    0.0000 [ 0.5878; 0.7426 ] 
#>   eta3 ~ eta1      0.4585      0.0877    5.2268    0.0000 [ 0.3091; 0.6004 ] 
#>   eta3 ~ eta2      0.3052      0.1003    3.0435    0.0023 [ 0.1412; 0.5057 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0320   20.6944    0.0000 [ 0.5971; 0.7213 ] 
#>   eta1 =~ y12      0.6493      0.0429   15.1177    0.0000 [ 0.5667; 0.7213 ] 
#>   eta1 =~ y13      0.7613      0.0363   20.9513    0.0000 [ 0.6996; 0.8329 ] 
#>   eta2 =~ y21      0.5165      0.0593    8.7088    0.0000 [ 0.3827; 0.5800 ] 
#>   eta2 =~ y22      0.7554      0.0366   20.6129    0.0000 [ 0.6655; 0.8077 ] 
#>   eta2 =~ y23      0.7997      0.0379   21.0769    0.0000 [ 0.7182; 0.8691 ] 
#>   eta3 =~ y31      0.8223      0.0239   34.4362    0.0000 [ 0.7717; 0.8512 ] 
#>   eta3 =~ y32      0.6581      0.0342   19.2564    0.0000 [ 0.6072; 0.7391 ] 
#>   eta3 =~ y33      0.7474      0.0349   21.4410    0.0000 [ 0.6810; 0.8198 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0206   19.1993    0.0000 [ 0.3523; 0.4265 ] 
#>   eta1 <~ y12      0.3873      0.0172   22.4933    0.0000 [ 0.3586; 0.4188 ] 
#>   eta1 <~ y13      0.4542      0.0204   22.2945    0.0000 [ 0.4241; 0.5007 ] 
#>   eta2 <~ y21      0.3058      0.0340    8.9870    0.0000 [ 0.2427; 0.3492 ] 
#>   eta2 <~ y22      0.4473      0.0224   19.9646    0.0000 [ 0.4122; 0.4914 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7844    0.0000 [ 0.4436; 0.5126 ] 
#>   eta3 <~ y31      0.4400      0.0143   30.7860    0.0000 [ 0.4088; 0.4599 ] 
#>   eta3 <~ y32      0.3521      0.0159   22.1220    0.0000 [ 0.3313; 0.3875 ] 
#>   eta3 <~ y33      0.3999      0.0163   24.5351    0.0000 [ 0.3633; 0.4232 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0496   13.5255    0.0000 [ 0.5878; 0.7426 ] 
#>   eta3 ~ eta1       0.6634      0.0374   17.7449    0.0000 [ 0.5994; 0.7327 ] 
#>   eta3 ~ eta2       0.3052      0.1003    3.0435    0.0023 [ 0.1412; 0.5057 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0680    3.0106    0.0026 [ 0.0943; 0.3159 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03204103 20.694401  3.890358e-95
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04294833 15.117652  1.238814e-51
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03633885 20.951288  1.826243e-97
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05930248  8.708824  3.070449e-18
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03664643 20.612857  2.104383e-94
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03794032 21.076886  1.296373e-98
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02387827 34.436218 7.241542e-260
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03417402 19.256411  1.247294e-82
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03485953 21.441024 5.538123e-102
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5970964          0.7213366
#> 2          0.5666581          0.7213452
#> 3          0.6995503          0.8328663
#> 4          0.3827471          0.5799981
#> 5          0.6654966          0.8076599
#> 6          0.7182153          0.8690920
#> 7          0.7717494          0.8511691
#> 8          0.6071819          0.7391393
#> 9          0.6810082          0.8198265

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
#>  Random seed                        = 1509378279
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
#>   eta2 ~ eta1      0.6713      0.0496   13.5255    0.0000 [ 0.5379; 0.7946 ] 
#>   eta3 ~ eta1      0.4585      0.0877    5.2268    0.0000 [ 0.2296; 0.6832 ] 
#>   eta3 ~ eta2      0.3052      0.1003    3.0435    0.0023 [ 0.0422; 0.5607 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0320   20.6944    0.0000 [ 0.5820; 0.7477 ] 
#>   eta1 =~ y12      0.6493      0.0429   15.1177    0.0000 [ 0.5380; 0.7601 ] 
#>   eta1 =~ y13      0.7613      0.0363   20.9513    0.0000 [ 0.6649; 0.8528 ] 
#>   eta2 =~ y21      0.5165      0.0593    8.7088    0.0000 [ 0.3801; 0.6868 ] 
#>   eta2 =~ y22      0.7554      0.0366   20.6129    0.0000 [ 0.6656; 0.8551 ] 
#>   eta2 =~ y23      0.7997      0.0379   21.0769    0.0000 [ 0.7119; 0.9081 ] 
#>   eta3 =~ y31      0.8223      0.0239   34.4362    0.0000 [ 0.7672; 0.8906 ] 
#>   eta3 =~ y32      0.6581      0.0342   19.2564    0.0000 [ 0.5552; 0.7319 ] 
#>   eta3 =~ y33      0.7474      0.0349   21.4410    0.0000 [ 0.6596; 0.8399 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0206   19.1993    0.0000 [ 0.3436; 0.4501 ] 
#>   eta1 <~ y12      0.3873      0.0172   22.4933    0.0000 [ 0.3434; 0.4325 ] 
#>   eta1 <~ y13      0.4542      0.0204   22.2945    0.0000 [ 0.4004; 0.5058 ] 
#>   eta2 <~ y21      0.3058      0.0340    8.9870    0.0000 [ 0.2220; 0.3980 ] 
#>   eta2 <~ y22      0.4473      0.0224   19.9646    0.0000 [ 0.3830; 0.4988 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7844    0.0000 [ 0.4137; 0.5261 ] 
#>   eta3 <~ y31      0.4400      0.0143   30.7860    0.0000 [ 0.4080; 0.4819 ] 
#>   eta3 <~ y32      0.3521      0.0159   22.1220    0.0000 [ 0.3046; 0.3869 ] 
#>   eta3 <~ y33      0.3999      0.0163   24.5351    0.0000 [ 0.3605; 0.4448 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0496   13.5255    0.0000 [ 0.5379; 0.7946 ] 
#>   eta3 ~ eta1       0.6634      0.0374   17.7449    0.0000 [ 0.5610; 0.7544 ] 
#>   eta3 ~ eta2       0.3052      0.1003    3.0435    0.0023 [ 0.0422; 0.5607 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0680    3.0106    0.0026 [ 0.0254; 0.3773 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04963466 13.525496 1.105954e-41
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08772161  5.226840 1.724320e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10026251  3.043522 2.338265e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.53791871          0.7946012          0.5687413          0.7637786
#> 2         0.22957483          0.6832215          0.2840489          0.6287474
#> 3         0.04216146          0.5606625          0.1044233          0.4984007
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5736388          0.8273269          0.5877785          0.7425994
#> 2          0.2377872          0.6295372          0.3090826          0.6004447
#> 3          0.1136520          0.5360343          0.1411568          0.5056832
```
