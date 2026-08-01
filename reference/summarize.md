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
#>  Random seed                        = 30236866
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
#>   eta2 ~ eta1      0.6713      0.0390   17.2247    0.0000 [ 0.5937; 0.7198 ] 
#>   eta3 ~ eta1      0.4585      0.0729    6.2880    0.0000 [ 0.3468; 0.5939 ] 
#>   eta3 ~ eta2      0.3052      0.0771    3.9576    0.0001 [ 0.1475; 0.4271 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0396   16.7237    0.0000 [ 0.6110; 0.7660 ] 
#>   eta1 =~ y12      0.6493      0.0366   17.7383    0.0000 [ 0.5744; 0.7078 ] 
#>   eta1 =~ y13      0.7613      0.0315   24.1727    0.0000 [ 0.7237; 0.8207 ] 
#>   eta2 =~ y21      0.5165      0.0500   10.3301    0.0000 [ 0.4287; 0.6112 ] 
#>   eta2 =~ y22      0.7554      0.0331   22.7919    0.0000 [ 0.7094; 0.8278 ] 
#>   eta2 =~ y23      0.7997      0.0262   30.5499    0.0000 [ 0.7664; 0.8577 ] 
#>   eta3 =~ y31      0.8223      0.0294   27.9657    0.0000 [ 0.7671; 0.8762 ] 
#>   eta3 =~ y32      0.6581      0.0380   17.3387    0.0000 [ 0.5994; 0.7198 ] 
#>   eta3 =~ y33      0.7474      0.0424   17.6129    0.0000 [ 0.6641; 0.8110 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0212   18.6353    0.0000 [ 0.3612; 0.4330 ] 
#>   eta1 <~ y12      0.3873      0.0179   21.6158    0.0000 [ 0.3478; 0.4082 ] 
#>   eta1 <~ y13      0.4542      0.0192   23.6524    0.0000 [ 0.4190; 0.4841 ] 
#>   eta2 <~ y21      0.3058      0.0217   14.0701    0.0000 [ 0.2533; 0.3490 ] 
#>   eta2 <~ y22      0.4473      0.0183   24.3996    0.0000 [ 0.4155; 0.4808 ] 
#>   eta2 <~ y23      0.4735      0.0195   24.2852    0.0000 [ 0.4399; 0.5051 ] 
#>   eta3 <~ y31      0.4400      0.0190   23.1173    0.0000 [ 0.4136; 0.4746 ] 
#>   eta3 <~ y32      0.3521      0.0171   20.6274    0.0000 [ 0.3223; 0.3863 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.8999    0.0000 [ 0.3567; 0.4299 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0390   17.2247    0.0000 [ 0.5937; 0.7198 ] 
#>   eta3 ~ eta1       0.6634      0.0372   17.8510    0.0000 [ 0.5784; 0.7133 ] 
#>   eta3 ~ eta2       0.3052      0.0771    3.9576    0.0001 [ 0.1475; 0.4271 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0523    3.9195    0.0001 [ 0.0935; 0.2959 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03964859 16.72367  8.812472e-63
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03660316 17.73830  2.122569e-70
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03149615 24.17266 4.314590e-129
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04999492 10.33015  5.148224e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03314283 22.79189 5.518333e-115
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02617566 30.54990 5.671577e-205
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02940309 27.96568 4.250449e-172
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03795366 17.33875  2.399094e-67
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04243623 17.61288  1.962152e-69
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6109638          0.7660181
#> 2          0.5743664          0.7077677
#> 3          0.7237120          0.8207160
#> 4          0.4287077          0.6112214
#> 5          0.7094361          0.8278288
#> 6          0.7664154          0.8576810
#> 7          0.7671133          0.8762200
#> 8          0.5994239          0.7197658
#> 9          0.6641111          0.8110286

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
#>  Random seed                        = 30236866
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
#>   eta2 ~ eta1      0.6713      0.0390   17.2247    0.0000 [ 0.5875; 0.7891 ] 
#>   eta3 ~ eta1      0.4585      0.0729    6.2880    0.0000 [ 0.2679; 0.6450 ] 
#>   eta3 ~ eta2      0.3052      0.0771    3.9576    0.0001 [ 0.1112; 0.5099 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0396   16.7237    0.0000 [ 0.5501; 0.7551 ] 
#>   eta1 =~ y12      0.6493      0.0366   17.7383    0.0000 [ 0.5569; 0.7462 ] 
#>   eta1 =~ y13      0.7613      0.0315   24.1727    0.0000 [ 0.6735; 0.8364 ] 
#>   eta2 =~ y21      0.5165      0.0500   10.3301    0.0000 [ 0.3919; 0.6505 ] 
#>   eta2 =~ y22      0.7554      0.0331   22.7919    0.0000 [ 0.6617; 0.8331 ] 
#>   eta2 =~ y23      0.7997      0.0262   30.5499    0.0000 [ 0.7237; 0.8591 ] 
#>   eta3 =~ y31      0.8223      0.0294   27.9657    0.0000 [ 0.7416; 0.8936 ] 
#>   eta3 =~ y32      0.6581      0.0380   17.3387    0.0000 [ 0.5580; 0.7543 ] 
#>   eta3 =~ y33      0.7474      0.0424   17.6129    0.0000 [ 0.6454; 0.8649 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0212   18.6353    0.0000 [ 0.3388; 0.4486 ] 
#>   eta1 <~ y12      0.3873      0.0179   21.6158    0.0000 [ 0.3467; 0.4394 ] 
#>   eta1 <~ y13      0.4542      0.0192   23.6524    0.0000 [ 0.4056; 0.5049 ] 
#>   eta2 <~ y21      0.3058      0.0217   14.0701    0.0000 [ 0.2560; 0.3684 ] 
#>   eta2 <~ y22      0.4473      0.0183   24.3996    0.0000 [ 0.3995; 0.4943 ] 
#>   eta2 <~ y23      0.4735      0.0195   24.2852    0.0000 [ 0.4226; 0.5234 ] 
#>   eta3 <~ y31      0.4400      0.0190   23.1173    0.0000 [ 0.3879; 0.4864 ] 
#>   eta3 <~ y32      0.3521      0.0171   20.6274    0.0000 [ 0.3069; 0.3952 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.8999    0.0000 [ 0.3521; 0.4560 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0390   17.2247    0.0000 [ 0.5875; 0.7891 ] 
#>   eta3 ~ eta1       0.6634      0.0372   17.8510    0.0000 [ 0.5740; 0.7662 ] 
#>   eta3 ~ eta2       0.3052      0.0771    3.9576    0.0001 [ 0.1112; 0.5099 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0523    3.9195    0.0001 [ 0.0785; 0.3488 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03897509 17.224680 1.733836e-66
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07291815  6.287965 3.216553e-10
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07710478  3.957616 7.570140e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5875437          0.7891009          0.6117468          0.7648978
#> 2          0.2678672          0.6449588          0.3131486          0.5996774
#> 3          0.1111976          0.5099400          0.1590788          0.4620588
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59155933          0.7307666          0.5937412          0.7198459
#> 2         0.34029074          0.6442016          0.3467919          0.5938904
#> 3         0.09640066          0.4311954          0.1474767          0.4271010
```
