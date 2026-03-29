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
#>  Random seed                        = -2034189219
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
#>   eta2 ~ eta1      0.6713      0.0350   19.1925    0.0000 [ 0.6094; 0.7230 ] 
#>   eta3 ~ eta1      0.4585      0.0805    5.6956    0.0000 [ 0.3134; 0.6047 ] 
#>   eta3 ~ eta2      0.3052      0.0758    4.0256    0.0001 [ 0.1650; 0.4471 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0426   15.5626    0.0000 [ 0.5827; 0.7423 ] 
#>   eta1 =~ y12      0.6493      0.0392   16.5629    0.0000 [ 0.5874; 0.7231 ] 
#>   eta1 =~ y13      0.7613      0.0345   22.0996    0.0000 [ 0.7090; 0.8411 ] 
#>   eta2 =~ y21      0.5165      0.0539    9.5800    0.0000 [ 0.3737; 0.5972 ] 
#>   eta2 =~ y22      0.7554      0.0376   20.0901    0.0000 [ 0.6855; 0.8145 ] 
#>   eta2 =~ y23      0.7997      0.0377   21.1876    0.0000 [ 0.7262; 0.8957 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.6905    0.0000 [ 0.7628; 0.8946 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.1992    0.0000 [ 0.5839; 0.7421 ] 
#>   eta3 =~ y33      0.7474      0.0460   16.2432    0.0000 [ 0.6554; 0.8081 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0227   17.4154    0.0000 [ 0.3540; 0.4329 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9452    0.0000 [ 0.3457; 0.4248 ] 
#>   eta1 <~ y13      0.4542      0.0197   23.0511    0.0000 [ 0.4234; 0.4888 ] 
#>   eta2 <~ y21      0.3058      0.0289   10.5881    0.0000 [ 0.2331; 0.3468 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6215    0.0000 [ 0.4183; 0.4987 ] 
#>   eta2 <~ y23      0.4735      0.0215   22.0659    0.0000 [ 0.4406; 0.5135 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4934    0.0000 [ 0.4019; 0.4957 ] 
#>   eta3 <~ y32      0.3521      0.0204   17.2392    0.0000 [ 0.3102; 0.3875 ] 
#>   eta3 <~ y33      0.3999      0.0206   19.4181    0.0000 [ 0.3592; 0.4339 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0350   19.1925    0.0000 [ 0.6094; 0.7230 ] 
#>   eta3 ~ eta1       0.6634      0.0424   15.6275    0.0000 [ 0.5838; 0.7360 ] 
#>   eta3 ~ eta2       0.3052      0.0758    4.0256    0.0001 [ 0.1650; 0.4471 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0503    4.0713    0.0000 [ 0.1097; 0.2831 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04260658 15.562618  1.306446e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03920085 16.562854  1.293121e-61
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03445068 22.099585 3.189825e-108
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05390958  9.580019  9.702722e-22
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03759999 20.090104  9.007312e-90
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03774212 21.187565  1.243655e-99
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03623882 22.690511 5.558762e-114
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04329620 15.199228  3.578112e-52
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04601456 16.243209  2.495236e-59
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5827215          0.7422772
#> 2          0.5873805          0.7230812
#> 3          0.7089831          0.8411482
#> 4          0.3736780          0.5971852
#> 5          0.6855421          0.8145364
#> 6          0.7261585          0.8957335
#> 7          0.7627652          0.8946297
#> 8          0.5839323          0.7420560
#> 9          0.6553724          0.8080873

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
#>  Random seed                        = -2034189219
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
#>   eta2 ~ eta1      0.6713      0.0350   19.1925    0.0000 [ 0.5918; 0.7727 ] 
#>   eta3 ~ eta1      0.4585      0.0805    5.6956    0.0000 [ 0.2384; 0.6547 ] 
#>   eta3 ~ eta2      0.3052      0.0758    4.0256    0.0001 [ 0.1172; 0.5092 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0426   15.5626    0.0000 [ 0.5511; 0.7714 ] 
#>   eta1 =~ y12      0.6493      0.0392   16.5629    0.0000 [ 0.5433; 0.7460 ] 
#>   eta1 =~ y13      0.7613      0.0345   22.0996    0.0000 [ 0.6652; 0.8433 ] 
#>   eta2 =~ y21      0.5165      0.0539    9.5800    0.0000 [ 0.3892; 0.6679 ] 
#>   eta2 =~ y22      0.7554      0.0376   20.0901    0.0000 [ 0.6538; 0.8483 ] 
#>   eta2 =~ y23      0.7997      0.0377   21.1876    0.0000 [ 0.6977; 0.8929 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.6905    0.0000 [ 0.7226; 0.9100 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.1992    0.0000 [ 0.5395; 0.7634 ] 
#>   eta3 =~ y33      0.7474      0.0460   16.2432    0.0000 [ 0.6371; 0.8751 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0227   17.4154    0.0000 [ 0.3399; 0.4574 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9452    0.0000 [ 0.3328; 0.4444 ] 
#>   eta1 <~ y13      0.4542      0.0197   23.0511    0.0000 [ 0.4037; 0.5056 ] 
#>   eta2 <~ y21      0.3058      0.0289   10.5881    0.0000 [ 0.2386; 0.3880 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6215    0.0000 [ 0.3857; 0.5036 ] 
#>   eta2 <~ y23      0.4735      0.0215   22.0659    0.0000 [ 0.4155; 0.5264 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4934    0.0000 [ 0.3796; 0.4963 ] 
#>   eta3 <~ y32      0.3521      0.0204   17.2392    0.0000 [ 0.2970; 0.4026 ] 
#>   eta3 <~ y33      0.3999      0.0206   19.4181    0.0000 [ 0.3527; 0.4592 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0350   19.1925    0.0000 [ 0.5918; 0.7727 ] 
#>   eta3 ~ eta1       0.6634      0.0424   15.6275    0.0000 [ 0.5504; 0.7700 ] 
#>   eta3 ~ eta2       0.3052      0.0758    4.0256    0.0001 [ 0.1172; 0.5092 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0503    4.0713    0.0000 [ 0.0835; 0.3437 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03497896 19.192492 4.276298e-82
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08050164  5.695620 1.229242e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07580281  4.025591 5.683238e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5917765          0.7726679          0.6134980          0.7509464
#> 2          0.2383973          0.6547064          0.2883879          0.6047158
#> 3          0.1171815          0.5091909          0.1642542          0.4621182
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5936077          0.7404834          0.6094205          0.7230052
#> 2          0.2950803          0.6293020          0.3134424          0.6047177
#> 3          0.1514754          0.4586915          0.1649770          0.4470648
```
