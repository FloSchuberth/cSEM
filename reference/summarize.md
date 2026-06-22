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
#>  Random seed                        = -71858927
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
#>   eta2 ~ eta1      0.6713      0.0459   14.6265    0.0000 [ 0.6007; 0.7457 ] 
#>   eta3 ~ eta1      0.4585      0.0829    5.5316    0.0000 [ 0.3499; 0.6360 ] 
#>   eta3 ~ eta2      0.3052      0.0919    3.3187    0.0009 [ 0.1531; 0.4571 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0443   14.9771    0.0000 [ 0.5761; 0.7240 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9633    0.0000 [ 0.5889; 0.7082 ] 
#>   eta1 =~ y13      0.7613      0.0318   23.9237    0.0000 [ 0.7095; 0.8263 ] 
#>   eta2 =~ y21      0.5165      0.0465   11.0979    0.0000 [ 0.4253; 0.5964 ] 
#>   eta2 =~ y22      0.7554      0.0372   20.3264    0.0000 [ 0.6856; 0.8197 ] 
#>   eta2 =~ y23      0.7997      0.0324   24.6967    0.0000 [ 0.7466; 0.8530 ] 
#>   eta3 =~ y31      0.8223      0.0364   22.5828    0.0000 [ 0.7747; 0.8978 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.7037    0.0000 [ 0.5772; 0.7051 ] 
#>   eta3 =~ y33      0.7474      0.0430   17.3749    0.0000 [ 0.6809; 0.8241 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0248   15.9613    0.0000 [ 0.3423; 0.4301 ] 
#>   eta1 <~ y12      0.3873      0.0200   19.3373    0.0000 [ 0.3551; 0.4300 ] 
#>   eta1 <~ y13      0.4542      0.0183   24.7664    0.0000 [ 0.4196; 0.4894 ] 
#>   eta2 <~ y21      0.3058      0.0213   14.3275    0.0000 [ 0.2571; 0.3401 ] 
#>   eta2 <~ y22      0.4473      0.0219   20.4219    0.0000 [ 0.4042; 0.4872 ] 
#>   eta2 <~ y23      0.4735      0.0207   22.9299    0.0000 [ 0.4342; 0.5067 ] 
#>   eta3 <~ y31      0.4400      0.0209   21.0619    0.0000 [ 0.4084; 0.4708 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9919    0.0000 [ 0.3071; 0.3815 ] 
#>   eta3 <~ y33      0.3999      0.0207   19.3197    0.0000 [ 0.3679; 0.4411 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0459   14.6265    0.0000 [ 0.6007; 0.7457 ] 
#>   eta3 ~ eta1       0.6634      0.0364   18.2282    0.0000 [ 0.6203; 0.7540 ] 
#>   eta3 ~ eta2       0.3052      0.0919    3.3187    0.0009 [ 0.1531; 0.4571 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0598    3.4257    0.0006 [ 0.1154; 0.3262 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04427222 14.97711  1.036257e-50
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03827539 16.96333  1.534106e-64
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03182392 23.92370 1.736119e-126
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04653612 11.09793  1.283763e-28
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03716296 20.32636  7.517496e-92
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03237942 24.69667 1.161392e-134
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03641166 22.58281 6.395683e-113
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04190529 15.70372  1.426274e-55
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04301751 17.37488  1.278729e-67
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5761016          0.7240195
#> 2          0.5889483          0.7082479
#> 3          0.7094598          0.8262884
#> 4          0.4253332          0.5964346
#> 5          0.6856454          0.8197394
#> 6          0.7465890          0.8530490
#> 7          0.7747312          0.8977658
#> 8          0.5771683          0.7050999
#> 9          0.6809308          0.8241147

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
#>  Random seed                        = -71858927
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
#>   eta2 ~ eta1      0.6713      0.0459   14.6265    0.0000 [ 0.5468; 0.7841 ] 
#>   eta3 ~ eta1      0.4585      0.0829    5.5316    0.0000 [ 0.2400; 0.6686 ] 
#>   eta3 ~ eta2      0.3052      0.0919    3.3187    0.0009 [ 0.0618; 0.5373 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0443   14.9771    0.0000 [ 0.5607; 0.7897 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9633    0.0000 [ 0.5481; 0.7461 ] 
#>   eta1 =~ y13      0.7613      0.0318   23.9237    0.0000 [ 0.6794; 0.8440 ] 
#>   eta2 =~ y21      0.5165      0.0465   11.0979    0.0000 [ 0.3991; 0.6398 ] 
#>   eta2 =~ y22      0.7554      0.0372   20.3264    0.0000 [ 0.6546; 0.8468 ] 
#>   eta2 =~ y23      0.7997      0.0324   24.6967    0.0000 [ 0.7160; 0.8834 ] 
#>   eta3 =~ y31      0.8223      0.0364   22.5828    0.0000 [ 0.7321; 0.9204 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.7037    0.0000 [ 0.5586; 0.7753 ] 
#>   eta3 =~ y33      0.7474      0.0430   17.3749    0.0000 [ 0.6294; 0.8519 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0248   15.9613    0.0000 [ 0.3361; 0.4643 ] 
#>   eta1 <~ y12      0.3873      0.0200   19.3373    0.0000 [ 0.3317; 0.4353 ] 
#>   eta1 <~ y13      0.4542      0.0183   24.7664    0.0000 [ 0.4038; 0.4987 ] 
#>   eta2 <~ y21      0.3058      0.0213   14.3275    0.0000 [ 0.2535; 0.3639 ] 
#>   eta2 <~ y22      0.4473      0.0219   20.4219    0.0000 [ 0.3888; 0.5021 ] 
#>   eta2 <~ y23      0.4735      0.0207   22.9299    0.0000 [ 0.4210; 0.5278 ] 
#>   eta3 <~ y31      0.4400      0.0209   21.0619    0.0000 [ 0.3865; 0.4945 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9919    0.0000 [ 0.3052; 0.4064 ] 
#>   eta3 <~ y33      0.3999      0.0207   19.3197    0.0000 [ 0.3415; 0.4485 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0459   14.6265    0.0000 [ 0.5468; 0.7841 ] 
#>   eta3 ~ eta1       0.6634      0.0364   18.2282    0.0000 [ 0.5610; 0.7492 ] 
#>   eta3 ~ eta2       0.3052      0.0919    3.3187    0.0009 [ 0.0618; 0.5373 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0598    3.4257    0.0006 [ 0.0461; 0.3554 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04589845 14.626494 1.903458e-48
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08288880  5.531589 3.173427e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09194981  3.318671 9.044708e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5467867          0.7841476          0.5752891          0.7556452
#> 2          0.2399766          0.6686306          0.2914496          0.6171577
#> 3          0.0617785          0.5372910          0.1188783          0.4801912
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5743785          0.7601829          0.6007121          0.7457458
#> 2          0.2976290          0.6553749          0.3498890          0.6360334
#> 3          0.1121695          0.5122804          0.1530648          0.4570622
```
