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
#>  Random seed                        = 1271966217
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
#>   eta2 ~ eta1      0.6713      0.0498   13.4782    0.0000 [ 0.5772; 0.7487 ] 
#>   eta3 ~ eta1      0.4585      0.0704    6.5161    0.0000 [ 0.3287; 0.5954 ] 
#>   eta3 ~ eta2      0.3052      0.0733    4.1644    0.0000 [ 0.1565; 0.4172 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0433   15.3225    0.0000 [ 0.5694; 0.7170 ] 
#>   eta1 =~ y12      0.6493      0.0371   17.5091    0.0000 [ 0.5941; 0.7247 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.2032    0.0000 [ 0.7104; 0.8453 ] 
#>   eta2 =~ y21      0.5165      0.0659    7.8407    0.0000 [ 0.3770; 0.6148 ] 
#>   eta2 =~ y22      0.7554      0.0386   19.5616    0.0000 [ 0.6715; 0.8145 ] 
#>   eta2 =~ y23      0.7997      0.0382   20.9544    0.0000 [ 0.7152; 0.8607 ] 
#>   eta3 =~ y31      0.8223      0.0327   25.1505    0.0000 [ 0.7490; 0.8759 ] 
#>   eta3 =~ y32      0.6581      0.0383   17.1777    0.0000 [ 0.5694; 0.7226 ] 
#>   eta3 =~ y33      0.7474      0.0450   16.6083    0.0000 [ 0.6425; 0.8251 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0220   17.9970    0.0000 [ 0.3372; 0.4156 ] 
#>   eta1 <~ y12      0.3873      0.0169   22.9680    0.0000 [ 0.3587; 0.4224 ] 
#>   eta1 <~ y13      0.4542      0.0236   19.2301    0.0000 [ 0.4261; 0.5069 ] 
#>   eta2 <~ y21      0.3058      0.0327    9.3651    0.0000 [ 0.2367; 0.3597 ] 
#>   eta2 <~ y22      0.4473      0.0282   15.8657    0.0000 [ 0.4058; 0.5080 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1847    0.0000 [ 0.4427; 0.5152 ] 
#>   eta3 <~ y31      0.4400      0.0199   22.0872    0.0000 [ 0.3969; 0.4738 ] 
#>   eta3 <~ y32      0.3521      0.0205   17.1912    0.0000 [ 0.3107; 0.3864 ] 
#>   eta3 <~ y33      0.3999      0.0191   20.9206    0.0000 [ 0.3642; 0.4285 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0498   13.4782    0.0000 [ 0.5772; 0.7487 ] 
#>   eta3 ~ eta1       0.6634      0.0444   14.9396    0.0000 [ 0.5838; 0.7308 ] 
#>   eta3 ~ eta2       0.3052      0.0733    4.1644    0.0000 [ 0.1565; 0.4172 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0516    3.9669    0.0001 [ 0.0999; 0.3001 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04327430 15.322487  5.410669e-53
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03708232 17.509096  1.221184e-68
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03590706 21.203235 8.915238e-100
#> 4 eta2 =~ y21  Common factor 0.5164548 0.06586808  7.840744  4.478835e-15
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03861585 19.561598  3.286149e-85
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03816204 20.954430  1.709643e-97
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03269433 25.150459 1.396939e-139
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03830946 17.177713  3.899847e-66
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04500295 16.608336  6.065231e-62
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5693617          0.7170226
#> 2          0.5940922          0.7246737
#> 3          0.7104231          0.8452851
#> 4          0.3769803          0.6147833
#> 5          0.6714926          0.8144782
#> 6          0.7151836          0.8606911
#> 7          0.7490346          0.8759361
#> 8          0.5693653          0.7226447
#> 9          0.6425066          0.8250668

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
#>  Random seed                        = 1271966217
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
#>   eta2 ~ eta1      0.6713      0.0498   13.4782    0.0000 [ 0.5448; 0.8023 ] 
#>   eta3 ~ eta1      0.4585      0.0704    6.5161    0.0000 [ 0.2899; 0.6537 ] 
#>   eta3 ~ eta2      0.3052      0.0733    4.1644    0.0000 [ 0.1012; 0.4801 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0433   15.3225    0.0000 [ 0.5591; 0.7829 ] 
#>   eta1 =~ y12      0.6493      0.0371   17.5091    0.0000 [ 0.5477; 0.7394 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.2032    0.0000 [ 0.6574; 0.8431 ] 
#>   eta2 =~ y21      0.5165      0.0659    7.8407    0.0000 [ 0.3436; 0.6842 ] 
#>   eta2 =~ y22      0.7554      0.0386   19.5616    0.0000 [ 0.6569; 0.8566 ] 
#>   eta2 =~ y23      0.7997      0.0382   20.9544    0.0000 [ 0.7000; 0.8974 ] 
#>   eta3 =~ y31      0.8223      0.0327   25.1505    0.0000 [ 0.7411; 0.9102 ] 
#>   eta3 =~ y32      0.6581      0.0383   17.1777    0.0000 [ 0.5576; 0.7557 ] 
#>   eta3 =~ y33      0.7474      0.0450   16.6083    0.0000 [ 0.6334; 0.8661 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0220   17.9970    0.0000 [ 0.3464; 0.4600 ] 
#>   eta1 <~ y12      0.3873      0.0169   22.9680    0.0000 [ 0.3432; 0.4304 ] 
#>   eta1 <~ y13      0.4542      0.0236   19.2301    0.0000 [ 0.3895; 0.5117 ] 
#>   eta2 <~ y21      0.3058      0.0327    9.3651    0.0000 [ 0.2214; 0.3903 ] 
#>   eta2 <~ y22      0.4473      0.0282   15.8657    0.0000 [ 0.3760; 0.5218 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1847    0.0000 [ 0.4214; 0.5270 ] 
#>   eta3 <~ y31      0.4400      0.0199   22.0872    0.0000 [ 0.3890; 0.4920 ] 
#>   eta3 <~ y32      0.3521      0.0205   17.1912    0.0000 [ 0.2975; 0.4034 ] 
#>   eta3 <~ y33      0.3999      0.0191   20.9206    0.0000 [ 0.3509; 0.4498 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0498   13.4782    0.0000 [ 0.5448; 0.8023 ] 
#>   eta3 ~ eta1       0.6634      0.0444   14.9396    0.0000 [ 0.5531; 0.7828 ] 
#>   eta3 ~ eta2       0.3052      0.0733    4.1644    0.0000 [ 0.1012; 0.4801 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0516    3.9669    0.0001 [ 0.0626; 0.3297 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04980879 13.478213 2.101345e-41
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07036550  6.516073 7.217173e-11
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07327684  4.164360 3.122273e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5447595          0.8023425          0.5756902          0.7714118
#> 2          0.2898519          0.6537425          0.3335480          0.6100463
#> 3          0.1011828          0.4801293          0.1466869          0.4346252
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5755270          0.7787497          0.5772022          0.7487006
#> 2          0.3185748          0.6127485          0.3287109          0.5953677
#> 3          0.1285221          0.4616874          0.1564957          0.4172136
```
