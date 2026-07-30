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
#>  Random seed                        = 1418278302
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
#>   eta2 ~ eta1      0.6713      0.0486   13.8243    0.0000 [ 0.5691; 0.7730 ] 
#>   eta3 ~ eta1      0.4585      0.0668    6.8676    0.0000 [ 0.3689; 0.5910 ] 
#>   eta3 ~ eta2      0.3052      0.0701    4.3550    0.0000 [ 0.1755; 0.3937 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0505   13.1180    0.0000 [ 0.5635; 0.7339 ] 
#>   eta1 =~ y12      0.6493      0.0339   19.1651    0.0000 [ 0.5878; 0.7046 ] 
#>   eta1 =~ y13      0.7613      0.0302   25.2362    0.0000 [ 0.7096; 0.8114 ] 
#>   eta2 =~ y21      0.5165      0.0604    8.5448    0.0000 [ 0.3958; 0.5967 ] 
#>   eta2 =~ y22      0.7554      0.0330   22.8802    0.0000 [ 0.7083; 0.8236 ] 
#>   eta2 =~ y23      0.7997      0.0440   18.1734    0.0000 [ 0.7151; 0.8647 ] 
#>   eta3 =~ y31      0.8223      0.0318   25.8711    0.0000 [ 0.7486; 0.8603 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1272    0.0000 [ 0.5995; 0.7425 ] 
#>   eta3 =~ y33      0.7474      0.0498   15.0072    0.0000 [ 0.6438; 0.8298 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0214   18.5185    0.0000 [ 0.3536; 0.4311 ] 
#>   eta1 <~ y12      0.3873      0.0181   21.4176    0.0000 [ 0.3606; 0.4299 ] 
#>   eta1 <~ y13      0.4542      0.0239   18.9713    0.0000 [ 0.4213; 0.5036 ] 
#>   eta2 <~ y21      0.3058      0.0331    9.2317    0.0000 [ 0.2496; 0.3511 ] 
#>   eta2 <~ y22      0.4473      0.0199   22.5005    0.0000 [ 0.4236; 0.4902 ] 
#>   eta2 <~ y23      0.4735      0.0249   18.9955    0.0000 [ 0.4293; 0.5119 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5332    0.0000 [ 0.3949; 0.4674 ] 
#>   eta3 <~ y32      0.3521      0.0206   17.0810    0.0000 [ 0.3236; 0.4046 ] 
#>   eta3 <~ y33      0.3999      0.0225   17.7859    0.0000 [ 0.3548; 0.4418 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0486   13.8243    0.0000 [ 0.5691; 0.7730 ] 
#>   eta3 ~ eta1       0.6634      0.0364   18.2308    0.0000 [ 0.6148; 0.7545 ] 
#>   eta3 ~ eta2       0.3052      0.0701    4.3550    0.0000 [ 0.1755; 0.3937 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0466    4.3981    0.0000 [ 0.1191; 0.2645 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.05054671 13.117964  2.598213e-39
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03387807 19.165140  7.236117e-82
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03016880 25.236194 1.605655e-140
#> 4 eta2 =~ y21  Common factor 0.5164548 0.06044077  8.544809  1.287488e-17
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03301484 22.880243 7.309214e-116
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04400182 18.173423  8.380355e-74
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03178357 25.871145 1.407089e-147
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04080493 16.127190  1.643294e-58
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04980441 15.007186  6.588359e-51
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5634909          0.7338544
#> 2          0.5877594          0.7046459
#> 3          0.7096343          0.8113896
#> 4          0.3958258          0.5966927
#> 5          0.7082681          0.8235918
#> 6          0.7150540          0.8647415
#> 7          0.7485966          0.8603363
#> 8          0.5995129          0.7425319
#> 9          0.6438045          0.8298484

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
#>  Random seed                        = 1418278302
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
#>   eta2 ~ eta1      0.6713      0.0486   13.8243    0.0000 [ 0.5499; 0.8010 ] 
#>   eta3 ~ eta1      0.4585      0.0668    6.8676    0.0000 [ 0.2692; 0.6145 ] 
#>   eta3 ~ eta2      0.3052      0.0701    4.3550    0.0000 [ 0.1353; 0.4976 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0505   13.1180    0.0000 [ 0.5417; 0.8031 ] 
#>   eta1 =~ y12      0.6493      0.0339   19.1651    0.0000 [ 0.5645; 0.7397 ] 
#>   eta1 =~ y13      0.7613      0.0302   25.2362    0.0000 [ 0.6834; 0.8394 ] 
#>   eta2 =~ y21      0.5165      0.0604    8.5448    0.0000 [ 0.3735; 0.6860 ] 
#>   eta2 =~ y22      0.7554      0.0330   22.8802    0.0000 [ 0.6704; 0.8411 ] 
#>   eta2 =~ y23      0.7997      0.0440   18.1734    0.0000 [ 0.6977; 0.9252 ] 
#>   eta3 =~ y31      0.8223      0.0318   25.8711    0.0000 [ 0.7544; 0.9188 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1272    0.0000 [ 0.5406; 0.7516 ] 
#>   eta3 =~ y33      0.7474      0.0498   15.0072    0.0000 [ 0.6190; 0.8766 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0214   18.5185    0.0000 [ 0.3432; 0.4536 ] 
#>   eta1 <~ y12      0.3873      0.0181   21.4176    0.0000 [ 0.3391; 0.4326 ] 
#>   eta1 <~ y13      0.4542      0.0239   18.9713    0.0000 [ 0.3882; 0.5120 ] 
#>   eta2 <~ y21      0.3058      0.0331    9.2317    0.0000 [ 0.2235; 0.3948 ] 
#>   eta2 <~ y22      0.4473      0.0199   22.5005    0.0000 [ 0.3887; 0.4915 ] 
#>   eta2 <~ y23      0.4735      0.0249   18.9955    0.0000 [ 0.4085; 0.5374 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5332    0.0000 [ 0.3962; 0.4971 ] 
#>   eta3 <~ y32      0.3521      0.0206   17.0810    0.0000 [ 0.2917; 0.3983 ] 
#>   eta3 <~ y33      0.3999      0.0225   17.7859    0.0000 [ 0.3414; 0.4577 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0486   13.8243    0.0000 [ 0.5499; 0.8010 ] 
#>   eta3 ~ eta1       0.6634      0.0364   18.2308    0.0000 [ 0.5620; 0.7502 ] 
#>   eta3 ~ eta2       0.3052      0.0701    4.3550    0.0000 [ 0.1353; 0.4976 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0466    4.3981    0.0000 [ 0.0938; 0.3347 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04856192 13.824276 1.819327e-43
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06676340  6.867637 6.527411e-12
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07006846  4.355042 1.330412e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5499118          0.8010467          0.5800682          0.7708903
#> 2          0.2692149          0.6144775          0.3106742          0.5730182
#> 3          0.1352638          0.4976184          0.1787755          0.4541067
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5683287          0.7819700          0.5691205          0.7730136
#> 2          0.3581208          0.5928304          0.3689443          0.5910221
#> 3          0.1747139          0.4034676          0.1755484          0.3936769
```
