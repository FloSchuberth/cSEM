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
#>  Random seed                        = -565493332
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
#>   eta2 ~ eta1      0.6713      0.0490   13.7030    0.0000 [ 0.5860; 0.7444 ] 
#>   eta3 ~ eta1      0.4585      0.0788    5.8194    0.0000 [ 0.3141; 0.5953 ] 
#>   eta3 ~ eta2      0.3052      0.0744    4.1021    0.0000 [ 0.1922; 0.4322 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0425   15.5928    0.0000 [ 0.5901; 0.7404 ] 
#>   eta1 =~ y12      0.6493      0.0378   17.1949    0.0000 [ 0.5628; 0.6938 ] 
#>   eta1 =~ y13      0.7613      0.0319   23.8413    0.0000 [ 0.7195; 0.8351 ] 
#>   eta2 =~ y21      0.5165      0.0502   10.2855    0.0000 [ 0.4553; 0.6171 ] 
#>   eta2 =~ y22      0.7554      0.0403   18.7562    0.0000 [ 0.6851; 0.8301 ] 
#>   eta2 =~ y23      0.7997      0.0419   19.1059    0.0000 [ 0.7210; 0.8693 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.7154    0.0000 [ 0.7326; 0.8719 ] 
#>   eta3 =~ y32      0.6581      0.0385   17.1045    0.0000 [ 0.5722; 0.7189 ] 
#>   eta3 =~ y33      0.7474      0.0360   20.7739    0.0000 [ 0.6836; 0.8137 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0228   17.3629    0.0000 [ 0.3504; 0.4316 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.8044    0.0000 [ 0.3468; 0.4335 ] 
#>   eta1 <~ y13      0.4542      0.0175   26.0223    0.0000 [ 0.4240; 0.4882 ] 
#>   eta2 <~ y21      0.3058      0.0267   11.4340    0.0000 [ 0.2677; 0.3584 ] 
#>   eta2 <~ y22      0.4473      0.0235   18.9994    0.0000 [ 0.4079; 0.4928 ] 
#>   eta2 <~ y23      0.4735      0.0218   21.6727    0.0000 [ 0.4334; 0.5180 ] 
#>   eta3 <~ y31      0.4400      0.0179   24.5546    0.0000 [ 0.4071; 0.4653 ] 
#>   eta3 <~ y32      0.3521      0.0168   20.9882    0.0000 [ 0.3159; 0.3819 ] 
#>   eta3 <~ y33      0.3999      0.0211   18.9150    0.0000 [ 0.3709; 0.4452 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0490   13.7030    0.0000 [ 0.5860; 0.7444 ] 
#>   eta3 ~ eta1       0.6634      0.0395   16.7995    0.0000 [ 0.5751; 0.7453 ] 
#>   eta3 ~ eta2       0.3052      0.0744    4.1021    0.0000 [ 0.1922; 0.4322 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0524    3.9128    0.0001 [ 0.1324; 0.2965 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04252413 15.59279  8.149530e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03775981 17.19495  2.897233e-66
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03193385 23.84134 1.245328e-125
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05021180 10.28553  8.189341e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04027400 18.75621  1.722407e-78
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04185434 19.10588  2.256130e-81
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03619915 22.71537 3.157635e-114
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03847348 17.10448  1.374187e-65
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03597900 20.77390  7.455222e-96
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5901045          0.7404429
#> 2          0.5628051          0.6937956
#> 3          0.7195227          0.8351205
#> 4          0.4552694          0.6170906
#> 5          0.6851124          0.8300545
#> 6          0.7209759          0.8693144
#> 7          0.7325563          0.8718536
#> 8          0.5722481          0.7189005
#> 9          0.6836148          0.8137411

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
#>  Random seed                        = -565493332
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
#>   eta2 ~ eta1      0.6713      0.0490   13.7030    0.0000 [ 0.5410; 0.7944 ] 
#>   eta3 ~ eta1      0.4585      0.0788    5.8194    0.0000 [ 0.2589; 0.6664 ] 
#>   eta3 ~ eta2      0.3052      0.0744    4.1021    0.0000 [ 0.1121; 0.4968 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0425   15.5928    0.0000 [ 0.5565; 0.7764 ] 
#>   eta1 =~ y12      0.6493      0.0378   17.1949    0.0000 [ 0.5498; 0.7451 ] 
#>   eta1 =~ y13      0.7613      0.0319   23.8413    0.0000 [ 0.6666; 0.8317 ] 
#>   eta2 =~ y21      0.5165      0.0502   10.2855    0.0000 [ 0.3762; 0.6358 ] 
#>   eta2 =~ y22      0.7554      0.0403   18.7562    0.0000 [ 0.6574; 0.8657 ] 
#>   eta2 =~ y23      0.7997      0.0419   19.1059    0.0000 [ 0.6949; 0.9113 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.7154    0.0000 [ 0.7355; 0.9227 ] 
#>   eta3 =~ y32      0.6581      0.0385   17.1045    0.0000 [ 0.5611; 0.7600 ] 
#>   eta3 =~ y33      0.7474      0.0360   20.7739    0.0000 [ 0.6562; 0.8423 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0228   17.3629    0.0000 [ 0.3421; 0.4599 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.8044    0.0000 [ 0.3333; 0.4458 ] 
#>   eta1 <~ y13      0.4542      0.0175   26.0223    0.0000 [ 0.4056; 0.4959 ] 
#>   eta2 <~ y21      0.3058      0.0267   11.4340    0.0000 [ 0.2307; 0.3690 ] 
#>   eta2 <~ y22      0.4473      0.0235   18.9994    0.0000 [ 0.3901; 0.5118 ] 
#>   eta2 <~ y23      0.4735      0.0218   21.6727    0.0000 [ 0.4191; 0.5321 ] 
#>   eta3 <~ y31      0.4400      0.0179   24.5546    0.0000 [ 0.3941; 0.4868 ] 
#>   eta3 <~ y32      0.3521      0.0168   20.9882    0.0000 [ 0.3077; 0.3944 ] 
#>   eta3 <~ y33      0.3999      0.0211   18.9150    0.0000 [ 0.3432; 0.4525 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0490   13.7030    0.0000 [ 0.5410; 0.7944 ] 
#>   eta3 ~ eta1       0.6634      0.0395   16.7995    0.0000 [ 0.5636; 0.7678 ] 
#>   eta3 ~ eta2       0.3052      0.0744    4.1021    0.0000 [ 0.1121; 0.4968 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0524    3.9128    0.0001 [ 0.0677; 0.3385 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04899161 13.70303 9.737505e-43
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07878895  5.81943 5.904860e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07438900  4.10210 4.094168e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5410079          0.7943649          0.5714311          0.7639416
#> 2          0.2589001          0.6663520          0.3078271          0.6174250
#> 3          0.1121268          0.4968247          0.1583215          0.4506300
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5416007          0.7482282          0.5860081          0.7443534
#> 2          0.2795634          0.6212660          0.3140517          0.5952637
#> 3          0.1467875          0.4547769          0.1922014          0.4322451
```
