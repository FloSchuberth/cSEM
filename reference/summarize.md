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
#>  Random seed                        = -1340897223
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
#>   eta2 ~ eta1      0.6713      0.0439   15.3091    0.0000 [ 0.5855; 0.7484 ] 
#>   eta3 ~ eta1      0.4585      0.0971    4.7213    0.0000 [ 0.2692; 0.6032 ] 
#>   eta3 ~ eta2      0.3052      0.0974    3.1335    0.0017 [ 0.1211; 0.4862 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0435   15.2376    0.0000 [ 0.5977; 0.7463 ] 
#>   eta1 =~ y12      0.6493      0.0385   16.8496    0.0000 [ 0.5731; 0.7012 ] 
#>   eta1 =~ y13      0.7613      0.0383   19.8906    0.0000 [ 0.6973; 0.8229 ] 
#>   eta2 =~ y21      0.5165      0.0422   12.2518    0.0000 [ 0.4312; 0.5757 ] 
#>   eta2 =~ y22      0.7554      0.0392   19.2557    0.0000 [ 0.6840; 0.8141 ] 
#>   eta2 =~ y23      0.7997      0.0401   19.9365    0.0000 [ 0.7148; 0.8499 ] 
#>   eta3 =~ y31      0.8223      0.0309   26.6202    0.0000 [ 0.7705; 0.8823 ] 
#>   eta3 =~ y32      0.6581      0.0398   16.5287    0.0000 [ 0.5964; 0.7263 ] 
#>   eta3 =~ y33      0.7474      0.0350   21.3623    0.0000 [ 0.6843; 0.8061 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0196   20.2158    0.0000 [ 0.3566; 0.4254 ] 
#>   eta1 <~ y12      0.3873      0.0198   19.5820    0.0000 [ 0.3517; 0.4223 ] 
#>   eta1 <~ y13      0.4542      0.0234   19.4388    0.0000 [ 0.4105; 0.5084 ] 
#>   eta2 <~ y21      0.3058      0.0263   11.6175    0.0000 [ 0.2567; 0.3484 ] 
#>   eta2 <~ y22      0.4473      0.0186   24.1069    0.0000 [ 0.4174; 0.4865 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3467    0.0000 [ 0.4331; 0.5121 ] 
#>   eta3 <~ y31      0.4400      0.0169   26.0075    0.0000 [ 0.4086; 0.4697 ] 
#>   eta3 <~ y32      0.3521      0.0164   21.4166    0.0000 [ 0.3188; 0.3763 ] 
#>   eta3 <~ y33      0.3999      0.0178   22.4570    0.0000 [ 0.3775; 0.4526 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0439   15.3091    0.0000 [ 0.5855; 0.7484 ] 
#>   eta3 ~ eta1       0.6634      0.0462   14.3674    0.0000 [ 0.5791; 0.7241 ] 
#>   eta3 ~ eta2       0.3052      0.0974    3.1335    0.0017 [ 0.1211; 0.4862 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0657    3.1178    0.0018 [ 0.0842; 0.3425 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04351538 15.23760  1.990539e-52
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03853384 16.84955  1.057107e-63
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03827669 19.89059  4.909785e-88
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04215329 12.25183  1.642546e-34
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03922926 19.25572  1.264041e-82
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04011058 19.93648  1.964251e-88
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03088924 26.62019 3.963950e-156
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03981373 16.52869  2.280337e-61
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03498804 21.36228 2.998337e-101
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5976742          0.7463465
#> 2          0.5731403          0.7012142
#> 3          0.6973373          0.8229348
#> 4          0.4311610          0.5757297
#> 5          0.6839529          0.8140934
#> 6          0.7148066          0.8498666
#> 7          0.7704985          0.8823328
#> 8          0.5964261          0.7263419
#> 9          0.6842663          0.8060559

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
#>  Random seed                        = -1340897223
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
#>   eta2 ~ eta1      0.6713      0.0439   15.3091    0.0000 [ 0.5504; 0.7772 ] 
#>   eta3 ~ eta1      0.4585      0.0971    4.7213    0.0000 [ 0.2167; 0.7190 ] 
#>   eta3 ~ eta2      0.3052      0.0974    3.1335    0.0017 [ 0.0498; 0.5534 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0435   15.2376    0.0000 [ 0.5437; 0.7687 ] 
#>   eta1 =~ y12      0.6493      0.0385   16.8496    0.0000 [ 0.5512; 0.7504 ] 
#>   eta1 =~ y13      0.7613      0.0383   19.8906    0.0000 [ 0.6628; 0.8608 ] 
#>   eta2 =~ y21      0.5165      0.0422   12.2518    0.0000 [ 0.4054; 0.6234 ] 
#>   eta2 =~ y22      0.7554      0.0392   19.2557    0.0000 [ 0.6534; 0.8563 ] 
#>   eta2 =~ y23      0.7997      0.0401   19.9365    0.0000 [ 0.7025; 0.9099 ] 
#>   eta3 =~ y31      0.8223      0.0309   26.6202    0.0000 [ 0.7405; 0.9002 ] 
#>   eta3 =~ y32      0.6581      0.0398   16.5287    0.0000 [ 0.5603; 0.7662 ] 
#>   eta3 =~ y33      0.7474      0.0350   21.3623    0.0000 [ 0.6510; 0.8319 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0196   20.2158    0.0000 [ 0.3426; 0.4438 ] 
#>   eta1 <~ y12      0.3873      0.0198   19.5820    0.0000 [ 0.3386; 0.4409 ] 
#>   eta1 <~ y13      0.4542      0.0234   19.4388    0.0000 [ 0.3956; 0.5164 ] 
#>   eta2 <~ y21      0.3058      0.0263   11.6175    0.0000 [ 0.2356; 0.3717 ] 
#>   eta2 <~ y22      0.4473      0.0186   24.1069    0.0000 [ 0.3980; 0.4940 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3467    0.0000 [ 0.4215; 0.5311 ] 
#>   eta3 <~ y31      0.4400      0.0169   26.0075    0.0000 [ 0.3961; 0.4836 ] 
#>   eta3 <~ y32      0.3521      0.0164   21.4166    0.0000 [ 0.3134; 0.3984 ] 
#>   eta3 <~ y33      0.3999      0.0178   22.4570    0.0000 [ 0.3516; 0.4437 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0439   15.3091    0.0000 [ 0.5504; 0.7772 ] 
#>   eta3 ~ eta1       0.6634      0.0462   14.3674    0.0000 [ 0.5494; 0.7882 ] 
#>   eta3 ~ eta2       0.3052      0.0974    3.1335    0.0017 [ 0.0498; 0.5534 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0657    3.1178    0.0018 [ 0.0310; 0.3708 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04385206 15.309051 6.652772e-53
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09711440  4.721306 2.343352e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09738473  3.133460 1.727587e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55043129          0.7772094          0.5776629          0.7499778
#> 2         0.21674602          0.7189668          0.2770529          0.6586599
#> 3         0.04980941          0.5534283          0.1102842          0.4929535
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5618515          0.7515140          0.5854629          0.7484187
#> 2          0.2320661          0.6175332          0.2692137          0.6031861
#> 3          0.1143368          0.5086852          0.1210909          0.4861837
```
