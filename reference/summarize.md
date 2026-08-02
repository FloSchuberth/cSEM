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
#>  Random seed                        = -766233315
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
#>   eta2 ~ eta1      0.6713      0.0406   16.5445    0.0000 [ 0.5921; 0.7371 ] 
#>   eta3 ~ eta1      0.4585      0.0951    4.8218    0.0000 [ 0.3313; 0.6808 ] 
#>   eta3 ~ eta2      0.3052      0.1008    3.0273    0.0025 [ 0.0703; 0.4607 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0407   16.2744    0.0000 [ 0.5913; 0.7267 ] 
#>   eta1 =~ y12      0.6493      0.0461   14.0779    0.0000 [ 0.5765; 0.7244 ] 
#>   eta1 =~ y13      0.7613      0.0303   25.1634    0.0000 [ 0.7002; 0.8172 ] 
#>   eta2 =~ y21      0.5165      0.0484   10.6761    0.0000 [ 0.4414; 0.6023 ] 
#>   eta2 =~ y22      0.7554      0.0359   21.0220    0.0000 [ 0.6915; 0.8374 ] 
#>   eta2 =~ y23      0.7997      0.0395   20.2376    0.0000 [ 0.7324; 0.8601 ] 
#>   eta3 =~ y31      0.8223      0.0289   28.4322    0.0000 [ 0.7847; 0.8797 ] 
#>   eta3 =~ y32      0.6581      0.0440   14.9469    0.0000 [ 0.5639; 0.7189 ] 
#>   eta3 =~ y33      0.7474      0.0307   24.3485    0.0000 [ 0.6939; 0.7923 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0250   15.8427    0.0000 [ 0.3538; 0.4472 ] 
#>   eta1 <~ y12      0.3873      0.0201   19.2304    0.0000 [ 0.3513; 0.4226 ] 
#>   eta1 <~ y13      0.4542      0.0193   23.5710    0.0000 [ 0.4205; 0.4843 ] 
#>   eta2 <~ y21      0.3058      0.0260   11.7641    0.0000 [ 0.2604; 0.3396 ] 
#>   eta2 <~ y22      0.4473      0.0214   20.8834    0.0000 [ 0.4183; 0.4855 ] 
#>   eta2 <~ y23      0.4735      0.0213   22.1861    0.0000 [ 0.4314; 0.5078 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1118    0.0000 [ 0.4168; 0.4830 ] 
#>   eta3 <~ y32      0.3521      0.0167   21.1117    0.0000 [ 0.3191; 0.3719 ] 
#>   eta3 <~ y33      0.3999      0.0189   21.1687    0.0000 [ 0.3601; 0.4327 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0406   16.5445    0.0000 [ 0.5921; 0.7371 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6880    0.0000 [ 0.6104; 0.7602 ] 
#>   eta3 ~ eta2       0.3052      0.1008    3.0273    0.0025 [ 0.0703; 0.4607 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0690    2.9669    0.0030 [ 0.0471; 0.3199 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04074318 16.27438  1.500412e-59
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04612021 14.07795  5.189449e-45
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03025602 25.16344 1.007120e-139
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04837502 10.67606  1.317392e-26
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03593311 21.02205  4.122415e-98
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03951381 20.23758  4.570946e-91
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02892063 28.43221 8.088092e-178
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04402712 14.94690  1.631704e-50
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03069694 24.34849 6.015433e-131
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5912627          0.7266564
#> 2          0.5765162          0.7244263
#> 3          0.7002083          0.8171646
#> 4          0.4414491          0.6023413
#> 5          0.6915010          0.8374415
#> 6          0.7323555          0.8601154
#> 7          0.7846663          0.8796931
#> 8          0.5638892          0.7189388
#> 9          0.6939432          0.7922533

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
#>  Random seed                        = -766233315
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
#>   eta2 ~ eta1      0.6713      0.0406   16.5445    0.0000 [ 0.5553; 0.7651 ] 
#>   eta3 ~ eta1      0.4585      0.0951    4.8218    0.0000 [ 0.2120; 0.7038 ] 
#>   eta3 ~ eta2      0.3052      0.1008    3.0273    0.0025 [ 0.0404; 0.5617 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0407   16.2744    0.0000 [ 0.5623; 0.7730 ] 
#>   eta1 =~ y12      0.6493      0.0461   14.0779    0.0000 [ 0.5262; 0.7647 ] 
#>   eta1 =~ y13      0.7613      0.0303   25.1634    0.0000 [ 0.6874; 0.8438 ] 
#>   eta2 =~ y21      0.5165      0.0484   10.6761    0.0000 [ 0.4010; 0.6511 ] 
#>   eta2 =~ y22      0.7554      0.0359   21.0220    0.0000 [ 0.6601; 0.8459 ] 
#>   eta2 =~ y23      0.7997      0.0395   20.2376    0.0000 [ 0.6927; 0.8971 ] 
#>   eta3 =~ y31      0.8223      0.0289   28.4322    0.0000 [ 0.7326; 0.8822 ] 
#>   eta3 =~ y32      0.6581      0.0440   14.9469    0.0000 [ 0.5594; 0.7871 ] 
#>   eta3 =~ y33      0.7474      0.0307   24.3485    0.0000 [ 0.6770; 0.8357 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0250   15.8427    0.0000 [ 0.3324; 0.4615 ] 
#>   eta1 <~ y12      0.3873      0.0201   19.2304    0.0000 [ 0.3321; 0.4363 ] 
#>   eta1 <~ y13      0.4542      0.0193   23.5710    0.0000 [ 0.4053; 0.5050 ] 
#>   eta2 <~ y21      0.3058      0.0260   11.7641    0.0000 [ 0.2446; 0.3790 ] 
#>   eta2 <~ y22      0.4473      0.0214   20.8834    0.0000 [ 0.3905; 0.5013 ] 
#>   eta2 <~ y23      0.4735      0.0213   22.1861    0.0000 [ 0.4156; 0.5260 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1118    0.0000 [ 0.3860; 0.4732 ] 
#>   eta3 <~ y32      0.3521      0.0167   21.1117    0.0000 [ 0.3157; 0.4019 ] 
#>   eta3 <~ y33      0.3999      0.0189   21.1687    0.0000 [ 0.3536; 0.4513 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0406   16.5445    0.0000 [ 0.5553; 0.7651 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6880    0.0000 [ 0.5541; 0.7597 ] 
#>   eta3 ~ eta2       0.3052      0.1008    3.0273    0.0025 [ 0.0404; 0.5617 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0690    2.9669    0.0030 [ 0.0205; 0.3776 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04057748 16.544482 1.754680e-61
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09508967  4.821836 1.422431e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10080117  3.027258 2.467834e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55527905          0.7651229          0.5804772          0.7399247
#> 2         0.21200455          0.7037546          0.2710541          0.6447051
#> 3         0.04044389          0.5617306          0.1030402          0.4991343
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5917876          0.7422265         0.59213978          0.7371374
#> 2          0.3213136          0.7104236         0.33132078          0.6808351
#> 3          0.0669674          0.4721194         0.07031943          0.4607440
```
