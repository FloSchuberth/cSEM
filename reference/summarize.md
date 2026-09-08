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
#>  Random seed                        = 1594960081
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
#>   eta2 ~ eta1      0.6713      0.0423   15.8716    0.0000 [ 0.5795; 0.7347 ] 
#>   eta3 ~ eta1      0.4585      0.0661    6.9399    0.0000 [ 0.3865; 0.6464 ] 
#>   eta3 ~ eta2      0.3052      0.0726    4.2036    0.0000 [ 0.1022; 0.4013 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0395   16.7700    0.0000 [ 0.5779; 0.7203 ] 
#>   eta1 =~ y12      0.6493      0.0378   17.1578    0.0000 [ 0.5587; 0.7129 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.1787    0.0000 [ 0.6695; 0.8020 ] 
#>   eta2 =~ y21      0.5165      0.0437   11.8112    0.0000 [ 0.4317; 0.5860 ] 
#>   eta2 =~ y22      0.7554      0.0344   21.9847    0.0000 [ 0.6839; 0.8067 ] 
#>   eta2 =~ y23      0.7997      0.0371   21.5731    0.0000 [ 0.7403; 0.8739 ] 
#>   eta3 =~ y31      0.8223      0.0350   23.4694    0.0000 [ 0.7570; 0.8751 ] 
#>   eta3 =~ y32      0.6581      0.0478   13.7602    0.0000 [ 0.5501; 0.7376 ] 
#>   eta3 =~ y33      0.7474      0.0487   15.3463    0.0000 [ 0.6518; 0.8350 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0200   19.8066    0.0000 [ 0.3631; 0.4289 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.4883    0.0000 [ 0.3508; 0.4314 ] 
#>   eta1 <~ y13      0.4542      0.0225   20.2013    0.0000 [ 0.4073; 0.4996 ] 
#>   eta2 <~ y21      0.3058      0.0218   14.0595    0.0000 [ 0.2628; 0.3376 ] 
#>   eta2 <~ y22      0.4473      0.0195   22.9175    0.0000 [ 0.4105; 0.4832 ] 
#>   eta2 <~ y23      0.4735      0.0194   24.3601    0.0000 [ 0.4459; 0.5114 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1256    0.0000 [ 0.4124; 0.4639 ] 
#>   eta3 <~ y32      0.3521      0.0234   15.0635    0.0000 [ 0.3165; 0.3874 ] 
#>   eta3 <~ y33      0.3999      0.0242   16.4999    0.0000 [ 0.3578; 0.4455 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0423   15.8716    0.0000 [ 0.5795; 0.7347 ] 
#>   eta3 ~ eta1       0.6634      0.0336   19.7495    0.0000 [ 0.6269; 0.7413 ] 
#>   eta3 ~ eta2       0.3052      0.0726    4.2036    0.0000 [ 0.1022; 0.4013 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0483    4.2418    0.0000 [ 0.0718; 0.2673 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03953900 16.77002  4.043630e-63
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03784146 17.15785  5.491357e-66
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03594868 21.17868  1.501721e-99
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04372573 11.81123  3.415126e-32
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03435968 21.98471 4.033555e-107
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03706762 21.57311 3.213628e-103
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03503621 23.46936 8.386858e-122
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04782416 13.76018  4.424440e-43
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04870379 15.34632  3.748362e-53
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5778684          0.7202608
#> 2          0.5586623          0.7128565
#> 3          0.6694749          0.8020344
#> 4          0.4317495          0.5859896
#> 5          0.6839151          0.8066618
#> 6          0.7402870          0.8738720
#> 7          0.7570500          0.8750572
#> 8          0.5500559          0.7375631
#> 9          0.6518265          0.8350059

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
#>  Random seed                        = 1594960081
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
#>   eta2 ~ eta1      0.6713      0.0423   15.8716    0.0000 [ 0.5544; 0.7731 ] 
#>   eta3 ~ eta1      0.4585      0.0661    6.9399    0.0000 [ 0.2603; 0.6019 ] 
#>   eta3 ~ eta2      0.3052      0.0726    4.2036    0.0000 [ 0.1353; 0.5108 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0395   16.7700    0.0000 [ 0.5597; 0.7641 ] 
#>   eta1 =~ y12      0.6493      0.0378   17.1578    0.0000 [ 0.5616; 0.7573 ] 
#>   eta1 =~ y13      0.7613      0.0359   21.1787    0.0000 [ 0.6716; 0.8575 ] 
#>   eta2 =~ y21      0.5165      0.0437   11.8112    0.0000 [ 0.4062; 0.6323 ] 
#>   eta2 =~ y22      0.7554      0.0344   21.9847    0.0000 [ 0.6689; 0.8466 ] 
#>   eta2 =~ y23      0.7997      0.0371   21.5731    0.0000 [ 0.7045; 0.8962 ] 
#>   eta3 =~ y31      0.8223      0.0350   23.4694    0.0000 [ 0.7385; 0.9197 ] 
#>   eta3 =~ y32      0.6581      0.0478   13.7602    0.0000 [ 0.5295; 0.7768 ] 
#>   eta3 =~ y33      0.7474      0.0487   15.3463    0.0000 [ 0.6291; 0.8809 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0200   19.8066    0.0000 [ 0.3401; 0.4434 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.4883    0.0000 [ 0.3389; 0.4417 ] 
#>   eta1 <~ y13      0.4542      0.0225   20.2013    0.0000 [ 0.3941; 0.5104 ] 
#>   eta2 <~ y21      0.3058      0.0218   14.0595    0.0000 [ 0.2504; 0.3629 ] 
#>   eta2 <~ y22      0.4473      0.0195   22.9175    0.0000 [ 0.3966; 0.4975 ] 
#>   eta2 <~ y23      0.4735      0.0194   24.3601    0.0000 [ 0.4220; 0.5225 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1256    0.0000 [ 0.3973; 0.4844 ] 
#>   eta3 <~ y32      0.3521      0.0234   15.0635    0.0000 [ 0.2869; 0.4078 ] 
#>   eta3 <~ y33      0.3999      0.0242   16.4999    0.0000 [ 0.3389; 0.4642 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0423   15.8716    0.0000 [ 0.5544; 0.7731 ] 
#>   eta3 ~ eta1       0.6634      0.0336   19.7495    0.0000 [ 0.5596; 0.7333 ] 
#>   eta3 ~ eta2       0.3052      0.0726    4.2036    0.0000 [ 0.1353; 0.5108 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0483    4.2418    0.0000 [ 0.0905; 0.3403 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04229781 15.871590 9.968386e-57
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06606839  6.939881 3.924315e-12
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07259346  4.203562 2.627466e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5543980          0.7731383          0.5806644          0.7468719
#> 2          0.2602598          0.6019282          0.3012875          0.5609005
#> 3          0.1353414          0.5107538          0.1804211          0.4656741
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.57443183          0.7578779          0.5794938          0.7347298
#> 2         0.37934548          0.6727342          0.3864760          0.6464350
#> 3         0.07155002          0.4056928          0.1021553          0.4013073
```
