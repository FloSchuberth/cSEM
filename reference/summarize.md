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
#>  Random seed                        = -347248402
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
#>   eta2 ~ eta1      0.6713      0.0384   17.4859    0.0000 [ 0.6162; 0.7428 ] 
#>   eta3 ~ eta1      0.4585      0.0987    4.6431    0.0000 [ 0.2739; 0.6385 ] 
#>   eta3 ~ eta2      0.3052      0.0996    3.0643    0.0022 [ 0.1431; 0.4760 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0427   15.5328    0.0000 [ 0.5742; 0.7136 ] 
#>   eta1 =~ y12      0.6493      0.0406   15.9904    0.0000 [ 0.5791; 0.7049 ] 
#>   eta1 =~ y13      0.7613      0.0354   21.5300    0.0000 [ 0.7028; 0.8360 ] 
#>   eta2 =~ y21      0.5165      0.0543    9.5170    0.0000 [ 0.4460; 0.6275 ] 
#>   eta2 =~ y22      0.7554      0.0430   17.5594    0.0000 [ 0.6733; 0.8281 ] 
#>   eta2 =~ y23      0.7997      0.0312   25.6305    0.0000 [ 0.7484; 0.8405 ] 
#>   eta3 =~ y31      0.8223      0.0281   29.2481    0.0000 [ 0.7644; 0.8545 ] 
#>   eta3 =~ y32      0.6581      0.0392   16.7943    0.0000 [ 0.5898; 0.7292 ] 
#>   eta3 =~ y33      0.7474      0.0387   19.3377    0.0000 [ 0.6789; 0.8164 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0232   17.0649    0.0000 [ 0.3474; 0.4303 ] 
#>   eta1 <~ y12      0.3873      0.0228   16.9729    0.0000 [ 0.3460; 0.4219 ] 
#>   eta1 <~ y13      0.4542      0.0173   26.2756    0.0000 [ 0.4339; 0.4947 ] 
#>   eta2 <~ y21      0.3058      0.0285   10.7453    0.0000 [ 0.2557; 0.3536 ] 
#>   eta2 <~ y22      0.4473      0.0203   22.0022    0.0000 [ 0.4125; 0.4837 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1354    0.0000 [ 0.4305; 0.5229 ] 
#>   eta3 <~ y31      0.4400      0.0184   23.9293    0.0000 [ 0.4068; 0.4699 ] 
#>   eta3 <~ y32      0.3521      0.0199   17.6514    0.0000 [ 0.3155; 0.3842 ] 
#>   eta3 <~ y33      0.3999      0.0142   28.2459    0.0000 [ 0.3727; 0.4244 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0384   17.4859    0.0000 [ 0.6162; 0.7428 ] 
#>   eta3 ~ eta1       0.6634      0.0469   14.1563    0.0000 [ 0.5912; 0.7361 ] 
#>   eta3 ~ eta2       0.3052      0.0996    3.0643    0.0022 [ 0.1431; 0.4760 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0685    2.9905    0.0028 [ 0.0932; 0.3319 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04268824 15.532848  2.079373e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04060414 15.990437  1.489823e-57
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03536203 21.530038 8.146641e-103
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05426643  9.517023  1.782117e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04301894 17.559421  5.038698e-69
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03119964 25.630542 6.968597e-145
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02811389 29.248084 4.749331e-188
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03918407 16.794295  2.686775e-63
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03865114 19.337698  2.587714e-83
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5741546          0.7136323
#> 2          0.5790528          0.7049417
#> 3          0.7028252          0.8359716
#> 4          0.4460122          0.6274578
#> 5          0.6732684          0.8280945
#> 6          0.7484014          0.8404883
#> 7          0.7643887          0.8545406
#> 8          0.5898221          0.7291947
#> 9          0.6789133          0.8164405

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
#>  Random seed                        = -347248402
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
#>   eta2 ~ eta1      0.6713      0.0384   17.4859    0.0000 [ 0.5688; 0.7674 ] 
#>   eta3 ~ eta1      0.4585      0.0987    4.6431    0.0000 [ 0.2036; 0.7142 ] 
#>   eta3 ~ eta2      0.3052      0.0996    3.0643    0.0022 [ 0.0418; 0.5568 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0427   15.5328    0.0000 [ 0.5622; 0.7830 ] 
#>   eta1 =~ y12      0.6493      0.0406   15.9904    0.0000 [ 0.5513; 0.7613 ] 
#>   eta1 =~ y13      0.7613      0.0354   21.5300    0.0000 [ 0.6616; 0.8445 ] 
#>   eta2 =~ y21      0.5165      0.0543    9.5170    0.0000 [ 0.3679; 0.6485 ] 
#>   eta2 =~ y22      0.7554      0.0430   17.5594    0.0000 [ 0.6413; 0.8638 ] 
#>   eta2 =~ y23      0.7997      0.0312   25.6305    0.0000 [ 0.7224; 0.8838 ] 
#>   eta3 =~ y31      0.8223      0.0281   29.2481    0.0000 [ 0.7600; 0.9054 ] 
#>   eta3 =~ y32      0.6581      0.0392   16.7943    0.0000 [ 0.5507; 0.7534 ] 
#>   eta3 =~ y33      0.7474      0.0387   19.3377    0.0000 [ 0.6424; 0.8422 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0232   17.0649    0.0000 [ 0.3394; 0.4593 ] 
#>   eta1 <~ y12      0.3873      0.0228   16.9729    0.0000 [ 0.3306; 0.4486 ] 
#>   eta1 <~ y13      0.4542      0.0173   26.2756    0.0000 [ 0.4023; 0.4917 ] 
#>   eta2 <~ y21      0.3058      0.0285   10.7453    0.0000 [ 0.2294; 0.3766 ] 
#>   eta2 <~ y22      0.4473      0.0203   22.0022    0.0000 [ 0.3957; 0.5009 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1354    0.0000 [ 0.4171; 0.5387 ] 
#>   eta3 <~ y31      0.4400      0.0184   23.9293    0.0000 [ 0.3981; 0.4931 ] 
#>   eta3 <~ y32      0.3521      0.0199   17.6514    0.0000 [ 0.2975; 0.4007 ] 
#>   eta3 <~ y33      0.3999      0.0142   28.2459    0.0000 [ 0.3610; 0.4342 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0384   17.4859    0.0000 [ 0.5688; 0.7674 ] 
#>   eta3 ~ eta1       0.6634      0.0469   14.1563    0.0000 [ 0.5380; 0.7804 ] 
#>   eta3 ~ eta2       0.3052      0.0996    3.0643    0.0022 [ 0.0418; 0.5568 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0685    2.9905    0.0028 [ 0.0232; 0.3774 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03839285 17.485896 1.835047e-68
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09874940  4.643135 3.431627e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09958352  3.064273 2.181995e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56882800          0.7673742          0.5926695          0.7435326
#> 2         0.20355995          0.7142361          0.2648822          0.6529139
#> 3         0.04178622          0.5567760          0.1036264          0.4949358
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.61589350          0.7505902          0.6162267          0.7427760
#> 2         0.25377742          0.6709168          0.2738655          0.6385144
#> 3         0.09229256          0.4958052          0.1430889          0.4759622
```
