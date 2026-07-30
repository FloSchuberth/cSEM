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
#>  Random seed                        = 552982520
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
#>   eta2 ~ eta1      0.6713      0.0360   18.6550    0.0000 [ 0.6180; 0.7290 ] 
#>   eta3 ~ eta1      0.4585      0.0751    6.1033    0.0000 [ 0.3090; 0.5754 ] 
#>   eta3 ~ eta2      0.3052      0.0748    4.0787    0.0000 [ 0.1939; 0.4350 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0363   18.2878    0.0000 [ 0.5967; 0.7390 ] 
#>   eta1 =~ y12      0.6493      0.0338   19.2236    0.0000 [ 0.5790; 0.7104 ] 
#>   eta1 =~ y13      0.7613      0.0352   21.6220    0.0000 [ 0.7060; 0.8304 ] 
#>   eta2 =~ y21      0.5165      0.0468   11.0457    0.0000 [ 0.4376; 0.6015 ] 
#>   eta2 =~ y22      0.7554      0.0349   21.6714    0.0000 [ 0.6701; 0.8059 ] 
#>   eta2 =~ y23      0.7997      0.0394   20.3120    0.0000 [ 0.7197; 0.8599 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5584    0.0000 [ 0.7743; 0.8757 ] 
#>   eta3 =~ y32      0.6581      0.0389   16.9064    0.0000 [ 0.5819; 0.7179 ] 
#>   eta3 =~ y33      0.7474      0.0344   21.7572    0.0000 [ 0.6931; 0.8069 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0202   19.5588    0.0000 [ 0.3608; 0.4355 ] 
#>   eta1 <~ y12      0.3873      0.0160   24.2784    0.0000 [ 0.3587; 0.4093 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9418    0.0000 [ 0.4150; 0.4977 ] 
#>   eta2 <~ y21      0.3058      0.0237   12.8967    0.0000 [ 0.2736; 0.3489 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6039    0.0000 [ 0.4126; 0.4913 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5407    0.0000 [ 0.4332; 0.5130 ] 
#>   eta3 <~ y31      0.4400      0.0179   24.6292    0.0000 [ 0.4147; 0.4730 ] 
#>   eta3 <~ y32      0.3521      0.0170   20.6759    0.0000 [ 0.3167; 0.3814 ] 
#>   eta3 <~ y33      0.3999      0.0168   23.7984    0.0000 [ 0.3725; 0.4272 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0360   18.6550    0.0000 [ 0.6180; 0.7290 ] 
#>   eta3 ~ eta1       0.6634      0.0360   18.4184    0.0000 [ 0.6057; 0.7400 ] 
#>   eta3 ~ eta2       0.3052      0.0748    4.0787    0.0000 [ 0.1939; 0.4350 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0505    4.0582    0.0000 [ 0.1304; 0.2956 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03625750 18.28780  1.035101e-74
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03377500 19.22363  2.347721e-82
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03521170 21.62196 1.116410e-103
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04675627 11.04568  2.300216e-28
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03485641 21.67141 3.818433e-104
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03936911 20.31196  1.008057e-91
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03096109 26.55841 2.053277e-155
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03892428 16.90639  4.037110e-64
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03435295 21.75720 5.904378e-105
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5966517          0.7389552
#> 2          0.5790304          0.7104049
#> 3          0.7059972          0.8303937
#> 4          0.4375897          0.6014621
#> 5          0.6700884          0.8058660
#> 6          0.7197052          0.8599190
#> 7          0.7742620          0.8756975
#> 8          0.5818655          0.7178877
#> 9          0.6931383          0.8069198

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
#>  Random seed                        = 552982520
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
#>   eta2 ~ eta1      0.6713      0.0360   18.6550    0.0000 [ 0.5749; 0.7610 ] 
#>   eta3 ~ eta1      0.4585      0.0751    6.1033    0.0000 [ 0.2611; 0.6496 ] 
#>   eta3 ~ eta2      0.3052      0.0748    4.0787    0.0000 [ 0.1122; 0.4991 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0363   18.2878    0.0000 [ 0.5645; 0.7520 ] 
#>   eta1 =~ y12      0.6493      0.0338   19.2236    0.0000 [ 0.5629; 0.7375 ] 
#>   eta1 =~ y13      0.7613      0.0352   21.6220    0.0000 [ 0.6719; 0.8540 ] 
#>   eta2 =~ y21      0.5165      0.0468   11.0457    0.0000 [ 0.3900; 0.6318 ] 
#>   eta2 =~ y22      0.7554      0.0349   21.6714    0.0000 [ 0.6694; 0.8497 ] 
#>   eta2 =~ y23      0.7997      0.0394   20.3120    0.0000 [ 0.7032; 0.9068 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5584    0.0000 [ 0.7445; 0.9047 ] 
#>   eta3 =~ y32      0.6581      0.0389   16.9064    0.0000 [ 0.5645; 0.7658 ] 
#>   eta3 =~ y33      0.7474      0.0344   21.7572    0.0000 [ 0.6597; 0.8374 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0202   19.5588    0.0000 [ 0.3412; 0.4458 ] 
#>   eta1 <~ y12      0.3873      0.0160   24.2784    0.0000 [ 0.3475; 0.4300 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9418    0.0000 [ 0.3999; 0.5120 ] 
#>   eta2 <~ y21      0.3058      0.0237   12.8967    0.0000 [ 0.2406; 0.3632 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6039    0.0000 [ 0.3893; 0.5073 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5407    0.0000 [ 0.4210; 0.5297 ] 
#>   eta3 <~ y31      0.4400      0.0179   24.6292    0.0000 [ 0.3920; 0.4844 ] 
#>   eta3 <~ y32      0.3521      0.0170   20.6759    0.0000 [ 0.3097; 0.3978 ] 
#>   eta3 <~ y33      0.3999      0.0168   23.7984    0.0000 [ 0.3544; 0.4413 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0360   18.6550    0.0000 [ 0.5749; 0.7610 ] 
#>   eta3 ~ eta1       0.6634      0.0360   18.4184    0.0000 [ 0.5667; 0.7529 ] 
#>   eta3 ~ eta2       0.3052      0.0748    4.0787    0.0000 [ 0.1122; 0.4991 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0505    4.0582    0.0000 [ 0.0739; 0.3349 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03598671 18.655039 1.149178e-77
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07512483  6.103265 1.039233e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07481642  4.078665 4.529507e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5749308          0.7610337          0.5972781          0.7386864
#> 2          0.2611439          0.6496471          0.3077956          0.6029955
#> 3          0.1121544          0.4990627          0.1586145          0.4526026
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5868765          0.7436498          0.6179634          0.7289768
#> 2          0.3070522          0.5965196          0.3089847          0.5753866
#> 3          0.1717462          0.4449479          0.1938656          0.4349919
```
