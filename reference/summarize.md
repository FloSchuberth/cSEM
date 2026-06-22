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
#>  Random seed                        = 10403
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
#>   eta2 ~ eta1      0.6713      0.0424   15.8310    0.0000 [ 0.5800; 0.7347 ] 
#>   eta3 ~ eta1      0.4585      0.0815    5.6235    0.0000 [ 0.3241; 0.6081 ] 
#>   eta3 ~ eta2      0.3052      0.0838    3.6415    0.0003 [ 0.1764; 0.4233 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0375   17.6926    0.0000 [ 0.6094; 0.7360 ] 
#>   eta1 =~ y12      0.6493      0.0443   14.6426    0.0000 [ 0.5639; 0.7082 ] 
#>   eta1 =~ y13      0.7613      0.0335   22.7590    0.0000 [ 0.7067; 0.8204 ] 
#>   eta2 =~ y21      0.5165      0.0503   10.2589    0.0000 [ 0.4046; 0.5865 ] 
#>   eta2 =~ y22      0.7554      0.0448   16.8710    0.0000 [ 0.6601; 0.8259 ] 
#>   eta2 =~ y23      0.7997      0.0371   21.5379    0.0000 [ 0.7114; 0.8326 ] 
#>   eta3 =~ y31      0.8223      0.0293   28.0172    0.0000 [ 0.7718; 0.8718 ] 
#>   eta3 =~ y32      0.6581      0.0403   16.3276    0.0000 [ 0.5909; 0.7343 ] 
#>   eta3 =~ y33      0.7474      0.0398   18.7636    0.0000 [ 0.6614; 0.8083 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0204   19.4260    0.0000 [ 0.3620; 0.4377 ] 
#>   eta1 <~ y12      0.3873      0.0228   16.9549    0.0000 [ 0.3387; 0.4247 ] 
#>   eta1 <~ y13      0.4542      0.0222   20.4215    0.0000 [ 0.4214; 0.4910 ] 
#>   eta2 <~ y21      0.3058      0.0302   10.1121    0.0000 [ 0.2555; 0.3575 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6497    0.0000 [ 0.4191; 0.4929 ] 
#>   eta2 <~ y23      0.4735      0.0220   21.5057    0.0000 [ 0.4387; 0.5093 ] 
#>   eta3 <~ y31      0.4400      0.0167   26.3346    0.0000 [ 0.4143; 0.4733 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.5834    0.0000 [ 0.3087; 0.3787 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4233    0.0000 [ 0.3599; 0.4309 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0424   15.8310    0.0000 [ 0.5800; 0.7347 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6835    0.0000 [ 0.5957; 0.7282 ] 
#>   eta3 ~ eta2       0.3052      0.0838    3.6415    0.0003 [ 0.1764; 0.4233 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0584    3.5095    0.0004 [ 0.1196; 0.2965 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03747726 17.69259  4.782225e-70
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04434175 14.64259  1.502445e-48
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03345251 22.75901 1.168637e-114
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05034194 10.25894  1.078863e-24
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04477426 16.87102  7.351034e-64
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03712818 21.53792 6.872250e-103
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02934897 28.01725 1.001695e-172
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04030416 16.32757  6.284196e-60
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03983367 18.76363  1.498138e-78
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6094465          0.7360277
#> 2          0.5639472          0.7082055
#> 3          0.7066634          0.8203975
#> 4          0.4045539          0.5865389
#> 5          0.6601402          0.8258894
#> 6          0.7113685          0.8325669
#> 7          0.7717507          0.8717903
#> 8          0.5909261          0.7342961
#> 9          0.6614041          0.8083259

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
#>  Random seed                        = 10403
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
#>   eta2 ~ eta1      0.6713      0.0424   15.8310    0.0000 [ 0.5669; 0.7862 ] 
#>   eta3 ~ eta1      0.4585      0.0815    5.6235    0.0000 [ 0.2571; 0.6788 ] 
#>   eta3 ~ eta2      0.3052      0.0838    3.6415    0.0003 [ 0.0793; 0.5127 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0375   17.6926    0.0000 [ 0.5641; 0.7579 ] 
#>   eta1 =~ y12      0.6493      0.0443   14.6426    0.0000 [ 0.5443; 0.7737 ] 
#>   eta1 =~ y13      0.7613      0.0335   22.7590    0.0000 [ 0.6751; 0.8481 ] 
#>   eta2 =~ y21      0.5165      0.0503   10.2589    0.0000 [ 0.3872; 0.6475 ] 
#>   eta2 =~ y22      0.7554      0.0448   16.8710    0.0000 [ 0.6409; 0.8725 ] 
#>   eta2 =~ y23      0.7997      0.0371   21.5379    0.0000 [ 0.7191; 0.9111 ] 
#>   eta3 =~ y31      0.8223      0.0293   28.0172    0.0000 [ 0.7495; 0.9013 ] 
#>   eta3 =~ y32      0.6581      0.0403   16.3276    0.0000 [ 0.5639; 0.7723 ] 
#>   eta3 =~ y33      0.7474      0.0398   18.7636    0.0000 [ 0.6367; 0.8427 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0204   19.4260    0.0000 [ 0.3397; 0.4450 ] 
#>   eta1 <~ y12      0.3873      0.0228   16.9549    0.0000 [ 0.3324; 0.4505 ] 
#>   eta1 <~ y13      0.4542      0.0222   20.4215    0.0000 [ 0.3944; 0.5094 ] 
#>   eta2 <~ y21      0.3058      0.0302   10.1121    0.0000 [ 0.2244; 0.3808 ] 
#>   eta2 <~ y22      0.4473      0.0228   19.6497    0.0000 [ 0.3840; 0.5017 ] 
#>   eta2 <~ y23      0.4735      0.0220   21.5057    0.0000 [ 0.4200; 0.5339 ] 
#>   eta3 <~ y31      0.4400      0.0167   26.3346    0.0000 [ 0.3970; 0.4834 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.5834    0.0000 [ 0.3076; 0.4056 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4233    0.0000 [ 0.3440; 0.4453 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0424   15.8310    0.0000 [ 0.5669; 0.7862 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6835    0.0000 [ 0.5657; 0.7713 ] 
#>   eta3 ~ eta2       0.3052      0.0838    3.6415    0.0003 [ 0.0793; 0.5127 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0584    3.5095    0.0004 [ 0.0496; 0.3515 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04240636 15.830961 1.902931e-56
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08153471  5.623455 1.871753e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08379892  3.641469 2.710871e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5668645          0.7861663          0.5931984          0.7598324
#> 2          0.2571079          0.6787594          0.3077400          0.6281273
#> 3          0.0793031          0.5126638          0.1313412          0.4606256
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.57001358          0.7406270          0.5799989          0.7346613
#> 2         0.31317118          0.7040887          0.3240626          0.6081437
#> 3         0.02821837          0.4271713          0.1763610          0.4233146
```
