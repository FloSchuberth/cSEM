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
#>  Random seed                        = -1953414189
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
#>   eta2 ~ eta1      0.6713      0.0447   15.0152    0.0000 [ 0.5806; 0.7381 ] 
#>   eta3 ~ eta1      0.4585      0.0791    5.7983    0.0000 [ 0.3345; 0.6006 ] 
#>   eta3 ~ eta2      0.3052      0.0825    3.6996    0.0002 [ 0.1568; 0.4270 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0357   18.5795    0.0000 [ 0.5954; 0.7203 ] 
#>   eta1 =~ y12      0.6493      0.0349   18.6071    0.0000 [ 0.5511; 0.6853 ] 
#>   eta1 =~ y13      0.7613      0.0292   26.0498    0.0000 [ 0.7231; 0.8261 ] 
#>   eta2 =~ y21      0.5165      0.0482   10.7242    0.0000 [ 0.4079; 0.5869 ] 
#>   eta2 =~ y22      0.7554      0.0334   22.5871    0.0000 [ 0.6903; 0.8035 ] 
#>   eta2 =~ y23      0.7997      0.0403   19.8244    0.0000 [ 0.7305; 0.8695 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.7205    0.0000 [ 0.7772; 0.8818 ] 
#>   eta3 =~ y32      0.6581      0.0379   17.3600    0.0000 [ 0.5757; 0.7183 ] 
#>   eta3 =~ y33      0.7474      0.0400   18.6742    0.0000 [ 0.6728; 0.8188 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0191   20.7205    0.0000 [ 0.3673; 0.4284 ] 
#>   eta1 <~ y12      0.3873      0.0170   22.7775    0.0000 [ 0.3449; 0.4128 ] 
#>   eta1 <~ y13      0.4542      0.0158   28.6640    0.0000 [ 0.4336; 0.4930 ] 
#>   eta2 <~ y21      0.3058      0.0268   11.4115    0.0000 [ 0.2517; 0.3453 ] 
#>   eta2 <~ y22      0.4473      0.0167   26.7216    0.0000 [ 0.4222; 0.4796 ] 
#>   eta2 <~ y23      0.4735      0.0247   19.2052    0.0000 [ 0.4428; 0.5284 ] 
#>   eta3 <~ y31      0.4400      0.0188   23.3615    0.0000 [ 0.4148; 0.4791 ] 
#>   eta3 <~ y32      0.3521      0.0183   19.2452    0.0000 [ 0.3195; 0.3822 ] 
#>   eta3 <~ y33      0.3999      0.0202   19.7541    0.0000 [ 0.3655; 0.4378 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0447   15.0152    0.0000 [ 0.5806; 0.7381 ] 
#>   eta3 ~ eta1       0.6634      0.0372   17.8320    0.0000 [ 0.6095; 0.7351 ] 
#>   eta3 ~ eta2       0.3052      0.0825    3.6996    0.0002 [ 0.1568; 0.4270 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0556    3.6849    0.0002 [ 0.1068; 0.2956 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03568831 18.57947  4.711177e-77
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03489407 18.60711  2.813843e-77
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02922651 26.04984 1.351029e-149
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04815807 10.72416  7.839695e-27
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03344338 22.58706 5.808200e-113
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04033727 19.82444  1.832170e-87
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03326297 24.72051 6.438094e-135
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03790715 17.36002  1.656593e-67
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04002450 18.67416  8.033901e-78
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5953629          0.7202853
#> 2          0.5511344          0.6852992
#> 3          0.7230949          0.8260652
#> 4          0.4079309          0.5869394
#> 5          0.6903279          0.8034751
#> 6          0.7305474          0.8694777
#> 7          0.7772268          0.8818470
#> 8          0.5756953          0.7183392
#> 9          0.6727704          0.8188450

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
#>  Random seed                        = -1953414189
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
#>   eta2 ~ eta1      0.6713      0.0447   15.0152    0.0000 [ 0.5631; 0.7943 ] 
#>   eta3 ~ eta1      0.4585      0.0791    5.7983    0.0000 [ 0.2492; 0.6582 ] 
#>   eta3 ~ eta2      0.3052      0.0825    3.6996    0.0002 [ 0.0959; 0.5225 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0357   18.5795    0.0000 [ 0.5705; 0.7551 ] 
#>   eta1 =~ y12      0.6493      0.0349   18.6071    0.0000 [ 0.5680; 0.7485 ] 
#>   eta1 =~ y13      0.7613      0.0292   26.0498    0.0000 [ 0.6765; 0.8276 ] 
#>   eta2 =~ y21      0.5165      0.0482   10.7242    0.0000 [ 0.4073; 0.6564 ] 
#>   eta2 =~ y22      0.7554      0.0334   22.5871    0.0000 [ 0.6776; 0.8506 ] 
#>   eta2 =~ y23      0.7997      0.0403   19.8244    0.0000 [ 0.6900; 0.8986 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.7205    0.0000 [ 0.7255; 0.8975 ] 
#>   eta3 =~ y32      0.6581      0.0379   17.3600    0.0000 [ 0.5692; 0.7652 ] 
#>   eta3 =~ y33      0.7474      0.0400   18.6742    0.0000 [ 0.6456; 0.8526 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0191   20.7205    0.0000 [ 0.3465; 0.4452 ] 
#>   eta1 <~ y12      0.3873      0.0170   22.7775    0.0000 [ 0.3492; 0.4372 ] 
#>   eta1 <~ y13      0.4542      0.0158   28.6640    0.0000 [ 0.4081; 0.4900 ] 
#>   eta2 <~ y21      0.3058      0.0268   11.4115    0.0000 [ 0.2425; 0.3811 ] 
#>   eta2 <~ y22      0.4473      0.0167   26.7216    0.0000 [ 0.4043; 0.4909 ] 
#>   eta2 <~ y23      0.4735      0.0247   19.2052    0.0000 [ 0.4012; 0.5287 ] 
#>   eta3 <~ y31      0.4400      0.0188   23.3615    0.0000 [ 0.3858; 0.4832 ] 
#>   eta3 <~ y32      0.3521      0.0183   19.2452    0.0000 [ 0.3101; 0.4047 ] 
#>   eta3 <~ y33      0.3999      0.0202   19.7541    0.0000 [ 0.3489; 0.4536 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0447   15.0152    0.0000 [ 0.5631; 0.7943 ] 
#>   eta3 ~ eta1       0.6634      0.0372   17.8320    0.0000 [ 0.5676; 0.7600 ] 
#>   eta3 ~ eta2       0.3052      0.0825    3.6996    0.0002 [ 0.0959; 0.5225 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0556    3.6849    0.0002 [ 0.0663; 0.3538 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04471020 15.015218 5.836975e-51
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07907647  5.798271 6.700231e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08248198  3.699610 2.159312e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56307155          0.7942875          0.5908361          0.7665230
#> 2         0.24922303          0.6581619          0.2983286          0.6090563
#> 3         0.09593283          0.5224830          0.1471532          0.4712627
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5700227          0.7660965          0.5805873          0.7381223
#> 2          0.3236887          0.6191947          0.3345138          0.6006128
#> 3          0.1292012          0.4308051          0.1567785          0.4269993
```
