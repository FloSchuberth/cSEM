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
#>  Random seed                        = 1532025801
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
#>   eta2 ~ eta1      0.6713      0.0399   16.8205    0.0000 [ 0.5691; 0.7324 ] 
#>   eta3 ~ eta1      0.4585      0.0847    5.4131    0.0000 [ 0.3326; 0.6158 ] 
#>   eta3 ~ eta2      0.3052      0.0881    3.4622    0.0005 [ 0.1641; 0.4968 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0338   19.6227    0.0000 [ 0.6029; 0.7371 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9707    0.0000 [ 0.5849; 0.7159 ] 
#>   eta1 =~ y13      0.7613      0.0366   20.8302    0.0000 [ 0.6882; 0.8178 ] 
#>   eta2 =~ y21      0.5165      0.0465   11.1110    0.0000 [ 0.4351; 0.5967 ] 
#>   eta2 =~ y22      0.7554      0.0343   21.9952    0.0000 [ 0.7022; 0.8207 ] 
#>   eta2 =~ y23      0.7997      0.0356   22.4839    0.0000 [ 0.7325; 0.8532 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.6848    0.0000 [ 0.7653; 0.8735 ] 
#>   eta3 =~ y32      0.6581      0.0381   17.2657    0.0000 [ 0.6247; 0.7499 ] 
#>   eta3 =~ y33      0.7474      0.0460   16.2369    0.0000 [ 0.6541; 0.8201 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0193   20.5250    0.0000 [ 0.3715; 0.4357 ] 
#>   eta1 <~ y12      0.3873      0.0188   20.5658    0.0000 [ 0.3543; 0.4159 ] 
#>   eta1 <~ y13      0.4542      0.0201   22.5505    0.0000 [ 0.4220; 0.4938 ] 
#>   eta2 <~ y21      0.3058      0.0244   12.5461    0.0000 [ 0.2597; 0.3414 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0217    0.0000 [ 0.4107; 0.4974 ] 
#>   eta2 <~ y23      0.4735      0.0176   26.9153    0.0000 [ 0.4446; 0.5023 ] 
#>   eta3 <~ y31      0.4400      0.0196   22.4448    0.0000 [ 0.3941; 0.4670 ] 
#>   eta3 <~ y32      0.3521      0.0183   19.2022    0.0000 [ 0.3273; 0.3909 ] 
#>   eta3 <~ y33      0.3999      0.0212   18.8475    0.0000 [ 0.3518; 0.4301 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0399   16.8205    0.0000 [ 0.5691; 0.7324 ] 
#>   eta3 ~ eta1       0.6634      0.0435   15.2532    0.0000 [ 0.5851; 0.7330 ] 
#>   eta3 ~ eta2       0.3052      0.0881    3.4622    0.0005 [ 0.1641; 0.4968 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0546    3.7541    0.0002 [ 0.1095; 0.2841 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03379089 19.62274  9.887482e-86
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03825878 16.97069  1.353315e-64
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03655009 20.83020  2.304769e-96
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04648130 11.11102  1.108801e-28
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03434336 21.99516 3.203907e-107
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03556600 22.48394 5.961357e-112
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03331104 24.68483 1.556348e-134
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03811413 17.26575  8.519647e-67
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04603239 16.23692  2.764716e-59
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6029000          0.7370602
#> 2          0.5849439          0.7159462
#> 3          0.6881623          0.8177673
#> 4          0.4351167          0.5967080
#> 5          0.7021627          0.8207394
#> 6          0.7324527          0.8531890
#> 7          0.7653283          0.8734761
#> 8          0.6247465          0.7499256
#> 9          0.6540989          0.8201268

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
#>  Random seed                        = 1532025801
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
#>   eta2 ~ eta1      0.6713      0.0399   16.8205    0.0000 [ 0.5845; 0.7909 ] 
#>   eta3 ~ eta1      0.4585      0.0847    5.4131    0.0000 [ 0.2188; 0.6568 ] 
#>   eta3 ~ eta2      0.3052      0.0881    3.4622    0.0005 [ 0.0924; 0.5482 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0338   19.6227    0.0000 [ 0.5690; 0.7438 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9707    0.0000 [ 0.5531; 0.7509 ] 
#>   eta1 =~ y13      0.7613      0.0366   20.8302    0.0000 [ 0.6744; 0.8634 ] 
#>   eta2 =~ y21      0.5165      0.0465   11.1110    0.0000 [ 0.3958; 0.6362 ] 
#>   eta2 =~ y22      0.7554      0.0343   21.9952    0.0000 [ 0.6615; 0.8391 ] 
#>   eta2 =~ y23      0.7997      0.0356   22.4839    0.0000 [ 0.7052; 0.8891 ] 
#>   eta3 =~ y31      0.8223      0.0333   24.6848    0.0000 [ 0.7399; 0.9122 ] 
#>   eta3 =~ y32      0.6581      0.0381   17.2657    0.0000 [ 0.5414; 0.7385 ] 
#>   eta3 =~ y33      0.7474      0.0460   16.2369    0.0000 [ 0.6290; 0.8671 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0193   20.5250    0.0000 [ 0.3408; 0.4404 ] 
#>   eta1 <~ y12      0.3873      0.0188   20.5658    0.0000 [ 0.3395; 0.4369 ] 
#>   eta1 <~ y13      0.4542      0.0201   22.5505    0.0000 [ 0.4055; 0.5097 ] 
#>   eta2 <~ y21      0.3058      0.0244   12.5461    0.0000 [ 0.2447; 0.3708 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0217    0.0000 [ 0.3921; 0.5022 ] 
#>   eta2 <~ y23      0.4735      0.0176   26.9153    0.0000 [ 0.4297; 0.5207 ] 
#>   eta3 <~ y31      0.4400      0.0196   22.4448    0.0000 [ 0.3952; 0.4966 ] 
#>   eta3 <~ y32      0.3521      0.0183   19.2022    0.0000 [ 0.2984; 0.3932 ] 
#>   eta3 <~ y33      0.3999      0.0212   18.8475    0.0000 [ 0.3492; 0.4590 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0399   16.8205    0.0000 [ 0.5845; 0.7909 ] 
#>   eta3 ~ eta1       0.6634      0.0435   15.2532    0.0000 [ 0.5465; 0.7714 ] 
#>   eta3 ~ eta2       0.3052      0.0881    3.4622    0.0005 [ 0.0924; 0.5482 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0546    3.7541    0.0002 [ 0.0800; 0.3622 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03991158 16.820518 1.726437e-63
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08470341  5.413085 6.194807e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08813775  3.462207 5.357648e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5844980          0.7908981          0.6092826          0.7661135
#> 2          0.2187533          0.6567915          0.2713531          0.6041917
#> 3          0.0923742          0.5481729          0.1471067          0.4934403
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5685333          0.7361955          0.5691432          0.7324334
#> 2          0.2607821          0.6207545          0.3326326          0.6157712
#> 3          0.1324889          0.5370657          0.1640590          0.4968287
```
