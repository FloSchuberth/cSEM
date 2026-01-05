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
#>  Random seed                        = -1398550464
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
#>   eta2 ~ eta1      0.6713      0.0447   15.0327    0.0000 [ 0.5975; 0.7564 ] 
#>   eta3 ~ eta1      0.4585      0.0803    5.7070    0.0000 [ 0.3358; 0.5753 ] 
#>   eta3 ~ eta2      0.3052      0.0840    3.6319    0.0003 [ 0.1695; 0.4232 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0430   15.4361    0.0000 [ 0.5834; 0.7334 ] 
#>   eta1 =~ y12      0.6493      0.0323   20.1323    0.0000 [ 0.5881; 0.6944 ] 
#>   eta1 =~ y13      0.7613      0.0351   21.6903    0.0000 [ 0.6967; 0.8252 ] 
#>   eta2 =~ y21      0.5165      0.0540    9.5609    0.0000 [ 0.3995; 0.5962 ] 
#>   eta2 =~ y22      0.7554      0.0360   20.9690    0.0000 [ 0.6896; 0.8206 ] 
#>   eta2 =~ y23      0.7997      0.0396   20.1817    0.0000 [ 0.7256; 0.8539 ] 
#>   eta3 =~ y31      0.8223      0.0469   17.5446    0.0000 [ 0.7339; 0.8898 ] 
#>   eta3 =~ y32      0.6581      0.0381   17.2603    0.0000 [ 0.5920; 0.7314 ] 
#>   eta3 =~ y33      0.7474      0.0286   26.1266    0.0000 [ 0.6962; 0.8001 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0185   21.4335    0.0000 [ 0.3554; 0.4306 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5597    0.0000 [ 0.3567; 0.4305 ] 
#>   eta1 <~ y13      0.4542      0.0206   21.9959    0.0000 [ 0.4218; 0.4953 ] 
#>   eta2 <~ y21      0.3058      0.0280   10.9218    0.0000 [ 0.2556; 0.3563 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0016    0.0000 [ 0.4234; 0.4923 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1604    0.0000 [ 0.4392; 0.5181 ] 
#>   eta3 <~ y31      0.4400      0.0221   19.8833    0.0000 [ 0.3987; 0.4759 ] 
#>   eta3 <~ y32      0.3521      0.0162   21.7559    0.0000 [ 0.3155; 0.3801 ] 
#>   eta3 <~ y33      0.3999      0.0199   20.0766    0.0000 [ 0.3695; 0.4375 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0447   15.0327    0.0000 [ 0.5975; 0.7564 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6850    0.0000 [ 0.5920; 0.7361 ] 
#>   eta3 ~ eta2       0.3052      0.0840    3.6319    0.0003 [ 0.1695; 0.4232 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0556    3.6863    0.0002 [ 0.1064; 0.2924 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04295571 15.436128  9.355057e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03225052 20.132324  3.845331e-90
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03510074 21.690304 2.533080e-104
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05401732  9.560911  1.167247e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03602396 20.969036  1.257864e-97
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03962327 20.181668  1.418783e-90
#> 7 eta3 =~ y31  Common factor 0.8222773 0.04686771 17.544646  6.535887e-69
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03812626 17.260254  9.369833e-67
#> 9 eta3 =~ y33  Common factor 0.7474241 0.02860782 26.126563 1.820159e-150
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5834019          0.7334142
#> 2          0.5880596          0.6944044
#> 3          0.6967407          0.8252439
#> 4          0.3994614          0.5962338
#> 5          0.6896425          0.8205686
#> 6          0.7256037          0.8538837
#> 7          0.7338546          0.8897513
#> 8          0.5920483          0.7313699
#> 9          0.6962115          0.8000632

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
#>  Random seed                        = -1398550464
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
#>   eta2 ~ eta1      0.6713      0.0447   15.0327    0.0000 [ 0.5592; 0.7901 ] 
#>   eta3 ~ eta1      0.4585      0.0803    5.7070    0.0000 [ 0.2406; 0.6561 ] 
#>   eta3 ~ eta2      0.3052      0.0840    3.6319    0.0003 [ 0.0940; 0.5285 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0430   15.4361    0.0000 [ 0.5486; 0.7707 ] 
#>   eta1 =~ y12      0.6493      0.0323   20.1323    0.0000 [ 0.5692; 0.7360 ] 
#>   eta1 =~ y13      0.7613      0.0351   21.6903    0.0000 [ 0.6700; 0.8516 ] 
#>   eta2 =~ y21      0.5165      0.0540    9.5609    0.0000 [ 0.3950; 0.6743 ] 
#>   eta2 =~ y22      0.7554      0.0360   20.9690    0.0000 [ 0.6600; 0.8463 ] 
#>   eta2 =~ y23      0.7997      0.0396   20.1817    0.0000 [ 0.7028; 0.9077 ] 
#>   eta3 =~ y31      0.8223      0.0469   17.5446    0.0000 [ 0.7067; 0.9491 ] 
#>   eta3 =~ y32      0.6581      0.0381   17.2603    0.0000 [ 0.5545; 0.7517 ] 
#>   eta3 =~ y33      0.7474      0.0286   26.1266    0.0000 [ 0.6690; 0.8170 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0185   21.4335    0.0000 [ 0.3465; 0.4419 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5597    0.0000 [ 0.3356; 0.4435 ] 
#>   eta1 <~ y13      0.4542      0.0206   21.9959    0.0000 [ 0.4008; 0.5076 ] 
#>   eta2 <~ y21      0.3058      0.0280   10.9218    0.0000 [ 0.2408; 0.3856 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0016    0.0000 [ 0.3850; 0.4952 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1604    0.0000 [ 0.4099; 0.5314 ] 
#>   eta3 <~ y31      0.4400      0.0221   19.8833    0.0000 [ 0.3871; 0.5015 ] 
#>   eta3 <~ y32      0.3521      0.0162   21.7559    0.0000 [ 0.3087; 0.3924 ] 
#>   eta3 <~ y33      0.3999      0.0199   20.0766    0.0000 [ 0.3469; 0.4499 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0447   15.0327    0.0000 [ 0.5592; 0.7901 ] 
#>   eta3 ~ eta1       0.6634      0.0398   16.6850    0.0000 [ 0.5562; 0.7618 ] 
#>   eta3 ~ eta2       0.3052      0.0840    3.6319    0.0003 [ 0.0940; 0.5285 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0556    3.6863    0.0002 [ 0.0669; 0.3543 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04465826 15.032682 4.484756e-51
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08034133  5.706985 1.149949e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08402039  3.631870 2.813749e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55918275          0.7901301          0.5869150          0.7623978
#> 2         0.24061461          0.6560946          0.2905056          0.6062036
#> 3         0.09395541          0.5284614          0.1461311          0.4762857
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5553894          0.7637034          0.5975037          0.7564060
#> 2          0.2544898          0.6047921          0.3358256          0.5752793
#> 3          0.1658780          0.4944772          0.1695091          0.4231867
```
