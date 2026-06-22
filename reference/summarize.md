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
#>  Random seed                        = 2009094496
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
#>   eta2 ~ eta1      0.6713      0.0353   18.9916    0.0000 [ 0.6069; 0.7205 ] 
#>   eta3 ~ eta1      0.4585      0.0784    5.8461    0.0000 [ 0.3223; 0.5821 ] 
#>   eta3 ~ eta2      0.3052      0.0873    3.4957    0.0005 [ 0.1745; 0.4516 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0466   14.2342    0.0000 [ 0.5954; 0.7561 ] 
#>   eta1 =~ y12      0.6493      0.0376   17.2884    0.0000 [ 0.5893; 0.7151 ] 
#>   eta1 =~ y13      0.7613      0.0358   21.2841    0.0000 [ 0.6890; 0.8133 ] 
#>   eta2 =~ y21      0.5165      0.0556    9.2822    0.0000 [ 0.4241; 0.6074 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1433    0.0000 [ 0.6836; 0.8198 ] 
#>   eta2 =~ y23      0.7997      0.0332   24.0922    0.0000 [ 0.7223; 0.8449 ] 
#>   eta3 =~ y31      0.8223      0.0409   20.1090    0.0000 [ 0.7573; 0.8945 ] 
#>   eta3 =~ y32      0.6581      0.0378   17.4078    0.0000 [ 0.5847; 0.7177 ] 
#>   eta3 =~ y33      0.7474      0.0363   20.5978    0.0000 [ 0.6806; 0.7995 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0209   18.8956    0.0000 [ 0.3656; 0.4338 ] 
#>   eta1 <~ y12      0.3873      0.0175   22.1377    0.0000 [ 0.3552; 0.4154 ] 
#>   eta1 <~ y13      0.4542      0.0240   18.9112    0.0000 [ 0.4139; 0.4954 ] 
#>   eta2 <~ y21      0.3058      0.0266   11.4966    0.0000 [ 0.2608; 0.3516 ] 
#>   eta2 <~ y22      0.4473      0.0225   19.8426    0.0000 [ 0.4104; 0.4878 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3549    0.0000 [ 0.4205; 0.5059 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.8543    0.0000 [ 0.4134; 0.4732 ] 
#>   eta3 <~ y32      0.3521      0.0195   18.0726    0.0000 [ 0.3121; 0.3867 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4447    0.0000 [ 0.3679; 0.4365 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0353   18.9916    0.0000 [ 0.6069; 0.7205 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7354    0.0000 [ 0.5950; 0.7142 ] 
#>   eta3 ~ eta2       0.3052      0.0873    3.4957    0.0005 [ 0.1745; 0.4516 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0588    3.4855    0.0005 [ 0.1200; 0.2978 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04658300 14.234160  5.623439e-46
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03755568 17.288409  5.751790e-67
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03577056 21.284146 1.592223e-100
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05563922  9.282208  1.660027e-20
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03572701 21.143324  3.178844e-99
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03319178 24.092220 3.016115e-128
#> 7 eta3 =~ y31  Common factor 0.8222773 0.04089111 20.108954  6.160984e-90
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03780303 17.407834  7.195473e-68
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03628662 20.597790  2.872561e-94
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5953830          0.7560797
#> 2          0.5892515          0.7150632
#> 3          0.6889715          0.8133090
#> 4          0.4240660          0.6074093
#> 5          0.6835600          0.8198458
#> 6          0.7223420          0.8448911
#> 7          0.7573344          0.8944525
#> 8          0.5847392          0.7176640
#> 9          0.6806197          0.7995315

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
#>  Random seed                        = 2009094496
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
#>   eta2 ~ eta1      0.6713      0.0353   18.9916    0.0000 [ 0.5830; 0.7658 ] 
#>   eta3 ~ eta1      0.4585      0.0784    5.8461    0.0000 [ 0.2632; 0.6688 ] 
#>   eta3 ~ eta2      0.3052      0.0873    3.4957    0.0005 [ 0.0736; 0.5251 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0466   14.2342    0.0000 [ 0.5394; 0.7803 ] 
#>   eta1 =~ y12      0.6493      0.0376   17.2884    0.0000 [ 0.5513; 0.7455 ] 
#>   eta1 =~ y13      0.7613      0.0358   21.2841    0.0000 [ 0.6669; 0.8519 ] 
#>   eta2 =~ y21      0.5165      0.0556    9.2822    0.0000 [ 0.3705; 0.6583 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1433    0.0000 [ 0.6612; 0.8460 ] 
#>   eta2 =~ y23      0.7997      0.0332   24.0922    0.0000 [ 0.7147; 0.8864 ] 
#>   eta3 =~ y31      0.8223      0.0409   20.1090    0.0000 [ 0.7149; 0.9264 ] 
#>   eta3 =~ y32      0.6581      0.0378   17.4078    0.0000 [ 0.5645; 0.7600 ] 
#>   eta3 =~ y33      0.7474      0.0363   20.5978    0.0000 [ 0.6599; 0.8476 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0209   18.8956    0.0000 [ 0.3416; 0.4499 ] 
#>   eta1 <~ y12      0.3873      0.0175   22.1377    0.0000 [ 0.3435; 0.4340 ] 
#>   eta1 <~ y13      0.4542      0.0240   18.9112    0.0000 [ 0.3928; 0.5170 ] 
#>   eta2 <~ y21      0.3058      0.0266   11.4966    0.0000 [ 0.2373; 0.3748 ] 
#>   eta2 <~ y22      0.4473      0.0225   19.8426    0.0000 [ 0.3891; 0.5057 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3549    0.0000 [ 0.4179; 0.5325 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.8543    0.0000 [ 0.3910; 0.4826 ] 
#>   eta3 <~ y32      0.3521      0.0195   18.0726    0.0000 [ 0.3020; 0.4028 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4447    0.0000 [ 0.3505; 0.4516 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0353   18.9916    0.0000 [ 0.5830; 0.7658 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7354    0.0000 [ 0.5766; 0.7597 ] 
#>   eta3 ~ eta2       0.3052      0.0873    3.4957    0.0005 [ 0.0736; 0.5251 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0588    3.4855    0.0005 [ 0.0501; 0.3541 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03534896 18.991604 2.001206e-80
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07842945  5.846104 5.032181e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08729399  3.495672 4.728701e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.58304188          0.7658467          0.6049932          0.7438954
#> 2         0.26323030          0.6688231          0.3119341          0.6201194
#> 3         0.07362163          0.5250569          0.1278302          0.4708483
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5980143          0.7290205          0.6069120          0.7204667
#> 2          0.3220523          0.6657716          0.3222788          0.5820684
#> 3          0.1129429          0.4580834          0.1744578          0.4516051
```
