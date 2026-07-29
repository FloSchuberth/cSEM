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
#>  Random seed                        = -810224398
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
#>   eta2 ~ eta1      0.6713      0.0426   15.7587    0.0000 [ 0.5943; 0.7492 ] 
#>   eta3 ~ eta1      0.4585      0.0953    4.8089    0.0000 [ 0.3104; 0.6539 ] 
#>   eta3 ~ eta2      0.3052      0.1017    3.0007    0.0027 [ 0.1046; 0.4820 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0426   15.5746    0.0000 [ 0.6026; 0.7289 ] 
#>   eta1 =~ y12      0.6493      0.0361   17.9901    0.0000 [ 0.5868; 0.7044 ] 
#>   eta1 =~ y13      0.7613      0.0308   24.6946    0.0000 [ 0.7025; 0.8183 ] 
#>   eta2 =~ y21      0.5165      0.0642    8.0427    0.0000 [ 0.4008; 0.6111 ] 
#>   eta2 =~ y22      0.7554      0.0377   20.0148    0.0000 [ 0.7093; 0.8342 ] 
#>   eta2 =~ y23      0.7997      0.0330   24.1976    0.0000 [ 0.7258; 0.8437 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5540    0.0000 [ 0.7588; 0.8704 ] 
#>   eta3 =~ y32      0.6581      0.0362   18.1889    0.0000 [ 0.5863; 0.7052 ] 
#>   eta3 =~ y33      0.7474      0.0384   19.4855    0.0000 [ 0.6714; 0.8028 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0210   18.8019    0.0000 [ 0.3679; 0.4296 ] 
#>   eta1 <~ y12      0.3873      0.0183   21.1473    0.0000 [ 0.3467; 0.4161 ] 
#>   eta1 <~ y13      0.4542      0.0219   20.7139    0.0000 [ 0.4173; 0.4888 ] 
#>   eta2 <~ y21      0.3058      0.0317    9.6512    0.0000 [ 0.2475; 0.3548 ] 
#>   eta2 <~ y22      0.4473      0.0221   20.2547    0.0000 [ 0.4093; 0.5000 ] 
#>   eta2 <~ y23      0.4735      0.0231   20.5360    0.0000 [ 0.4334; 0.5109 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5957    0.0000 [ 0.4089; 0.4681 ] 
#>   eta3 <~ y32      0.3521      0.0181   19.4836    0.0000 [ 0.3108; 0.3864 ] 
#>   eta3 <~ y33      0.3999      0.0144   27.7759    0.0000 [ 0.3694; 0.4265 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0426   15.7587    0.0000 [ 0.5943; 0.7492 ] 
#>   eta3 ~ eta1       0.6634      0.0381   17.3931    0.0000 [ 0.5952; 0.7268 ] 
#>   eta3 ~ eta2       0.3052      0.1017    3.0007    0.0027 [ 0.1046; 0.4820 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0655    3.1284    0.0018 [ 0.0682; 0.3112 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04257392 15.574556  1.084034e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03609082 17.990114  2.328793e-72
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03083048 24.694581 1.222913e-134
#> 4 eta2 =~ y21  Common factor 0.5164548 0.06421405  8.042707  8.787511e-16
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03774149 20.014781  4.094311e-89
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03304724 24.197593 2.358374e-129
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03096621 26.554023 2.307640e-155
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03617974 18.188878  6.322165e-74
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03835804 19.485464  1.458525e-84
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6026462          0.7289272
#> 2          0.5868114          0.7044049
#> 3          0.7024611          0.8182601
#> 4          0.4008349          0.6111213
#> 5          0.7092720          0.8342248
#> 6          0.7258267          0.8436656
#> 7          0.7587953          0.8704267
#> 8          0.5863311          0.7052039
#> 9          0.6713788          0.8028250

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
#>  Random seed                        = -810224398
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
#>   eta2 ~ eta1      0.6713      0.0426   15.7587    0.0000 [ 0.5638; 0.7841 ] 
#>   eta3 ~ eta1      0.4585      0.0953    4.8089    0.0000 [ 0.2116; 0.7047 ] 
#>   eta3 ~ eta2      0.3052      0.1017    3.0007    0.0027 [ 0.0410; 0.5669 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0426   15.5746    0.0000 [ 0.5513; 0.7714 ] 
#>   eta1 =~ y12      0.6493      0.0361   17.9901    0.0000 [ 0.5618; 0.7484 ] 
#>   eta1 =~ y13      0.7613      0.0308   24.6946    0.0000 [ 0.6851; 0.8445 ] 
#>   eta2 =~ y21      0.5165      0.0642    8.0427    0.0000 [ 0.3580; 0.6901 ] 
#>   eta2 =~ y22      0.7554      0.0377   20.0148    0.0000 [ 0.6502; 0.8453 ] 
#>   eta2 =~ y23      0.7997      0.0330   24.1976    0.0000 [ 0.7172; 0.8881 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5540    0.0000 [ 0.7492; 0.9093 ] 
#>   eta3 =~ y32      0.6581      0.0362   18.1889    0.0000 [ 0.5667; 0.7538 ] 
#>   eta3 =~ y33      0.7474      0.0384   19.4855    0.0000 [ 0.6517; 0.8501 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0210   18.8019    0.0000 [ 0.3384; 0.4472 ] 
#>   eta1 <~ y12      0.3873      0.0183   21.1473    0.0000 [ 0.3416; 0.4364 ] 
#>   eta1 <~ y13      0.4542      0.0219   20.7139    0.0000 [ 0.3971; 0.5105 ] 
#>   eta2 <~ y21      0.3058      0.0317    9.6512    0.0000 [ 0.2290; 0.3929 ] 
#>   eta2 <~ y22      0.4473      0.0221   20.2547    0.0000 [ 0.3855; 0.4997 ] 
#>   eta2 <~ y23      0.4735      0.0231   20.5360    0.0000 [ 0.4154; 0.5346 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5957    0.0000 [ 0.3895; 0.4902 ] 
#>   eta3 <~ y32      0.3521      0.0181   19.4836    0.0000 [ 0.3036; 0.3971 ] 
#>   eta3 <~ y33      0.3999      0.0144   27.7759    0.0000 [ 0.3614; 0.4359 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0426   15.7587    0.0000 [ 0.5638; 0.7841 ] 
#>   eta3 ~ eta1       0.6634      0.0381   17.3931    0.0000 [ 0.5655; 0.7627 ] 
#>   eta3 ~ eta2       0.3052      0.1017    3.0007    0.0027 [ 0.0410; 0.5669 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0655    3.1284    0.0018 [ 0.0366; 0.3752 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04260069 15.758745 5.980934e-56
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09534466  4.808940 1.517327e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10169442  3.000667 2.693887e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56377708          0.7840838          0.5902316          0.7576293
#> 2         0.21164350          0.7047123          0.2708514          0.6455043
#> 3         0.04100787          0.5669140          0.1041589          0.5037630
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.56165537          0.7790995          0.5942619          0.7492062
#> 2         0.28899441          0.6852258          0.3103687          0.6538973
#> 3         0.08522016          0.4927913          0.1045784          0.4820052
```
