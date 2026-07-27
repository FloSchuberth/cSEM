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
#>  Random seed                        = -232404463
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
#>   eta2 ~ eta1      0.6713      0.0393   17.0805    0.0000 [ 0.6041; 0.7375 ] 
#>   eta3 ~ eta1      0.4585      0.0863    5.3127    0.0000 [ 0.3077; 0.6126 ] 
#>   eta3 ~ eta2      0.3052      0.0945    3.2286    0.0012 [ 0.1449; 0.4727 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0381   17.4069    0.0000 [ 0.5928; 0.7291 ] 
#>   eta1 =~ y12      0.6493      0.0431   15.0702    0.0000 [ 0.5647; 0.7242 ] 
#>   eta1 =~ y13      0.7613      0.0318   23.9319    0.0000 [ 0.6784; 0.8058 ] 
#>   eta2 =~ y21      0.5165      0.0509   10.1517    0.0000 [ 0.4184; 0.6048 ] 
#>   eta2 =~ y22      0.7554      0.0342   22.1104    0.0000 [ 0.6916; 0.8096 ] 
#>   eta2 =~ y23      0.7997      0.0313   25.5612    0.0000 [ 0.7498; 0.8526 ] 
#>   eta3 =~ y31      0.8223      0.0311   26.4634    0.0000 [ 0.7546; 0.8787 ] 
#>   eta3 =~ y32      0.6581      0.0316   20.8199    0.0000 [ 0.6019; 0.7025 ] 
#>   eta3 =~ y33      0.7474      0.0353   21.1479    0.0000 [ 0.6903; 0.8209 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0237   16.6589    0.0000 [ 0.3660; 0.4531 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5733    0.0000 [ 0.3526; 0.4238 ] 
#>   eta1 <~ y13      0.4542      0.0199   22.7976    0.0000 [ 0.4166; 0.4867 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3795    0.0000 [ 0.2687; 0.3434 ] 
#>   eta2 <~ y22      0.4473      0.0191   23.3845    0.0000 [ 0.4137; 0.4792 ] 
#>   eta2 <~ y23      0.4735      0.0184   25.7284    0.0000 [ 0.4463; 0.5051 ] 
#>   eta3 <~ y31      0.4400      0.0186   23.6064    0.0000 [ 0.3985; 0.4680 ] 
#>   eta3 <~ y32      0.3521      0.0150   23.4180    0.0000 [ 0.3275; 0.3781 ] 
#>   eta3 <~ y33      0.3999      0.0153   26.1651    0.0000 [ 0.3765; 0.4277 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0393   17.0805    0.0000 [ 0.6041; 0.7375 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5413    0.0000 [ 0.5875; 0.7237 ] 
#>   eta3 ~ eta2       0.3052      0.0945    3.2286    0.0012 [ 0.1449; 0.4727 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0615    3.3336    0.0009 [ 0.0972; 0.3088 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03809240 17.40688  7.316465e-68
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04308348 15.07023  2.542263e-51
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03181298 23.93192 1.425518e-126
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05087374 10.15170  3.256516e-24
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03416439 22.11038 2.511267e-108
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03128431 25.56118 4.124830e-144
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03107224 26.46341 2.557757e-154
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03160774 20.81986  2.859901e-96
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03534279 21.14785  2.887991e-99
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5928418          0.7290964
#> 2          0.5647402          0.7241854
#> 3          0.6784332          0.8057939
#> 4          0.4184052          0.6047535
#> 5          0.6915508          0.8096188
#> 6          0.7497910          0.8525737
#> 7          0.7545941          0.8787074
#> 8          0.6019462          0.7024804
#> 9          0.6902639          0.8208732

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
#>  Random seed                        = -232404463
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
#>   eta2 ~ eta1      0.6713      0.0393   17.0805    0.0000 [ 0.5732; 0.7765 ] 
#>   eta3 ~ eta1      0.4585      0.0863    5.3127    0.0000 [ 0.2281; 0.6744 ] 
#>   eta3 ~ eta2      0.3052      0.0945    3.2286    0.0012 [ 0.0637; 0.5525 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0381   17.4069    0.0000 [ 0.5544; 0.7514 ] 
#>   eta1 =~ y12      0.6493      0.0431   15.0702    0.0000 [ 0.5370; 0.7598 ] 
#>   eta1 =~ y13      0.7613      0.0318   23.9319    0.0000 [ 0.6866; 0.8511 ] 
#>   eta2 =~ y21      0.5165      0.0509   10.1517    0.0000 [ 0.3826; 0.6457 ] 
#>   eta2 =~ y22      0.7554      0.0342   22.1104    0.0000 [ 0.6728; 0.8495 ] 
#>   eta2 =~ y23      0.7997      0.0313   25.5612    0.0000 [ 0.7132; 0.8750 ] 
#>   eta3 =~ y31      0.8223      0.0311   26.4634    0.0000 [ 0.7461; 0.9068 ] 
#>   eta3 =~ y32      0.6581      0.0316   20.8199    0.0000 [ 0.5770; 0.7405 ] 
#>   eta3 =~ y33      0.7474      0.0353   21.1479    0.0000 [ 0.6470; 0.8298 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0237   16.6589    0.0000 [ 0.3292; 0.4520 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5733    0.0000 [ 0.3343; 0.4421 ] 
#>   eta1 <~ y13      0.4542      0.0199   22.7976    0.0000 [ 0.4084; 0.5114 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3795    0.0000 [ 0.2416; 0.3693 ] 
#>   eta2 <~ y22      0.4473      0.0191   23.3845    0.0000 [ 0.4020; 0.5010 ] 
#>   eta2 <~ y23      0.4735      0.0184   25.7284    0.0000 [ 0.4235; 0.5186 ] 
#>   eta3 <~ y31      0.4400      0.0186   23.6064    0.0000 [ 0.3953; 0.4917 ] 
#>   eta3 <~ y32      0.3521      0.0150   23.4180    0.0000 [ 0.3148; 0.3925 ] 
#>   eta3 <~ y33      0.3999      0.0153   26.1651    0.0000 [ 0.3570; 0.4360 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0393   17.0805    0.0000 [ 0.5732; 0.7765 ] 
#>   eta3 ~ eta1       0.6634      0.0358   18.5413    0.0000 [ 0.5675; 0.7525 ] 
#>   eta3 ~ eta2       0.3052      0.0945    3.2286    0.0012 [ 0.0637; 0.5525 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0615    3.3336    0.0009 [ 0.0499; 0.3677 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03930404 17.080518 2.072723e-65
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08630329  5.312738 1.079902e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09451583  3.228571 1.244102e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.57321198          0.7764703          0.5976193          0.7520630
#> 2         0.22808853          0.6744004          0.2816819          0.6208071
#> 3         0.06371628          0.5524988          0.1224095          0.4938056
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5836398          0.7456938          0.6041176          0.7374793
#> 2          0.2747954          0.6162851          0.3076839          0.6126183
#> 3          0.1223480          0.4953736          0.1448694          0.4727208
```
