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
#>  Random seed                        = 2130309188
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
#>   eta2 ~ eta1      0.6713      0.0408   16.4481    0.0000 [ 0.5730; 0.7210 ] 
#>   eta3 ~ eta1      0.4585      0.0793    5.7796    0.0000 [ 0.3538; 0.6421 ] 
#>   eta3 ~ eta2      0.3052      0.0868    3.5157    0.0004 [ 0.1364; 0.4568 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0365   18.1454    0.0000 [ 0.5889; 0.7248 ] 
#>   eta1 =~ y12      0.6493      0.0350   18.5742    0.0000 [ 0.5842; 0.7005 ] 
#>   eta1 =~ y13      0.7613      0.0267   28.5117    0.0000 [ 0.7029; 0.7928 ] 
#>   eta2 =~ y21      0.5165      0.0518    9.9664    0.0000 [ 0.4368; 0.6095 ] 
#>   eta2 =~ y22      0.7554      0.0348   21.7022    0.0000 [ 0.6952; 0.8068 ] 
#>   eta2 =~ y23      0.7997      0.0388   20.6253    0.0000 [ 0.7096; 0.8437 ] 
#>   eta3 =~ y31      0.8223      0.0236   34.8817    0.0000 [ 0.7777; 0.8643 ] 
#>   eta3 =~ y32      0.6581      0.0470   13.9948    0.0000 [ 0.5587; 0.7360 ] 
#>   eta3 =~ y33      0.7474      0.0357   20.9654    0.0000 [ 0.6821; 0.8260 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0155   25.4523    0.0000 [ 0.3692; 0.4208 ] 
#>   eta1 <~ y12      0.3873      0.0158   24.5672    0.0000 [ 0.3645; 0.4215 ] 
#>   eta1 <~ y13      0.4542      0.0219   20.7168    0.0000 [ 0.4201; 0.4954 ] 
#>   eta2 <~ y21      0.3058      0.0292   10.4873    0.0000 [ 0.2585; 0.3556 ] 
#>   eta2 <~ y22      0.4473      0.0207   21.5775    0.0000 [ 0.4121; 0.4855 ] 
#>   eta2 <~ y23      0.4735      0.0218   21.6949    0.0000 [ 0.4420; 0.5193 ] 
#>   eta3 <~ y31      0.4400      0.0161   27.3224    0.0000 [ 0.4082; 0.4626 ] 
#>   eta3 <~ y32      0.3521      0.0195   18.0297    0.0000 [ 0.3153; 0.3814 ] 
#>   eta3 <~ y33      0.3999      0.0186   21.4986    0.0000 [ 0.3750; 0.4459 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0408   16.4481    0.0000 [ 0.5730; 0.7210 ] 
#>   eta3 ~ eta1       0.6634      0.0429   15.4497    0.0000 [ 0.5925; 0.7485 ] 
#>   eta3 ~ eta2       0.3052      0.0868    3.5157    0.0004 [ 0.1364; 0.4568 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0550    3.7223    0.0002 [ 0.0853; 0.2763 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03654214 18.145350  1.397420e-73
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03495586 18.574221  5.195314e-77
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02670295 28.511672 8.395353e-179
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05181935  9.966446  2.137397e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03480701 21.702173 1.956949e-104
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03877104 20.625285  1.627679e-94
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02357329 34.881734 1.406859e-266
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04702224 13.994844  1.675980e-44
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03565038 20.965390  1.358046e-97
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5888571          0.7248196
#> 2          0.5841606          0.7004981
#> 3          0.7028842          0.7928161
#> 4          0.4368357          0.6095396
#> 5          0.6951724          0.8067690
#> 6          0.7095774          0.8436540
#> 7          0.7776933          0.8642796
#> 8          0.5587215          0.7360335
#> 9          0.6821382          0.8260232

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
#>  Random seed                        = 2130309188
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
#>   eta2 ~ eta1      0.6713      0.0408   16.4481    0.0000 [ 0.5716; 0.7827 ] 
#>   eta3 ~ eta1      0.4585      0.0793    5.7796    0.0000 [ 0.2428; 0.6531 ] 
#>   eta3 ~ eta2      0.3052      0.0868    3.5157    0.0004 [ 0.0887; 0.5375 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0365   18.1454    0.0000 [ 0.5722; 0.7612 ] 
#>   eta1 =~ y12      0.6493      0.0350   18.5742    0.0000 [ 0.5536; 0.7344 ] 
#>   eta1 =~ y13      0.7613      0.0267   28.5117    0.0000 [ 0.6989; 0.8370 ] 
#>   eta2 =~ y21      0.5165      0.0518    9.9664    0.0000 [ 0.3904; 0.6584 ] 
#>   eta2 =~ y22      0.7554      0.0348   21.7022    0.0000 [ 0.6700; 0.8500 ] 
#>   eta2 =~ y23      0.7997      0.0388   20.6253    0.0000 [ 0.7059; 0.9064 ] 
#>   eta3 =~ y31      0.8223      0.0236   34.8817    0.0000 [ 0.7643; 0.8862 ] 
#>   eta3 =~ y32      0.6581      0.0470   13.9948    0.0000 [ 0.5276; 0.7708 ] 
#>   eta3 =~ y33      0.7474      0.0357   20.9654    0.0000 [ 0.6545; 0.8389 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0155   25.4523    0.0000 [ 0.3564; 0.4368 ] 
#>   eta1 <~ y12      0.3873      0.0158   24.5672    0.0000 [ 0.3422; 0.4238 ] 
#>   eta1 <~ y13      0.4542      0.0219   20.7168    0.0000 [ 0.3996; 0.5130 ] 
#>   eta2 <~ y21      0.3058      0.0292   10.4873    0.0000 [ 0.2315; 0.3823 ] 
#>   eta2 <~ y22      0.4473      0.0207   21.5775    0.0000 [ 0.3909; 0.4981 ] 
#>   eta2 <~ y23      0.4735      0.0218   21.6949    0.0000 [ 0.4151; 0.5280 ] 
#>   eta3 <~ y31      0.4400      0.0161   27.3224    0.0000 [ 0.4018; 0.4850 ] 
#>   eta3 <~ y32      0.3521      0.0195   18.0297    0.0000 [ 0.2988; 0.3998 ] 
#>   eta3 <~ y33      0.3999      0.0186   21.4986    0.0000 [ 0.3533; 0.4495 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0408   16.4481    0.0000 [ 0.5716; 0.7827 ] 
#>   eta3 ~ eta1       0.6634      0.0429   15.4497    0.0000 [ 0.5499; 0.7719 ] 
#>   eta3 ~ eta2       0.3052      0.0868    3.5157    0.0004 [ 0.0887; 0.5375 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0550    3.7223    0.0002 [ 0.0707; 0.3553 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04081527 16.448096 8.654583e-61
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07933135  5.779642 7.485993e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08679789  3.515651 4.386764e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.57162518          0.7826987          0.5969710          0.7573529
#> 2         0.24282153          0.6530785          0.2920854          0.6038146
#> 3         0.08866589          0.5375356          0.1425664          0.4836351
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5579148          0.7213909          0.5730154          0.7209943
#> 2          0.3194253          0.6497264          0.3537907          0.6420939
#> 3          0.1238834          0.4590790          0.1363514          0.4567956
```
