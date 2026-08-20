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
#>  Random seed                        = -1810241874
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
#>   eta2 ~ eta1      0.6713      0.0367   18.3001    0.0000 [ 0.6148; 0.7300 ] 
#>   eta3 ~ eta1      0.4585      0.0912    5.0284    0.0000 [ 0.2635; 0.6164 ] 
#>   eta3 ~ eta2      0.3052      0.0874    3.4902    0.0005 [ 0.1450; 0.4666 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0400   16.5920    0.0000 [ 0.5754; 0.7371 ] 
#>   eta1 =~ y12      0.6493      0.0369   17.5807    0.0000 [ 0.5572; 0.7024 ] 
#>   eta1 =~ y13      0.7613      0.0372   20.4708    0.0000 [ 0.6851; 0.8093 ] 
#>   eta2 =~ y21      0.5165      0.0580    8.9041    0.0000 [ 0.3693; 0.6031 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6714    0.0000 [ 0.6896; 0.8032 ] 
#>   eta2 =~ y23      0.7997      0.0363   22.0119    0.0000 [ 0.6995; 0.8457 ] 
#>   eta3 =~ y31      0.8223      0.0338   24.3152    0.0000 [ 0.7664; 0.8853 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.2033    0.0000 [ 0.5567; 0.7109 ] 
#>   eta3 =~ y33      0.7474      0.0427   17.5130    0.0000 [ 0.6590; 0.8124 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0226   17.5048    0.0000 [ 0.3702; 0.4506 ] 
#>   eta1 <~ y12      0.3873      0.0177   21.8357    0.0000 [ 0.3574; 0.4200 ] 
#>   eta1 <~ y13      0.4542      0.0226   20.1096    0.0000 [ 0.4191; 0.5024 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5576    0.0000 [ 0.2271; 0.3502 ] 
#>   eta2 <~ y22      0.4473      0.0204   21.9419    0.0000 [ 0.4154; 0.4819 ] 
#>   eta2 <~ y23      0.4735      0.0201   23.5501    0.0000 [ 0.4365; 0.5097 ] 
#>   eta3 <~ y31      0.4400      0.0209   21.0153    0.0000 [ 0.4167; 0.4810 ] 
#>   eta3 <~ y32      0.3521      0.0180   19.5750    0.0000 [ 0.3199; 0.3852 ] 
#>   eta3 <~ y33      0.3999      0.0227   17.6340    0.0000 [ 0.3657; 0.4433 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0367   18.3001    0.0000 [ 0.6148; 0.7300 ] 
#>   eta3 ~ eta1       0.6634      0.0449   14.7801    0.0000 [ 0.5652; 0.7257 ] 
#>   eta3 ~ eta2       0.3052      0.0874    3.4902    0.0005 [ 0.1450; 0.4666 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0572    3.5823    0.0003 [ 0.0952; 0.2884 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03996325 16.591989  7.963804e-62
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03693135 17.580674  3.464347e-69
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03719172 20.470844  3.917775e-93
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05800221  8.904054  5.384339e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03331896 22.671404 8.581346e-114
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03632875 22.011867 2.216735e-107
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03381746 24.315171 1.354878e-130
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04328450 15.203339  3.360452e-52
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04267816 17.513033  1.139563e-68
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5753603          0.7371291
#> 2          0.5571560          0.7024144
#> 3          0.6850923          0.8092606
#> 4          0.3693289          0.6031293
#> 5          0.6896266          0.8032310
#> 6          0.6995470          0.8457138
#> 7          0.7664428          0.8853263
#> 8          0.5566858          0.7108673
#> 9          0.6590239          0.8124069

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
#>  Random seed                        = -1810241874
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
#>   eta2 ~ eta1      0.6713      0.0367   18.3001    0.0000 [ 0.5814; 0.7711 ] 
#>   eta3 ~ eta1      0.4585      0.0912    5.0284    0.0000 [ 0.2230; 0.6945 ] 
#>   eta3 ~ eta2      0.3052      0.0874    3.4902    0.0005 [ 0.0811; 0.5332 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0400   16.5920    0.0000 [ 0.5618; 0.7685 ] 
#>   eta1 =~ y12      0.6493      0.0369   17.5807    0.0000 [ 0.5696; 0.7606 ] 
#>   eta1 =~ y13      0.7613      0.0372   20.4708    0.0000 [ 0.6714; 0.8637 ] 
#>   eta2 =~ y21      0.5165      0.0580    8.9041    0.0000 [ 0.3679; 0.6679 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6714    0.0000 [ 0.6731; 0.8455 ] 
#>   eta2 =~ y23      0.7997      0.0363   22.0119    0.0000 [ 0.7152; 0.9031 ] 
#>   eta3 =~ y31      0.8223      0.0338   24.3152    0.0000 [ 0.7386; 0.9135 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.2033    0.0000 [ 0.5512; 0.7750 ] 
#>   eta3 =~ y33      0.7474      0.0427   17.5130    0.0000 [ 0.6443; 0.8650 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0226   17.5048    0.0000 [ 0.3317; 0.4486 ] 
#>   eta1 <~ y12      0.3873      0.0177   21.8357    0.0000 [ 0.3447; 0.4364 ] 
#>   eta1 <~ y13      0.4542      0.0226   20.1096    0.0000 [ 0.3918; 0.5086 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5576    0.0000 [ 0.2212; 0.3867 ] 
#>   eta2 <~ y22      0.4473      0.0204   21.9419    0.0000 [ 0.3925; 0.4979 ] 
#>   eta2 <~ y23      0.4735      0.0201   23.5501    0.0000 [ 0.4226; 0.5266 ] 
#>   eta3 <~ y31      0.4400      0.0209   21.0153    0.0000 [ 0.3831; 0.4914 ] 
#>   eta3 <~ y32      0.3521      0.0180   19.5750    0.0000 [ 0.3049; 0.3979 ] 
#>   eta3 <~ y33      0.3999      0.0227   17.6340    0.0000 [ 0.3410; 0.4583 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0367   18.3001    0.0000 [ 0.5814; 0.7711 ] 
#>   eta3 ~ eta1       0.6634      0.0449   14.7801    0.0000 [ 0.5508; 0.7829 ] 
#>   eta3 ~ eta2       0.3052      0.0874    3.4902    0.0005 [ 0.0811; 0.5332 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0572    3.5823    0.0003 [ 0.0603; 0.3560 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03668474 18.300073 8.263811e-75
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09118354  5.028394 4.946053e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08743058  3.490210 4.826404e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.58142568          0.7711384          0.6042065          0.7483577
#> 2         0.22298842          0.6945382          0.2796123          0.6379143
#> 3         0.08108427          0.5332259          0.1353776          0.4789325
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5858481          0.7432412          0.6147853          0.7300022
#> 2          0.2615849          0.6291233          0.2634964          0.6163903
#> 3          0.1345788          0.5065499          0.1450246          0.4665815
```
