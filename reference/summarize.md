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
#>  Random seed                        = 332756229
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
#>   eta2 ~ eta1      0.6713      0.0389   17.2682    0.0000 [ 0.6010; 0.7500 ] 
#>   eta3 ~ eta1      0.4585      0.0704    6.5155    0.0000 [ 0.3117; 0.5781 ] 
#>   eta3 ~ eta2      0.3052      0.0829    3.6799    0.0002 [ 0.1612; 0.4719 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0380   17.4489    0.0000 [ 0.5998; 0.7284 ] 
#>   eta1 =~ y12      0.6493      0.0418   15.5259    0.0000 [ 0.5817; 0.7177 ] 
#>   eta1 =~ y13      0.7613      0.0313   24.3271    0.0000 [ 0.7154; 0.8311 ] 
#>   eta2 =~ y21      0.5165      0.0584    8.8440    0.0000 [ 0.3770; 0.6085 ] 
#>   eta2 =~ y22      0.7554      0.0379   19.9290    0.0000 [ 0.6801; 0.8199 ] 
#>   eta2 =~ y23      0.7997      0.0353   22.6530    0.0000 [ 0.7253; 0.8389 ] 
#>   eta3 =~ y31      0.8223      0.0309   26.6526    0.0000 [ 0.7780; 0.8840 ] 
#>   eta3 =~ y32      0.6581      0.0377   17.4502    0.0000 [ 0.5987; 0.7150 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2332    0.0000 [ 0.6691; 0.8093 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0201   19.7055    0.0000 [ 0.3611; 0.4274 ] 
#>   eta1 <~ y12      0.3873      0.0201   19.2378    0.0000 [ 0.3522; 0.4184 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.7877    0.0000 [ 0.4243; 0.4995 ] 
#>   eta2 <~ y21      0.3058      0.0297   10.3094    0.0000 [ 0.2313; 0.3454 ] 
#>   eta2 <~ y22      0.4473      0.0203   22.0874    0.0000 [ 0.4066; 0.4830 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3622    0.0000 [ 0.4305; 0.5066 ] 
#>   eta3 <~ y31      0.4400      0.0182   24.1634    0.0000 [ 0.4065; 0.4780 ] 
#>   eta3 <~ y32      0.3521      0.0186   18.9349    0.0000 [ 0.3174; 0.3820 ] 
#>   eta3 <~ y33      0.3999      0.0200   19.9620    0.0000 [ 0.3611; 0.4256 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0389   17.2682    0.0000 [ 0.6010; 0.7500 ] 
#>   eta3 ~ eta1       0.6634      0.0338   19.6256    0.0000 [ 0.5988; 0.7179 ] 
#>   eta3 ~ eta2       0.3052      0.0829    3.6799    0.0002 [ 0.1612; 0.4719 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0565    3.6241    0.0003 [ 0.1189; 0.3224 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03800070 17.448884  3.510227e-68
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04181910 15.525869  2.318411e-54
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03129622 24.327085 1.013554e-130
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05839612  8.843993  9.235779e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03790394 19.929001  2.280709e-88
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03530061 22.652972 1.304114e-113
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03085167 26.652600 1.669688e-156
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03771134 17.450157  3.432851e-68
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04099238 18.233247  2.811224e-74
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5998177          0.7283871
#> 2          0.5817235          0.7176747
#> 3          0.7154215          0.8310837
#> 4          0.3769581          0.6084877
#> 5          0.6800973          0.8198539
#> 6          0.7253381          0.8388848
#> 7          0.7780048          0.8839608
#> 8          0.5986839          0.7149727
#> 9          0.6691169          0.8093254

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
#>  Random seed                        = 332756229
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
#>   eta2 ~ eta1      0.6713      0.0389   17.2682    0.0000 [ 0.5630; 0.7640 ] 
#>   eta3 ~ eta1      0.4585      0.0704    6.5155    0.0000 [ 0.2801; 0.6441 ] 
#>   eta3 ~ eta2      0.3052      0.0829    3.6799    0.0002 [ 0.0954; 0.5242 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0380   17.4489    0.0000 [ 0.5690; 0.7655 ] 
#>   eta1 =~ y12      0.6493      0.0418   15.5259    0.0000 [ 0.5414; 0.7576 ] 
#>   eta1 =~ y13      0.7613      0.0313   24.3271    0.0000 [ 0.6792; 0.8411 ] 
#>   eta2 =~ y21      0.5165      0.0584    8.8440    0.0000 [ 0.3679; 0.6699 ] 
#>   eta2 =~ y22      0.7554      0.0379   19.9290    0.0000 [ 0.6507; 0.8467 ] 
#>   eta2 =~ y23      0.7997      0.0353   22.6530    0.0000 [ 0.7143; 0.8969 ] 
#>   eta3 =~ y31      0.8223      0.0309   26.6526    0.0000 [ 0.7380; 0.8975 ] 
#>   eta3 =~ y32      0.6581      0.0377   17.4502    0.0000 [ 0.5609; 0.7560 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2332    0.0000 [ 0.6543; 0.8663 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0201   19.7055    0.0000 [ 0.3455; 0.4493 ] 
#>   eta1 <~ y12      0.3873      0.0201   19.2378    0.0000 [ 0.3349; 0.4390 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.7877    0.0000 [ 0.3986; 0.5064 ] 
#>   eta2 <~ y21      0.3058      0.0297   10.3094    0.0000 [ 0.2310; 0.3844 ] 
#>   eta2 <~ y22      0.4473      0.0203   22.0874    0.0000 [ 0.3909; 0.4957 ] 
#>   eta2 <~ y23      0.4735      0.0212   22.3622    0.0000 [ 0.4221; 0.5316 ] 
#>   eta3 <~ y31      0.4400      0.0182   24.1634    0.0000 [ 0.3880; 0.4821 ] 
#>   eta3 <~ y32      0.3521      0.0186   18.9349    0.0000 [ 0.3024; 0.3985 ] 
#>   eta3 <~ y33      0.3999      0.0200   19.9620    0.0000 [ 0.3529; 0.4566 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0389   17.2682    0.0000 [ 0.5630; 0.7640 ] 
#>   eta3 ~ eta1       0.6634      0.0338   19.6256    0.0000 [ 0.5808; 0.7556 ] 
#>   eta3 ~ eta2       0.3052      0.0829    3.6799    0.0002 [ 0.0954; 0.5242 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0565    3.6241    0.0003 [ 0.0599; 0.3522 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03887683 17.268213 8.163173e-67
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07037131  6.515536 7.243065e-11
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08292340  3.679916 2.333111e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56299057          0.7640396          0.5871326          0.7398975
#> 2         0.28014321          0.6440639          0.3238430          0.6003641
#> 3         0.09539759          0.5242306          0.1468920          0.4727361
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5962560          0.7761199          0.6009832          0.7499921
#> 2          0.3106602          0.5960557          0.3117140          0.5781191
#> 3          0.1579491          0.4893256          0.1611630          0.4719317
```
