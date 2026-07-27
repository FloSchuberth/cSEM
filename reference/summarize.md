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
#>  Random seed                        = -969749304
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
#>   eta2 ~ eta1      0.6713      0.0434   15.4779    0.0000 [ 0.5852; 0.7409 ] 
#>   eta3 ~ eta1      0.4585      0.0637    7.1933    0.0000 [ 0.3237; 0.5591 ] 
#>   eta3 ~ eta2      0.3052      0.0681    4.4797    0.0000 [ 0.1991; 0.4550 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0404   16.4083    0.0000 [ 0.5850; 0.7337 ] 
#>   eta1 =~ y12      0.6493      0.0422   15.3765    0.0000 [ 0.5748; 0.7357 ] 
#>   eta1 =~ y13      0.7613      0.0352   21.6102    0.0000 [ 0.6818; 0.8149 ] 
#>   eta2 =~ y21      0.5165      0.0482   10.7070    0.0000 [ 0.4279; 0.5939 ] 
#>   eta2 =~ y22      0.7554      0.0327   23.0987    0.0000 [ 0.7071; 0.8184 ] 
#>   eta2 =~ y23      0.7997      0.0343   23.3221    0.0000 [ 0.7386; 0.8572 ] 
#>   eta3 =~ y31      0.8223      0.0316   26.0583    0.0000 [ 0.7652; 0.8732 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1406    0.0000 [ 0.5870; 0.7321 ] 
#>   eta3 =~ y33      0.7474      0.0433   17.2680    0.0000 [ 0.6731; 0.8216 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0215   18.4120    0.0000 [ 0.3595; 0.4282 ] 
#>   eta1 <~ y12      0.3873      0.0228   17.0186    0.0000 [ 0.3437; 0.4249 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0920    0.0000 [ 0.4144; 0.4942 ] 
#>   eta2 <~ y21      0.3058      0.0250   12.2225    0.0000 [ 0.2624; 0.3512 ] 
#>   eta2 <~ y22      0.4473      0.0179   25.0161    0.0000 [ 0.4126; 0.4772 ] 
#>   eta2 <~ y23      0.4735      0.0203   23.3501    0.0000 [ 0.4374; 0.5005 ] 
#>   eta3 <~ y31      0.4400      0.0192   22.9408    0.0000 [ 0.4109; 0.4792 ] 
#>   eta3 <~ y32      0.3521      0.0198   17.7720    0.0000 [ 0.3198; 0.3868 ] 
#>   eta3 <~ y33      0.3999      0.0205   19.5140    0.0000 [ 0.3666; 0.4434 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0434   15.4779    0.0000 [ 0.5852; 0.7409 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7502    0.0000 [ 0.5994; 0.7014 ] 
#>   eta3 ~ eta2       0.3052      0.0681    4.4797    0.0000 [ 0.1991; 0.4550 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0447    4.5876    0.0000 [ 0.1413; 0.2991 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04041055 16.40833  1.667181e-60
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04222537 15.37649  2.353658e-53
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03523093 21.61016 1.441601e-103
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04823512 10.70703  9.434239e-27
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03270261 23.09870 4.771728e-118
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03428784 23.32208 2.647067e-120
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03155528 26.05832 1.082921e-149
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04077108 16.14058  1.322889e-58
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04328367 17.26804  8.187903e-67
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5849861          0.7336582
#> 2          0.5747834          0.7357452
#> 3          0.6817897          0.8149327
#> 4          0.4279270          0.5939404
#> 5          0.7070877          0.8184240
#> 6          0.7386411          0.8572324
#> 7          0.7652060          0.8731680
#> 8          0.5869927          0.7321430
#> 9          0.6731086          0.8215967

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
#>  Random seed                        = -969749304
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
#>   eta2 ~ eta1      0.6713      0.0434   15.4779    0.0000 [ 0.5633; 0.7876 ] 
#>   eta3 ~ eta1      0.4585      0.0637    7.1933    0.0000 [ 0.3119; 0.6415 ] 
#>   eta3 ~ eta2      0.3052      0.0681    4.4797    0.0000 [ 0.1052; 0.4575 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0404   16.4083    0.0000 [ 0.5616; 0.7705 ] 
#>   eta1 =~ y12      0.6493      0.0422   15.3765    0.0000 [ 0.5320; 0.7504 ] 
#>   eta1 =~ y13      0.7613      0.0352   21.6102    0.0000 [ 0.6729; 0.8551 ] 
#>   eta2 =~ y21      0.5165      0.0482   10.7070    0.0000 [ 0.3785; 0.6280 ] 
#>   eta2 =~ y22      0.7554      0.0327   23.0987    0.0000 [ 0.6668; 0.8359 ] 
#>   eta2 =~ y23      0.7997      0.0343   23.3221    0.0000 [ 0.7114; 0.8887 ] 
#>   eta3 =~ y31      0.8223      0.0316   26.0583    0.0000 [ 0.7438; 0.9070 ] 
#>   eta3 =~ y32      0.6581      0.0408   16.1406    0.0000 [ 0.5496; 0.7604 ] 
#>   eta3 =~ y33      0.7474      0.0433   17.2680    0.0000 [ 0.6336; 0.8574 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0215   18.4120    0.0000 [ 0.3427; 0.4538 ] 
#>   eta1 <~ y12      0.3873      0.0228   17.0186    0.0000 [ 0.3246; 0.4423 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0920    0.0000 [ 0.4035; 0.5099 ] 
#>   eta2 <~ y21      0.3058      0.0250   12.2225    0.0000 [ 0.2372; 0.3666 ] 
#>   eta2 <~ y22      0.4473      0.0179   25.0161    0.0000 [ 0.4038; 0.4963 ] 
#>   eta2 <~ y23      0.4735      0.0203   23.3501    0.0000 [ 0.4267; 0.5316 ] 
#>   eta3 <~ y31      0.4400      0.0192   22.9408    0.0000 [ 0.3926; 0.4918 ] 
#>   eta3 <~ y32      0.3521      0.0198   17.7720    0.0000 [ 0.2999; 0.4024 ] 
#>   eta3 <~ y33      0.3999      0.0205   19.5140    0.0000 [ 0.3467; 0.4526 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0434   15.4779    0.0000 [ 0.5633; 0.7876 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7502    0.0000 [ 0.5761; 0.7590 ] 
#>   eta3 ~ eta2       0.3052      0.0681    4.4797    0.0000 [ 0.1052; 0.4575 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0447    4.5876    0.0000 [ 0.0754; 0.3063 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04337367 15.477901 4.891755e-54
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06374039  7.193347 6.322183e-13
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.06811908  4.479672 7.475799e-06
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5633139          0.7876180          0.5902484          0.7606835
#> 2          0.3118817          0.6415110          0.3514638          0.6019290
#> 3          0.1052262          0.4574997          0.1475274          0.4151985
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5795856          0.7628139          0.5851605          0.7408574
#> 2          0.3089496          0.5979726          0.3237176          0.5591471
#> 3          0.1882120          0.4569946          0.1990880          0.4549640
```
