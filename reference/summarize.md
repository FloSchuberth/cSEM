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
#>  Random seed                        = -1726277882
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
#>   eta2 ~ eta1      0.6713      0.0413   16.2671    0.0000 [ 0.6131; 0.7573 ] 
#>   eta3 ~ eta1      0.4585      0.0888    5.1615    0.0000 [ 0.2692; 0.6115 ] 
#>   eta3 ~ eta2      0.3052      0.0864    3.5339    0.0004 [ 0.1708; 0.4532 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0426   15.5746    0.0000 [ 0.5618; 0.7416 ] 
#>   eta1 =~ y12      0.6493      0.0390   16.6686    0.0000 [ 0.5674; 0.7159 ] 
#>   eta1 =~ y13      0.7613      0.0326   23.3848    0.0000 [ 0.7076; 0.8265 ] 
#>   eta2 =~ y21      0.5165      0.0500   10.3374    0.0000 [ 0.4354; 0.5907 ] 
#>   eta2 =~ y22      0.7554      0.0322   23.4844    0.0000 [ 0.7092; 0.8307 ] 
#>   eta2 =~ y23      0.7997      0.0361   22.1279    0.0000 [ 0.7204; 0.8566 ] 
#>   eta3 =~ y31      0.8223      0.0302   27.2342    0.0000 [ 0.7603; 0.8666 ] 
#>   eta3 =~ y32      0.6581      0.0442   14.8729    0.0000 [ 0.5817; 0.7278 ] 
#>   eta3 =~ y33      0.7474      0.0357   20.9172    0.0000 [ 0.6702; 0.8148 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0201   19.6765    0.0000 [ 0.3567; 0.4367 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9263    0.0000 [ 0.3532; 0.4282 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0513    0.0000 [ 0.4277; 0.5064 ] 
#>   eta2 <~ y21      0.3058      0.0295   10.3833    0.0000 [ 0.2606; 0.3518 ] 
#>   eta2 <~ y22      0.4473      0.0176   25.4676    0.0000 [ 0.4171; 0.4759 ] 
#>   eta2 <~ y23      0.4735      0.0174   27.1964    0.0000 [ 0.4401; 0.5025 ] 
#>   eta3 <~ y31      0.4400      0.0169   25.9592    0.0000 [ 0.4103; 0.4770 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5438    0.0000 [ 0.3183; 0.3891 ] 
#>   eta3 <~ y33      0.3999      0.0179   22.2851    0.0000 [ 0.3703; 0.4304 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0413   16.2671    0.0000 [ 0.6131; 0.7573 ] 
#>   eta3 ~ eta1       0.6634      0.0463   14.3423    0.0000 [ 0.5833; 0.7577 ] 
#>   eta3 ~ eta2       0.3052      0.0864    3.5339    0.0004 [ 0.1708; 0.4532 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0620    3.3051    0.0009 [ 0.1195; 0.3225 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04257375 15.57462  1.082980e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03895216 16.66860  2.217267e-62
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03255728 23.38481 6.100747e-121
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04995991 10.33739  4.773859e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03216554 23.48437 5.891649e-122
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03613829 22.12788 1.703972e-108
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03019281 27.23421 2.556338e-163
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04424618 14.87290  4.942579e-50
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03573254 20.91718  3.735235e-97
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5618321          0.7416211
#> 2          0.5673701          0.7158822
#> 3          0.7076496          0.8264904
#> 4          0.4354018          0.5907319
#> 5          0.7091780          0.8306916
#> 6          0.7203547          0.8566078
#> 7          0.7603195          0.8665804
#> 8          0.5816735          0.7277860
#> 9          0.6701617          0.8148235

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
#>  Random seed                        = -1726277882
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
#>   eta2 ~ eta1      0.6713      0.0413   16.2671    0.0000 [ 0.5551; 0.7685 ] 
#>   eta3 ~ eta1      0.4585      0.0888    5.1615    0.0000 [ 0.2239; 0.6833 ] 
#>   eta3 ~ eta2      0.3052      0.0864    3.5339    0.0004 [ 0.0826; 0.5292 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0426   15.5746    0.0000 [ 0.5618; 0.7819 ] 
#>   eta1 =~ y12      0.6493      0.0390   16.6686    0.0000 [ 0.5487; 0.7502 ] 
#>   eta1 =~ y13      0.7613      0.0326   23.3848    0.0000 [ 0.6756; 0.8439 ] 
#>   eta2 =~ y21      0.5165      0.0500   10.3374    0.0000 [ 0.3800; 0.6384 ] 
#>   eta2 =~ y22      0.7554      0.0322   23.4844    0.0000 [ 0.6709; 0.8373 ] 
#>   eta2 =~ y23      0.7997      0.0361   22.1279    0.0000 [ 0.7096; 0.8965 ] 
#>   eta3 =~ y31      0.8223      0.0302   27.2342    0.0000 [ 0.7463; 0.9025 ] 
#>   eta3 =~ y32      0.6581      0.0442   14.8729    0.0000 [ 0.5406; 0.7694 ] 
#>   eta3 =~ y33      0.7474      0.0357   20.9172    0.0000 [ 0.6558; 0.8406 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0201   19.6765    0.0000 [ 0.3472; 0.4512 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9263    0.0000 [ 0.3298; 0.4415 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.0513    0.0000 [ 0.3977; 0.5042 ] 
#>   eta2 <~ y21      0.3058      0.0295   10.3833    0.0000 [ 0.2265; 0.3788 ] 
#>   eta2 <~ y22      0.4473      0.0176   25.4676    0.0000 [ 0.4029; 0.4937 ] 
#>   eta2 <~ y23      0.4735      0.0174   27.1964    0.0000 [ 0.4324; 0.5225 ] 
#>   eta3 <~ y31      0.4400      0.0169   25.9592    0.0000 [ 0.3973; 0.4850 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5438    0.0000 [ 0.3017; 0.3999 ] 
#>   eta3 <~ y33      0.3999      0.0179   22.2851    0.0000 [ 0.3541; 0.4469 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0413   16.2671    0.0000 [ 0.5551; 0.7685 ] 
#>   eta3 ~ eta1       0.6634      0.0463   14.3423    0.0000 [ 0.5363; 0.7755 ] 
#>   eta3 ~ eta2       0.3052      0.0864    3.5339    0.0004 [ 0.0826; 0.5292 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0620    3.3051    0.0009 [ 0.0421; 0.3626 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04126945 16.267079 1.690323e-59
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08883197  5.161506 2.449710e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08635057  3.533864 4.095324e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55508418          0.7685065          0.5807120          0.7428786
#> 2         0.22386790          0.6832567          0.2790315          0.6280931
#> 3         0.08261007          0.5291664          0.1362328          0.4755438
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5945244          0.7652777          0.6130504          0.7572833
#> 2          0.2396568          0.6385464          0.2692339          0.6115492
#> 3          0.1414910          0.5020673          0.1707817          0.4531757
```
