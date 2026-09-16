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
#>  Random seed                        = 1442818945
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
#>   eta2 ~ eta1      0.6713      0.0432   15.5277    0.0000 [ 0.5897; 0.7353 ] 
#>   eta3 ~ eta1      0.4585      0.1034    4.4352    0.0000 [ 0.2694; 0.6199 ] 
#>   eta3 ~ eta2      0.3052      0.1060    2.8781    0.0040 [ 0.1606; 0.5361 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0385   17.2087    0.0000 [ 0.5897; 0.7224 ] 
#>   eta1 =~ y12      0.6493      0.0459   14.1567    0.0000 [ 0.5745; 0.7338 ] 
#>   eta1 =~ y13      0.7613      0.0341   22.3348    0.0000 [ 0.6976; 0.8144 ] 
#>   eta2 =~ y21      0.5165      0.0579    8.9218    0.0000 [ 0.3970; 0.6285 ] 
#>   eta2 =~ y22      0.7554      0.0282   26.7565    0.0000 [ 0.7002; 0.7941 ] 
#>   eta2 =~ y23      0.7997      0.0292   27.3455    0.0000 [ 0.7277; 0.8436 ] 
#>   eta3 =~ y31      0.8223      0.0332   24.7938    0.0000 [ 0.7768; 0.8786 ] 
#>   eta3 =~ y32      0.6581      0.0343   19.1820    0.0000 [ 0.6049; 0.7235 ] 
#>   eta3 =~ y33      0.7474      0.0344   21.7020    0.0000 [ 0.6900; 0.8134 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0207   19.0686    0.0000 [ 0.3532; 0.4229 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.5004    0.0000 [ 0.3575; 0.4292 ] 
#>   eta1 <~ y13      0.4542      0.0232   19.5890    0.0000 [ 0.4114; 0.4956 ] 
#>   eta2 <~ y21      0.3058      0.0279   10.9476    0.0000 [ 0.2549; 0.3596 ] 
#>   eta2 <~ y22      0.4473      0.0190   23.5622    0.0000 [ 0.4112; 0.4871 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7862    0.0000 [ 0.4361; 0.5240 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.8067    0.0000 [ 0.4089; 0.4685 ] 
#>   eta3 <~ y32      0.3521      0.0160   22.0158    0.0000 [ 0.3234; 0.3809 ] 
#>   eta3 <~ y33      0.3999      0.0173   23.1563    0.0000 [ 0.3727; 0.4387 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0432   15.5277    0.0000 [ 0.5897; 0.7353 ] 
#>   eta3 ~ eta1       0.6634      0.0436   15.2117    0.0000 [ 0.5838; 0.7290 ] 
#>   eta3 ~ eta2       0.3052      0.1060    2.8781    0.0040 [ 0.1606; 0.5361 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0756    2.7082    0.0068 [ 0.1055; 0.3771 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03853099 17.208743  2.283326e-66
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04586355 14.156732  1.697056e-45
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03408792 22.334769 1.698054e-110
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05788653  8.921848  4.585727e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.02823190 26.756528 1.036667e-157
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02924302 27.345454 1.223064e-164
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03316465 24.793787 1.046065e-135
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03430667 19.181952  5.237651e-82
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03444032 21.702006 1.964063e-104
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5897361          0.7224446
#> 2          0.5744676          0.7337505
#> 3          0.6975510          0.8144445
#> 4          0.3970253          0.6284920
#> 5          0.7002344          0.7941320
#> 6          0.7276508          0.8435651
#> 7          0.7768045          0.8785843
#> 8          0.6048722          0.7235324
#> 9          0.6899673          0.8134106

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
#>  Random seed                        = 1442818945
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
#>   eta2 ~ eta1      0.6713      0.0432   15.5277    0.0000 [ 0.5586; 0.7822 ] 
#>   eta3 ~ eta1      0.4585      0.1034    4.4352    0.0000 [ 0.1949; 0.7295 ] 
#>   eta3 ~ eta2      0.3052      0.1060    2.8781    0.0040 [ 0.0276; 0.5759 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0385   17.2087    0.0000 [ 0.5751; 0.7743 ] 
#>   eta1 =~ y12      0.6493      0.0459   14.1567    0.0000 [ 0.5310; 0.7682 ] 
#>   eta1 =~ y13      0.7613      0.0341   22.3348    0.0000 [ 0.6720; 0.8482 ] 
#>   eta2 =~ y21      0.5165      0.0579    8.9218    0.0000 [ 0.3707; 0.6701 ] 
#>   eta2 =~ y22      0.7554      0.0282   26.7565    0.0000 [ 0.6936; 0.8396 ] 
#>   eta2 =~ y23      0.7997      0.0292   27.3455    0.0000 [ 0.7244; 0.8757 ] 
#>   eta3 =~ y31      0.8223      0.0332   24.7938    0.0000 [ 0.7348; 0.9063 ] 
#>   eta3 =~ y32      0.6581      0.0343   19.1820    0.0000 [ 0.5753; 0.7527 ] 
#>   eta3 =~ y33      0.7474      0.0344   21.7020    0.0000 [ 0.6551; 0.8332 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0207   19.0686    0.0000 [ 0.3461; 0.4534 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.5004    0.0000 [ 0.3338; 0.4365 ] 
#>   eta1 <~ y13      0.4542      0.0232   19.5890    0.0000 [ 0.3900; 0.5099 ] 
#>   eta2 <~ y21      0.3058      0.0279   10.9476    0.0000 [ 0.2335; 0.3780 ] 
#>   eta2 <~ y22      0.4473      0.0190   23.5622    0.0000 [ 0.4002; 0.4984 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7862    0.0000 [ 0.4125; 0.5249 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.8067    0.0000 [ 0.3912; 0.4867 ] 
#>   eta3 <~ y32      0.3521      0.0160   22.0158    0.0000 [ 0.3140; 0.3967 ] 
#>   eta3 <~ y33      0.3999      0.0173   23.1563    0.0000 [ 0.3535; 0.4428 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0432   15.5277    0.0000 [ 0.5586; 0.7822 ] 
#>   eta3 ~ eta1       0.6634      0.0436   15.2117    0.0000 [ 0.5511; 0.7766 ] 
#>   eta3 ~ eta2       0.3052      0.1060    2.8781    0.0040 [ 0.0276; 0.5759 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0756    2.7082    0.0068 [ 0.0061; 0.3973 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04323458 15.527696 2.253323e-54
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.10337892  4.435205 9.198452e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10602662  2.878062 4.001271e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55860220          0.7821870          0.5854504          0.7553389
#> 2         0.19485359          0.7294710          0.2590507          0.6652739
#> 3         0.02760801          0.5759178          0.0934493          0.5100765
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5684671          0.7452703          0.5897015          0.7352855
#> 2          0.2181043          0.6316498          0.2694241          0.6198725
#> 3          0.1313051          0.5680062          0.1606341          0.5360529
```
