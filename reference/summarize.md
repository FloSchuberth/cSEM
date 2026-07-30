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
#>  Random seed                        = 1556851927
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
#>   eta2 ~ eta1      0.6713      0.0464   14.4577    0.0000 [ 0.5823; 0.7525 ] 
#>   eta3 ~ eta1      0.4585      0.0842    5.4436    0.0000 [ 0.3419; 0.6480 ] 
#>   eta3 ~ eta2      0.3052      0.0938    3.2537    0.0011 [ 0.1134; 0.4441 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0426   15.5816    0.0000 [ 0.5983; 0.7464 ] 
#>   eta1 =~ y12      0.6493      0.0402   16.1625    0.0000 [ 0.5704; 0.7154 ] 
#>   eta1 =~ y13      0.7613      0.0327   23.2591    0.0000 [ 0.6965; 0.8194 ] 
#>   eta2 =~ y21      0.5165      0.0474   10.8916    0.0000 [ 0.4262; 0.5783 ] 
#>   eta2 =~ y22      0.7554      0.0345   21.9102    0.0000 [ 0.6866; 0.7990 ] 
#>   eta2 =~ y23      0.7997      0.0352   22.7193    0.0000 [ 0.7395; 0.8927 ] 
#>   eta3 =~ y31      0.8223      0.0316   25.9963    0.0000 [ 0.7680; 0.8800 ] 
#>   eta3 =~ y32      0.6581      0.0421   15.6182    0.0000 [ 0.5624; 0.7243 ] 
#>   eta3 =~ y33      0.7474      0.0395   18.9302    0.0000 [ 0.6910; 0.8330 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0197   20.0687    0.0000 [ 0.3630; 0.4287 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.7752    0.0000 [ 0.3439; 0.4235 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.8487    0.0000 [ 0.4173; 0.4869 ] 
#>   eta2 <~ y21      0.3058      0.0244   12.5360    0.0000 [ 0.2466; 0.3400 ] 
#>   eta2 <~ y22      0.4473      0.0188   23.8098    0.0000 [ 0.4169; 0.4806 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5056    0.0000 [ 0.4456; 0.5232 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5750    0.0000 [ 0.4060; 0.4711 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1897    0.0000 [ 0.3129; 0.3758 ] 
#>   eta3 <~ y33      0.3999      0.0198   20.1839    0.0000 [ 0.3742; 0.4460 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0464   14.4577    0.0000 [ 0.5823; 0.7525 ] 
#>   eta3 ~ eta1       0.6634      0.0335   19.7971    0.0000 [ 0.6044; 0.7156 ] 
#>   eta3 ~ eta2       0.3052      0.0938    3.2537    0.0011 [ 0.1134; 0.4441 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0628    3.2632    0.0011 [ 0.0852; 0.3030 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04255478 15.58156  9.715583e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04017188 16.16250  9.272755e-59
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03273321 23.25912 1.150141e-119
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04741780 10.89158  1.264254e-27
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03447650 21.91022 2.075778e-106
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03519759 22.71928 2.889218e-114
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03163051 25.99634 5.447691e-149
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04213466 15.61823  5.470120e-55
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03948317 18.93019  6.432633e-80
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5982529          0.7463581
#> 2          0.5704209          0.7154147
#> 3          0.6964987          0.8194307
#> 4          0.4262154          0.5782878
#> 5          0.6865882          0.7989533
#> 6          0.7395411          0.8926966
#> 7          0.7679993          0.8800132
#> 8          0.5624098          0.7242733
#> 9          0.6910445          0.8330285

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
#>  Random seed                        = 1556851927
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
#>   eta2 ~ eta1      0.6713      0.0464   14.4577    0.0000 [ 0.5477; 0.7878 ] 
#>   eta3 ~ eta1      0.4585      0.0842    5.4436    0.0000 [ 0.2357; 0.6712 ] 
#>   eta3 ~ eta2      0.3052      0.0938    3.2537    0.0011 [ 0.0747; 0.5597 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0426   15.5816    0.0000 [ 0.5495; 0.7696 ] 
#>   eta1 =~ y12      0.6493      0.0402   16.1625    0.0000 [ 0.5434; 0.7512 ] 
#>   eta1 =~ y13      0.7613      0.0327   23.2591    0.0000 [ 0.6770; 0.8462 ] 
#>   eta2 =~ y21      0.5165      0.0474   10.8916    0.0000 [ 0.4035; 0.6487 ] 
#>   eta2 =~ y22      0.7554      0.0345   21.9102    0.0000 [ 0.6731; 0.8514 ] 
#>   eta2 =~ y23      0.7997      0.0352   22.7193    0.0000 [ 0.7065; 0.8885 ] 
#>   eta3 =~ y31      0.8223      0.0316   25.9963    0.0000 [ 0.7381; 0.9017 ] 
#>   eta3 =~ y32      0.6581      0.0421   15.6182    0.0000 [ 0.5567; 0.7746 ] 
#>   eta3 =~ y33      0.7474      0.0395   18.9302    0.0000 [ 0.6365; 0.8407 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0197   20.0687    0.0000 [ 0.3444; 0.4463 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.7752    0.0000 [ 0.3315; 0.4441 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.8487    0.0000 [ 0.4023; 0.5098 ] 
#>   eta2 <~ y21      0.3058      0.0244   12.5360    0.0000 [ 0.2461; 0.3722 ] 
#>   eta2 <~ y22      0.4473      0.0188   23.8098    0.0000 [ 0.3989; 0.4960 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5056    0.0000 [ 0.4136; 0.5224 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5750    0.0000 [ 0.3895; 0.4903 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1897    0.0000 [ 0.3145; 0.4005 ] 
#>   eta3 <~ y33      0.3999      0.0198   20.1839    0.0000 [ 0.3453; 0.4477 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0464   14.4577    0.0000 [ 0.5477; 0.7878 ] 
#>   eta3 ~ eta1       0.6634      0.0335   19.7971    0.0000 [ 0.5796; 0.7529 ] 
#>   eta3 ~ eta2       0.3052      0.0938    3.2537    0.0011 [ 0.0747; 0.5597 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0628    3.2632    0.0011 [ 0.0505; 0.3752 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04643422 14.457732 2.240608e-47
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08422839  5.443613 5.221053e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09378480  3.253738 1.138975e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.54769564          0.7878272          0.5765307          0.7589921
#> 2         0.23566699          0.6712487          0.2879718          0.6189438
#> 3         0.07473714          0.5597392          0.1329764          0.5014999
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.57369975          0.7783204          0.5823188          0.7525222
#> 2         0.32580079          0.6633126          0.3419354          0.6480372
#> 3         0.07512644          0.4560733          0.1133771          0.4441015
```
