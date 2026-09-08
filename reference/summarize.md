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
#>  Random seed                        = -1380591191
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
#>   eta2 ~ eta1      0.6713      0.0365   18.3937    0.0000 [ 0.5962; 0.7398 ] 
#>   eta3 ~ eta1      0.4585      0.0836    5.4845    0.0000 [ 0.3220; 0.6285 ] 
#>   eta3 ~ eta2      0.3052      0.0885    3.4496    0.0006 [ 0.1158; 0.4129 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0421   15.7560    0.0000 [ 0.5712; 0.7388 ] 
#>   eta1 =~ y12      0.6493      0.0379   17.1494    0.0000 [ 0.5690; 0.7100 ] 
#>   eta1 =~ y13      0.7613      0.0336   22.6382    0.0000 [ 0.6940; 0.8087 ] 
#>   eta2 =~ y21      0.5165      0.0544    9.4988    0.0000 [ 0.4328; 0.6288 ] 
#>   eta2 =~ y22      0.7554      0.0304   24.8161    0.0000 [ 0.6910; 0.7897 ] 
#>   eta2 =~ y23      0.7997      0.0397   20.1207    0.0000 [ 0.7302; 0.8668 ] 
#>   eta3 =~ y31      0.8223      0.0324   25.3869    0.0000 [ 0.7731; 0.8896 ] 
#>   eta3 =~ y32      0.6581      0.0396   16.6112    0.0000 [ 0.5709; 0.7209 ] 
#>   eta3 =~ y33      0.7474      0.0376   19.8661    0.0000 [ 0.6579; 0.8000 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0220   17.9730    0.0000 [ 0.3553; 0.4357 ] 
#>   eta1 <~ y12      0.3873      0.0213   18.1539    0.0000 [ 0.3488; 0.4269 ] 
#>   eta1 <~ y13      0.4542      0.0199   22.8317    0.0000 [ 0.4209; 0.5046 ] 
#>   eta2 <~ y21      0.3058      0.0289   10.5877    0.0000 [ 0.2641; 0.3582 ] 
#>   eta2 <~ y22      0.4473      0.0221   20.1953    0.0000 [ 0.3930; 0.4804 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.4708    0.0000 [ 0.4446; 0.5104 ] 
#>   eta3 <~ y31      0.4400      0.0200   21.9905    0.0000 [ 0.4116; 0.4820 ] 
#>   eta3 <~ y32      0.3521      0.0179   19.6782    0.0000 [ 0.3172; 0.3812 ] 
#>   eta3 <~ y33      0.3999      0.0177   22.5878    0.0000 [ 0.3613; 0.4277 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0365   18.3937    0.0000 [ 0.5962; 0.7398 ] 
#>   eta3 ~ eta1       0.6634      0.0402   16.5135    0.0000 [ 0.5964; 0.7451 ] 
#>   eta3 ~ eta2       0.3052      0.0885    3.4496    0.0006 [ 0.1158; 0.4129 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0571    3.5893    0.0003 [ 0.0857; 0.2861 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04208363 15.756006  6.245727e-56
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03786004 17.149427  6.347512e-66
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03363103 22.638192 1.823694e-113
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05437051  9.498803  2.123161e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03043936 24.816147 6.002024e-136
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03974343 20.120654  4.866198e-90
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03238977 25.386949 3.514613e-142
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03961601 16.611187  5.783753e-62
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03762313 19.866082  8.000882e-88
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5711815          0.7387950
#> 2          0.5689967          0.7099592
#> 3          0.6939828          0.8087085
#> 4          0.4328015          0.6287529
#> 5          0.6910162          0.7897135
#> 6          0.7302051          0.8668243
#> 7          0.7731037          0.8895546
#> 8          0.5708984          0.7208604
#> 9          0.6578688          0.7999795

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
#>  Random seed                        = -1380591191
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
#>   eta2 ~ eta1      0.6713      0.0365   18.3937    0.0000 [ 0.5710; 0.7597 ] 
#>   eta3 ~ eta1      0.4585      0.0836    5.4845    0.0000 [ 0.2235; 0.6558 ] 
#>   eta3 ~ eta2      0.3052      0.0885    3.4496    0.0006 [ 0.0955; 0.5530 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0421   15.7560    0.0000 [ 0.5560; 0.7736 ] 
#>   eta1 =~ y12      0.6493      0.0379   17.1494    0.0000 [ 0.5540; 0.7497 ] 
#>   eta1 =~ y13      0.7613      0.0336   22.6382    0.0000 [ 0.6715; 0.8454 ] 
#>   eta2 =~ y21      0.5165      0.0544    9.4988    0.0000 [ 0.3692; 0.6503 ] 
#>   eta2 =~ y22      0.7554      0.0304   24.8161    0.0000 [ 0.6861; 0.8436 ] 
#>   eta2 =~ y23      0.7997      0.0397   20.1207    0.0000 [ 0.7035; 0.9091 ] 
#>   eta3 =~ y31      0.8223      0.0324   25.3869    0.0000 [ 0.7394; 0.9069 ] 
#>   eta3 =~ y32      0.6581      0.0396   16.6112    0.0000 [ 0.5621; 0.7669 ] 
#>   eta3 =~ y33      0.7474      0.0376   19.8661    0.0000 [ 0.6522; 0.8467 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0220   17.9730    0.0000 [ 0.3396; 0.4535 ] 
#>   eta1 <~ y12      0.3873      0.0213   18.1539    0.0000 [ 0.3336; 0.4439 ] 
#>   eta1 <~ y13      0.4542      0.0199   22.8317    0.0000 [ 0.4008; 0.5037 ] 
#>   eta2 <~ y21      0.3058      0.0289   10.5877    0.0000 [ 0.2253; 0.3747 ] 
#>   eta2 <~ y22      0.4473      0.0221   20.1953    0.0000 [ 0.3923; 0.5069 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.4708    0.0000 [ 0.4243; 0.5243 ] 
#>   eta3 <~ y31      0.4400      0.0200   21.9905    0.0000 [ 0.3860; 0.4895 ] 
#>   eta3 <~ y32      0.3521      0.0179   19.6782    0.0000 [ 0.3074; 0.3999 ] 
#>   eta3 <~ y33      0.3999      0.0177   22.5878    0.0000 [ 0.3530; 0.4445 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0365   18.3937    0.0000 [ 0.5710; 0.7597 ] 
#>   eta3 ~ eta1       0.6634      0.0402   16.5135    0.0000 [ 0.5533; 0.7610 ] 
#>   eta3 ~ eta2       0.3052      0.0885    3.4496    0.0006 [ 0.0955; 0.5530 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0571    3.5893    0.0003 [ 0.0699; 0.3651 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03649794 18.393734 1.474646e-75
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08360100  5.484465 4.147221e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08846019  3.449587 5.614449e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.57099035          0.7597371          0.5936551          0.7370723
#> 2         0.22350773          0.6558449          0.2754230          0.6039296
#> 3         0.09549186          0.5529580          0.1504246          0.4980253
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59379719          0.7402954          0.5962094          0.7397751
#> 2         0.31124247          0.6523399          0.3219938          0.6285401
#> 3         0.09946342          0.4364090          0.1158162          0.4128555
```
