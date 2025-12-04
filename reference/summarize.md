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
#>  Random seed                        = 698173000
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
#>   eta2 ~ eta1      0.6713      0.0438   15.3437    0.0000 [ 0.5971; 0.7461 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8142    0.0000 [ 0.3325; 0.6245 ] 
#>   eta3 ~ eta2      0.3052      0.0853    3.5760    0.0003 [ 0.1175; 0.4152 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0447   14.8473    0.0000 [ 0.5799; 0.7353 ] 
#>   eta1 =~ y12      0.6493      0.0370   17.5460    0.0000 [ 0.6037; 0.7290 ] 
#>   eta1 =~ y13      0.7613      0.0279   27.2476    0.0000 [ 0.7072; 0.8132 ] 
#>   eta2 =~ y21      0.5165      0.0594    8.6993    0.0000 [ 0.3851; 0.6190 ] 
#>   eta2 =~ y22      0.7554      0.0292   25.8427    0.0000 [ 0.6862; 0.8013 ] 
#>   eta2 =~ y23      0.7997      0.0401   19.9245    0.0000 [ 0.7214; 0.8593 ] 
#>   eta3 =~ y31      0.8223      0.0273   30.1558    0.0000 [ 0.7730; 0.8787 ] 
#>   eta3 =~ y32      0.6581      0.0462   14.2510    0.0000 [ 0.5729; 0.7543 ] 
#>   eta3 =~ y33      0.7474      0.0416   17.9455    0.0000 [ 0.6573; 0.8092 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0195   20.3046    0.0000 [ 0.3579; 0.4225 ] 
#>   eta1 <~ y12      0.3873      0.0179   21.6938    0.0000 [ 0.3465; 0.4146 ] 
#>   eta1 <~ y13      0.4542      0.0220   20.6899    0.0000 [ 0.4160; 0.4884 ] 
#>   eta2 <~ y21      0.3058      0.0303   10.0856    0.0000 [ 0.2381; 0.3475 ] 
#>   eta2 <~ y22      0.4473      0.0232   19.2675    0.0000 [ 0.4056; 0.4817 ] 
#>   eta2 <~ y23      0.4735      0.0209   22.6128    0.0000 [ 0.4334; 0.5090 ] 
#>   eta3 <~ y31      0.4400      0.0127   34.6667    0.0000 [ 0.4186; 0.4640 ] 
#>   eta3 <~ y32      0.3521      0.0224   15.7353    0.0000 [ 0.3169; 0.3964 ] 
#>   eta3 <~ y33      0.3999      0.0225   17.7798    0.0000 [ 0.3581; 0.4344 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0438   15.3437    0.0000 [ 0.5971; 0.7461 ] 
#>   eta3 ~ eta1       0.6634      0.0363   18.2907    0.0000 [ 0.6016; 0.7282 ] 
#>   eta3 ~ eta2       0.3052      0.0853    3.5760    0.0003 [ 0.1175; 0.4152 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0587    3.4919    0.0005 [ 0.0806; 0.2712 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04465927 14.847307  7.242041e-50
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03700443 17.545952  6.387311e-69
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02794175 27.247607 1.773986e-163
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05936742  8.699297  3.339457e-18
#> 5 eta2 =~ y22  Common factor 0.7553877 0.02923025 25.842664 2.941870e-147
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04013460 19.924547  2.492911e-88
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02726764 30.155791 9.006089e-200
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04617704 14.250999  4.419151e-46
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04164971 17.945484  5.205540e-72
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5799158          0.7352652
#> 2          0.6036634          0.7289581
#> 3          0.7072042          0.8132110
#> 4          0.3851496          0.6189821
#> 5          0.6862486          0.8012740
#> 6          0.7213763          0.8592970
#> 7          0.7730287          0.8787350
#> 8          0.5728548          0.7543002
#> 9          0.6573226          0.8092066

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
#>  Random seed                        = 698173000
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
#>   eta2 ~ eta1      0.6713      0.0438   15.3437    0.0000 [ 0.5580; 0.7843 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8142    0.0000 [ 0.2467; 0.6545 ] 
#>   eta3 ~ eta2      0.3052      0.0853    3.5760    0.0003 [ 0.0994; 0.5407 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0447   14.8473    0.0000 [ 0.5575; 0.7884 ] 
#>   eta1 =~ y12      0.6493      0.0370   17.5460    0.0000 [ 0.5547; 0.7461 ] 
#>   eta1 =~ y13      0.7613      0.0279   27.2476    0.0000 [ 0.6862; 0.8307 ] 
#>   eta2 =~ y21      0.5165      0.0594    8.6993    0.0000 [ 0.3661; 0.6731 ] 
#>   eta2 =~ y22      0.7554      0.0292   25.8427    0.0000 [ 0.6856; 0.8368 ] 
#>   eta2 =~ y23      0.7997      0.0401   19.9245    0.0000 [ 0.7031; 0.9107 ] 
#>   eta3 =~ y31      0.8223      0.0273   30.1558    0.0000 [ 0.7486; 0.8896 ] 
#>   eta3 =~ y32      0.6581      0.0462   14.2510    0.0000 [ 0.5351; 0.7740 ] 
#>   eta3 =~ y33      0.7474      0.0416   17.9455    0.0000 [ 0.6430; 0.8584 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0195   20.3046    0.0000 [ 0.3494; 0.4501 ] 
#>   eta1 <~ y12      0.3873      0.0179   21.6938    0.0000 [ 0.3399; 0.4323 ] 
#>   eta1 <~ y13      0.4542      0.0220   20.6899    0.0000 [ 0.3930; 0.5065 ] 
#>   eta2 <~ y21      0.3058      0.0303   10.0856    0.0000 [ 0.2265; 0.3833 ] 
#>   eta2 <~ y22      0.4473      0.0232   19.2675    0.0000 [ 0.3856; 0.5057 ] 
#>   eta2 <~ y23      0.4735      0.0209   22.6128    0.0000 [ 0.4187; 0.5270 ] 
#>   eta3 <~ y31      0.4400      0.0127   34.6667    0.0000 [ 0.4067; 0.4724 ] 
#>   eta3 <~ y32      0.3521      0.0224   15.7353    0.0000 [ 0.2935; 0.4092 ] 
#>   eta3 <~ y33      0.3999      0.0225   17.7798    0.0000 [ 0.3446; 0.4609 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0438   15.3437    0.0000 [ 0.5580; 0.7843 ] 
#>   eta3 ~ eta1       0.6634      0.0363   18.2907    0.0000 [ 0.5717; 0.7593 ] 
#>   eta3 ~ eta2       0.3052      0.0853    3.5760    0.0003 [ 0.0994; 0.5407 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0587    3.4919    0.0005 [ 0.0632; 0.3666 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04375296 15.343725 3.901485e-53
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07885992  5.814193 6.092727e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08533421  3.575953 3.489543e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55804509          0.7843107          0.5852152          0.7571406
#> 2         0.24667590          0.6544949          0.2956470          0.6055238
#> 3         0.09937726          0.5406776          0.1523688          0.4876861
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59292284          0.7512488          0.5970572          0.7461329
#> 2         0.31583304          0.6411891          0.3324953          0.6244741
#> 3         0.06222874          0.4714591          0.1175475          0.4152268
```
