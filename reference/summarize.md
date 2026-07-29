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
#>  Random seed                        = 586124825
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
#>   eta2 ~ eta1      0.6713      0.0446   15.0631    0.0000 [ 0.6061; 0.7883 ] 
#>   eta3 ~ eta1      0.4585      0.1065    4.3048    0.0000 [ 0.3027; 0.6393 ] 
#>   eta3 ~ eta2      0.3052      0.1118    2.7298    0.0063 [ 0.0946; 0.4906 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0410   16.1741    0.0000 [ 0.5587; 0.7125 ] 
#>   eta1 =~ y12      0.6493      0.0392   16.5508    0.0000 [ 0.5793; 0.7307 ] 
#>   eta1 =~ y13      0.7613      0.0328   23.1849    0.0000 [ 0.7144; 0.8111 ] 
#>   eta2 =~ y21      0.5165      0.0654    7.9028    0.0000 [ 0.3623; 0.6278 ] 
#>   eta2 =~ y22      0.7554      0.0432   17.4997    0.0000 [ 0.6678; 0.8488 ] 
#>   eta2 =~ y23      0.7997      0.0378   21.1402    0.0000 [ 0.7259; 0.8622 ] 
#>   eta3 =~ y31      0.8223      0.0423   19.4530    0.0000 [ 0.7485; 0.8973 ] 
#>   eta3 =~ y32      0.6581      0.0449   14.6441    0.0000 [ 0.5788; 0.7278 ] 
#>   eta3 =~ y33      0.7474      0.0403   18.5464    0.0000 [ 0.6628; 0.8142 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0182   21.6946    0.0000 [ 0.3498; 0.4240 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.4234    0.0000 [ 0.3624; 0.4253 ] 
#>   eta1 <~ y13      0.4542      0.0209   21.7196    0.0000 [ 0.4275; 0.4989 ] 
#>   eta2 <~ y21      0.3058      0.0347    8.8048    0.0000 [ 0.2298; 0.3611 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.6291    0.0000 [ 0.4098; 0.4990 ] 
#>   eta2 <~ y23      0.4735      0.0232   20.4335    0.0000 [ 0.4337; 0.5158 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.8366    0.0000 [ 0.4105; 0.4729 ] 
#>   eta3 <~ y32      0.3521      0.0213   16.5528    0.0000 [ 0.3183; 0.3858 ] 
#>   eta3 <~ y33      0.3999      0.0248   16.1459    0.0000 [ 0.3456; 0.4372 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0446   15.0631    0.0000 [ 0.6061; 0.7883 ] 
#>   eta3 ~ eta1       0.6634      0.0453   14.6374    0.0000 [ 0.5991; 0.7514 ] 
#>   eta3 ~ eta2       0.3052      0.1118    2.7298    0.0063 [ 0.0946; 0.4906 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0766    2.6743    0.0075 [ 0.0651; 0.3359 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04099566 16.17415  7.675268e-59
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03922944 16.55078  1.580364e-61
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03283794 23.18495 6.460023e-119
#> 4 eta2 =~ y21  Common factor 0.5164548 0.06535070  7.90282  2.726627e-15
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04316585 17.49966  1.441329e-68
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03782665 21.14022  3.394841e-99
#> 7 eta3 =~ y31  Common factor 0.8222773 0.04226997 19.45299  2.749233e-84
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04493746 14.64411  1.469213e-48
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04030013 18.54644  8.712785e-77
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5587314          0.7124935
#> 2          0.5793409          0.7306949
#> 3          0.7144001          0.8110851
#> 4          0.3622926          0.6278149
#> 5          0.6677968          0.8488033
#> 6          0.7258929          0.8621918
#> 7          0.7485403          0.8973149
#> 8          0.5787847          0.7277942
#> 9          0.6627649          0.8141800

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
#>  Random seed                        = 586124825
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
#>   eta2 ~ eta1      0.6713      0.0446   15.0631    0.0000 [ 0.5480; 0.7785 ] 
#>   eta3 ~ eta1      0.4585      0.1065    4.3048    0.0000 [ 0.1751; 0.7260 ] 
#>   eta3 ~ eta2      0.3052      0.1118    2.7298    0.0063 [ 0.0234; 0.6014 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0410   16.1741    0.0000 [ 0.5690; 0.7810 ] 
#>   eta1 =~ y12      0.6493      0.0392   16.5508    0.0000 [ 0.5489; 0.7517 ] 
#>   eta1 =~ y13      0.7613      0.0328   23.1849    0.0000 [ 0.6732; 0.8430 ] 
#>   eta2 =~ y21      0.5165      0.0654    7.9028    0.0000 [ 0.3414; 0.6794 ] 
#>   eta2 =~ y22      0.7554      0.0432   17.4997    0.0000 [ 0.6486; 0.8719 ] 
#>   eta2 =~ y23      0.7997      0.0378   21.1402    0.0000 [ 0.7039; 0.8995 ] 
#>   eta3 =~ y31      0.8223      0.0423   19.4530    0.0000 [ 0.7075; 0.9261 ] 
#>   eta3 =~ y32      0.6581      0.0449   14.6441    0.0000 [ 0.5390; 0.7714 ] 
#>   eta3 =~ y33      0.7474      0.0403   18.5464    0.0000 [ 0.6531; 0.8615 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0182   21.6946    0.0000 [ 0.3532; 0.4475 ] 
#>   eta1 <~ y12      0.3873      0.0199   19.4234    0.0000 [ 0.3340; 0.4371 ] 
#>   eta1 <~ y13      0.4542      0.0209   21.7196    0.0000 [ 0.3951; 0.5032 ] 
#>   eta2 <~ y21      0.3058      0.0347    8.8048    0.0000 [ 0.2129; 0.3925 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.6291    0.0000 [ 0.3881; 0.5122 ] 
#>   eta2 <~ y23      0.4735      0.0232   20.4335    0.0000 [ 0.4146; 0.5345 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.8366    0.0000 [ 0.3892; 0.4847 ] 
#>   eta3 <~ y32      0.3521      0.0213   16.5528    0.0000 [ 0.2955; 0.4056 ] 
#>   eta3 <~ y33      0.3999      0.0248   16.1459    0.0000 [ 0.3408; 0.4689 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0446   15.0631    0.0000 [ 0.5480; 0.7785 ] 
#>   eta3 ~ eta1       0.6634      0.0453   14.6374    0.0000 [ 0.5411; 0.7754 ] 
#>   eta3 ~ eta2       0.3052      0.1118    2.7298    0.0063 [ 0.0234; 0.6014 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0766    2.6743    0.0075 [ 0.0096; 0.4058 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04456819 15.063061 2.833562e-51
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.10651064  4.304798 1.671382e-05
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.11178477  2.729810 6.337092e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.54803932          0.7785208         0.57571564          0.7508445
#> 2         0.17514225          0.7259551         0.24128411          0.6598133
#> 3         0.02335936          0.6014471         0.09277639          0.5320300
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59877709          0.7910033          0.6061241          0.7883243
#> 2         0.29462146          0.6580497          0.3027141          0.6393190
#> 3         0.08022551          0.5072961          0.0946089          0.4906467
```
