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
#>  Random seed                        = 2093351596
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
#>   eta2 ~ eta1      0.6713      0.0361   18.5862    0.0000 [ 0.6078; 0.7374 ] 
#>   eta3 ~ eta1      0.4585      0.0997    4.5994    0.0000 [ 0.2549; 0.6177 ] 
#>   eta3 ~ eta2      0.3052      0.0972    3.1391    0.0017 [ 0.1524; 0.4702 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0379   17.4752    0.0000 [ 0.5877; 0.7254 ] 
#>   eta1 =~ y12      0.6493      0.0446   14.5634    0.0000 [ 0.5444; 0.7110 ] 
#>   eta1 =~ y13      0.7613      0.0253   30.0589    0.0000 [ 0.7071; 0.7951 ] 
#>   eta2 =~ y21      0.5165      0.0496   10.4213    0.0000 [ 0.4228; 0.6122 ] 
#>   eta2 =~ y22      0.7554      0.0458   16.4802    0.0000 [ 0.6818; 0.8260 ] 
#>   eta2 =~ y23      0.7997      0.0369   21.6781    0.0000 [ 0.7280; 0.8647 ] 
#>   eta3 =~ y31      0.8223      0.0287   28.6379    0.0000 [ 0.7536; 0.8629 ] 
#>   eta3 =~ y32      0.6581      0.0365   18.0520    0.0000 [ 0.5932; 0.7201 ] 
#>   eta3 =~ y33      0.7474      0.0387   19.3338    0.0000 [ 0.6779; 0.8257 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0183   21.6282    0.0000 [ 0.3658; 0.4329 ] 
#>   eta1 <~ y12      0.3873      0.0219   17.6616    0.0000 [ 0.3306; 0.4188 ] 
#>   eta1 <~ y13      0.4542      0.0209   21.7478    0.0000 [ 0.4119; 0.4978 ] 
#>   eta2 <~ y21      0.3058      0.0265   11.5292    0.0000 [ 0.2612; 0.3496 ] 
#>   eta2 <~ y22      0.4473      0.0184   24.3358    0.0000 [ 0.4122; 0.4778 ] 
#>   eta2 <~ y23      0.4735      0.0248   19.1076    0.0000 [ 0.4285; 0.5183 ] 
#>   eta3 <~ y31      0.4400      0.0171   25.7526    0.0000 [ 0.4040; 0.4640 ] 
#>   eta3 <~ y32      0.3521      0.0176   19.9901    0.0000 [ 0.3220; 0.3891 ] 
#>   eta3 <~ y33      0.3999      0.0177   22.6041    0.0000 [ 0.3689; 0.4302 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0361   18.5862    0.0000 [ 0.6078; 0.7374 ] 
#>   eta3 ~ eta1       0.6634      0.0455   14.5739    0.0000 [ 0.5753; 0.7440 ] 
#>   eta3 ~ eta2       0.3052      0.0972    3.1391    0.0017 [ 0.1524; 0.4702 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0668    3.0676    0.0022 [ 0.1091; 0.3305 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03794343 17.47522  2.212792e-68
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04458274 14.56344  4.798234e-48
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02532849 30.05887 1.671732e-198
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04955740 10.42135  1.981279e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04583610 16.48019  5.092341e-61
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03688809 21.67810 3.302247e-104
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02871291 28.63790 2.267953e-180
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03645402 18.05202  7.605168e-73
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03865886 19.33384  2.788842e-83
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5877145          0.7254114
#> 2          0.5444311          0.7109685
#> 3          0.7070720          0.7951252
#> 4          0.4228469          0.6121672
#> 5          0.6817748          0.8260465
#> 6          0.7279769          0.8647414
#> 7          0.7535902          0.8629457
#> 8          0.5931501          0.7200547
#> 9          0.6778842          0.8256975

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
#>  Random seed                        = 2093351596
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
#>   eta2 ~ eta1      0.6713      0.0361   18.5862    0.0000 [ 0.5658; 0.7526 ] 
#>   eta3 ~ eta1      0.4585      0.0997    4.5994    0.0000 [ 0.2044; 0.7200 ] 
#>   eta3 ~ eta2      0.3052      0.0972    3.1391    0.0017 [ 0.0553; 0.5580 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0379   17.4752    0.0000 [ 0.5694; 0.7656 ] 
#>   eta1 =~ y12      0.6493      0.0446   14.5634    0.0000 [ 0.5423; 0.7728 ] 
#>   eta1 =~ y13      0.7613      0.0253   30.0589    0.0000 [ 0.7026; 0.8336 ] 
#>   eta2 =~ y21      0.5165      0.0496   10.4213    0.0000 [ 0.3909; 0.6472 ] 
#>   eta2 =~ y22      0.7554      0.0458   16.4802    0.0000 [ 0.6379; 0.8749 ] 
#>   eta2 =~ y23      0.7997      0.0369   21.6781    0.0000 [ 0.7073; 0.8981 ] 
#>   eta3 =~ y31      0.8223      0.0287   28.6379    0.0000 [ 0.7526; 0.9011 ] 
#>   eta3 =~ y32      0.6581      0.0365   18.0520    0.0000 [ 0.5610; 0.7495 ] 
#>   eta3 =~ y33      0.7474      0.0387   19.3338    0.0000 [ 0.6428; 0.8427 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0183   21.6282    0.0000 [ 0.3458; 0.4403 ] 
#>   eta1 <~ y12      0.3873      0.0219   17.6616    0.0000 [ 0.3306; 0.4440 ] 
#>   eta1 <~ y13      0.4542      0.0209   21.7478    0.0000 [ 0.3978; 0.5058 ] 
#>   eta2 <~ y21      0.3058      0.0265   11.5292    0.0000 [ 0.2377; 0.3749 ] 
#>   eta2 <~ y22      0.4473      0.0184   24.3358    0.0000 [ 0.3989; 0.4940 ] 
#>   eta2 <~ y23      0.4735      0.0248   19.1076    0.0000 [ 0.4091; 0.5372 ] 
#>   eta3 <~ y31      0.4400      0.0171   25.7526    0.0000 [ 0.3991; 0.4875 ] 
#>   eta3 <~ y32      0.3521      0.0176   19.9901    0.0000 [ 0.3059; 0.3970 ] 
#>   eta3 <~ y33      0.3999      0.0177   22.6041    0.0000 [ 0.3527; 0.4442 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0361   18.5862    0.0000 [ 0.5658; 0.7526 ] 
#>   eta3 ~ eta1       0.6634      0.0455   14.5739    0.0000 [ 0.5470; 0.7824 ] 
#>   eta3 ~ eta2       0.3052      0.0972    3.1391    0.0017 [ 0.0553; 0.5580 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0668    3.0676    0.0022 [ 0.0298; 0.3752 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03612007 18.586162 4.158948e-77
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09968850  4.599395 4.237203e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09721099  3.139060 1.694907e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56581223          0.7526048          0.5882424          0.7301747
#> 2         0.20442859          0.7199612          0.2663340          0.6580558
#> 3         0.05526515          0.5579855          0.1156320          0.4976186
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.6040202          0.7492837          0.6078277          0.7373714
#> 2          0.2436767          0.6785346          0.2548985          0.6176922
#> 3          0.1046972          0.4929247          0.1523539          0.4702009
```
