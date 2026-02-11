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
#>  Random seed                        = -1918252376
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
#>   eta2 ~ eta1      0.6713      0.0509   13.1765    0.0000 [ 0.5719; 0.7635 ] 
#>   eta3 ~ eta1      0.4585      0.0838    5.4709    0.0000 [ 0.2548; 0.5819 ] 
#>   eta3 ~ eta2      0.3052      0.0821    3.7169    0.0002 [ 0.2108; 0.5283 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0429   15.4422    0.0000 [ 0.5592; 0.7254 ] 
#>   eta1 =~ y12      0.6493      0.0377   17.2182    0.0000 [ 0.5710; 0.7118 ] 
#>   eta1 =~ y13      0.7613      0.0323   23.5837    0.0000 [ 0.7005; 0.8266 ] 
#>   eta2 =~ y21      0.5165      0.0553    9.3370    0.0000 [ 0.3667; 0.5939 ] 
#>   eta2 =~ y22      0.7554      0.0378   19.9797    0.0000 [ 0.6916; 0.8296 ] 
#>   eta2 =~ y23      0.7997      0.0291   27.5148    0.0000 [ 0.7401; 0.8534 ] 
#>   eta3 =~ y31      0.8223      0.0254   32.3761    0.0000 [ 0.7687; 0.8631 ] 
#>   eta3 =~ y32      0.6581      0.0354   18.5677    0.0000 [ 0.5855; 0.7213 ] 
#>   eta3 =~ y33      0.7474      0.0364   20.5203    0.0000 [ 0.6796; 0.8034 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0210   18.8071    0.0000 [ 0.3559; 0.4370 ] 
#>   eta1 <~ y12      0.3873      0.0217   17.8199    0.0000 [ 0.3431; 0.4220 ] 
#>   eta1 <~ y13      0.4542      0.0202   22.4536    0.0000 [ 0.4262; 0.5001 ] 
#>   eta2 <~ y21      0.3058      0.0277   11.0229    0.0000 [ 0.2380; 0.3415 ] 
#>   eta2 <~ y22      0.4473      0.0194   23.0717    0.0000 [ 0.4222; 0.4916 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5428    0.0000 [ 0.4425; 0.5157 ] 
#>   eta3 <~ y31      0.4400      0.0175   25.1602    0.0000 [ 0.4136; 0.4726 ] 
#>   eta3 <~ y32      0.3521      0.0153   23.0826    0.0000 [ 0.3216; 0.3799 ] 
#>   eta3 <~ y33      0.3999      0.0152   26.3929    0.0000 [ 0.3643; 0.4323 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0509   13.1765    0.0000 [ 0.5719; 0.7635 ] 
#>   eta3 ~ eta1       0.6634      0.0394   16.8288    0.0000 [ 0.6039; 0.7465 ] 
#>   eta3 ~ eta2       0.3052      0.0821    3.7169    0.0002 [ 0.2108; 0.5283 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0609    3.3621    0.0008 [ 0.1374; 0.3863 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04293882 15.442202  8.514331e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03770879 17.218212  1.938874e-66
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03228277 23.583660 5.670770e-123
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05531255  9.337028  9.907561e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03780777 19.979692  8.273222e-89
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02906305 27.514785 1.168343e-166
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02539763 32.376144 5.947957e-230
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03544158 18.567705  5.865666e-77
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03642358 20.520336  1.417296e-93
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5591969          0.7254355
#> 2          0.5709714          0.7118402
#> 3          0.7004817          0.8266129
#> 4          0.3666848          0.5938803
#> 5          0.6915566          0.8296345
#> 6          0.7401248          0.8534206
#> 7          0.7686908          0.8631131
#> 8          0.5854796          0.7212753
#> 9          0.6795723          0.8034046

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
#>  Random seed                        = -1918252376
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
#>   eta2 ~ eta1      0.6713      0.0509   13.1765    0.0000 [ 0.5371; 0.8006 ] 
#>   eta3 ~ eta1      0.4585      0.0838    5.4709    0.0000 [ 0.2441; 0.6775 ] 
#>   eta3 ~ eta2      0.3052      0.0821    3.7169    0.0002 [ 0.0833; 0.5079 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0429   15.4422    0.0000 [ 0.5490; 0.7711 ] 
#>   eta1 =~ y12      0.6493      0.0377   17.2182    0.0000 [ 0.5522; 0.7472 ] 
#>   eta1 =~ y13      0.7613      0.0323   23.5837    0.0000 [ 0.6780; 0.8449 ] 
#>   eta2 =~ y21      0.5165      0.0553    9.3370    0.0000 [ 0.3789; 0.6649 ] 
#>   eta2 =~ y22      0.7554      0.0378   19.9797    0.0000 [ 0.6520; 0.8476 ] 
#>   eta2 =~ y23      0.7997      0.0291   27.5148    0.0000 [ 0.7293; 0.8796 ] 
#>   eta3 =~ y31      0.8223      0.0254   32.3761    0.0000 [ 0.7586; 0.8899 ] 
#>   eta3 =~ y32      0.6581      0.0354   18.5677    0.0000 [ 0.5724; 0.7557 ] 
#>   eta3 =~ y33      0.7474      0.0364   20.5203    0.0000 [ 0.6534; 0.8418 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0210   18.8071    0.0000 [ 0.3405; 0.4493 ] 
#>   eta1 <~ y12      0.3873      0.0217   17.8199    0.0000 [ 0.3322; 0.4446 ] 
#>   eta1 <~ y13      0.4542      0.0202   22.4536    0.0000 [ 0.4029; 0.5075 ] 
#>   eta2 <~ y21      0.3058      0.0277   11.0229    0.0000 [ 0.2372; 0.3807 ] 
#>   eta2 <~ y22      0.4473      0.0194   23.0717    0.0000 [ 0.3931; 0.4933 ] 
#>   eta2 <~ y23      0.4735      0.0210   22.5428    0.0000 [ 0.4209; 0.5295 ] 
#>   eta3 <~ y31      0.4400      0.0175   25.1602    0.0000 [ 0.3933; 0.4838 ] 
#>   eta3 <~ y32      0.3521      0.0153   23.0826    0.0000 [ 0.3142; 0.3931 ] 
#>   eta3 <~ y33      0.3999      0.0152   26.3929    0.0000 [ 0.3589; 0.4373 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0509   13.1765    0.0000 [ 0.5371; 0.8006 ] 
#>   eta3 ~ eta1       0.6634      0.0394   16.8288    0.0000 [ 0.5561; 0.7600 ] 
#>   eta3 ~ eta2       0.3052      0.0821    3.7169    0.0002 [ 0.0833; 0.5079 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0609    3.3621    0.0008 [ 0.0397; 0.3548 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.05094922 13.176520 1.197910e-39
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08380862  5.470878 4.478108e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08209790  3.716918 2.016682e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.53707273          0.8005533          0.5687116          0.7689145
#> 2         0.24408818          0.6774991          0.2961324          0.6254549
#> 3         0.08332241          0.5078864          0.1343042          0.4569046
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5383378          0.7679064          0.5719101          0.7634733
#> 2          0.2447831          0.6338303          0.2547572          0.5819221
#> 3          0.1421308          0.5496254          0.2107850          0.5283282
```
