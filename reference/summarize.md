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
#>  Random seed                        = 1166480474
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
#>   eta2 ~ eta1      0.6713      0.0467   14.3864    0.0000 [ 0.5973; 0.7546 ] 
#>   eta3 ~ eta1      0.4585      0.0717    6.3919    0.0000 [ 0.3215; 0.5824 ] 
#>   eta3 ~ eta2      0.3052      0.0716    4.2592    0.0000 [ 0.2123; 0.4340 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0410   16.1609    0.0000 [ 0.5849; 0.7357 ] 
#>   eta1 =~ y12      0.6493      0.0340   19.0903    0.0000 [ 0.5761; 0.7014 ] 
#>   eta1 =~ y13      0.7613      0.0342   22.2605    0.0000 [ 0.7049; 0.8217 ] 
#>   eta2 =~ y21      0.5165      0.0460   11.2383    0.0000 [ 0.4329; 0.5999 ] 
#>   eta2 =~ y22      0.7554      0.0337   22.4389    0.0000 [ 0.6769; 0.7959 ] 
#>   eta2 =~ y23      0.7997      0.0410   19.4983    0.0000 [ 0.7365; 0.8558 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5534    0.0000 [ 0.7482; 0.8660 ] 
#>   eta3 =~ y32      0.6581      0.0435   15.1353    0.0000 [ 0.5871; 0.7340 ] 
#>   eta3 =~ y33      0.7474      0.0355   21.0423    0.0000 [ 0.6874; 0.7983 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0209   18.8832    0.0000 [ 0.3531; 0.4393 ] 
#>   eta1 <~ y12      0.3873      0.0168   23.0740    0.0000 [ 0.3557; 0.4139 ] 
#>   eta1 <~ y13      0.4542      0.0211   21.5042    0.0000 [ 0.4213; 0.4977 ] 
#>   eta2 <~ y21      0.3058      0.0228   13.4288    0.0000 [ 0.2567; 0.3424 ] 
#>   eta2 <~ y22      0.4473      0.0214   20.9203    0.0000 [ 0.4089; 0.4829 ] 
#>   eta2 <~ y23      0.4735      0.0190   24.8723    0.0000 [ 0.4426; 0.5076 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.6192    0.0000 [ 0.4006; 0.4666 ] 
#>   eta3 <~ y32      0.3521      0.0160   21.9882    0.0000 [ 0.3219; 0.3752 ] 
#>   eta3 <~ y33      0.3999      0.0193   20.6904    0.0000 [ 0.3731; 0.4403 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0467   14.3864    0.0000 [ 0.5973; 0.7546 ] 
#>   eta3 ~ eta1       0.6634      0.0408   16.2588    0.0000 [ 0.6032; 0.7620 ] 
#>   eta3 ~ eta2       0.3052      0.0716    4.2592    0.0000 [ 0.2123; 0.4340 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0464    4.4154    0.0000 [ 0.1493; 0.2953 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04102937 16.16086  9.522299e-59
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03401089 19.09030  3.040414e-81
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03420172 22.26045 8.934023e-110
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04595475 11.23833  2.643226e-29
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03366420 22.43890 1.642450e-111
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04101192 19.49832  1.134376e-84
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03096695 26.55339 2.346908e-155
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04347904 15.13531  9.472925e-52
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03552007 21.04231  2.689695e-98
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5848839          0.7356909
#> 2          0.5761087          0.7014447
#> 3          0.7049331          0.8216912
#> 4          0.4329106          0.5999276
#> 5          0.6768503          0.7958620
#> 6          0.7365339          0.8558219
#> 7          0.7481622          0.8660153
#> 8          0.5870562          0.7340084
#> 9          0.6874050          0.7983443

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
#>  Random seed                        = 1166480474
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
#>   eta2 ~ eta1      0.6713      0.0467   14.3864    0.0000 [ 0.5421; 0.7835 ] 
#>   eta3 ~ eta1      0.4585      0.0717    6.3919    0.0000 [ 0.2926; 0.6635 ] 
#>   eta3 ~ eta2      0.3052      0.0716    4.2592    0.0000 [ 0.0928; 0.4633 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0410   16.1609    0.0000 [ 0.5634; 0.7756 ] 
#>   eta1 =~ y12      0.6493      0.0340   19.0903    0.0000 [ 0.5648; 0.7407 ] 
#>   eta1 =~ y13      0.7613      0.0342   22.2605    0.0000 [ 0.6727; 0.8496 ] 
#>   eta2 =~ y21      0.5165      0.0460   11.2383    0.0000 [ 0.3862; 0.6239 ] 
#>   eta2 =~ y22      0.7554      0.0337   22.4389    0.0000 [ 0.6755; 0.8496 ] 
#>   eta2 =~ y23      0.7997      0.0410   19.4983    0.0000 [ 0.6930; 0.9051 ] 
#>   eta3 =~ y31      0.8223      0.0310   26.5534    0.0000 [ 0.7493; 0.9094 ] 
#>   eta3 =~ y32      0.6581      0.0435   15.1353    0.0000 [ 0.5505; 0.7753 ] 
#>   eta3 =~ y33      0.7474      0.0355   21.0423    0.0000 [ 0.6517; 0.8354 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0209   18.8832    0.0000 [ 0.3428; 0.4511 ] 
#>   eta1 <~ y12      0.3873      0.0168   23.0740    0.0000 [ 0.3436; 0.4304 ] 
#>   eta1 <~ y13      0.4542      0.0211   21.5042    0.0000 [ 0.3964; 0.5056 ] 
#>   eta2 <~ y21      0.3058      0.0228   13.4288    0.0000 [ 0.2413; 0.3590 ] 
#>   eta2 <~ y22      0.4473      0.0214   20.9203    0.0000 [ 0.3973; 0.5079 ] 
#>   eta2 <~ y23      0.4735      0.0190   24.8723    0.0000 [ 0.4254; 0.5238 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.6192    0.0000 [ 0.3910; 0.4916 ] 
#>   eta3 <~ y32      0.3521      0.0160   21.9882    0.0000 [ 0.3118; 0.3947 ] 
#>   eta3 <~ y33      0.3999      0.0193   20.6904    0.0000 [ 0.3457; 0.4457 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0467   14.3864    0.0000 [ 0.5421; 0.7835 ] 
#>   eta3 ~ eta1       0.6634      0.0408   16.2588    0.0000 [ 0.5575; 0.7685 ] 
#>   eta3 ~ eta2       0.3052      0.0716    4.2592    0.0000 [ 0.0928; 0.4633 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0464    4.4154    0.0000 [ 0.0650; 0.3049 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04666460 14.386352 6.303396e-47
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07173213  6.391930 1.638048e-10
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07164562  4.259174 2.051841e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.54213872          0.7834617          0.5711169          0.7544835
#> 2         0.29255194          0.6635100          0.3370968          0.6189652
#> 3         0.09281842          0.4633291          0.1373095          0.4188380
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5947518          0.7694056          0.5972533          0.7546205
#> 2          0.2883053          0.5872651          0.3214818          0.5824216
#> 3          0.1920832          0.4353357          0.2123420          0.4339856
```
