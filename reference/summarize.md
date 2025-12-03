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
#>  Random seed                        = 52038992
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
#>   eta2 ~ eta1      0.6713      0.0369   18.1794    0.0000 [ 0.5925; 0.7220 ] 
#>   eta3 ~ eta1      0.4585      0.0814    5.6298    0.0000 [ 0.3215; 0.6152 ] 
#>   eta3 ~ eta2      0.3052      0.0817    3.7365    0.0002 [ 0.1365; 0.4251 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0398   16.6757    0.0000 [ 0.5972; 0.7408 ] 
#>   eta1 =~ y12      0.6493      0.0363   17.8632    0.0000 [ 0.5838; 0.7212 ] 
#>   eta1 =~ y13      0.7613      0.0350   21.7790    0.0000 [ 0.7027; 0.8280 ] 
#>   eta2 =~ y21      0.5165      0.0583    8.8639    0.0000 [ 0.3891; 0.5816 ] 
#>   eta2 =~ y22      0.7554      0.0385   19.6185    0.0000 [ 0.6805; 0.8378 ] 
#>   eta2 =~ y23      0.7997      0.0417   19.1542    0.0000 [ 0.7188; 0.8671 ] 
#>   eta3 =~ y31      0.8223      0.0295   27.8673    0.0000 [ 0.7699; 0.8721 ] 
#>   eta3 =~ y32      0.6581      0.0463   14.2162    0.0000 [ 0.5768; 0.7482 ] 
#>   eta3 =~ y33      0.7474      0.0363   20.6049    0.0000 [ 0.6819; 0.8157 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0217   18.1888    0.0000 [ 0.3578; 0.4379 ] 
#>   eta1 <~ y12      0.3873      0.0194   19.9240    0.0000 [ 0.3556; 0.4173 ] 
#>   eta1 <~ y13      0.4542      0.0203   22.3251    0.0000 [ 0.4126; 0.4911 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5703    0.0000 [ 0.2448; 0.3484 ] 
#>   eta2 <~ y22      0.4473      0.0196   22.8446    0.0000 [ 0.4128; 0.4788 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1777    0.0000 [ 0.4323; 0.5196 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5802    0.0000 [ 0.4073; 0.4661 ] 
#>   eta3 <~ y32      0.3521      0.0207   17.0306    0.0000 [ 0.3132; 0.3867 ] 
#>   eta3 <~ y33      0.3999      0.0198   20.2243    0.0000 [ 0.3522; 0.4260 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0369   18.1794    0.0000 [ 0.5925; 0.7220 ] 
#>   eta3 ~ eta1       0.6634      0.0395   16.8109    0.0000 [ 0.5894; 0.7165 ] 
#>   eta3 ~ eta2       0.3052      0.0817    3.7365    0.0002 [ 0.1365; 0.4251 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0528    3.8822    0.0001 [ 0.0865; 0.2916 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03976269 16.675679  1.969605e-62
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03634720 17.863217  2.281032e-71
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03495781 21.778992 3.670575e-105
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05826474  8.863934  7.723981e-19
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03850376 19.618542  1.073891e-85
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04174880 19.154171  8.933670e-82
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02950688 27.867314 6.645455e-171
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04629014 14.216180  7.271561e-46
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03627417 20.604856  2.482585e-94
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5972021          0.7407677
#> 2          0.5838440          0.7212359
#> 3          0.7026864          0.8279739
#> 4          0.3890736          0.5816052
#> 5          0.6805072          0.8377654
#> 6          0.7187934          0.8671207
#> 7          0.7698705          0.8721190
#> 8          0.5767855          0.7482126
#> 9          0.6819407          0.8157079

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
#>  Random seed                        = 52038992
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
#>   eta2 ~ eta1      0.6713      0.0369   18.1794    0.0000 [ 0.5857; 0.7767 ] 
#>   eta3 ~ eta1      0.4585      0.0814    5.6298    0.0000 [ 0.2493; 0.6705 ] 
#>   eta3 ~ eta2      0.3052      0.0817    3.7365    0.0002 [ 0.0838; 0.5061 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0398   16.6757    0.0000 [ 0.5584; 0.7640 ] 
#>   eta1 =~ y12      0.6493      0.0363   17.8632    0.0000 [ 0.5579; 0.7459 ] 
#>   eta1 =~ y13      0.7613      0.0350   21.7790    0.0000 [ 0.6692; 0.8499 ] 
#>   eta2 =~ y21      0.5165      0.0583    8.8639    0.0000 [ 0.3756; 0.6769 ] 
#>   eta2 =~ y22      0.7554      0.0385   19.6185    0.0000 [ 0.6589; 0.8580 ] 
#>   eta2 =~ y23      0.7997      0.0417   19.1542    0.0000 [ 0.6910; 0.9069 ] 
#>   eta3 =~ y31      0.8223      0.0295   27.8673    0.0000 [ 0.7470; 0.8996 ] 
#>   eta3 =~ y32      0.6581      0.0463   14.2162    0.0000 [ 0.5360; 0.7754 ] 
#>   eta3 =~ y33      0.7474      0.0363   20.6049    0.0000 [ 0.6435; 0.8311 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0217   18.1888    0.0000 [ 0.3388; 0.4512 ] 
#>   eta1 <~ y12      0.3873      0.0194   19.9240    0.0000 [ 0.3392; 0.4397 ] 
#>   eta1 <~ y13      0.4542      0.0203   22.3251    0.0000 [ 0.4010; 0.5063 ] 
#>   eta2 <~ y21      0.3058      0.0320    9.5703    0.0000 [ 0.2273; 0.3925 ] 
#>   eta2 <~ y22      0.4473      0.0196   22.8446    0.0000 [ 0.3955; 0.4967 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1777    0.0000 [ 0.4092; 0.5305 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5802    0.0000 [ 0.3996; 0.4885 ] 
#>   eta3 <~ y32      0.3521      0.0207   17.0306    0.0000 [ 0.3005; 0.4074 ] 
#>   eta3 <~ y33      0.3999      0.0198   20.2243    0.0000 [ 0.3467; 0.4489 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0369   18.1794    0.0000 [ 0.5857; 0.7767 ] 
#>   eta3 ~ eta1       0.6634      0.0395   16.8109    0.0000 [ 0.5593; 0.7634 ] 
#>   eta3 ~ eta2       0.3052      0.0817    3.7365    0.0002 [ 0.0838; 0.5061 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0528    3.8822    0.0001 [ 0.0651; 0.3380 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03692820 18.179424 7.511838e-74
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08144260  5.629815 1.804028e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08166776  3.736494 1.866038e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.58569472          0.7766665          0.6086267          0.7537345
#> 2         0.24928512          0.6704602          0.2998600          0.6198853
#> 3         0.08379688          0.5061364          0.1345116          0.4554217
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5850493          0.7223758          0.5925360          0.7219955
#> 2          0.3106273          0.6194183          0.3215135          0.6152075
#> 3          0.1197643          0.4513607          0.1365358          0.4250663
```
