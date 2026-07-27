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
#>  Random seed                        = -1426784409
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
#>   eta2 ~ eta1      0.6713      0.0442   15.1731    0.0000 [ 0.6031; 0.7427 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8113    0.0000 [ 0.3048; 0.6071 ] 
#>   eta3 ~ eta2      0.3052      0.0866    3.5217    0.0004 [ 0.1497; 0.4749 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0311   21.3357    0.0000 [ 0.6118; 0.7174 ] 
#>   eta1 =~ y12      0.6493      0.0337   19.2434    0.0000 [ 0.5887; 0.7064 ] 
#>   eta1 =~ y13      0.7613      0.0282   27.0401    0.0000 [ 0.7212; 0.8106 ] 
#>   eta2 =~ y21      0.5165      0.0543    9.5082    0.0000 [ 0.4155; 0.6124 ] 
#>   eta2 =~ y22      0.7554      0.0286   26.4501    0.0000 [ 0.7093; 0.8028 ] 
#>   eta2 =~ y23      0.7997      0.0328   24.3722    0.0000 [ 0.7340; 0.8625 ] 
#>   eta3 =~ y31      0.8223      0.0379   21.6782    0.0000 [ 0.7579; 0.8969 ] 
#>   eta3 =~ y32      0.6581      0.0376   17.5053    0.0000 [ 0.5858; 0.7331 ] 
#>   eta3 =~ y33      0.7474      0.0404   18.5044    0.0000 [ 0.6747; 0.8083 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0178   22.2523    0.0000 [ 0.3534; 0.4246 ] 
#>   eta1 <~ y12      0.3873      0.0164   23.6632    0.0000 [ 0.3569; 0.4183 ] 
#>   eta1 <~ y13      0.4542      0.0159   28.5168    0.0000 [ 0.4277; 0.4827 ] 
#>   eta2 <~ y21      0.3058      0.0250   12.2313    0.0000 [ 0.2640; 0.3512 ] 
#>   eta2 <~ y22      0.4473      0.0217   20.5967    0.0000 [ 0.4166; 0.4899 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.4894    0.0000 [ 0.4380; 0.5047 ] 
#>   eta3 <~ y31      0.4400      0.0197   22.3868    0.0000 [ 0.4062; 0.4724 ] 
#>   eta3 <~ y32      0.3521      0.0209   16.8536    0.0000 [ 0.3219; 0.3932 ] 
#>   eta3 <~ y33      0.3999      0.0178   22.4392    0.0000 [ 0.3767; 0.4344 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0442   15.1731    0.0000 [ 0.6031; 0.7427 ] 
#>   eta3 ~ eta1       0.6634      0.0412   16.1164    0.0000 [ 0.5737; 0.7430 ] 
#>   eta3 ~ eta2       0.3052      0.0866    3.5217    0.0004 [ 0.1497; 0.4749 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0567    3.6160    0.0003 [ 0.0946; 0.2984 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03107790 21.335734 5.290768e-101
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03374026 19.243418  1.602816e-82
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02815613 27.040139 4.988647e-161
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05431705  9.508153  1.940777e-21
#> 5 eta2 =~ y22  Common factor 0.7553877 0.02855899 26.450084 3.640625e-154
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03281046 24.372219 3.370981e-131
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03793109 21.678188 3.295991e-104
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03759266 17.505250  1.306516e-68
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04039170 18.504400  1.902847e-76
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6118103          0.7173887
#> 2          0.5886583          0.7064199
#> 3          0.7212187          0.8106309
#> 4          0.4155365          0.6123876
#> 5          0.7092853          0.8027919
#> 6          0.7339690          0.8625369
#> 7          0.7579345          0.8968966
#> 8          0.5858303          0.7331485
#> 9          0.6747120          0.8083028

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
#>  Random seed                        = -1426784409
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
#>   eta2 ~ eta1      0.6713      0.0442   15.1731    0.0000 [ 0.5594; 0.7882 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8113    0.0000 [ 0.2499; 0.6579 ] 
#>   eta3 ~ eta2      0.3052      0.0866    3.5217    0.0004 [ 0.0781; 0.5262 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0311   21.3357    0.0000 [ 0.5795; 0.7402 ] 
#>   eta1 =~ y12      0.6493      0.0337   19.2434    0.0000 [ 0.5612; 0.7356 ] 
#>   eta1 =~ y13      0.7613      0.0282   27.0401    0.0000 [ 0.6835; 0.8291 ] 
#>   eta2 =~ y21      0.5165      0.0543    9.5082    0.0000 [ 0.3749; 0.6558 ] 
#>   eta2 =~ y22      0.7554      0.0286   26.4501    0.0000 [ 0.6819; 0.8296 ] 
#>   eta2 =~ y23      0.7997      0.0328   24.3722    0.0000 [ 0.7199; 0.8895 ] 
#>   eta3 =~ y31      0.8223      0.0379   21.6782    0.0000 [ 0.7281; 0.9243 ] 
#>   eta3 =~ y32      0.6581      0.0376   17.5053    0.0000 [ 0.5606; 0.7550 ] 
#>   eta3 =~ y33      0.7474      0.0404   18.5044    0.0000 [ 0.6473; 0.8562 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0178   22.2523    0.0000 [ 0.3504; 0.4423 ] 
#>   eta1 <~ y12      0.3873      0.0164   23.6632    0.0000 [ 0.3473; 0.4319 ] 
#>   eta1 <~ y13      0.4542      0.0159   28.5168    0.0000 [ 0.4131; 0.4955 ] 
#>   eta2 <~ y21      0.3058      0.0250   12.2313    0.0000 [ 0.2403; 0.3696 ] 
#>   eta2 <~ y22      0.4473      0.0217   20.5967    0.0000 [ 0.3899; 0.5022 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.4894    0.0000 [ 0.4252; 0.5252 ] 
#>   eta3 <~ y31      0.4400      0.0197   22.3868    0.0000 [ 0.3890; 0.4906 ] 
#>   eta3 <~ y32      0.3521      0.0209   16.8536    0.0000 [ 0.2961; 0.4041 ] 
#>   eta3 <~ y33      0.3999      0.0178   22.4392    0.0000 [ 0.3543; 0.4464 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0442   15.1731    0.0000 [ 0.5594; 0.7882 ] 
#>   eta3 ~ eta1       0.6634      0.0412   16.1164    0.0000 [ 0.5517; 0.7645 ] 
#>   eta3 ~ eta2       0.3052      0.0866    3.5217    0.0004 [ 0.0781; 0.5262 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0567    3.6160    0.0003 [ 0.0577; 0.3507 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04424498 15.173099 5.329931e-52
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07889879  5.811328 6.197922e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08664915  3.521686 4.288115e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55938336          0.7881934          0.5868590          0.7607178
#> 2         0.24989210          0.6579121          0.2988873          0.6089169
#> 3         0.07810909          0.5262096          0.1319172          0.4724015
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5570760          0.7613111          0.6030870          0.7426664
#> 2          0.2879003          0.6243140          0.3048196          0.6070739
#> 3          0.1181642          0.5009258          0.1496769          0.4748726
```
