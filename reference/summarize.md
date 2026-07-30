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
#>  Random seed                        = -1624992484
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
#>   eta2 ~ eta1      0.6713      0.0472   14.2368    0.0000 [ 0.6002; 0.7728 ] 
#>   eta3 ~ eta1      0.4585      0.0857    5.3486    0.0000 [ 0.3163; 0.6115 ] 
#>   eta3 ~ eta2      0.3052      0.0871    3.5023    0.0005 [ 0.1430; 0.4530 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0442   14.9961    0.0000 [ 0.6010; 0.7540 ] 
#>   eta1 =~ y12      0.6493      0.0428   15.1606    0.0000 [ 0.5696; 0.7334 ] 
#>   eta1 =~ y13      0.7613      0.0374   20.3413    0.0000 [ 0.7002; 0.8417 ] 
#>   eta2 =~ y21      0.5165      0.0509   10.1544    0.0000 [ 0.4192; 0.5998 ] 
#>   eta2 =~ y22      0.7554      0.0450   16.7863    0.0000 [ 0.6920; 0.8541 ] 
#>   eta2 =~ y23      0.7997      0.0355   22.5216    0.0000 [ 0.7213; 0.8555 ] 
#>   eta3 =~ y31      0.8223      0.0317   25.9688    0.0000 [ 0.7744; 0.8838 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.6965    0.0000 [ 0.5750; 0.7093 ] 
#>   eta3 =~ y33      0.7474      0.0380   19.6629    0.0000 [ 0.6785; 0.8080 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0265   14.9181    0.0000 [ 0.3433; 0.4345 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9582    0.0000 [ 0.3401; 0.4254 ] 
#>   eta1 <~ y13      0.4542      0.0226   20.1317    0.0000 [ 0.4089; 0.4979 ] 
#>   eta2 <~ y21      0.3058      0.0254   12.0597    0.0000 [ 0.2574; 0.3515 ] 
#>   eta2 <~ y22      0.4473      0.0271   16.4799    0.0000 [ 0.3962; 0.5056 ] 
#>   eta2 <~ y23      0.4735      0.0211   22.4612    0.0000 [ 0.4333; 0.5118 ] 
#>   eta3 <~ y31      0.4400      0.0174   25.2939    0.0000 [ 0.4103; 0.4776 ] 
#>   eta3 <~ y32      0.3521      0.0199   17.7362    0.0000 [ 0.3105; 0.3818 ] 
#>   eta3 <~ y33      0.3999      0.0199   20.1380    0.0000 [ 0.3605; 0.4302 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0472   14.2368    0.0000 [ 0.6002; 0.7728 ] 
#>   eta3 ~ eta1       0.6634      0.0386   17.1800    0.0000 [ 0.5859; 0.7391 ] 
#>   eta3 ~ eta2       0.3052      0.0871    3.5023    0.0005 [ 0.1430; 0.4530 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0568    3.6073    0.0003 [ 0.0982; 0.2998 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04421612 14.99611  7.785024e-51
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04282668 15.16059  6.448364e-52
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03742850 20.34134  5.539567e-92
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05086005 10.15443  3.166557e-24
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04500037 16.78625  3.076593e-63
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03550646 22.52164 2.547466e-112
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03166399 25.96884 1.114026e-148
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04192446 15.69654  1.597187e-55
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03801188 19.66291  4.482774e-86
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6010499          0.7540092
#> 2          0.5696365          0.7334210
#> 3          0.7002101          0.8417000
#> 4          0.4192109          0.5998424
#> 5          0.6919693          0.8540826
#> 6          0.7212963          0.8554697
#> 7          0.7743917          0.8838324
#> 8          0.5750004          0.7093353
#> 9          0.6785484          0.8080432

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
#>  Random seed                        = -1624992484
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
#>   eta2 ~ eta1      0.6713      0.0472   14.2368    0.0000 [ 0.5410; 0.7848 ] 
#>   eta3 ~ eta1      0.4585      0.0857    5.3486    0.0000 [ 0.2229; 0.6662 ] 
#>   eta3 ~ eta2      0.3052      0.0871    3.5023    0.0005 [ 0.0972; 0.5478 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0442   14.9961    0.0000 [ 0.5506; 0.7792 ] 
#>   eta1 =~ y12      0.6493      0.0428   15.1606    0.0000 [ 0.5404; 0.7618 ] 
#>   eta1 =~ y13      0.7613      0.0374   20.3413    0.0000 [ 0.6565; 0.8500 ] 
#>   eta2 =~ y21      0.5165      0.0509   10.1544    0.0000 [ 0.3867; 0.6497 ] 
#>   eta2 =~ y22      0.7554      0.0450   16.7863    0.0000 [ 0.6400; 0.8727 ] 
#>   eta2 =~ y23      0.7997      0.0355   22.5216    0.0000 [ 0.7171; 0.9007 ] 
#>   eta3 =~ y31      0.8223      0.0317   25.9688    0.0000 [ 0.7356; 0.8993 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.6965    0.0000 [ 0.5521; 0.7689 ] 
#>   eta3 =~ y33      0.7474      0.0380   19.6629    0.0000 [ 0.6535; 0.8500 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0265   14.9181    0.0000 [ 0.3296; 0.4668 ] 
#>   eta1 <~ y12      0.3873      0.0216   17.9582    0.0000 [ 0.3344; 0.4460 ] 
#>   eta1 <~ y13      0.4542      0.0226   20.1317    0.0000 [ 0.3929; 0.5095 ] 
#>   eta2 <~ y21      0.3058      0.0254   12.0597    0.0000 [ 0.2392; 0.3704 ] 
#>   eta2 <~ y22      0.4473      0.0271   16.4799    0.0000 [ 0.3740; 0.5144 ] 
#>   eta2 <~ y23      0.4735      0.0211   22.4612    0.0000 [ 0.4207; 0.5297 ] 
#>   eta3 <~ y31      0.4400      0.0174   25.2939    0.0000 [ 0.3921; 0.4820 ] 
#>   eta3 <~ y32      0.3521      0.0199   17.7362    0.0000 [ 0.3020; 0.4046 ] 
#>   eta3 <~ y33      0.3999      0.0199   20.1380    0.0000 [ 0.3506; 0.4533 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0472   14.2368    0.0000 [ 0.5410; 0.7848 ] 
#>   eta3 ~ eta1       0.6634      0.0386   17.1800    0.0000 [ 0.5604; 0.7601 ] 
#>   eta3 ~ eta2       0.3052      0.0871    3.5023    0.0005 [ 0.0972; 0.5478 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0568    3.6073    0.0003 [ 0.0688; 0.3625 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04715464 14.236849 5.411264e-46
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08572449  5.348609 8.863270e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08712822  3.502323 4.612207e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.54099016          0.7848473          0.5702726          0.7555648
#> 2         0.22288563          0.6662043          0.2761195          0.6129704
#> 3         0.09718482          0.5477628          0.1512904          0.4936572
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58423042          0.7786282          0.6001733          0.7727705
#> 2         0.28755354          0.6759235          0.3163103          0.6115360
#> 3         0.08652417          0.4608133          0.1430343          0.4529830
```
