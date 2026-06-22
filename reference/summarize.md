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
#>  Random seed                        = 323052898
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
#>   eta2 ~ eta1      0.6713      0.0419   16.0053    0.0000 [ 0.5930; 0.7580 ] 
#>   eta3 ~ eta1      0.4585      0.0798    5.7484    0.0000 [ 0.3217; 0.6404 ] 
#>   eta3 ~ eta2      0.3052      0.0912    3.3461    0.0008 [ 0.1250; 0.4639 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0391   16.9414    0.0000 [ 0.5924; 0.7289 ] 
#>   eta1 =~ y12      0.6493      0.0344   18.8674    0.0000 [ 0.5790; 0.7061 ] 
#>   eta1 =~ y13      0.7613      0.0293   26.0163    0.0000 [ 0.7172; 0.8314 ] 
#>   eta2 =~ y21      0.5165      0.0589    8.7734    0.0000 [ 0.3930; 0.6047 ] 
#>   eta2 =~ y22      0.7554      0.0304   24.8227    0.0000 [ 0.6987; 0.8062 ] 
#>   eta2 =~ y23      0.7997      0.0332   24.0902    0.0000 [ 0.7229; 0.8459 ] 
#>   eta3 =~ y31      0.8223      0.0298   27.5692    0.0000 [ 0.7654; 0.8594 ] 
#>   eta3 =~ y32      0.6581      0.0398   16.5176    0.0000 [ 0.5887; 0.7368 ] 
#>   eta3 =~ y33      0.7474      0.0346   21.6220    0.0000 [ 0.6905; 0.7905 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0166   23.8052    0.0000 [ 0.3601; 0.4186 ] 
#>   eta1 <~ y12      0.3873      0.0169   22.9186    0.0000 [ 0.3531; 0.4187 ] 
#>   eta1 <~ y13      0.4542      0.0215   21.1544    0.0000 [ 0.4184; 0.4952 ] 
#>   eta2 <~ y21      0.3058      0.0296   10.3428    0.0000 [ 0.2540; 0.3608 ] 
#>   eta2 <~ y22      0.4473      0.0219   20.4600    0.0000 [ 0.4114; 0.4854 ] 
#>   eta2 <~ y23      0.4735      0.0201   23.5598    0.0000 [ 0.4436; 0.5059 ] 
#>   eta3 <~ y31      0.4400      0.0176   25.0244    0.0000 [ 0.4048; 0.4649 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1594    0.0000 [ 0.3297; 0.3807 ] 
#>   eta3 <~ y33      0.3999      0.0168   23.8024    0.0000 [ 0.3581; 0.4263 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0419   16.0053    0.0000 [ 0.5930; 0.7580 ] 
#>   eta3 ~ eta1       0.6634      0.0370   17.9129    0.0000 [ 0.6131; 0.7391 ] 
#>   eta3 ~ eta2       0.3052      0.0912    3.3461    0.0008 [ 0.1250; 0.4639 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0622    3.2945    0.0010 [ 0.0804; 0.3073 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03913894 16.941438  2.226174e-64
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03441269 18.867398  2.114570e-79
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02926421 26.016276 3.241013e-149
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05886629  8.773354  1.734224e-18
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03043129 24.822726 5.096416e-136
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03319458 24.090187 3.167771e-128
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02982590 27.569240 2.602315e-167
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03984051 16.517584  2.741617e-61
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03456772 21.622023 1.114821e-103
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5923591          0.7288913
#> 2          0.5790468          0.7061371
#> 3          0.7171849          0.8313635
#> 4          0.3930477          0.6047024
#> 5          0.6987001          0.8062475
#> 6          0.7229347          0.8459304
#> 7          0.7654178          0.8593557
#> 8          0.5886731          0.7368346
#> 9          0.6905475          0.7905109

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
#>  Random seed                        = 323052898
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
#>   eta2 ~ eta1      0.6713      0.0419   16.0053    0.0000 [ 0.5642; 0.7812 ] 
#>   eta3 ~ eta1      0.4585      0.0798    5.7484    0.0000 [ 0.2465; 0.6590 ] 
#>   eta3 ~ eta2      0.3052      0.0912    3.3461    0.0008 [ 0.0725; 0.5442 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0391   16.9414    0.0000 [ 0.5722; 0.7746 ] 
#>   eta1 =~ y12      0.6493      0.0344   18.8674    0.0000 [ 0.5667; 0.7447 ] 
#>   eta1 =~ y13      0.7613      0.0293   26.0163    0.0000 [ 0.6750; 0.8264 ] 
#>   eta2 =~ y21      0.5165      0.0589    8.7734    0.0000 [ 0.3760; 0.6805 ] 
#>   eta2 =~ y22      0.7554      0.0304   24.8227    0.0000 [ 0.6862; 0.8436 ] 
#>   eta2 =~ y23      0.7997      0.0332   24.0902    0.0000 [ 0.7223; 0.8940 ] 
#>   eta3 =~ y31      0.8223      0.0298   27.5692    0.0000 [ 0.7486; 0.9028 ] 
#>   eta3 =~ y32      0.6581      0.0398   16.5176    0.0000 [ 0.5469; 0.7530 ] 
#>   eta3 =~ y33      0.7474      0.0346   21.6220    0.0000 [ 0.6626; 0.8414 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0166   23.8052    0.0000 [ 0.3576; 0.4435 ] 
#>   eta1 <~ y12      0.3873      0.0169   22.9186    0.0000 [ 0.3462; 0.4336 ] 
#>   eta1 <~ y13      0.4542      0.0215   21.1544    0.0000 [ 0.3904; 0.5014 ] 
#>   eta2 <~ y21      0.3058      0.0296   10.3428    0.0000 [ 0.2309; 0.3839 ] 
#>   eta2 <~ y22      0.4473      0.0219   20.4600    0.0000 [ 0.3874; 0.5004 ] 
#>   eta2 <~ y23      0.4735      0.0201   23.5598    0.0000 [ 0.4171; 0.5211 ] 
#>   eta3 <~ y31      0.4400      0.0176   25.0244    0.0000 [ 0.3962; 0.4872 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1594    0.0000 [ 0.3050; 0.3910 ] 
#>   eta3 <~ y33      0.3999      0.0168   23.8024    0.0000 [ 0.3590; 0.4459 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0419   16.0053    0.0000 [ 0.5642; 0.7812 ] 
#>   eta3 ~ eta1       0.6634      0.0370   17.9129    0.0000 [ 0.5647; 0.7562 ] 
#>   eta3 ~ eta2       0.3052      0.0912    3.3461    0.0008 [ 0.0725; 0.5442 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0622    3.2945    0.0010 [ 0.0469; 0.3685 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04194439 16.005321 1.173058e-57
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07976238  5.748409 9.008735e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09119701  3.346065 8.196712e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5642410          0.7811537          0.5902880          0.7551067
#> 2          0.2464736          0.6589596          0.2960051          0.6094281
#> 3          0.0725469          0.5441663          0.1291792          0.4875340
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58146269          0.7690577          0.5929515          0.7579651
#> 2         0.30989795          0.6745311          0.3217450          0.6404081
#> 3         0.08961286          0.4723827          0.1250430          0.4638914
```
