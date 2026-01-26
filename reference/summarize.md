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
#>  Random seed                        = -1022885568
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
#>   eta2 ~ eta1      0.6713      0.0412   16.3113    0.0000 [ 0.5970; 0.7435 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8121    0.0000 [ 0.3634; 0.6200 ] 
#>   eta3 ~ eta2      0.3052      0.0872    3.4993    0.0005 [ 0.1045; 0.4364 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0475   13.9677    0.0000 [ 0.5877; 0.7350 ] 
#>   eta1 =~ y12      0.6493      0.0408   15.9080    0.0000 [ 0.5722; 0.7084 ] 
#>   eta1 =~ y13      0.7613      0.0243   31.2954    0.0000 [ 0.7201; 0.8020 ] 
#>   eta2 =~ y21      0.5165      0.0603    8.5663    0.0000 [ 0.4041; 0.5990 ] 
#>   eta2 =~ y22      0.7554      0.0380   19.9010    0.0000 [ 0.6869; 0.8254 ] 
#>   eta2 =~ y23      0.7997      0.0370   21.5881    0.0000 [ 0.7330; 0.8564 ] 
#>   eta3 =~ y31      0.8223      0.0291   28.2406    0.0000 [ 0.7567; 0.8622 ] 
#>   eta3 =~ y32      0.6581      0.0369   17.8484    0.0000 [ 0.5762; 0.7302 ] 
#>   eta3 =~ y33      0.7474      0.0417   17.9275    0.0000 [ 0.6696; 0.8272 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0242   16.3483    0.0000 [ 0.3570; 0.4309 ] 
#>   eta1 <~ y12      0.3873      0.0236   16.4248    0.0000 [ 0.3429; 0.4298 ] 
#>   eta1 <~ y13      0.4542      0.0167   27.1224    0.0000 [ 0.4295; 0.4919 ] 
#>   eta2 <~ y21      0.3058      0.0314    9.7317    0.0000 [ 0.2513; 0.3481 ] 
#>   eta2 <~ y22      0.4473      0.0230   19.4885    0.0000 [ 0.3985; 0.4801 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7895    0.0000 [ 0.4361; 0.5128 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5786    0.0000 [ 0.4076; 0.4692 ] 
#>   eta3 <~ y32      0.3521      0.0172   20.4567    0.0000 [ 0.3247; 0.3836 ] 
#>   eta3 <~ y33      0.3999      0.0197   20.2841    0.0000 [ 0.3693; 0.4374 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0412   16.3113    0.0000 [ 0.5970; 0.7435 ] 
#>   eta3 ~ eta1       0.6634      0.0339   19.5414    0.0000 [ 0.6073; 0.7520 ] 
#>   eta3 ~ eta2       0.3052      0.0872    3.4993    0.0005 [ 0.1045; 0.4364 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0594    3.4462    0.0006 [ 0.0732; 0.2803 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04747181 13.967656  2.455771e-44
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04081456 15.907999  5.576761e-57
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02432769 31.295442 5.382617e-215
#> 4 eta2 =~ y21  Common factor 0.5164548 0.06028891  8.566333  1.068331e-17
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03795734 19.900965  3.991656e-88
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03704187 21.588105 2.323447e-103
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02911687 28.240586 1.857532e-175
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03686984 17.848433  2.972612e-71
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04169147 17.927506  7.193514e-72
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5877206          0.7349740
#> 2          0.5721692          0.7084450
#> 3          0.7201243          0.8020405
#> 4          0.4040739          0.5989763
#> 5          0.6869211          0.8254403
#> 6          0.7330244          0.8564005
#> 7          0.7567071          0.8621813
#> 8          0.5762205          0.7301662
#> 9          0.6695769          0.8271990

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
#>  Random seed                        = -1022885568
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
#>   eta2 ~ eta1      0.6713      0.0412   16.3113    0.0000 [ 0.5572; 0.7700 ] 
#>   eta3 ~ eta1      0.4585      0.0789    5.8121    0.0000 [ 0.2370; 0.6450 ] 
#>   eta3 ~ eta2      0.3052      0.0872    3.4993    0.0005 [ 0.0971; 0.5481 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0475   13.9677    0.0000 [ 0.5492; 0.7947 ] 
#>   eta1 =~ y12      0.6493      0.0408   15.9080    0.0000 [ 0.5471; 0.7582 ] 
#>   eta1 =~ y13      0.7613      0.0243   31.2954    0.0000 [ 0.6996; 0.8254 ] 
#>   eta2 =~ y21      0.5165      0.0603    8.5663    0.0000 [ 0.3598; 0.6716 ] 
#>   eta2 =~ y22      0.7554      0.0380   19.9010    0.0000 [ 0.6649; 0.8612 ] 
#>   eta2 =~ y23      0.7997      0.0370   21.5881    0.0000 [ 0.7036; 0.8951 ] 
#>   eta3 =~ y31      0.8223      0.0291   28.2406    0.0000 [ 0.7509; 0.9015 ] 
#>   eta3 =~ y32      0.6581      0.0369   17.8484    0.0000 [ 0.5582; 0.7489 ] 
#>   eta3 =~ y33      0.7474      0.0417   17.9275    0.0000 [ 0.6419; 0.8575 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0242   16.3483    0.0000 [ 0.3350; 0.4602 ] 
#>   eta1 <~ y12      0.3873      0.0236   16.4248    0.0000 [ 0.3250; 0.4469 ] 
#>   eta1 <~ y13      0.4542      0.0167   27.1224    0.0000 [ 0.4074; 0.4940 ] 
#>   eta2 <~ y21      0.3058      0.0314    9.7317    0.0000 [ 0.2234; 0.3859 ] 
#>   eta2 <~ y22      0.4473      0.0230   19.4885    0.0000 [ 0.3907; 0.5094 ] 
#>   eta2 <~ y23      0.4735      0.0217   21.7895    0.0000 [ 0.4153; 0.5276 ] 
#>   eta3 <~ y31      0.4400      0.0172   25.5786    0.0000 [ 0.3971; 0.4860 ] 
#>   eta3 <~ y32      0.3521      0.0172   20.4567    0.0000 [ 0.3049; 0.3939 ] 
#>   eta3 <~ y33      0.3999      0.0197   20.2841    0.0000 [ 0.3499; 0.4518 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0412   16.3113    0.0000 [ 0.5572; 0.7700 ] 
#>   eta3 ~ eta1       0.6634      0.0339   19.5414    0.0000 [ 0.5679; 0.7434 ] 
#>   eta3 ~ eta2       0.3052      0.0872    3.4993    0.0005 [ 0.0971; 0.5481 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0594    3.4462    0.0006 [ 0.0609; 0.3683 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04115768 16.311254 8.209041e-60
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07888806  5.812119 6.168710e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08720320  3.499311 4.664620e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55717293          0.7700172          0.5827314          0.7444588
#> 2         0.23704884          0.6450133          0.2860374          0.5960248
#> 3         0.09713142          0.5480971          0.1512836          0.4939450
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5482299          0.7455560          0.5970357          0.7435157
#> 2          0.3627213          0.6869485          0.3634105          0.6199617
#> 3          0.1038958          0.4459362          0.1045096          0.4364100
```
