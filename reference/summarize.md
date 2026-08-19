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
#>  Random seed                        = 548689761
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
#>   eta2 ~ eta1      0.6713      0.0400   16.7898    0.0000 [ 0.5876; 0.7304 ] 
#>   eta3 ~ eta1      0.4585      0.0675    6.7945    0.0000 [ 0.3233; 0.5438 ] 
#>   eta3 ~ eta2      0.3052      0.0759    4.0220    0.0001 [ 0.2008; 0.4226 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0321   20.6619    0.0000 [ 0.6164; 0.7261 ] 
#>   eta1 =~ y12      0.6493      0.0419   15.4920    0.0000 [ 0.5974; 0.7252 ] 
#>   eta1 =~ y13      0.7613      0.0361   21.1019    0.0000 [ 0.7001; 0.8455 ] 
#>   eta2 =~ y21      0.5165      0.0497   10.3831    0.0000 [ 0.4173; 0.5828 ] 
#>   eta2 =~ y22      0.7554      0.0358   21.0798    0.0000 [ 0.7024; 0.8232 ] 
#>   eta2 =~ y23      0.7997      0.0326   24.5403    0.0000 [ 0.7495; 0.8590 ] 
#>   eta3 =~ y31      0.8223      0.0351   23.4524    0.0000 [ 0.7698; 0.8702 ] 
#>   eta3 =~ y32      0.6581      0.0445   14.7859    0.0000 [ 0.5667; 0.7603 ] 
#>   eta3 =~ y33      0.7474      0.0400   18.6736    0.0000 [ 0.6709; 0.8040 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0170   23.3270    0.0000 [ 0.3653; 0.4230 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5334    0.0000 [ 0.3521; 0.4258 ] 
#>   eta1 <~ y13      0.4542      0.0211   21.5197    0.0000 [ 0.4166; 0.4893 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3756    0.0000 [ 0.2561; 0.3386 ] 
#>   eta2 <~ y22      0.4473      0.0227   19.7209    0.0000 [ 0.4162; 0.4923 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.5149    0.0000 [ 0.4467; 0.5201 ] 
#>   eta3 <~ y31      0.4400      0.0210   20.9051    0.0000 [ 0.4075; 0.4874 ] 
#>   eta3 <~ y32      0.3521      0.0205   17.2025    0.0000 [ 0.3189; 0.4000 ] 
#>   eta3 <~ y33      0.3999      0.0182   21.9917    0.0000 [ 0.3668; 0.4351 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0400   16.7898    0.0000 [ 0.5876; 0.7304 ] 
#>   eta3 ~ eta1       0.6634      0.0365   18.1900    0.0000 [ 0.5832; 0.7289 ] 
#>   eta3 ~ eta2       0.3052      0.0759    4.0220    0.0001 [ 0.2008; 0.4226 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0495    4.1423    0.0000 [ 0.1389; 0.2944 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03209141 20.66191  7.628365e-95
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04191044 15.49203  3.926649e-54
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03607948 21.10191  7.639068e-99
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04973974 10.38314  2.958721e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03583465 21.07981  1.218688e-98
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03258578 24.54027 5.493905e-133
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03506153 23.45241 1.249213e-121
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04450641 14.78593  1.805330e-49
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04002567 18.67362  8.116221e-78
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6164446          0.7260687
#> 2          0.5973674          0.7251592
#> 3          0.7000707          0.8455435
#> 4          0.4172541          0.5828380
#> 5          0.7024202          0.8232454
#> 6          0.7495133          0.8589941
#> 7          0.7697708          0.8701712
#> 8          0.5666872          0.7603468
#> 9          0.6709249          0.8040454

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
#>  Random seed                        = 548689761
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
#>   eta2 ~ eta1      0.6713      0.0400   16.7898    0.0000 [ 0.5671; 0.7739 ] 
#>   eta3 ~ eta1      0.4585      0.0675    6.7945    0.0000 [ 0.2940; 0.6429 ] 
#>   eta3 ~ eta2      0.3052      0.0759    4.0220    0.0001 [ 0.1004; 0.4927 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0321   20.6619    0.0000 [ 0.5809; 0.7468 ] 
#>   eta1 =~ y12      0.6493      0.0419   15.4920    0.0000 [ 0.5317; 0.7484 ] 
#>   eta1 =~ y13      0.7613      0.0361   21.1019    0.0000 [ 0.6653; 0.8519 ] 
#>   eta2 =~ y21      0.5165      0.0497   10.3831    0.0000 [ 0.3981; 0.6553 ] 
#>   eta2 =~ y22      0.7554      0.0358   21.0798    0.0000 [ 0.6567; 0.8420 ] 
#>   eta2 =~ y23      0.7997      0.0326   24.5403    0.0000 [ 0.7128; 0.8813 ] 
#>   eta3 =~ y31      0.8223      0.0351   23.4524    0.0000 [ 0.7299; 0.9112 ] 
#>   eta3 =~ y32      0.6581      0.0445   14.7859    0.0000 [ 0.5366; 0.7667 ] 
#>   eta3 =~ y33      0.7474      0.0400   18.6736    0.0000 [ 0.6510; 0.8580 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0170   23.3270    0.0000 [ 0.3555; 0.4432 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5334    0.0000 [ 0.3312; 0.4393 ] 
#>   eta1 <~ y13      0.4542      0.0211   21.5197    0.0000 [ 0.4017; 0.5108 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3756    0.0000 [ 0.2487; 0.3764 ] 
#>   eta2 <~ y22      0.4473      0.0227   19.7209    0.0000 [ 0.3854; 0.5027 ] 
#>   eta2 <~ y23      0.4735      0.0193   24.5149    0.0000 [ 0.4224; 0.5223 ] 
#>   eta3 <~ y31      0.4400      0.0210   20.9051    0.0000 [ 0.3849; 0.4937 ] 
#>   eta3 <~ y32      0.3521      0.0205   17.2025    0.0000 [ 0.2963; 0.4021 ] 
#>   eta3 <~ y33      0.3999      0.0182   21.9917    0.0000 [ 0.3572; 0.4512 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0400   16.7898    0.0000 [ 0.5671; 0.7739 ] 
#>   eta3 ~ eta1       0.6634      0.0365   18.1900    0.0000 [ 0.5738; 0.7624 ] 
#>   eta3 ~ eta2       0.3052      0.0759    4.0220    0.0001 [ 0.1004; 0.4927 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0495    4.1423    0.0000 [ 0.0718; 0.3276 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03998450 16.789843 2.896115e-63
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06748250  6.794454 1.087234e-11
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07587049  4.022000 5.770597e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5671166          0.7738939          0.5919465          0.7490640
#> 2          0.2939592          0.6429406          0.3358650          0.6010348
#> 3          0.1003811          0.4927404          0.1474958          0.4456258
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5818234          0.7428559          0.5875978          0.7303533
#> 2          0.2931461          0.5521410          0.3233047          0.5438442
#> 3          0.1542021          0.4375582          0.2008045          0.4225855
```
