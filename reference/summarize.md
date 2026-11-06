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
#>  Random seed                        = 650698612
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
#>   eta2 ~ eta1      0.6713      0.0376   17.8605    0.0000 [ 0.6170; 0.7519 ] 
#>   eta3 ~ eta1      0.4585      0.0884    5.1845    0.0000 [ 0.3234; 0.6513 ] 
#>   eta3 ~ eta2      0.3052      0.0938    3.2528    0.0011 [ 0.1139; 0.4556 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0387   17.1131    0.0000 [ 0.5913; 0.7315 ] 
#>   eta1 =~ y12      0.6493      0.0403   16.1040    0.0000 [ 0.5661; 0.7084 ] 
#>   eta1 =~ y13      0.7613      0.0313   24.3322    0.0000 [ 0.6939; 0.8102 ] 
#>   eta2 =~ y21      0.5165      0.0525    9.8313    0.0000 [ 0.4092; 0.6020 ] 
#>   eta2 =~ y22      0.7554      0.0364   20.7628    0.0000 [ 0.6673; 0.8074 ] 
#>   eta2 =~ y23      0.7997      0.0419   19.0943    0.0000 [ 0.7250; 0.8698 ] 
#>   eta3 =~ y31      0.8223      0.0332   24.7606    0.0000 [ 0.7565; 0.8706 ] 
#>   eta3 =~ y32      0.6581      0.0345   19.0557    0.0000 [ 0.6023; 0.7119 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2280    0.0000 [ 0.6827; 0.8444 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0230   17.2256    0.0000 [ 0.3569; 0.4454 ] 
#>   eta1 <~ y12      0.3873      0.0176   22.0687    0.0000 [ 0.3536; 0.4147 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.7911    0.0000 [ 0.4215; 0.4945 ] 
#>   eta2 <~ y21      0.3058      0.0288   10.6367    0.0000 [ 0.2475; 0.3554 ] 
#>   eta2 <~ y22      0.4473      0.0189   23.6648    0.0000 [ 0.4138; 0.4793 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1907    0.0000 [ 0.4364; 0.5161 ] 
#>   eta3 <~ y31      0.4400      0.0169   25.9936    0.0000 [ 0.3995; 0.4643 ] 
#>   eta3 <~ y32      0.3521      0.0170   20.7488    0.0000 [ 0.3348; 0.3875 ] 
#>   eta3 <~ y33      0.3999      0.0189   21.1329    0.0000 [ 0.3692; 0.4371 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0376   17.8605    0.0000 [ 0.6170; 0.7519 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7278    0.0000 [ 0.6005; 0.7270 ] 
#>   eta3 ~ eta2       0.3052      0.0938    3.2528    0.0011 [ 0.1139; 0.4556 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0618    3.3173    0.0009 [ 0.0742; 0.2953 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03874642 17.113060  1.186017e-65
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04031783 16.103988  2.391785e-58
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03128963 24.332210 8.945444e-131
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05253172  9.831295  8.255074e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03638178 20.762799  9.393248e-96
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04187974 19.094286  2.816847e-81
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03320910 24.760606 2.383359e-135
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03453400 19.055682  5.894301e-81
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04100417 18.228004  3.094076e-74
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5913065          0.7314933
#> 2          0.5661136          0.7084161
#> 3          0.6938654          0.8102472
#> 4          0.4091701          0.6019524
#> 5          0.6673131          0.8074468
#> 6          0.7249511          0.8697544
#> 7          0.7565330          0.8706117
#> 8          0.6022716          0.7118808
#> 9          0.6827193          0.8443663

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
#>  Random seed                        = 650698612
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
#>   eta2 ~ eta1      0.6713      0.0376   17.8605    0.0000 [ 0.5614; 0.7558 ] 
#>   eta3 ~ eta1      0.4585      0.0884    5.1845    0.0000 [ 0.2090; 0.6663 ] 
#>   eta3 ~ eta2      0.3052      0.0938    3.2528    0.0011 [ 0.0850; 0.5701 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0387   17.1131    0.0000 [ 0.5677; 0.7680 ] 
#>   eta1 =~ y12      0.6493      0.0403   16.1040    0.0000 [ 0.5482; 0.7567 ] 
#>   eta1 =~ y13      0.7613      0.0313   24.3322    0.0000 [ 0.6818; 0.8436 ] 
#>   eta2 =~ y21      0.5165      0.0525    9.8313    0.0000 [ 0.3724; 0.6440 ] 
#>   eta2 =~ y22      0.7554      0.0364   20.7628    0.0000 [ 0.6676; 0.8557 ] 
#>   eta2 =~ y23      0.7997      0.0419   19.0943    0.0000 [ 0.6975; 0.9140 ] 
#>   eta3 =~ y31      0.8223      0.0332   24.7606    0.0000 [ 0.7437; 0.9155 ] 
#>   eta3 =~ y32      0.6581      0.0345   19.0557    0.0000 [ 0.5615; 0.7401 ] 
#>   eta3 =~ y33      0.7474      0.0410   18.2280    0.0000 [ 0.6422; 0.8543 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0230   17.2256    0.0000 [ 0.3366; 0.4553 ] 
#>   eta1 <~ y12      0.3873      0.0176   22.0687    0.0000 [ 0.3418; 0.4325 ] 
#>   eta1 <~ y13      0.4542      0.0208   21.7911    0.0000 [ 0.3981; 0.5059 ] 
#>   eta2 <~ y21      0.3058      0.0288   10.6367    0.0000 [ 0.2258; 0.3745 ] 
#>   eta2 <~ y22      0.4473      0.0189   23.6648    0.0000 [ 0.4007; 0.4985 ] 
#>   eta2 <~ y23      0.4735      0.0235   20.1907    0.0000 [ 0.4149; 0.5362 ] 
#>   eta3 <~ y31      0.4400      0.0169   25.9936    0.0000 [ 0.3998; 0.4874 ] 
#>   eta3 <~ y32      0.3521      0.0170   20.7488    0.0000 [ 0.3041; 0.3919 ] 
#>   eta3 <~ y33      0.3999      0.0189   21.1329    0.0000 [ 0.3513; 0.4491 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0376   17.8605    0.0000 [ 0.5614; 0.7558 ] 
#>   eta3 ~ eta1       0.6634      0.0354   18.7278    0.0000 [ 0.5632; 0.7464 ] 
#>   eta3 ~ eta2       0.3052      0.0938    3.2528    0.0011 [ 0.0850; 0.5701 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0618    3.3173    0.0009 [ 0.0575; 0.3769 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03758756 17.860522 2.393895e-71
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08843756  5.184525 2.165657e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09381166  3.252806 1.142715e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56137543          0.7557571          0.5847169          0.7324156
#> 2         0.20896604          0.6663152          0.2638847          0.6113965
#> 3         0.08495901          0.5700999          0.1432150          0.5118440
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59666807          0.7695673          0.6169803          0.7518558
#> 2         0.30484089          0.6579970          0.3233986          0.6513390
#> 3         0.08660898          0.4597294          0.1139061          0.4555865
```
