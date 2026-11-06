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
#>  Random seed                        = -1863265758
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
#>   eta2 ~ eta1      0.6713      0.0467   14.3837    0.0000 [ 0.5884; 0.7583 ] 
#>   eta3 ~ eta1      0.4585      0.0707    6.4858    0.0000 [ 0.3459; 0.5784 ] 
#>   eta3 ~ eta2      0.3052      0.0779    3.9183    0.0001 [ 0.1772; 0.4364 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0365   18.1501    0.0000 [ 0.6070; 0.7281 ] 
#>   eta1 =~ y12      0.6493      0.0467   13.9103    0.0000 [ 0.5704; 0.7297 ] 
#>   eta1 =~ y13      0.7613      0.0280   27.1489    0.0000 [ 0.7231; 0.8185 ] 
#>   eta2 =~ y21      0.5165      0.0561    9.2113    0.0000 [ 0.4137; 0.6077 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1636    0.0000 [ 0.6977; 0.8312 ] 
#>   eta2 =~ y23      0.7997      0.0294   27.1998    0.0000 [ 0.7508; 0.8572 ] 
#>   eta3 =~ y31      0.8223      0.0307   26.7928    0.0000 [ 0.7790; 0.8900 ] 
#>   eta3 =~ y32      0.6581      0.0460   14.3060    0.0000 [ 0.5572; 0.7342 ] 
#>   eta3 =~ y33      0.7474      0.0386   19.3594    0.0000 [ 0.6899; 0.8197 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0194   20.4264    0.0000 [ 0.3507; 0.4247 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5766    0.0000 [ 0.3532; 0.4227 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8399    0.0000 [ 0.4243; 0.4970 ] 
#>   eta2 <~ y21      0.3058      0.0273   11.2070    0.0000 [ 0.2459; 0.3432 ] 
#>   eta2 <~ y22      0.4473      0.0254   17.5951    0.0000 [ 0.4097; 0.5038 ] 
#>   eta2 <~ y23      0.4735      0.0169   28.0380    0.0000 [ 0.4438; 0.4987 ] 
#>   eta3 <~ y31      0.4400      0.0189   23.3326    0.0000 [ 0.4171; 0.4813 ] 
#>   eta3 <~ y32      0.3521      0.0202   17.4421    0.0000 [ 0.3067; 0.3775 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.3833    0.0000 [ 0.3633; 0.4359 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0467   14.3837    0.0000 [ 0.5884; 0.7583 ] 
#>   eta3 ~ eta1       0.6634      0.0330   20.1128    0.0000 [ 0.5964; 0.7221 ] 
#>   eta3 ~ eta2       0.3052      0.0779    3.9183    0.0001 [ 0.1772; 0.4364 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0545    3.7610    0.0002 [ 0.1048; 0.3083 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03653257 18.150103  1.281614e-73
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04667609 13.910288  5.485789e-44
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02804335 27.148886 2.609839e-162
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05606750  9.211304  3.221876e-20
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03569280 21.163588  2.068663e-99
#> 6 eta2 =~ y23  Common factor 0.7996637 0.02939957 27.199843 6.522735e-163
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03069025 26.792787 3.921279e-158
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04599956 14.305984  2.007787e-46
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03860775 19.359434  1.697388e-83
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6070481          0.7281384
#> 2          0.5704034          0.7297282
#> 3          0.7230827          0.8184823
#> 4          0.4136947          0.6076773
#> 5          0.6977311          0.8312112
#> 6          0.7508255          0.8572356
#> 7          0.7790333          0.8900352
#> 8          0.5572032          0.7341834
#> 9          0.6899423          0.8197123

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
#>  Random seed                        = -1863265758
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
#>   eta2 ~ eta1      0.6713      0.0467   14.3837    0.0000 [ 0.5504; 0.7918 ] 
#>   eta3 ~ eta1      0.4585      0.0707    6.4858    0.0000 [ 0.2680; 0.6336 ] 
#>   eta3 ~ eta2      0.3052      0.0779    3.9183    0.0001 [ 0.1089; 0.5117 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0365   18.1501    0.0000 [ 0.5673; 0.7562 ] 
#>   eta1 =~ y12      0.6493      0.0467   13.9103    0.0000 [ 0.5254; 0.7668 ] 
#>   eta1 =~ y13      0.7613      0.0280   27.1489    0.0000 [ 0.6826; 0.8276 ] 
#>   eta2 =~ y21      0.5165      0.0561    9.2113    0.0000 [ 0.3715; 0.6614 ] 
#>   eta2 =~ y22      0.7554      0.0357   21.1636    0.0000 [ 0.6541; 0.8386 ] 
#>   eta2 =~ y23      0.7997      0.0294   27.1998    0.0000 [ 0.7207; 0.8728 ] 
#>   eta3 =~ y31      0.8223      0.0307   26.7928    0.0000 [ 0.7328; 0.8915 ] 
#>   eta3 =~ y32      0.6581      0.0460   14.3060    0.0000 [ 0.5437; 0.7816 ] 
#>   eta3 =~ y33      0.7474      0.0386   19.3594    0.0000 [ 0.6487; 0.8483 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0194   20.4264    0.0000 [ 0.3480; 0.4481 ] 
#>   eta1 <~ y12      0.3873      0.0209   18.5766    0.0000 [ 0.3351; 0.4429 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8399    0.0000 [ 0.3976; 0.5103 ] 
#>   eta2 <~ y21      0.3058      0.0273   11.2070    0.0000 [ 0.2388; 0.3799 ] 
#>   eta2 <~ y22      0.4473      0.0254   17.5951    0.0000 [ 0.3804; 0.5119 ] 
#>   eta2 <~ y23      0.4735      0.0169   28.0380    0.0000 [ 0.4328; 0.5201 ] 
#>   eta3 <~ y31      0.4400      0.0189   23.3326    0.0000 [ 0.3874; 0.4849 ] 
#>   eta3 <~ y32      0.3521      0.0202   17.4421    0.0000 [ 0.3040; 0.4084 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.3833    0.0000 [ 0.3513; 0.4528 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0467   14.3837    0.0000 [ 0.5504; 0.7918 ] 
#>   eta3 ~ eta1       0.6634      0.0330   20.1128    0.0000 [ 0.5739; 0.7444 ] 
#>   eta3 ~ eta2       0.3052      0.0779    3.9183    0.0001 [ 0.1089; 0.5117 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0545    3.7610    0.0002 [ 0.0675; 0.3492 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04667328 14.383676 6.551952e-47
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07069415  6.485780 8.827362e-11
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07787941  3.918252 8.919360e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5504083          0.7917762          0.5793918          0.7627926
#> 2          0.2679818          0.6335720          0.3118820          0.5896718
#> 3          0.1089445          0.5116929          0.1573067          0.4633307
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5831004          0.7636568          0.5883877          0.7582777
#> 2          0.3246617          0.5963867          0.3459152          0.5784272
#> 3          0.1388299          0.4515078          0.1772349          0.4363900
```
