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
#>  Random seed                        = -282670885
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
#>   eta2 ~ eta1      0.6713      0.0363   18.4978    0.0000 [ 0.5962; 0.7379 ] 
#>   eta3 ~ eta1      0.4585      0.0782    5.8632    0.0000 [ 0.3347; 0.5798 ] 
#>   eta3 ~ eta2      0.3052      0.0857    3.5593    0.0004 [ 0.1308; 0.4282 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0370   17.9385    0.0000 [ 0.5873; 0.7233 ] 
#>   eta1 =~ y12      0.6493      0.0439   14.7930    0.0000 [ 0.5404; 0.7238 ] 
#>   eta1 =~ y13      0.7613      0.0269   28.2754    0.0000 [ 0.6940; 0.7996 ] 
#>   eta2 =~ y21      0.5165      0.0516   10.0088    0.0000 [ 0.4212; 0.6082 ] 
#>   eta2 =~ y22      0.7554      0.0489   15.4589    0.0000 [ 0.6708; 0.8329 ] 
#>   eta2 =~ y23      0.7997      0.0392   20.4174    0.0000 [ 0.7426; 0.8709 ] 
#>   eta3 =~ y31      0.8223      0.0383   21.4732    0.0000 [ 0.7675; 0.9180 ] 
#>   eta3 =~ y32      0.6581      0.0355   18.5536    0.0000 [ 0.6003; 0.7321 ] 
#>   eta3 =~ y33      0.7474      0.0402   18.6078    0.0000 [ 0.6635; 0.7953 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0174   22.6914    0.0000 [ 0.3662; 0.4261 ] 
#>   eta1 <~ y12      0.3873      0.0215   18.0275    0.0000 [ 0.3441; 0.4200 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8668    0.0000 [ 0.4235; 0.5031 ] 
#>   eta2 <~ y21      0.3058      0.0321    9.5416    0.0000 [ 0.2642; 0.3606 ] 
#>   eta2 <~ y22      0.4473      0.0220   20.3035    0.0000 [ 0.4062; 0.4756 ] 
#>   eta2 <~ y23      0.4735      0.0185   25.6121    0.0000 [ 0.4411; 0.5135 ] 
#>   eta3 <~ y31      0.4400      0.0218   20.1523    0.0000 [ 0.4061; 0.4880 ] 
#>   eta3 <~ y32      0.3521      0.0163   21.6058    0.0000 [ 0.3215; 0.3873 ] 
#>   eta3 <~ y33      0.3999      0.0180   22.1765    0.0000 [ 0.3564; 0.4229 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0363   18.4978    0.0000 [ 0.5962; 0.7379 ] 
#>   eta3 ~ eta1       0.6634      0.0357   18.5928    0.0000 [ 0.5972; 0.7264 ] 
#>   eta3 ~ eta2       0.3052      0.0857    3.5593    0.0004 [ 0.1308; 0.4282 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4331    0.0006 [ 0.0906; 0.2980 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03696341 17.93855  5.897578e-72
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04389090 14.79299  1.625565e-49
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02692611 28.27537 6.942403e-176
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05160028 10.00876  1.394901e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04886441 15.45885  6.576004e-54
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03916576 20.41742  1.170898e-92
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03829312 21.47324 2.770282e-102
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03546853 18.55360  7.627218e-77
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04016717 18.60784  2.776003e-77
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5872688          0.7232551
#> 2          0.5403960          0.7238132
#> 3          0.6939970          0.7995512
#> 4          0.4211772          0.6082315
#> 5          0.6708220          0.8329407
#> 6          0.7426094          0.8709396
#> 7          0.7674814          0.9179822
#> 8          0.6002868          0.7321058
#> 9          0.6635428          0.7953431

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
#>  Random seed                        = -282670885
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
#>   eta2 ~ eta1      0.6713      0.0363   18.4978    0.0000 [ 0.5839; 0.7716 ] 
#>   eta3 ~ eta1      0.4585      0.0782    5.8632    0.0000 [ 0.2475; 0.6519 ] 
#>   eta3 ~ eta2      0.3052      0.0857    3.5593    0.0004 [ 0.0912; 0.5345 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0370   17.9385    0.0000 [ 0.5648; 0.7559 ] 
#>   eta1 =~ y12      0.6493      0.0439   14.7930    0.0000 [ 0.5445; 0.7714 ] 
#>   eta1 =~ y13      0.7613      0.0269   28.2754    0.0000 [ 0.6909; 0.8302 ] 
#>   eta2 =~ y21      0.5165      0.0516   10.0088    0.0000 [ 0.3739; 0.6407 ] 
#>   eta2 =~ y22      0.7554      0.0489   15.4589    0.0000 [ 0.6380; 0.8907 ] 
#>   eta2 =~ y23      0.7997      0.0392   20.4174    0.0000 [ 0.6998; 0.9024 ] 
#>   eta3 =~ y31      0.8223      0.0383   21.4732    0.0000 [ 0.7163; 0.9143 ] 
#>   eta3 =~ y32      0.6581      0.0355   18.5536    0.0000 [ 0.5582; 0.7416 ] 
#>   eta3 =~ y33      0.7474      0.0402   18.6078    0.0000 [ 0.6533; 0.8610 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0174   22.6914    0.0000 [ 0.3478; 0.4380 ] 
#>   eta1 <~ y12      0.3873      0.0215   18.0275    0.0000 [ 0.3360; 0.4471 ] 
#>   eta1 <~ y13      0.4542      0.0218   20.8668    0.0000 [ 0.3957; 0.5083 ] 
#>   eta2 <~ y21      0.3058      0.0321    9.5416    0.0000 [ 0.2170; 0.3827 ] 
#>   eta2 <~ y22      0.4473      0.0220   20.3035    0.0000 [ 0.3954; 0.5094 ] 
#>   eta2 <~ y23      0.4735      0.0185   25.6121    0.0000 [ 0.4261; 0.5217 ] 
#>   eta3 <~ y31      0.4400      0.0218   20.1523    0.0000 [ 0.3814; 0.4943 ] 
#>   eta3 <~ y32      0.3521      0.0163   21.6058    0.0000 [ 0.3071; 0.3914 ] 
#>   eta3 <~ y33      0.3999      0.0180   22.1765    0.0000 [ 0.3601; 0.4534 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0363   18.4978    0.0000 [ 0.5839; 0.7716 ] 
#>   eta3 ~ eta1       0.6634      0.0357   18.5928    0.0000 [ 0.5690; 0.7535 ] 
#>   eta3 ~ eta2       0.3052      0.0857    3.5593    0.0004 [ 0.0912; 0.5345 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4331    0.0006 [ 0.0573; 0.3659 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03629255 18.497833 2.149454e-76
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07820016  5.863246 4.539059e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08573297  3.559321 3.718153e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.58387597          0.7715605          0.6064132          0.7490233
#> 2         0.24746958          0.6518767          0.2960310          0.6033153
#> 3         0.09116102          0.5345235          0.1444002          0.4812843
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.58392444          0.7460989          0.5961625          0.7379122
#> 2         0.32846201          0.6478088          0.3347290          0.5797768
#> 3         0.09790394          0.4305031          0.1307639          0.4281507
```
