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
#>  Random seed                        = -641162702
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
#>   eta2 ~ eta1      0.6713      0.0407   16.4828    0.0000 [ 0.5996; 0.7425 ] 
#>   eta3 ~ eta1      0.4585      0.0976    4.7001    0.0000 [ 0.3348; 0.6737 ] 
#>   eta3 ~ eta2      0.3052      0.0946    3.2265    0.0013 [ 0.1352; 0.4491 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0334   19.8567    0.0000 [ 0.6057; 0.7333 ] 
#>   eta1 =~ y12      0.6493      0.0382   16.9846    0.0000 [ 0.5778; 0.7183 ] 
#>   eta1 =~ y13      0.7613      0.0370   20.5857    0.0000 [ 0.7114; 0.8300 ] 
#>   eta2 =~ y21      0.5165      0.0486   10.6236    0.0000 [ 0.4115; 0.5703 ] 
#>   eta2 =~ y22      0.7554      0.0257   29.3807    0.0000 [ 0.7132; 0.7986 ] 
#>   eta2 =~ y23      0.7997      0.0391   20.4714    0.0000 [ 0.7294; 0.8615 ] 
#>   eta3 =~ y31      0.8223      0.0284   28.9392    0.0000 [ 0.7694; 0.8760 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.6977    0.0000 [ 0.5947; 0.7456 ] 
#>   eta3 =~ y33      0.7474      0.0404   18.5143    0.0000 [ 0.6644; 0.7927 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0215   18.4341    0.0000 [ 0.3521; 0.4258 ] 
#>   eta1 <~ y12      0.3873      0.0184   21.0571    0.0000 [ 0.3488; 0.4093 ] 
#>   eta1 <~ y13      0.4542      0.0169   26.8647    0.0000 [ 0.4234; 0.4806 ] 
#>   eta2 <~ y21      0.3058      0.0230   13.2948    0.0000 [ 0.2586; 0.3419 ] 
#>   eta2 <~ y22      0.4473      0.0196   22.8151    0.0000 [ 0.4194; 0.4846 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3671    0.0000 [ 0.4488; 0.5171 ] 
#>   eta3 <~ y31      0.4400      0.0175   25.1568    0.0000 [ 0.3988; 0.4710 ] 
#>   eta3 <~ y32      0.3521      0.0203   17.3493    0.0000 [ 0.3170; 0.3913 ] 
#>   eta3 <~ y33      0.3999      0.0181   22.0758    0.0000 [ 0.3697; 0.4298 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0407   16.4828    0.0000 [ 0.5996; 0.7425 ] 
#>   eta3 ~ eta1       0.6634      0.0449   14.7634    0.0000 [ 0.5939; 0.7593 ] 
#>   eta3 ~ eta2       0.3052      0.0946    3.2265    0.0013 [ 0.1352; 0.4491 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0642    3.1894    0.0014 [ 0.0923; 0.3071 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03339281 19.85667  9.650758e-88
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03822740 16.98462  1.067426e-64
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03698422 20.58569  3.687328e-94
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04861392 10.62360  2.314614e-26
#> 5 eta2 =~ y22  Common factor 0.7553877 0.02571037 29.38067 9.699688e-190
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03906245 20.47142  3.872064e-93
#> 7 eta3 =~ y31  Common factor 0.8222773 0.02841399 28.93917 3.840660e-184
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04192138 15.69769  1.568469e-55
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04037001 18.51434  1.582201e-76
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6057014          0.7332905
#> 2          0.5778401          0.7183043
#> 3          0.7114112          0.8299595
#> 4          0.4114507          0.5702707
#> 5          0.7132180          0.7986416
#> 6          0.7293594          0.8614760
#> 7          0.7693730          0.8759891
#> 8          0.5946736          0.7456049
#> 9          0.6644386          0.7926635

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
#>  Random seed                        = -641162702
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
#>   eta2 ~ eta1      0.6713      0.0407   16.4828    0.0000 [ 0.5635; 0.7742 ] 
#>   eta3 ~ eta1      0.4585      0.0976    4.7001    0.0000 [ 0.1907; 0.6952 ] 
#>   eta3 ~ eta2      0.3052      0.0946    3.2265    0.0013 [ 0.0715; 0.5606 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0334   19.8567    0.0000 [ 0.5682; 0.7409 ] 
#>   eta1 =~ y12      0.6493      0.0382   16.9846    0.0000 [ 0.5547; 0.7524 ] 
#>   eta1 =~ y13      0.7613      0.0370   20.5857    0.0000 [ 0.6660; 0.8572 ] 
#>   eta2 =~ y21      0.5165      0.0486   10.6236    0.0000 [ 0.4077; 0.6591 ] 
#>   eta2 =~ y22      0.7554      0.0257   29.3807    0.0000 [ 0.6895; 0.8225 ] 
#>   eta2 =~ y23      0.7997      0.0391   20.4714    0.0000 [ 0.6989; 0.9009 ] 
#>   eta3 =~ y31      0.8223      0.0284   28.9392    0.0000 [ 0.7523; 0.8993 ] 
#>   eta3 =~ y32      0.6581      0.0419   15.6977    0.0000 [ 0.5415; 0.7583 ] 
#>   eta3 =~ y33      0.7474      0.0404   18.5143    0.0000 [ 0.6540; 0.8628 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0215   18.4341    0.0000 [ 0.3361; 0.4471 ] 
#>   eta1 <~ y12      0.3873      0.0184   21.0571    0.0000 [ 0.3437; 0.4388 ] 
#>   eta1 <~ y13      0.4542      0.0169   26.8647    0.0000 [ 0.4122; 0.4996 ] 
#>   eta2 <~ y21      0.3058      0.0230   13.2948    0.0000 [ 0.2537; 0.3727 ] 
#>   eta2 <~ y22      0.4473      0.0196   22.8151    0.0000 [ 0.3920; 0.4934 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3671    0.0000 [ 0.4114; 0.5260 ] 
#>   eta3 <~ y31      0.4400      0.0175   25.1568    0.0000 [ 0.3946; 0.4851 ] 
#>   eta3 <~ y32      0.3521      0.0203   17.3493    0.0000 [ 0.2938; 0.3988 ] 
#>   eta3 <~ y33      0.3999      0.0181   22.0758    0.0000 [ 0.3574; 0.4511 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0407   16.4828    0.0000 [ 0.5635; 0.7742 ] 
#>   eta3 ~ eta1       0.6634      0.0449   14.7634    0.0000 [ 0.5388; 0.7711 ] 
#>   eta3 ~ eta2       0.3052      0.0946    3.2265    0.0013 [ 0.0715; 0.5606 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0642    3.1894    0.0014 [ 0.0460; 0.3781 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04072938 16.482781 4.878814e-61
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09755179  4.700137 2.599873e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09457724  3.226475 1.253251e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5635332          0.7741625          0.5888257          0.7488701
#> 2          0.1906683          0.6951511          0.2512468          0.6345726
#> 3          0.0714863          0.5605864          0.1302177          0.5018550
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5691840          0.7425456          0.5995606          0.7425085
#> 2          0.2931203          0.6828063          0.3348161          0.6737250
#> 3          0.1047122          0.4981155          0.1352161          0.4491348
```
