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
#>  Random seed                        = 607542032
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
#>   eta2 ~ eta1      0.6713      0.0496   13.5237    0.0000 [ 0.5955; 0.7545 ] 
#>   eta3 ~ eta1      0.4585      0.0820    5.5938    0.0000 [ 0.3046; 0.6071 ] 
#>   eta3 ~ eta2      0.3052      0.0717    4.2550    0.0000 [ 0.1601; 0.4372 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0438   15.1363    0.0000 [ 0.5961; 0.7597 ] 
#>   eta1 =~ y12      0.6493      0.0368   17.6530    0.0000 [ 0.5971; 0.7094 ] 
#>   eta1 =~ y13      0.7613      0.0330   23.0651    0.0000 [ 0.7124; 0.8220 ] 
#>   eta2 =~ y21      0.5165      0.0503   10.2634    0.0000 [ 0.4578; 0.6118 ] 
#>   eta2 =~ y22      0.7554      0.0371   20.3833    0.0000 [ 0.7012; 0.8224 ] 
#>   eta2 =~ y23      0.7997      0.0313   25.5138    0.0000 [ 0.7350; 0.8386 ] 
#>   eta3 =~ y31      0.8223      0.0306   26.9078    0.0000 [ 0.7763; 0.8833 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.1845    0.0000 [ 0.5798; 0.7187 ] 
#>   eta3 =~ y33      0.7474      0.0342   21.8843    0.0000 [ 0.6981; 0.8156 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0232   17.0853    0.0000 [ 0.3526; 0.4427 ] 
#>   eta1 <~ y12      0.3873      0.0204   18.9461    0.0000 [ 0.3532; 0.4215 ] 
#>   eta1 <~ y13      0.4542      0.0193   23.5148    0.0000 [ 0.4244; 0.4923 ] 
#>   eta2 <~ y21      0.3058      0.0269   11.3813    0.0000 [ 0.2727; 0.3605 ] 
#>   eta2 <~ y22      0.4473      0.0179   24.9200    0.0000 [ 0.4110; 0.4688 ] 
#>   eta2 <~ y23      0.4735      0.0188   25.1722    0.0000 [ 0.4260; 0.5004 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1331    0.0000 [ 0.4097; 0.4702 ] 
#>   eta3 <~ y32      0.3521      0.0185   19.0505    0.0000 [ 0.3182; 0.3763 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.9032    0.0000 [ 0.3632; 0.4358 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0496   13.5237    0.0000 [ 0.5955; 0.7545 ] 
#>   eta3 ~ eta1       0.6634      0.0457   14.5067    0.0000 [ 0.5850; 0.7444 ] 
#>   eta3 ~ eta2       0.3052      0.0717    4.2550    0.0000 [ 0.1601; 0.4372 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0492    4.1667    0.0000 [ 0.1202; 0.3042 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04380647 15.13635  9.325441e-52
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03678002 17.65301  9.648257e-70
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03300861 23.06507 1.038574e-117
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05032011 10.26339  1.030275e-24
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03705914 20.38330  2.352454e-92
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03134240 25.51380 1.385673e-143
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03055905 26.90782 1.779042e-159
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04333825 15.18448  4.480932e-52
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03415347 21.88428 3.667650e-106
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5961121          0.7597281
#> 2          0.5970787          0.7094142
#> 3          0.7124406          0.8219589
#> 4          0.4577553          0.6117582
#> 5          0.7012033          0.8223616
#> 6          0.7350178          0.8386477
#> 7          0.7762560          0.8833458
#> 8          0.5797953          0.7187237
#> 9          0.6980563          0.8156166

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
#>  Random seed                        = 607542032
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
#>   eta2 ~ eta1      0.6713      0.0496   13.5237    0.0000 [ 0.5395; 0.7962 ] 
#>   eta3 ~ eta1      0.4585      0.0820    5.5938    0.0000 [ 0.2715; 0.6954 ] 
#>   eta3 ~ eta2      0.3052      0.0717    4.2550    0.0000 [ 0.0950; 0.4659 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0438   15.1363    0.0000 [ 0.5535; 0.7801 ] 
#>   eta1 =~ y12      0.6493      0.0368   17.6530    0.0000 [ 0.5492; 0.7394 ] 
#>   eta1 =~ y13      0.7613      0.0330   23.0651    0.0000 [ 0.6709; 0.8416 ] 
#>   eta2 =~ y21      0.5165      0.0503   10.2634    0.0000 [ 0.3658; 0.6260 ] 
#>   eta2 =~ y22      0.7554      0.0371   20.3833    0.0000 [ 0.6621; 0.8537 ] 
#>   eta2 =~ y23      0.7997      0.0313   25.5138    0.0000 [ 0.7294; 0.8915 ] 
#>   eta3 =~ y31      0.8223      0.0306   26.9078    0.0000 [ 0.7358; 0.8939 ] 
#>   eta3 =~ y32      0.6581      0.0433   15.1845    0.0000 [ 0.5564; 0.7805 ] 
#>   eta3 =~ y33      0.7474      0.0342   21.8843    0.0000 [ 0.6533; 0.8299 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0232   17.0853    0.0000 [ 0.3401; 0.4598 ] 
#>   eta1 <~ y12      0.3873      0.0204   18.9461    0.0000 [ 0.3335; 0.4392 ] 
#>   eta1 <~ y13      0.4542      0.0193   23.5148    0.0000 [ 0.4035; 0.5034 ] 
#>   eta2 <~ y21      0.3058      0.0269   11.3813    0.0000 [ 0.2255; 0.3645 ] 
#>   eta2 <~ y22      0.4473      0.0179   24.9200    0.0000 [ 0.4041; 0.4969 ] 
#>   eta2 <~ y23      0.4735      0.0188   25.1722    0.0000 [ 0.4329; 0.5302 ] 
#>   eta3 <~ y31      0.4400      0.0168   26.1331    0.0000 [ 0.3937; 0.4808 ] 
#>   eta3 <~ y32      0.3521      0.0185   19.0505    0.0000 [ 0.3111; 0.4067 ] 
#>   eta3 <~ y33      0.3999      0.0201   19.9032    0.0000 [ 0.3459; 0.4498 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0496   13.5237    0.0000 [ 0.5395; 0.7962 ] 
#>   eta3 ~ eta1       0.6634      0.0457   14.5067    0.0000 [ 0.5528; 0.7893 ] 
#>   eta3 ~ eta2       0.3052      0.0717    4.2550    0.0000 [ 0.0950; 0.4659 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0492    4.1667    0.0000 [ 0.0604; 0.3147 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04964117 13.523724 1.132929e-41
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08196761  5.593756 2.222096e-08
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07171588  4.255001 2.090482e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.53946503          0.7961811          0.5702916          0.7653545
#> 2         0.27154034          0.6954305          0.3224413          0.6445296
#> 3         0.09504043          0.4659145          0.1395752          0.4213797
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5933947          0.7697636          0.5954597          0.7545019
#> 2          0.3025274          0.6295894          0.3046435          0.6070843
#> 3          0.1378835          0.4528036          0.1601449          0.4372080
```
