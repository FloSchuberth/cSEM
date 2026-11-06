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
#>  Random seed                        = -673478419
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
#>   eta2 ~ eta1      0.6713      0.0439   15.2810    0.0000 [ 0.6059; 0.7439 ] 
#>   eta3 ~ eta1      0.4585      0.0628    7.2985    0.0000 [ 0.3297; 0.5679 ] 
#>   eta3 ~ eta2      0.3052      0.0714    4.2732    0.0000 [ 0.1657; 0.4525 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0415   15.9905    0.0000 [ 0.5989; 0.7316 ] 
#>   eta1 =~ y12      0.6493      0.0316   20.5306    0.0000 [ 0.6063; 0.7172 ] 
#>   eta1 =~ y13      0.7613      0.0325   23.4046    0.0000 [ 0.6863; 0.8124 ] 
#>   eta2 =~ y21      0.5165      0.0499   10.3589    0.0000 [ 0.4255; 0.6050 ] 
#>   eta2 =~ y22      0.7554      0.0301   25.0741    0.0000 [ 0.6934; 0.7999 ] 
#>   eta2 =~ y23      0.7997      0.0320   25.0145    0.0000 [ 0.7263; 0.8575 ] 
#>   eta3 =~ y31      0.8223      0.0374   22.0059    0.0000 [ 0.7656; 0.9085 ] 
#>   eta3 =~ y32      0.6581      0.0401   16.4010    0.0000 [ 0.5692; 0.7046 ] 
#>   eta3 =~ y33      0.7474      0.0416   17.9714    0.0000 [ 0.6665; 0.8186 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0232   17.0536    0.0000 [ 0.3598; 0.4441 ] 
#>   eta1 <~ y12      0.3873      0.0171   22.6000    0.0000 [ 0.3651; 0.4187 ] 
#>   eta1 <~ y13      0.4542      0.0194   23.4439    0.0000 [ 0.4160; 0.4825 ] 
#>   eta2 <~ y21      0.3058      0.0268   11.4250    0.0000 [ 0.2569; 0.3633 ] 
#>   eta2 <~ y22      0.4473      0.0154   29.0270    0.0000 [ 0.4252; 0.4780 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1879    0.0000 [ 0.4378; 0.5128 ] 
#>   eta3 <~ y31      0.4400      0.0198   22.2395    0.0000 [ 0.4105; 0.4894 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9785    0.0000 [ 0.3086; 0.3790 ] 
#>   eta3 <~ y33      0.3999      0.0204   19.5722    0.0000 [ 0.3614; 0.4404 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0439   15.2810    0.0000 [ 0.6059; 0.7439 ] 
#>   eta3 ~ eta1       0.6634      0.0369   17.9857    0.0000 [ 0.5925; 0.7196 ] 
#>   eta3 ~ eta2       0.3052      0.0714    4.2732    0.0000 [ 0.1657; 0.4525 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0441    4.6442    0.0000 [ 0.1262; 0.2784 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04146658 15.99047  1.489132e-57
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03162483 20.53064  1.146631e-93
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03252970 23.40463 3.833625e-121
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04985624 10.35888  3.814021e-25
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03012627 25.07405 9.545228e-139
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03196796 25.01454 4.247637e-138
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03736629 22.00586 2.530740e-107
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04012377 16.40097  1.881983e-60
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04158969 17.97138  3.264844e-72
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5989068          0.7315825
#> 2          0.6063191          0.7171970
#> 3          0.6862524          0.8124219
#> 4          0.4254816          0.6050362
#> 5          0.6933817          0.7999481
#> 6          0.7262527          0.8575049
#> 7          0.7655936          0.9084903
#> 8          0.5692165          0.7046324
#> 9          0.6664994          0.8186440

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
#>  Random seed                        = -673478419
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
#>   eta2 ~ eta1      0.6713      0.0439   15.2810    0.0000 [ 0.5578; 0.7850 ] 
#>   eta3 ~ eta1      0.4585      0.0628    7.2985    0.0000 [ 0.3167; 0.6416 ] 
#>   eta3 ~ eta2      0.3052      0.0714    4.2732    0.0000 [ 0.1071; 0.4763 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0415   15.9905    0.0000 [ 0.5608; 0.7753 ] 
#>   eta1 =~ y12      0.6493      0.0316   20.5306    0.0000 [ 0.5624; 0.7259 ] 
#>   eta1 =~ y13      0.7613      0.0325   23.4046    0.0000 [ 0.6857; 0.8540 ] 
#>   eta2 =~ y21      0.5165      0.0499   10.3589    0.0000 [ 0.3873; 0.6452 ] 
#>   eta2 =~ y22      0.7554      0.0301   25.0741    0.0000 [ 0.6760; 0.8318 ] 
#>   eta2 =~ y23      0.7997      0.0320   25.0145    0.0000 [ 0.7143; 0.8796 ] 
#>   eta3 =~ y31      0.8223      0.0374   22.0059    0.0000 [ 0.7247; 0.9179 ] 
#>   eta3 =~ y32      0.6581      0.0401   16.4010    0.0000 [ 0.5615; 0.7690 ] 
#>   eta3 =~ y33      0.7474      0.0416   17.9714    0.0000 [ 0.6434; 0.8585 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0232   17.0536    0.0000 [ 0.3364; 0.4563 ] 
#>   eta1 <~ y12      0.3873      0.0171   22.6000    0.0000 [ 0.3378; 0.4264 ] 
#>   eta1 <~ y13      0.4542      0.0194   23.4439    0.0000 [ 0.4066; 0.5067 ] 
#>   eta2 <~ y21      0.3058      0.0268   11.4250    0.0000 [ 0.2379; 0.3764 ] 
#>   eta2 <~ y22      0.4473      0.0154   29.0270    0.0000 [ 0.4084; 0.4881 ] 
#>   eta2 <~ y23      0.4735      0.0204   23.1879    0.0000 [ 0.4209; 0.5265 ] 
#>   eta3 <~ y31      0.4400      0.0198   22.2395    0.0000 [ 0.3856; 0.4879 ] 
#>   eta3 <~ y32      0.3521      0.0196   17.9785    0.0000 [ 0.3033; 0.4046 ] 
#>   eta3 <~ y33      0.3999      0.0204   19.5722    0.0000 [ 0.3466; 0.4523 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0439   15.2810    0.0000 [ 0.5578; 0.7850 ] 
#>   eta3 ~ eta1       0.6634      0.0369   17.9857    0.0000 [ 0.5810; 0.7717 ] 
#>   eta3 ~ eta2       0.3052      0.0714    4.2732    0.0000 [ 0.1071; 0.4763 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0441    4.6442    0.0000 [ 0.0831; 0.3113 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04393257 15.280995 1.023653e-52
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06282190  7.298518 2.909548e-13
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07141037  4.273205 1.926836e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5578059          0.7850004          0.5850876          0.7577188
#> 2          0.3167148          0.6415942          0.3557264          0.6025825
#> 3          0.1070514          0.4763455          0.1513964          0.4320005
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5883587          0.7591735          0.6058883          0.7438566
#> 2          0.3186324          0.5872075          0.3297196          0.5678706
#> 3          0.1529400          0.4941914          0.1656537          0.4525289
```
