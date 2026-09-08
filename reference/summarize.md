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
#>  Random seed                        = 597647898
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
#>   eta2 ~ eta1      0.6713      0.0456   14.7291    0.0000 [ 0.5983; 0.7500 ] 
#>   eta3 ~ eta1      0.4585      0.0746    6.1492    0.0000 [ 0.2864; 0.5716 ] 
#>   eta3 ~ eta2      0.3052      0.0690    4.4240    0.0000 [ 0.1981; 0.4545 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0351   18.8783    0.0000 [ 0.5906; 0.7138 ] 
#>   eta1 =~ y12      0.6493      0.0443   14.6586    0.0000 [ 0.5694; 0.7337 ] 
#>   eta1 =~ y13      0.7613      0.0304   25.0113    0.0000 [ 0.6934; 0.8158 ] 
#>   eta2 =~ y21      0.5165      0.0479   10.7743    0.0000 [ 0.4350; 0.5956 ] 
#>   eta2 =~ y22      0.7554      0.0350   21.5755    0.0000 [ 0.6839; 0.8098 ] 
#>   eta2 =~ y23      0.7997      0.0352   22.6917    0.0000 [ 0.7322; 0.8496 ] 
#>   eta3 =~ y31      0.8223      0.0328   25.0759    0.0000 [ 0.7653; 0.8736 ] 
#>   eta3 =~ y32      0.6581      0.0424   15.5326    0.0000 [ 0.5596; 0.7261 ] 
#>   eta3 =~ y33      0.7474      0.0382   19.5875    0.0000 [ 0.6720; 0.8117 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0189   20.9668    0.0000 [ 0.3609; 0.4237 ] 
#>   eta1 <~ y12      0.3873      0.0237   16.3106    0.0000 [ 0.3418; 0.4325 ] 
#>   eta1 <~ y13      0.4542      0.0195   23.2796    0.0000 [ 0.4238; 0.4857 ] 
#>   eta2 <~ y21      0.3058      0.0275   11.1156    0.0000 [ 0.2561; 0.3505 ] 
#>   eta2 <~ y22      0.4473      0.0167   26.7567    0.0000 [ 0.4180; 0.4795 ] 
#>   eta2 <~ y23      0.4735      0.0199   23.7623    0.0000 [ 0.4451; 0.5189 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.9097    0.0000 [ 0.4089; 0.4789 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.6040    0.0000 [ 0.3137; 0.3859 ] 
#>   eta3 <~ y33      0.3999      0.0217   18.4392    0.0000 [ 0.3561; 0.4416 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0456   14.7291    0.0000 [ 0.5983; 0.7500 ] 
#>   eta3 ~ eta1       0.6634      0.0400   16.5635    0.0000 [ 0.5889; 0.7242 ] 
#>   eta3 ~ eta2       0.3052      0.0690    4.4240    0.0000 [ 0.1981; 0.4545 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0465    4.4035    0.0000 [ 0.1317; 0.3043 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03512343 18.87828  1.721096e-79
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04429318 14.65864  1.186205e-48
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03044003 25.01134 4.602141e-138
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04793385 10.77432  4.551115e-27
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03501141 21.57547 3.053680e-103
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03524030 22.69174 5.405412e-114
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03279153 25.07591 9.110279e-139
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04236708 15.53255  2.088930e-54
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03815825 19.58749  1.977147e-85
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5905546          0.7137729
#> 2          0.5693709          0.7337194
#> 3          0.6934382          0.8158048
#> 4          0.4350408          0.5956191
#> 5          0.6839447          0.8098455
#> 6          0.7322217          0.8495508
#> 7          0.7653145          0.8736044
#> 8          0.5595574          0.7260866
#> 9          0.6719869          0.8117305

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
#>  Random seed                        = 597647898
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
#>   eta2 ~ eta1      0.6713      0.0456   14.7291    0.0000 [ 0.5460; 0.7817 ] 
#>   eta3 ~ eta1      0.4585      0.0746    6.1492    0.0000 [ 0.2744; 0.6600 ] 
#>   eta3 ~ eta2      0.3052      0.0690    4.4240    0.0000 [ 0.1120; 0.4687 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0351   18.8783    0.0000 [ 0.5795; 0.7612 ] 
#>   eta1 =~ y12      0.6493      0.0443   14.6586    0.0000 [ 0.5343; 0.7634 ] 
#>   eta1 =~ y13      0.7613      0.0304   25.0113    0.0000 [ 0.6836; 0.8410 ] 
#>   eta2 =~ y21      0.5165      0.0479   10.7743    0.0000 [ 0.4016; 0.6495 ] 
#>   eta2 =~ y22      0.7554      0.0350   21.5755    0.0000 [ 0.6667; 0.8477 ] 
#>   eta2 =~ y23      0.7997      0.0352   22.6917    0.0000 [ 0.7081; 0.8904 ] 
#>   eta3 =~ y31      0.8223      0.0328   25.0759    0.0000 [ 0.7365; 0.9061 ] 
#>   eta3 =~ y32      0.6581      0.0424   15.5326    0.0000 [ 0.5540; 0.7731 ] 
#>   eta3 =~ y33      0.7474      0.0382   19.5875    0.0000 [ 0.6608; 0.8581 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0189   20.9668    0.0000 [ 0.3492; 0.4468 ] 
#>   eta1 <~ y12      0.3873      0.0237   16.3106    0.0000 [ 0.3239; 0.4467 ] 
#>   eta1 <~ y13      0.4542      0.0195   23.2796    0.0000 [ 0.4019; 0.5028 ] 
#>   eta2 <~ y21      0.3058      0.0275   11.1156    0.0000 [ 0.2383; 0.3806 ] 
#>   eta2 <~ y22      0.4473      0.0167   26.7567    0.0000 [ 0.4026; 0.4891 ] 
#>   eta2 <~ y23      0.4735      0.0199   23.7623    0.0000 [ 0.4189; 0.5220 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.9097    0.0000 [ 0.3890; 0.4804 ] 
#>   eta3 <~ y32      0.3521      0.0189   18.6040    0.0000 [ 0.3026; 0.4005 ] 
#>   eta3 <~ y33      0.3999      0.0217   18.4392    0.0000 [ 0.3459; 0.4581 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0456   14.7291    0.0000 [ 0.5460; 0.7817 ] 
#>   eta3 ~ eta1       0.6634      0.0400   16.5635    0.0000 [ 0.5568; 0.7639 ] 
#>   eta3 ~ eta2       0.3052      0.0690    4.4240    0.0000 [ 0.1120; 0.4687 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0465    4.4035    0.0000 [ 0.0729; 0.3135 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04557878 14.729080 4.193771e-49
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07456420  6.149155 7.789699e-10
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.06897661  4.423980 9.689899e-06
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5459533          0.7816610          0.5742572          0.7533571
#> 2          0.2743550          0.6599589          0.3206585          0.6136554
#> 3          0.1119619          0.4686699          0.1547955          0.4258363
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5550999          0.7582827          0.5983286          0.7500304
#> 2          0.2861326          0.5945839          0.2863744          0.5715740
#> 3          0.1824729          0.4907267          0.1980541          0.4545318
```
