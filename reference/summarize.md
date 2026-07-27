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
#>  Random seed                        = -1233018149
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
#>   eta2 ~ eta1      0.6713      0.0487   13.7959    0.0000 [ 0.5827; 0.7449 ] 
#>   eta3 ~ eta1      0.4585      0.0735    6.2358    0.0000 [ 0.3476; 0.6057 ] 
#>   eta3 ~ eta2      0.3052      0.0920    3.3159    0.0009 [ 0.0655; 0.4161 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0389   17.0372    0.0000 [ 0.5675; 0.7002 ] 
#>   eta1 =~ y12      0.6493      0.0431   15.0669    0.0000 [ 0.5563; 0.7071 ] 
#>   eta1 =~ y13      0.7613      0.0343   22.1673    0.0000 [ 0.6749; 0.8080 ] 
#>   eta2 =~ y21      0.5165      0.0449   11.5027    0.0000 [ 0.4060; 0.5769 ] 
#>   eta2 =~ y22      0.7554      0.0326   23.1917    0.0000 [ 0.6999; 0.8069 ] 
#>   eta2 =~ y23      0.7997      0.0387   20.6829    0.0000 [ 0.7423; 0.8919 ] 
#>   eta3 =~ y31      0.8223      0.0358   22.9426    0.0000 [ 0.7538; 0.8858 ] 
#>   eta3 =~ y32      0.6581      0.0505   13.0371    0.0000 [ 0.5652; 0.7291 ] 
#>   eta3 =~ y33      0.7474      0.0478   15.6415    0.0000 [ 0.6491; 0.8148 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0241   16.3959    0.0000 [ 0.3438; 0.4325 ] 
#>   eta1 <~ y12      0.3873      0.0215   17.9752    0.0000 [ 0.3467; 0.4284 ] 
#>   eta1 <~ y13      0.4542      0.0200   22.6715    0.0000 [ 0.4245; 0.4941 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3956    0.0000 [ 0.2494; 0.3374 ] 
#>   eta2 <~ y22      0.4473      0.0185   24.2016    0.0000 [ 0.4129; 0.4809 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3526    0.0000 [ 0.4397; 0.5175 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4482    0.0000 [ 0.4029; 0.4841 ] 
#>   eta3 <~ y32      0.3521      0.0240   14.6612    0.0000 [ 0.3112; 0.3841 ] 
#>   eta3 <~ y33      0.3999      0.0232   17.2139    0.0000 [ 0.3542; 0.4343 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0487   13.7959    0.0000 [ 0.5827; 0.7449 ] 
#>   eta3 ~ eta1       0.6634      0.0339   19.5440    0.0000 [ 0.6131; 0.7243 ] 
#>   eta3 ~ eta2       0.3052      0.0920    3.3159    0.0009 [ 0.0655; 0.4161 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4310    0.0006 [ 0.0480; 0.2797 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03891896 17.03720  4.351019e-65
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04309292 15.06693  2.672575e-51
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03434551 22.16726 7.111307e-109
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04489864 11.50268  1.278791e-30
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03257149 23.19168 5.523828e-119
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03866302 20.68291  4.937483e-95
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03584058 22.94264 1.745317e-116
#> 8 eta3 =~ y32  Common factor 0.6580689 0.05047668 13.03709  7.527679e-39
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04778470 15.64149  3.797158e-55
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5674786          0.7002047
#> 2          0.5563433          0.7070976
#> 3          0.6748896          0.8079733
#> 4          0.4060232          0.5768525
#> 5          0.6999160          0.8069114
#> 6          0.7422595          0.8919155
#> 7          0.7538170          0.8857999
#> 8          0.5651826          0.7290633
#> 9          0.6490797          0.8147903

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
#>  Random seed                        = -1233018149
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
#>   eta2 ~ eta1      0.6713      0.0487   13.7959    0.0000 [ 0.5540; 0.8057 ] 
#>   eta3 ~ eta1      0.4585      0.0735    6.2358    0.0000 [ 0.2608; 0.6411 ] 
#>   eta3 ~ eta2      0.3052      0.0920    3.3159    0.0009 [ 0.0731; 0.5491 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0389   17.0372    0.0000 [ 0.5786; 0.7798 ] 
#>   eta1 =~ y12      0.6493      0.0431   15.0669    0.0000 [ 0.5349; 0.7578 ] 
#>   eta1 =~ y13      0.7613      0.0343   22.1673    0.0000 [ 0.6726; 0.8503 ] 
#>   eta2 =~ y21      0.5165      0.0449   11.5027    0.0000 [ 0.4035; 0.6357 ] 
#>   eta2 =~ y22      0.7554      0.0326   23.1917    0.0000 [ 0.6707; 0.8392 ] 
#>   eta2 =~ y23      0.7997      0.0387   20.6829    0.0000 [ 0.6938; 0.8938 ] 
#>   eta3 =~ y31      0.8223      0.0358   22.9426    0.0000 [ 0.7288; 0.9142 ] 
#>   eta3 =~ y32      0.6581      0.0505   13.0371    0.0000 [ 0.5235; 0.7846 ] 
#>   eta3 =~ y33      0.7474      0.0478   15.6415    0.0000 [ 0.6237; 0.8708 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0241   16.3959    0.0000 [ 0.3393; 0.4640 ] 
#>   eta1 <~ y12      0.3873      0.0215   17.9752    0.0000 [ 0.3267; 0.4381 ] 
#>   eta1 <~ y13      0.4542      0.0200   22.6715    0.0000 [ 0.3984; 0.5020 ] 
#>   eta2 <~ y21      0.3058      0.0247   12.3956    0.0000 [ 0.2451; 0.3727 ] 
#>   eta2 <~ y22      0.4473      0.0185   24.2016    0.0000 [ 0.4008; 0.4964 ] 
#>   eta2 <~ y23      0.4735      0.0222   21.3526    0.0000 [ 0.4143; 0.5290 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4482    0.0000 [ 0.3827; 0.4997 ] 
#>   eta3 <~ y32      0.3521      0.0240   14.6612    0.0000 [ 0.2895; 0.4137 ] 
#>   eta3 <~ y33      0.3999      0.0232   17.2139    0.0000 [ 0.3415; 0.4616 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0487   13.7959    0.0000 [ 0.5540; 0.8057 ] 
#>   eta3 ~ eta1       0.6634      0.0339   19.5440    0.0000 [ 0.5764; 0.7520 ] 
#>   eta3 ~ eta2       0.3052      0.0920    3.3159    0.0009 [ 0.0731; 0.5491 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0597    3.4310    0.0006 [ 0.0589; 0.3677 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04866169 13.795935 2.696355e-43
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07352799  6.235813 4.494385e-10
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09202619  3.315916 9.134325e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5540195          0.8056702          0.5842378          0.7754519
#> 2          0.2608060          0.6410512          0.3064660          0.5953912
#> 3          0.0731481          0.5490556          0.1302953          0.4919084
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.56296595          0.7653665         0.58274782          0.7448786
#> 2         0.34606397          0.6373089         0.34764784          0.6056818
#> 3         0.05951257          0.4284571         0.06552304          0.4161390
```
