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
#>  Random seed                        = 1844491508
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
#>   eta2 ~ eta1      0.6713      0.0466   14.4143    0.0000 [ 0.5738; 0.7411 ] 
#>   eta3 ~ eta1      0.4585      0.0686    6.6883    0.0000 [ 0.3417; 0.5772 ] 
#>   eta3 ~ eta2      0.3052      0.0731    4.1751    0.0000 [ 0.1730; 0.4258 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0381   17.3814    0.0000 [ 0.5977; 0.7461 ] 
#>   eta1 =~ y12      0.6493      0.0453   14.3387    0.0000 [ 0.5377; 0.7151 ] 
#>   eta1 =~ y13      0.7613      0.0303   25.1543    0.0000 [ 0.6943; 0.8062 ] 
#>   eta2 =~ y21      0.5165      0.0533    9.6957    0.0000 [ 0.4195; 0.6058 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6748    0.0000 [ 0.7009; 0.8199 ] 
#>   eta2 =~ y23      0.7997      0.0374   21.3541    0.0000 [ 0.7099; 0.8675 ] 
#>   eta3 =~ y31      0.8223      0.0406   20.2426    0.0000 [ 0.7463; 0.8728 ] 
#>   eta3 =~ y32      0.6581      0.0414   15.8782    0.0000 [ 0.5843; 0.7259 ] 
#>   eta3 =~ y33      0.7474      0.0433   17.2809    0.0000 [ 0.6639; 0.8208 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0250   15.7955    0.0000 [ 0.3681; 0.4533 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.7885    0.0000 [ 0.3343; 0.4168 ] 
#>   eta1 <~ y13      0.4542      0.0148   30.7229    0.0000 [ 0.4292; 0.4862 ] 
#>   eta2 <~ y21      0.3058      0.0278   10.9845    0.0000 [ 0.2481; 0.3495 ] 
#>   eta2 <~ y22      0.4473      0.0226   19.8037    0.0000 [ 0.4111; 0.4972 ] 
#>   eta2 <~ y23      0.4735      0.0183   25.8149    0.0000 [ 0.4343; 0.5038 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4975    0.0000 [ 0.3908; 0.4676 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5056    0.0000 [ 0.3227; 0.3882 ] 
#>   eta3 <~ y33      0.3999      0.0224   17.8851    0.0000 [ 0.3511; 0.4344 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0466   14.4143    0.0000 [ 0.5738; 0.7411 ] 
#>   eta3 ~ eta1       0.6634      0.0407   16.2867    0.0000 [ 0.5786; 0.7171 ] 
#>   eta3 ~ eta2       0.3052      0.0731    4.1751    0.0000 [ 0.1730; 0.4258 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0500    4.0942    0.0000 [ 0.1139; 0.2920 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03814822 17.381409  1.141148e-67
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04528157 14.338681  1.254165e-46
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03026707 25.154259 1.269409e-139
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05326618  9.695736  3.143615e-22
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03331390 22.674847 7.935660e-114
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03744778 21.354101 3.571744e-101
#> 7 eta3 =~ y31  Common factor 0.8222773 0.04062113 20.242600  4.128044e-91
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04144472 15.878234  8.966824e-57
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04325145 17.280904  6.551337e-67
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5977057          0.7460821
#> 2          0.5376602          0.7151399
#> 3          0.6942649          0.8061783
#> 4          0.4195180          0.6057908
#> 5          0.7009343          0.8198650
#> 6          0.7099216          0.8675393
#> 7          0.7462887          0.8727912
#> 8          0.5843400          0.7259091
#> 9          0.6638975          0.8207772

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
#>  Random seed                        = 1844491508
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
#>   eta2 ~ eta1      0.6713      0.0466   14.4143    0.0000 [ 0.5546; 0.7954 ] 
#>   eta3 ~ eta1      0.4585      0.0686    6.6883    0.0000 [ 0.2805; 0.6350 ] 
#>   eta3 ~ eta2      0.3052      0.0731    4.1751    0.0000 [ 0.1228; 0.5007 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0381   17.3814    0.0000 [ 0.5584; 0.7557 ] 
#>   eta1 =~ y12      0.6493      0.0453   14.3387    0.0000 [ 0.5405; 0.7747 ] 
#>   eta1 =~ y13      0.7613      0.0303   25.1543    0.0000 [ 0.6848; 0.8413 ] 
#>   eta2 =~ y21      0.5165      0.0533    9.6957    0.0000 [ 0.3839; 0.6594 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6748    0.0000 [ 0.6617; 0.8340 ] 
#>   eta2 =~ y23      0.7997      0.0374   21.3541    0.0000 [ 0.7029; 0.8965 ] 
#>   eta3 =~ y31      0.8223      0.0406   20.2426    0.0000 [ 0.7220; 0.9321 ] 
#>   eta3 =~ y32      0.6581      0.0414   15.8782    0.0000 [ 0.5441; 0.7584 ] 
#>   eta3 =~ y33      0.7474      0.0433   17.2809    0.0000 [ 0.6398; 0.8635 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0250   15.7955    0.0000 [ 0.3261; 0.4556 ] 
#>   eta1 <~ y12      0.3873      0.0218   17.7885    0.0000 [ 0.3354; 0.4480 ] 
#>   eta1 <~ y13      0.4542      0.0148   30.7229    0.0000 [ 0.4160; 0.4924 ] 
#>   eta2 <~ y21      0.3058      0.0278   10.9845    0.0000 [ 0.2382; 0.3822 ] 
#>   eta2 <~ y22      0.4473      0.0226   19.8037    0.0000 [ 0.3857; 0.5025 ] 
#>   eta2 <~ y23      0.4735      0.0183   25.8149    0.0000 [ 0.4278; 0.5226 ] 
#>   eta3 <~ y31      0.4400      0.0226   19.4975    0.0000 [ 0.3836; 0.5003 ] 
#>   eta3 <~ y32      0.3521      0.0190   18.5056    0.0000 [ 0.2990; 0.3974 ] 
#>   eta3 <~ y33      0.3999      0.0224   17.8851    0.0000 [ 0.3439; 0.4595 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0466   14.4143    0.0000 [ 0.5546; 0.7954 ] 
#>   eta3 ~ eta1       0.6634      0.0407   16.2867    0.0000 [ 0.5632; 0.7738 ] 
#>   eta3 ~ eta2       0.3052      0.0731    4.1751    0.0000 [ 0.1228; 0.5007 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0500    4.0942    0.0000 [ 0.0813; 0.3401 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04657408 14.414313 4.206029e-47
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06855376  6.688281 2.258077e-11
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07308856  4.175087 2.978716e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5545609          0.7954158          0.5834829          0.7664938
#> 2          0.2804963          0.6350176          0.3230674          0.5924465
#> 3          0.1227719          0.5007447          0.1681590          0.4553575
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5554239          0.7503827          0.5738252          0.7411109
#> 2          0.3032774          0.5792078          0.3417428          0.5772095
#> 3          0.1646972          0.4260797          0.1729629          0.4257645
```
