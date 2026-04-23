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
#>  Random seed                        = 573249391
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
#>   eta2 ~ eta1      0.6713      0.0523   12.8304    0.0000 [ 0.5741; 0.7535 ] 
#>   eta3 ~ eta1      0.4585      0.0667    6.8724    0.0000 [ 0.3637; 0.5876 ] 
#>   eta3 ~ eta2      0.3052      0.0666    4.5851    0.0000 [ 0.1770; 0.4082 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0423   15.6606    0.0000 [ 0.5902; 0.7364 ] 
#>   eta1 =~ y12      0.6493      0.0351   18.4919    0.0000 [ 0.5874; 0.7101 ] 
#>   eta1 =~ y13      0.7613      0.0358   21.2537    0.0000 [ 0.6969; 0.8265 ] 
#>   eta2 =~ y21      0.5165      0.0518    9.9633    0.0000 [ 0.3981; 0.5931 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6600    0.0000 [ 0.6816; 0.8191 ] 
#>   eta2 =~ y23      0.7997      0.0346   23.1357    0.0000 [ 0.7422; 0.8615 ] 
#>   eta3 =~ y31      0.8223      0.0321   25.6004    0.0000 [ 0.7700; 0.8780 ] 
#>   eta3 =~ y32      0.6581      0.0417   15.7759    0.0000 [ 0.5957; 0.7330 ] 
#>   eta3 =~ y33      0.7474      0.0385   19.3982    0.0000 [ 0.6551; 0.8027 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0209   18.9510    0.0000 [ 0.3582; 0.4362 ] 
#>   eta1 <~ y12      0.3873      0.0185   20.9883    0.0000 [ 0.3503; 0.4181 ] 
#>   eta1 <~ y13      0.4542      0.0204   22.2229    0.0000 [ 0.4208; 0.5036 ] 
#>   eta2 <~ y21      0.3058      0.0267   11.4585    0.0000 [ 0.2419; 0.3458 ] 
#>   eta2 <~ y22      0.4473      0.0225   19.8577    0.0000 [ 0.4024; 0.4936 ] 
#>   eta2 <~ y23      0.4735      0.0174   27.1643    0.0000 [ 0.4504; 0.5083 ] 
#>   eta3 <~ y31      0.4400      0.0232   18.9700    0.0000 [ 0.4037; 0.4938 ] 
#>   eta3 <~ y32      0.3521      0.0179   19.6815    0.0000 [ 0.3248; 0.3788 ] 
#>   eta3 <~ y33      0.3999      0.0165   24.2271    0.0000 [ 0.3766; 0.4369 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0523   12.8304    0.0000 [ 0.5741; 0.7535 ] 
#>   eta3 ~ eta1       0.6634      0.0363   18.2599    0.0000 [ 0.5998; 0.7193 ] 
#>   eta3 ~ eta2       0.3052      0.0666    4.5851    0.0000 [ 0.1770; 0.4082 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0421    4.8662    0.0000 [ 0.1240; 0.2649 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04234001 15.660597  2.812426e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03511149 18.491893  2.399799e-76
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03582176 21.253724 3.045286e-100
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05183588  9.963268  2.206869e-23
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03333573 22.660002 1.111763e-113
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03456404 23.135714 2.024630e-118
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03211974 25.600377 1.510855e-144
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04171346 15.775937  4.555896e-56
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03853063 19.398179  7.995332e-84
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5902250          0.7364123
#> 2          0.5873579          0.7101325
#> 3          0.6969412          0.8264826
#> 4          0.3981068          0.5931479
#> 5          0.6816055          0.8190842
#> 6          0.7422076          0.8614642
#> 7          0.7700391          0.8780120
#> 8          0.5957282          0.7330341
#> 9          0.6551408          0.8027403

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
#>  Random seed                        = 573249391
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
#>   eta2 ~ eta1      0.6713      0.0523   12.8304    0.0000 [ 0.5386; 0.8092 ] 
#>   eta3 ~ eta1      0.4585      0.0667    6.8724    0.0000 [ 0.2818; 0.6268 ] 
#>   eta3 ~ eta2      0.3052      0.0666    4.5851    0.0000 [ 0.1341; 0.4783 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0423   15.6606    0.0000 [ 0.5474; 0.7663 ] 
#>   eta1 =~ y12      0.6493      0.0351   18.4919    0.0000 [ 0.5553; 0.7369 ] 
#>   eta1 =~ y13      0.7613      0.0358   21.2537    0.0000 [ 0.6699; 0.8552 ] 
#>   eta2 =~ y21      0.5165      0.0518    9.9633    0.0000 [ 0.3901; 0.6581 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6600    0.0000 [ 0.6716; 0.8440 ] 
#>   eta2 =~ y23      0.7997      0.0346   23.1357    0.0000 [ 0.7085; 0.8872 ] 
#>   eta3 =~ y31      0.8223      0.0321   25.6004    0.0000 [ 0.7377; 0.9038 ] 
#>   eta3 =~ y32      0.6581      0.0417   15.7759    0.0000 [ 0.5607; 0.7764 ] 
#>   eta3 =~ y33      0.7474      0.0385   19.3982    0.0000 [ 0.6562; 0.8554 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0209   18.9510    0.0000 [ 0.3405; 0.4484 ] 
#>   eta1 <~ y12      0.3873      0.0185   20.9883    0.0000 [ 0.3401; 0.4356 ] 
#>   eta1 <~ y13      0.4542      0.0204   22.2229    0.0000 [ 0.4047; 0.5104 ] 
#>   eta2 <~ y21      0.3058      0.0267   11.4585    0.0000 [ 0.2404; 0.3784 ] 
#>   eta2 <~ y22      0.4473      0.0225   19.8577    0.0000 [ 0.3883; 0.5048 ] 
#>   eta2 <~ y23      0.4735      0.0174   27.1643    0.0000 [ 0.4254; 0.5155 ] 
#>   eta3 <~ y31      0.4400      0.0232   18.9700    0.0000 [ 0.3740; 0.4939 ] 
#>   eta3 <~ y32      0.3521      0.0179   19.6815    0.0000 [ 0.3078; 0.4003 ] 
#>   eta3 <~ y33      0.3999      0.0165   24.2271    0.0000 [ 0.3575; 0.4428 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0523   12.8304    0.0000 [ 0.5386; 0.8092 ] 
#>   eta3 ~ eta1       0.6634      0.0363   18.2599    0.0000 [ 0.5677; 0.7556 ] 
#>   eta3 ~ eta2       0.3052      0.0666    4.5851    0.0000 [ 0.1341; 0.4783 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0421    4.8662    0.0000 [ 0.0985; 0.3162 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.05232382 12.830359 1.108488e-37
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06671683  6.872431 6.311709e-12
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.06655334  4.585061 4.538542e-06
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5385680          0.8091572          0.5710605          0.7766648
#> 2          0.2817782          0.6268000          0.3232086          0.5853696
#> 3          0.1340757          0.4782521          0.1754046          0.4369232
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5702585          0.7578220          0.5740672          0.7535386
#> 2          0.3600283          0.6107284          0.3637014          0.5875710
#> 3          0.1597997          0.4132562          0.1769655          0.4082285
```
