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
#>  Random seed                        = -1214398069
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
#>   eta2 ~ eta1      0.6713      0.0408   16.4419    0.0000 [ 0.6216; 0.7327 ] 
#>   eta3 ~ eta1      0.4585      0.0939    4.8818    0.0000 [ 0.3022; 0.6117 ] 
#>   eta3 ~ eta2      0.3052      0.1016    3.0032    0.0027 [ 0.0850; 0.4709 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0425   15.6021    0.0000 [ 0.5771; 0.7330 ] 
#>   eta1 =~ y12      0.6493      0.0465   13.9591    0.0000 [ 0.5635; 0.7210 ] 
#>   eta1 =~ y13      0.7613      0.0307   24.8331    0.0000 [ 0.7146; 0.8184 ] 
#>   eta2 =~ y21      0.5165      0.0423   12.1968    0.0000 [ 0.4554; 0.5910 ] 
#>   eta2 =~ y22      0.7554      0.0406   18.6108    0.0000 [ 0.6936; 0.8190 ] 
#>   eta2 =~ y23      0.7997      0.0359   22.3001    0.0000 [ 0.7154; 0.8504 ] 
#>   eta3 =~ y31      0.8223      0.0345   23.8632    0.0000 [ 0.7663; 0.8902 ] 
#>   eta3 =~ y32      0.6581      0.0371   17.7539    0.0000 [ 0.6016; 0.7425 ] 
#>   eta3 =~ y33      0.7474      0.0395   18.9275    0.0000 [ 0.6592; 0.8224 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0202   19.5910    0.0000 [ 0.3596; 0.4398 ] 
#>   eta1 <~ y12      0.3873      0.0259   14.9797    0.0000 [ 0.3287; 0.4309 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9494    0.0000 [ 0.4227; 0.4916 ] 
#>   eta2 <~ y21      0.3058      0.0252   12.1316    0.0000 [ 0.2579; 0.3626 ] 
#>   eta2 <~ y22      0.4473      0.0194   23.0506    0.0000 [ 0.4157; 0.4828 ] 
#>   eta2 <~ y23      0.4735      0.0186   25.5215    0.0000 [ 0.4367; 0.5056 ] 
#>   eta3 <~ y31      0.4400      0.0207   21.2940    0.0000 [ 0.4046; 0.4857 ] 
#>   eta3 <~ y32      0.3521      0.0185   19.0776    0.0000 [ 0.3272; 0.3836 ] 
#>   eta3 <~ y33      0.3999      0.0159   25.1457    0.0000 [ 0.3582; 0.4220 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0408   16.4419    0.0000 [ 0.6216; 0.7327 ] 
#>   eta3 ~ eta1       0.6634      0.0385   17.2093    0.0000 [ 0.5967; 0.7296 ] 
#>   eta3 ~ eta2       0.3052      0.1016    3.0032    0.0027 [ 0.0850; 0.4709 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0691    2.9649    0.0030 [ 0.0667; 0.3272 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04249864 15.60214  7.038885e-55
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04651287 13.95910  2.768946e-44
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03065854 24.83307 3.940103e-136
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04234352 12.19678  3.233438e-34
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04058864 18.61082  2.625903e-77
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03585926 22.30006 3.690126e-110
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03445792 23.86323 7.381468e-126
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03706617 17.75390  1.607997e-70
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03948869 18.92755  6.763682e-80
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5770987          0.7330232
#> 2          0.5634825          0.7210350
#> 3          0.7146327          0.8184145
#> 4          0.4553899          0.5910429
#> 5          0.6935787          0.8190499
#> 6          0.7154034          0.8503829
#> 7          0.7662660          0.8901598
#> 8          0.6016316          0.7424673
#> 9          0.6591510          0.8223993

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
#>  Random seed                        = -1214398069
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
#>   eta2 ~ eta1      0.6713      0.0408   16.4419    0.0000 [ 0.5571; 0.7682 ] 
#>   eta3 ~ eta1      0.4585      0.0939    4.8818    0.0000 [ 0.2276; 0.7134 ] 
#>   eta3 ~ eta2      0.3052      0.1016    3.0032    0.0027 [ 0.0330; 0.5584 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0425   15.6021    0.0000 [ 0.5519; 0.7717 ] 
#>   eta1 =~ y12      0.6493      0.0465   13.9591    0.0000 [ 0.5357; 0.7762 ] 
#>   eta1 =~ y13      0.7613      0.0307   24.8331    0.0000 [ 0.6748; 0.8334 ] 
#>   eta2 =~ y21      0.5165      0.0423   12.1968    0.0000 [ 0.3974; 0.6164 ] 
#>   eta2 =~ y22      0.7554      0.0406   18.6108    0.0000 [ 0.6521; 0.8620 ] 
#>   eta2 =~ y23      0.7997      0.0359   22.3001    0.0000 [ 0.7157; 0.9012 ] 
#>   eta3 =~ y31      0.8223      0.0345   23.8632    0.0000 [ 0.7176; 0.8958 ] 
#>   eta3 =~ y32      0.6581      0.0371   17.7539    0.0000 [ 0.5636; 0.7553 ] 
#>   eta3 =~ y33      0.7474      0.0395   18.9275    0.0000 [ 0.6554; 0.8596 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0202   19.5910    0.0000 [ 0.3437; 0.4481 ] 
#>   eta1 <~ y12      0.3873      0.0259   14.9797    0.0000 [ 0.3254; 0.4591 ] 
#>   eta1 <~ y13      0.4542      0.0217   20.9494    0.0000 [ 0.3946; 0.5067 ] 
#>   eta2 <~ y21      0.3058      0.0252   12.1316    0.0000 [ 0.2345; 0.3649 ] 
#>   eta2 <~ y22      0.4473      0.0194   23.0506    0.0000 [ 0.3977; 0.4980 ] 
#>   eta2 <~ y23      0.4735      0.0186   25.5215    0.0000 [ 0.4302; 0.5261 ] 
#>   eta3 <~ y31      0.4400      0.0207   21.2940    0.0000 [ 0.3796; 0.4865 ] 
#>   eta3 <~ y32      0.3521      0.0185   19.0776    0.0000 [ 0.3064; 0.4018 ] 
#>   eta3 <~ y33      0.3999      0.0159   25.1457    0.0000 [ 0.3657; 0.4480 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0408   16.4419    0.0000 [ 0.5571; 0.7682 ] 
#>   eta3 ~ eta1       0.6634      0.0385   17.2093    0.0000 [ 0.5678; 0.7671 ] 
#>   eta3 ~ eta2       0.3052      0.1016    3.0032    0.0027 [ 0.0330; 0.5584 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0691    2.9649    0.0030 [ 0.0183; 0.3756 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04083064 16.441902 9.586283e-61
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.09392245  4.881759 1.051438e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10161006  3.003159 2.671931e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55707854          0.7682316         0.58243389          0.7428762
#> 2         0.22764247          0.7133564         0.28596721          0.6550316
#> 3         0.03296863          0.5584385         0.09606729          0.4953398
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.61620752          0.7746139         0.62162372          0.7326808
#> 2         0.23646558          0.7284386         0.30218349          0.6116960
#> 3         0.03651861          0.5027369         0.08495699          0.4708621
```
