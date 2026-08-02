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
#>  Random seed                        = 284123341
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
#>   eta2 ~ eta1      0.6713      0.0392   17.1193    0.0000 [ 0.6177; 0.7491 ] 
#>   eta3 ~ eta1      0.4585      0.0754    6.0791    0.0000 [ 0.3197; 0.6048 ] 
#>   eta3 ~ eta2      0.3052      0.0888    3.4346    0.0006 [ 0.1647; 0.4936 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0406   16.3221    0.0000 [ 0.5948; 0.7350 ] 
#>   eta1 =~ y12      0.6493      0.0474   13.7044    0.0000 [ 0.5822; 0.7399 ] 
#>   eta1 =~ y13      0.7613      0.0302   25.2142    0.0000 [ 0.7050; 0.8108 ] 
#>   eta2 =~ y21      0.5165      0.0438   11.7780    0.0000 [ 0.4326; 0.5852 ] 
#>   eta2 =~ y22      0.7554      0.0436   17.3086    0.0000 [ 0.6817; 0.8462 ] 
#>   eta2 =~ y23      0.7997      0.0346   23.1012    0.0000 [ 0.7378; 0.8452 ] 
#>   eta3 =~ y31      0.8223      0.0366   22.4745    0.0000 [ 0.7588; 0.8848 ] 
#>   eta3 =~ y32      0.6581      0.0388   16.9419    0.0000 [ 0.5718; 0.7210 ] 
#>   eta3 =~ y33      0.7474      0.0374   19.9580    0.0000 [ 0.6588; 0.7955 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0191   20.7373    0.0000 [ 0.3578; 0.4324 ] 
#>   eta1 <~ y12      0.3873      0.0226   17.1127    0.0000 [ 0.3525; 0.4268 ] 
#>   eta1 <~ y13      0.4542      0.0228   19.9361    0.0000 [ 0.4088; 0.4882 ] 
#>   eta2 <~ y21      0.3058      0.0224   13.6666    0.0000 [ 0.2587; 0.3412 ] 
#>   eta2 <~ y22      0.4473      0.0245   18.2754    0.0000 [ 0.4142; 0.5019 ] 
#>   eta2 <~ y23      0.4735      0.0206   22.9376    0.0000 [ 0.4353; 0.5161 ] 
#>   eta3 <~ y31      0.4400      0.0170   25.8646    0.0000 [ 0.4117; 0.4739 ] 
#>   eta3 <~ y32      0.3521      0.0180   19.5480    0.0000 [ 0.3259; 0.3852 ] 
#>   eta3 <~ y33      0.3999      0.0209   19.1415    0.0000 [ 0.3621; 0.4406 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0392   17.1193    0.0000 [ 0.6177; 0.7491 ] 
#>   eta3 ~ eta1       0.6634      0.0287   23.0816    0.0000 [ 0.6263; 0.7300 ] 
#>   eta3 ~ eta2       0.3052      0.0888    3.4346    0.0006 [ 0.1647; 0.4936 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0575    3.5618    0.0004 [ 0.1171; 0.3255 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04062411 16.32208  6.875617e-60
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04737746 13.70436  9.560193e-43
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03019511 25.21421 2.798335e-140
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04384915 11.77799  5.068861e-32
#> 5 eta2 =~ y22  Common factor 0.7553877 0.04364245 17.30855  4.054793e-67
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03461573 23.10117 4.506777e-118
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03658716 22.47448 7.376057e-112
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03884271 16.94189  2.209135e-64
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03744983 19.95801  1.276886e-88
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5947632          0.7350073
#> 2          0.5821924          0.7398987
#> 3          0.7050429          0.8108004
#> 4          0.4326201          0.5851772
#> 5          0.6816719          0.8461707
#> 6          0.7378004          0.8452254
#> 7          0.7587647          0.8847972
#> 8          0.5717608          0.7209966
#> 9          0.6588223          0.7955064

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
#>  Random seed                        = 284123341
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
#>   eta2 ~ eta1      0.6713      0.0392   17.1193    0.0000 [ 0.5664; 0.7692 ] 
#>   eta3 ~ eta1      0.4585      0.0754    6.0791    0.0000 [ 0.2398; 0.6299 ] 
#>   eta3 ~ eta2      0.3052      0.0888    3.4346    0.0006 [ 0.0904; 0.5498 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0406   16.3221    0.0000 [ 0.5445; 0.7546 ] 
#>   eta1 =~ y12      0.6493      0.0474   13.7044    0.0000 [ 0.5268; 0.7718 ] 
#>   eta1 =~ y13      0.7613      0.0302   25.2142    0.0000 [ 0.6857; 0.8419 ] 
#>   eta2 =~ y21      0.5165      0.0438   11.7780    0.0000 [ 0.4160; 0.6427 ] 
#>   eta2 =~ y22      0.7554      0.0436   17.3086    0.0000 [ 0.6325; 0.8582 ] 
#>   eta2 =~ y23      0.7997      0.0346   23.1012    0.0000 [ 0.7127; 0.8917 ] 
#>   eta3 =~ y31      0.8223      0.0366   22.4745    0.0000 [ 0.7320; 0.9212 ] 
#>   eta3 =~ y32      0.6581      0.0388   16.9419    0.0000 [ 0.5479; 0.7488 ] 
#>   eta3 =~ y33      0.7474      0.0374   19.9580    0.0000 [ 0.6550; 0.8487 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0191   20.7373    0.0000 [ 0.3416; 0.4403 ] 
#>   eta1 <~ y12      0.3873      0.0226   17.1127    0.0000 [ 0.3323; 0.4493 ] 
#>   eta1 <~ y13      0.4542      0.0228   19.9361    0.0000 [ 0.4002; 0.5180 ] 
#>   eta2 <~ y21      0.3058      0.0224   13.6666    0.0000 [ 0.2553; 0.3710 ] 
#>   eta2 <~ y22      0.4473      0.0245   18.2754    0.0000 [ 0.3773; 0.5039 ] 
#>   eta2 <~ y23      0.4735      0.0206   22.9376    0.0000 [ 0.4207; 0.5275 ] 
#>   eta3 <~ y31      0.4400      0.0170   25.8646    0.0000 [ 0.3986; 0.4866 ] 
#>   eta3 <~ y32      0.3521      0.0180   19.5480    0.0000 [ 0.3007; 0.3939 ] 
#>   eta3 <~ y33      0.3999      0.0209   19.1415    0.0000 [ 0.3485; 0.4565 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0392   17.1193    0.0000 [ 0.5664; 0.7692 ] 
#>   eta3 ~ eta1       0.6634      0.0287   23.0816    0.0000 [ 0.5755; 0.7241 ] 
#>   eta3 ~ eta2       0.3052      0.0888    3.4346    0.0006 [ 0.0904; 0.5498 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0575    3.5618    0.0004 [ 0.0663; 0.3637 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03921491 17.119343 1.064695e-65
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07542312  6.079127 1.208384e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08884508  3.434643 5.933348e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.56642999          0.7692273          0.5907820          0.7448753
#> 2         0.23980686          0.6298527          0.2866437          0.5830158
#> 3         0.09035924          0.5498158          0.1455310          0.4946441
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5954490          0.7717187          0.6177240          0.7491351
#> 2          0.3052927          0.6170484          0.3196781          0.6047706
#> 3          0.1402464          0.4940417          0.1646726          0.4935862
```
