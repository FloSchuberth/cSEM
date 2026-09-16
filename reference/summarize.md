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
#>  Random seed                        = 2057376503
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
#>   eta2 ~ eta1      0.6713      0.0480   13.9774    0.0000 [ 0.5808; 0.7418 ] 
#>   eta3 ~ eta1      0.4585      0.0888    5.1617    0.0000 [ 0.3365; 0.6507 ] 
#>   eta3 ~ eta2      0.3052      0.0936    3.2617    0.0011 [ 0.1203; 0.4610 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0359   18.4518    0.0000 [ 0.5958; 0.7183 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9552    0.0000 [ 0.5809; 0.7128 ] 
#>   eta1 =~ y13      0.7613      0.0278   27.3713    0.0000 [ 0.7112; 0.8217 ] 
#>   eta2 =~ y21      0.5165      0.0562    9.1926    0.0000 [ 0.4179; 0.6046 ] 
#>   eta2 =~ y22      0.7554      0.0371   20.3429    0.0000 [ 0.6893; 0.8206 ] 
#>   eta2 =~ y23      0.7997      0.0405   19.7361    0.0000 [ 0.7322; 0.8729 ] 
#>   eta3 =~ y31      0.8223      0.0313   26.2572    0.0000 [ 0.7653; 0.8698 ] 
#>   eta3 =~ y32      0.6581      0.0353   18.6468    0.0000 [ 0.5943; 0.7020 ] 
#>   eta3 =~ y33      0.7474      0.0475   15.7386    0.0000 [ 0.6630; 0.8469 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0207   19.1510    0.0000 [ 0.3541; 0.4377 ] 
#>   eta1 <~ y12      0.3873      0.0182   21.2989    0.0000 [ 0.3518; 0.4166 ] 
#>   eta1 <~ y13      0.4542      0.0184   24.7045    0.0000 [ 0.4281; 0.4980 ] 
#>   eta2 <~ y21      0.3058      0.0306    9.9881    0.0000 [ 0.2506; 0.3553 ] 
#>   eta2 <~ y22      0.4473      0.0218   20.5002    0.0000 [ 0.4080; 0.4806 ] 
#>   eta2 <~ y23      0.4735      0.0234   20.2560    0.0000 [ 0.4253; 0.5105 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.7591    0.0000 [ 0.4062; 0.4706 ] 
#>   eta3 <~ y32      0.3521      0.0201   17.4823    0.0000 [ 0.3196; 0.3823 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4556    0.0000 [ 0.3617; 0.4334 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0480   13.9774    0.0000 [ 0.5808; 0.7418 ] 
#>   eta3 ~ eta1       0.6634      0.0424   15.6605    0.0000 [ 0.6011; 0.7425 ] 
#>   eta3 ~ eta2       0.3052      0.0936    3.2617    0.0011 [ 0.1203; 0.4610 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0624    3.2817    0.0010 [ 0.0903; 0.2969 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03593528 18.451779  5.045639e-76
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03829369 16.955222  1.760949e-64
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02781546 27.371322 6.021152e-165
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05618188  9.192552  3.836300e-20
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03713268 20.342933  5.362576e-92
#> 6 eta2 =~ y23  Common factor 0.7996637 0.04051781 19.736104  1.056152e-86
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03131624 26.257214 5.912726e-152
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03529126 18.646793  1.340817e-77
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04748992 15.738586  8.226164e-56
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5958284          0.7182618
#> 2          0.5808833          0.7127574
#> 3          0.7111643          0.8216891
#> 4          0.4178782          0.6045763
#> 5          0.6893102          0.8206416
#> 6          0.7322059          0.8729007
#> 7          0.7652797          0.8697752
#> 8          0.5943486          0.7019848
#> 9          0.6630451          0.8468526

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
#>  Random seed                        = 2057376503
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
#>   eta2 ~ eta1      0.6713      0.0480   13.9774    0.0000 [ 0.5438; 0.7922 ] 
#>   eta3 ~ eta1      0.4585      0.0888    5.1617    0.0000 [ 0.2246; 0.6840 ] 
#>   eta3 ~ eta2      0.3052      0.0936    3.2617    0.0011 [ 0.0727; 0.5565 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0359   18.4518    0.0000 [ 0.5772; 0.7630 ] 
#>   eta1 =~ y12      0.6493      0.0383   16.9552    0.0000 [ 0.5494; 0.7474 ] 
#>   eta1 =~ y13      0.7613      0.0278   27.3713    0.0000 [ 0.6886; 0.8325 ] 
#>   eta2 =~ y21      0.5165      0.0562    9.1926    0.0000 [ 0.3738; 0.6644 ] 
#>   eta2 =~ y22      0.7554      0.0371   20.3429    0.0000 [ 0.6560; 0.8480 ] 
#>   eta2 =~ y23      0.7997      0.0405   19.7361    0.0000 [ 0.6960; 0.9055 ] 
#>   eta3 =~ y31      0.8223      0.0313   26.2572    0.0000 [ 0.7374; 0.8994 ] 
#>   eta3 =~ y32      0.6581      0.0353   18.6468    0.0000 [ 0.5710; 0.7535 ] 
#>   eta3 =~ y33      0.7474      0.0475   15.7386    0.0000 [ 0.6233; 0.8689 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0207   19.1510    0.0000 [ 0.3450; 0.4518 ] 
#>   eta1 <~ y12      0.3873      0.0182   21.2989    0.0000 [ 0.3387; 0.4327 ] 
#>   eta1 <~ y13      0.4542      0.0184   24.7045    0.0000 [ 0.4045; 0.4996 ] 
#>   eta2 <~ y21      0.3058      0.0306    9.9881    0.0000 [ 0.2288; 0.3872 ] 
#>   eta2 <~ y22      0.4473      0.0218   20.5002    0.0000 [ 0.3894; 0.5022 ] 
#>   eta2 <~ y23      0.4735      0.0234   20.2560    0.0000 [ 0.4143; 0.5352 ] 
#>   eta3 <~ y31      0.4400      0.0185   23.7591    0.0000 [ 0.3905; 0.4863 ] 
#>   eta3 <~ y32      0.3521      0.0201   17.4823    0.0000 [ 0.3027; 0.4068 ] 
#>   eta3 <~ y33      0.3999      0.0196   20.4556    0.0000 [ 0.3495; 0.4506 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0480   13.9774    0.0000 [ 0.5438; 0.7922 ] 
#>   eta3 ~ eta1       0.6634      0.0424   15.6605    0.0000 [ 0.5560; 0.7751 ] 
#>   eta3 ~ eta2       0.3052      0.0936    3.2617    0.0011 [ 0.0727; 0.5565 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0624    3.2817    0.0010 [ 0.0499; 0.3727 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04803008 13.977355 2.143052e-44
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08882942  5.161654 2.447770e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.09355545  3.261714 1.107408e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5437689          0.7921533          0.5735950          0.7623272
#> 2          0.2245960          0.6839716          0.2797581          0.6288096
#> 3          0.0727328          0.5565488          0.1308296          0.4984519
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5776790          0.7590756          0.5808235          0.7417806
#> 2          0.3202873          0.6579863          0.3364951          0.6507260
#> 3          0.1187830          0.4677341          0.1203251          0.4610030
```
