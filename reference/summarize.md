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
#>  Random seed                        = -676976479
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
#>   eta2 ~ eta1      0.6713      0.0392   17.1196    0.0000 [ 0.6021; 0.7347 ] 
#>   eta3 ~ eta1      0.4585      0.1028    4.4598    0.0000 [ 0.3157; 0.6499 ] 
#>   eta3 ~ eta2      0.3052      0.1089    2.8014    0.0051 [ 0.0548; 0.4361 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0374   17.7231    0.0000 [ 0.5877; 0.7135 ] 
#>   eta1 =~ y12      0.6493      0.0405   16.0177    0.0000 [ 0.5879; 0.7004 ] 
#>   eta1 =~ y13      0.7613      0.0330   23.0707    0.0000 [ 0.6912; 0.8221 ] 
#>   eta2 =~ y21      0.5165      0.0534    9.6648    0.0000 [ 0.4038; 0.5908 ] 
#>   eta2 =~ y22      0.7554      0.0370   20.4331    0.0000 [ 0.6943; 0.8115 ] 
#>   eta2 =~ y23      0.7997      0.0342   23.4087    0.0000 [ 0.7487; 0.8558 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.7171    0.0000 [ 0.7642; 0.8873 ] 
#>   eta3 =~ y32      0.6581      0.0345   19.0792    0.0000 [ 0.5912; 0.7216 ] 
#>   eta3 =~ y33      0.7474      0.0423   17.6901    0.0000 [ 0.6654; 0.8389 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0223   17.7133    0.0000 [ 0.3630; 0.4316 ] 
#>   eta1 <~ y12      0.3873      0.0184   20.9963    0.0000 [ 0.3498; 0.4313 ] 
#>   eta1 <~ y13      0.4542      0.0191   23.7306    0.0000 [ 0.4216; 0.4941 ] 
#>   eta2 <~ y21      0.3058      0.0285   10.7221    0.0000 [ 0.2428; 0.3449 ] 
#>   eta2 <~ y22      0.4473      0.0224   19.9677    0.0000 [ 0.4150; 0.4878 ] 
#>   eta2 <~ y23      0.4735      0.0184   25.7070    0.0000 [ 0.4461; 0.5123 ] 
#>   eta3 <~ y31      0.4400      0.0193   22.8033    0.0000 [ 0.3975; 0.4801 ] 
#>   eta3 <~ y32      0.3521      0.0158   22.2985    0.0000 [ 0.3222; 0.3760 ] 
#>   eta3 <~ y33      0.3999      0.0206   19.3919    0.0000 [ 0.3662; 0.4349 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0392   17.1196    0.0000 [ 0.6021; 0.7347 ] 
#>   eta3 ~ eta1       0.6634      0.0407   16.2894    0.0000 [ 0.6111; 0.7433 ] 
#>   eta3 ~ eta2       0.3052      0.1089    2.8014    0.0051 [ 0.0548; 0.4361 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0725    2.8271    0.0047 [ 0.0426; 0.3002 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03741276 17.723096  2.781712e-70
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04053509 16.017674  9.618098e-58
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03300060 23.070669 9.124249e-118
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05343688  9.664764  4.256002e-22
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03696879 20.433118  8.490308e-93
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03416100 23.408674 3.487143e-121
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03619645 22.717071 3.038039e-114
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03449141 19.079209  3.759036e-81
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04225104 17.690077  5.000721e-70
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5877300          0.7134863
#> 2          0.5879266          0.7003760
#> 3          0.6911752          0.8221460
#> 4          0.4038198          0.5907715
#> 5          0.6942745          0.8114860
#> 6          0.7486961          0.8558004
#> 7          0.7641510          0.8873165
#> 8          0.5912076          0.7216199
#> 9          0.6654304          0.8388717

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
#>  Random seed                        = -676976479
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
#>   eta2 ~ eta1      0.6713      0.0392   17.1196    0.0000 [ 0.5730; 0.7757 ] 
#>   eta3 ~ eta1      0.4585      0.1028    4.4598    0.0000 [ 0.1643; 0.6960 ] 
#>   eta3 ~ eta2      0.3052      0.1089    2.8014    0.0051 [ 0.0584; 0.6218 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0374   17.7231    0.0000 [ 0.5780; 0.7715 ] 
#>   eta1 =~ y12      0.6493      0.0405   16.0177    0.0000 [ 0.5475; 0.7571 ] 
#>   eta1 =~ y13      0.7613      0.0330   23.0707    0.0000 [ 0.6744; 0.8451 ] 
#>   eta2 =~ y21      0.5165      0.0534    9.6648    0.0000 [ 0.3927; 0.6691 ] 
#>   eta2 =~ y22      0.7554      0.0370   20.4331    0.0000 [ 0.6603; 0.8515 ] 
#>   eta2 =~ y23      0.7997      0.0342   23.4087    0.0000 [ 0.7095; 0.8861 ] 
#>   eta3 =~ y31      0.8223      0.0362   22.7171    0.0000 [ 0.7251; 0.9123 ] 
#>   eta3 =~ y32      0.6581      0.0345   19.0792    0.0000 [ 0.5793; 0.7577 ] 
#>   eta3 =~ y33      0.7474      0.0423   17.6901    0.0000 [ 0.6334; 0.8519 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0223   17.7133    0.0000 [ 0.3413; 0.4568 ] 
#>   eta1 <~ y12      0.3873      0.0184   20.9963    0.0000 [ 0.3382; 0.4336 ] 
#>   eta1 <~ y13      0.4542      0.0191   23.7306    0.0000 [ 0.3996; 0.4986 ] 
#>   eta2 <~ y21      0.3058      0.0285   10.7221    0.0000 [ 0.2388; 0.3863 ] 
#>   eta2 <~ y22      0.4473      0.0224   19.9677    0.0000 [ 0.3865; 0.5023 ] 
#>   eta2 <~ y23      0.4735      0.0184   25.7070    0.0000 [ 0.4215; 0.5167 ] 
#>   eta3 <~ y31      0.4400      0.0193   22.8033    0.0000 [ 0.3878; 0.4876 ] 
#>   eta3 <~ y32      0.3521      0.0158   22.2985    0.0000 [ 0.3167; 0.3984 ] 
#>   eta3 <~ y33      0.3999      0.0206   19.3919    0.0000 [ 0.3438; 0.4505 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0392   17.1196    0.0000 [ 0.5730; 0.7757 ] 
#>   eta3 ~ eta1       0.6634      0.0407   16.2894    0.0000 [ 0.5550; 0.7656 ] 
#>   eta3 ~ eta2       0.3052      0.1089    2.8014    0.0051 [ 0.0584; 0.6218 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0725    2.8271    0.0047 [ 0.0428; 0.4176 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03921434 17.119592 1.060167e-65
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.10280902  4.459791 8.203956e-06
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.10892874  2.801383 5.088406e-03
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5729537          0.7757481          0.5973053          0.7513964
#> 2          0.1643128          0.6959831          0.2281560          0.6321399
#> 3          0.0584328          0.6217507          0.1260763          0.5541072
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.59634069          0.7708727         0.60206209          0.7347203
#> 2         0.28651595          0.6602413         0.31566112          0.6499195
#> 3         0.05377212          0.4945724         0.05476141          0.4361417
```
