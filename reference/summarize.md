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
#>  Random seed                        = -573941024
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
#>   eta2 ~ eta1      0.6713      0.0481   13.9557    0.0000 [ 0.6072; 0.7673 ] 
#>   eta3 ~ eta1      0.4585      0.0871    5.2650    0.0000 [ 0.3380; 0.6743 ] 
#>   eta3 ~ eta2      0.3052      0.0836    3.6505    0.0003 [ 0.1110; 0.4205 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0404   16.3975    0.0000 [ 0.5808; 0.7420 ] 
#>   eta1 =~ y12      0.6493      0.0340   19.1107    0.0000 [ 0.5841; 0.7074 ] 
#>   eta1 =~ y13      0.7613      0.0329   23.1250    0.0000 [ 0.6899; 0.8111 ] 
#>   eta2 =~ y21      0.5165      0.0455   11.3598    0.0000 [ 0.4457; 0.6191 ] 
#>   eta2 =~ y22      0.7554      0.0389   19.4030    0.0000 [ 0.6826; 0.8196 ] 
#>   eta2 =~ y23      0.7997      0.0365   21.9230    0.0000 [ 0.7300; 0.8619 ] 
#>   eta3 =~ y31      0.8223      0.0343   23.9911    0.0000 [ 0.7666; 0.8770 ] 
#>   eta3 =~ y32      0.6581      0.0386   17.0472    0.0000 [ 0.5702; 0.7074 ] 
#>   eta3 =~ y33      0.7474      0.0427   17.5238    0.0000 [ 0.6588; 0.8339 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0190   20.8674    0.0000 [ 0.3675; 0.4357 ] 
#>   eta1 <~ y12      0.3873      0.0206   18.7708    0.0000 [ 0.3495; 0.4211 ] 
#>   eta1 <~ y13      0.4542      0.0191   23.7725    0.0000 [ 0.4213; 0.4895 ] 
#>   eta2 <~ y21      0.3058      0.0271   11.2962    0.0000 [ 0.2663; 0.3690 ] 
#>   eta2 <~ y22      0.4473      0.0178   25.1675    0.0000 [ 0.4155; 0.4730 ] 
#>   eta2 <~ y23      0.4735      0.0208   22.7268    0.0000 [ 0.4329; 0.5114 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5400    0.0000 [ 0.4131; 0.4791 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1516    0.0000 [ 0.3171; 0.3770 ] 
#>   eta3 <~ y33      0.3999      0.0230   17.3664    0.0000 [ 0.3632; 0.4552 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0481   13.9557    0.0000 [ 0.6072; 0.7673 ] 
#>   eta3 ~ eta1       0.6634      0.0411   16.1275    0.0000 [ 0.5998; 0.7390 ] 
#>   eta3 ~ eta2       0.3052      0.0836    3.6505    0.0003 [ 0.1110; 0.4205 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0569    3.6026    0.0003 [ 0.0703; 0.2932 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04043721 16.39752  1.992161e-60
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03397459 19.11069  2.057110e-81
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03292304 23.12502 2.594120e-118
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04546336 11.35980  6.629438e-30
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03893152 19.40298  7.282201e-84
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03647608 21.92296 1.569115e-106
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03427430 23.99108 3.445936e-127
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03860285 17.04716  3.669298e-65
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04265187 17.52383  9.426097e-69
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5807571          0.7420448
#> 2          0.5840504          0.7073712
#> 3          0.6899278          0.8111353
#> 4          0.4456902          0.6190980
#> 5          0.6826491          0.8196479
#> 6          0.7299955          0.8618828
#> 7          0.7666229          0.8770072
#> 8          0.5702097          0.7074168
#> 9          0.6588020          0.8339195

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
#>  Random seed                        = -573941024
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
#>   eta2 ~ eta1      0.6713      0.0481   13.9557    0.0000 [ 0.5378; 0.7865 ] 
#>   eta3 ~ eta1      0.4585      0.0871    5.2650    0.0000 [ 0.2158; 0.6661 ] 
#>   eta3 ~ eta2      0.3052      0.0836    3.6505    0.0003 [ 0.1023; 0.5346 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0404   16.3975    0.0000 [ 0.5614; 0.7705 ] 
#>   eta1 =~ y12      0.6493      0.0340   19.1107    0.0000 [ 0.5666; 0.7423 ] 
#>   eta1 =~ y13      0.7613      0.0329   23.1250    0.0000 [ 0.6828; 0.8531 ] 
#>   eta2 =~ y21      0.5165      0.0455   11.3598    0.0000 [ 0.3910; 0.6261 ] 
#>   eta2 =~ y22      0.7554      0.0389   19.4030    0.0000 [ 0.6583; 0.8596 ] 
#>   eta2 =~ y23      0.7997      0.0365   21.9230    0.0000 [ 0.7030; 0.8916 ] 
#>   eta3 =~ y31      0.8223      0.0343   23.9911    0.0000 [ 0.7366; 0.9139 ] 
#>   eta3 =~ y32      0.6581      0.0386   17.0472    0.0000 [ 0.5695; 0.7691 ] 
#>   eta3 =~ y33      0.7474      0.0427   17.5238    0.0000 [ 0.6365; 0.8571 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0190   20.8674    0.0000 [ 0.3445; 0.4425 ] 
#>   eta1 <~ y12      0.3873      0.0206   18.7708    0.0000 [ 0.3330; 0.4397 ] 
#>   eta1 <~ y13      0.4542      0.0191   23.7725    0.0000 [ 0.4041; 0.5029 ] 
#>   eta2 <~ y21      0.3058      0.0271   11.2962    0.0000 [ 0.2326; 0.3726 ] 
#>   eta2 <~ y22      0.4473      0.0178   25.1675    0.0000 [ 0.4058; 0.4977 ] 
#>   eta2 <~ y23      0.4735      0.0208   22.7268    0.0000 [ 0.4205; 0.5282 ] 
#>   eta3 <~ y31      0.4400      0.0195   22.5400    0.0000 [ 0.3872; 0.4882 ] 
#>   eta3 <~ y32      0.3521      0.0166   21.1516    0.0000 [ 0.3123; 0.3984 ] 
#>   eta3 <~ y33      0.3999      0.0230   17.3664    0.0000 [ 0.3366; 0.4557 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0481   13.9557    0.0000 [ 0.5378; 0.7865 ] 
#>   eta3 ~ eta1       0.6634      0.0411   16.1275    0.0000 [ 0.5462; 0.7589 ] 
#>   eta3 ~ eta2       0.3052      0.0836    3.6505    0.0003 [ 0.1023; 0.5346 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0569    3.6026    0.0003 [ 0.0646; 0.3586 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04810444 13.955748 2.902421e-44
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.08708515  5.265040 1.401592e-07
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08359065  3.650541 2.616884e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5377748          0.7865438          0.5676471          0.7566715
#> 2          0.2157800          0.6661352          0.2698588          0.6120564
#> 3          0.1022958          0.5345794          0.1542046          0.4826706
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1         0.57058377          0.7738415          0.6071671          0.7673078
#> 2         0.26917950          0.6836862          0.3380054          0.6742888
#> 3         0.07806668          0.4276433          0.1109551          0.4204998
```
