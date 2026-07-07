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
#>  Random seed                        = 1063704626
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
#>   eta2 ~ eta1      0.6713      0.0435   15.4327    0.0000 [ 0.5714; 0.7407 ] 
#>   eta3 ~ eta1      0.4585      0.0771    5.9507    0.0000 [ 0.3588; 0.6415 ] 
#>   eta3 ~ eta2      0.3052      0.0854    3.5721    0.0004 [ 0.1054; 0.4482 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0367   18.0586    0.0000 [ 0.6018; 0.7181 ] 
#>   eta1 =~ y12      0.6493      0.0430   15.0839    0.0000 [ 0.5845; 0.7374 ] 
#>   eta1 =~ y13      0.7613      0.0341   22.2943    0.0000 [ 0.6834; 0.8023 ] 
#>   eta2 =~ y21      0.5165      0.0471   10.9703    0.0000 [ 0.4530; 0.6008 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6644    0.0000 [ 0.6871; 0.8055 ] 
#>   eta2 =~ y23      0.7997      0.0348   22.9665    0.0000 [ 0.7510; 0.8858 ] 
#>   eta3 =~ y31      0.8223      0.0305   26.9926    0.0000 [ 0.7639; 0.8640 ] 
#>   eta3 =~ y32      0.6581      0.0396   16.6286    0.0000 [ 0.5937; 0.7330 ] 
#>   eta3 =~ y33      0.7474      0.0319   23.4649    0.0000 [ 0.6969; 0.8077 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0205   19.3173    0.0000 [ 0.3527; 0.4248 ] 
#>   eta1 <~ y12      0.3873      0.0205   18.9120    0.0000 [ 0.3525; 0.4315 ] 
#>   eta1 <~ y13      0.4542      0.0200   22.7019    0.0000 [ 0.4186; 0.4893 ] 
#>   eta2 <~ y21      0.3058      0.0243   12.6003    0.0000 [ 0.2722; 0.3564 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0074    0.0000 [ 0.4067; 0.4794 ] 
#>   eta2 <~ y23      0.4735      0.0173   27.4370    0.0000 [ 0.4463; 0.5086 ] 
#>   eta3 <~ y31      0.4400      0.0152   28.9317    0.0000 [ 0.4071; 0.4631 ] 
#>   eta3 <~ y32      0.3521      0.0176   19.9706    0.0000 [ 0.3216; 0.3889 ] 
#>   eta3 <~ y33      0.3999      0.0174   23.0232    0.0000 [ 0.3681; 0.4368 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0435   15.4327    0.0000 [ 0.5714; 0.7407 ] 
#>   eta3 ~ eta1       0.6634      0.0334   19.8714    0.0000 [ 0.6022; 0.7249 ] 
#>   eta3 ~ eta2       0.3052      0.0854    3.5721    0.0004 [ 0.1054; 0.4482 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0588    3.4859    0.0005 [ 0.0768; 0.3267 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03671762 18.05863  6.747554e-73
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04304445 15.08389  2.067128e-51
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03414982 22.29428 4.198445e-110
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04707740 10.97033  5.307554e-28
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03332923 22.66442 1.005666e-113
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03481873 22.96648 1.008620e-116
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03046304 26.99262 1.804262e-160
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03957446 16.62863  4.323992e-62
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03185282 23.46493 9.307953e-122
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6017950          0.7180702
#> 2          0.5844578          0.7373969
#> 3          0.6834320          0.8022553
#> 4          0.4529604          0.6007737
#> 5          0.6871457          0.8054967
#> 6          0.7510139          0.8858403
#> 7          0.7639083          0.8639952
#> 8          0.5936932          0.7330223
#> 9          0.6969210          0.8077251

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
#>  Random seed                        = 1063704626
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
#>   eta2 ~ eta1      0.6713      0.0435   15.4327    0.0000 [ 0.5553; 0.7802 ] 
#>   eta3 ~ eta1      0.4585      0.0771    5.9507    0.0000 [ 0.2499; 0.6484 ] 
#>   eta3 ~ eta2      0.3052      0.0854    3.5721    0.0004 [ 0.0909; 0.5327 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0367   18.0586    0.0000 [ 0.5712; 0.7611 ] 
#>   eta1 =~ y12      0.6493      0.0430   15.0839    0.0000 [ 0.5308; 0.7534 ] 
#>   eta1 =~ y13      0.7613      0.0341   22.2943    0.0000 [ 0.6779; 0.8545 ] 
#>   eta2 =~ y21      0.5165      0.0471   10.9703    0.0000 [ 0.3927; 0.6361 ] 
#>   eta2 =~ y22      0.7554      0.0333   22.6644    0.0000 [ 0.6717; 0.8440 ] 
#>   eta2 =~ y23      0.7997      0.0348   22.9665    0.0000 [ 0.7094; 0.8895 ] 
#>   eta3 =~ y31      0.8223      0.0305   26.9926    0.0000 [ 0.7503; 0.9079 ] 
#>   eta3 =~ y32      0.6581      0.0396   16.6286    0.0000 [ 0.5480; 0.7526 ] 
#>   eta3 =~ y33      0.7474      0.0319   23.4649    0.0000 [ 0.6712; 0.8359 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0205   19.3173    0.0000 [ 0.3443; 0.4502 ] 
#>   eta1 <~ y12      0.3873      0.0205   18.9120    0.0000 [ 0.3302; 0.4361 ] 
#>   eta1 <~ y13      0.4542      0.0200   22.7019    0.0000 [ 0.4051; 0.5086 ] 
#>   eta2 <~ y21      0.3058      0.0243   12.6003    0.0000 [ 0.2422; 0.3677 ] 
#>   eta2 <~ y22      0.4473      0.0213   21.0074    0.0000 [ 0.3937; 0.5038 ] 
#>   eta2 <~ y23      0.4735      0.0173   27.4370    0.0000 [ 0.4289; 0.5182 ] 
#>   eta3 <~ y31      0.4400      0.0152   28.9317    0.0000 [ 0.4027; 0.4814 ] 
#>   eta3 <~ y32      0.3521      0.0176   19.9706    0.0000 [ 0.3012; 0.3924 ] 
#>   eta3 <~ y33      0.3999      0.0174   23.0232    0.0000 [ 0.3568; 0.4466 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0435   15.4327    0.0000 [ 0.5553; 0.7802 ] 
#>   eta3 ~ eta1       0.6634      0.0334   19.8714    0.0000 [ 0.5715; 0.7441 ] 
#>   eta3 ~ eta2       0.3052      0.0854    3.5721    0.0004 [ 0.0909; 0.5327 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0588    3.4859    0.0005 [ 0.0567; 0.3606 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.04350075 15.432687 9.867516e-54
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07705122  5.950675 2.670389e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08542584  3.572117 3.541068e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1         0.55527490          0.7802362          0.5822884          0.7532227
#> 2         0.24993528          0.6484007          0.2977832          0.6005528
#> 3         0.09093383          0.5327081          0.1439823          0.4796596
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5622483          0.7486554          0.5714095          0.7406845
#> 2          0.3558145          0.6582279          0.3588101          0.6415384
#> 3          0.1047681          0.4645176          0.1054336          0.4481986
```
