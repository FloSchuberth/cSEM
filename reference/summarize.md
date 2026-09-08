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
#>  Random seed                        = 1520353910
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
#>   eta2 ~ eta1      0.6713      0.0384   17.4839    0.0000 [ 0.6080; 0.7351 ] 
#>   eta3 ~ eta1      0.4585      0.0669    6.8487    0.0000 [ 0.3432; 0.5748 ] 
#>   eta3 ~ eta2      0.3052      0.0723    4.2179    0.0000 [ 0.1900; 0.4410 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0376   17.6492    0.0000 [ 0.6014; 0.7392 ] 
#>   eta1 =~ y12      0.6493      0.0433   14.9829    0.0000 [ 0.5523; 0.7116 ] 
#>   eta1 =~ y13      0.7613      0.0289   26.3258    0.0000 [ 0.7140; 0.8074 ] 
#>   eta2 =~ y21      0.5165      0.0480   10.7562    0.0000 [ 0.4233; 0.5937 ] 
#>   eta2 =~ y22      0.7554      0.0349   21.6226    0.0000 [ 0.7049; 0.8243 ] 
#>   eta2 =~ y23      0.7997      0.0396   20.1954    0.0000 [ 0.7390; 0.8976 ] 
#>   eta3 =~ y31      0.8223      0.0337   24.4103    0.0000 [ 0.7650; 0.8766 ] 
#>   eta3 =~ y32      0.6581      0.0342   19.2570    0.0000 [ 0.6066; 0.7187 ] 
#>   eta3 =~ y33      0.7474      0.0318   23.5152    0.0000 [ 0.6598; 0.7913 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0211   18.7637    0.0000 [ 0.3673; 0.4395 ] 
#>   eta1 <~ y12      0.3873      0.0172   22.4625    0.0000 [ 0.3509; 0.4121 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.1008    0.0000 [ 0.4123; 0.4914 ] 
#>   eta2 <~ y21      0.3058      0.0258   11.8303    0.0000 [ 0.2523; 0.3477 ] 
#>   eta2 <~ y22      0.4473      0.0188   23.7445    0.0000 [ 0.4141; 0.4817 ] 
#>   eta2 <~ y23      0.4735      0.0225   21.0054    0.0000 [ 0.4266; 0.5043 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.8491    0.0000 [ 0.4102; 0.4723 ] 
#>   eta3 <~ y32      0.3521      0.0169   20.8814    0.0000 [ 0.3316; 0.3852 ] 
#>   eta3 <~ y33      0.3999      0.0161   24.8661    0.0000 [ 0.3620; 0.4180 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0384   17.4839    0.0000 [ 0.6080; 0.7351 ] 
#>   eta3 ~ eta1       0.6634      0.0381   17.4187    0.0000 [ 0.6013; 0.7376 ] 
#>   eta3 ~ eta2       0.3052      0.0723    4.2179    0.0000 [ 0.1900; 0.4410 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0487    4.2080    0.0000 [ 0.1263; 0.2993 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err   t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.03756940 17.64920  1.032087e-69
#> 2 eta1 =~ y12  Common factor 0.6492779 0.04333448 14.98294  9.492203e-51
#> 3 eta1 =~ y13  Common factor 0.7613458 0.02892015 26.32579 9.718976e-153
#> 4 eta2 =~ y21  Common factor 0.5164548 0.04801473 10.75617  5.542355e-27
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03493514 21.62257 1.101619e-103
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03959640 20.19536  1.075363e-90
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03368561 24.41034 1.328082e-131
#> 8 eta3 =~ y32  Common factor 0.6580689 0.03417302 19.25697  1.233801e-82
#> 9 eta3 =~ y33  Common factor 0.7474241 0.03178474 23.51518 2.852500e-122
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.6013776          0.7391574
#> 2          0.5523497          0.7115789
#> 3          0.7140494          0.8074379
#> 4          0.4233302          0.5937212
#> 5          0.7049092          0.8243037
#> 6          0.7389974          0.8976095
#> 7          0.7650435          0.8766429
#> 8          0.6066127          0.7187244
#> 9          0.6597590          0.7912754

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
#>  Random seed                        = 1520353910
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
#>   eta2 ~ eta1      0.6713      0.0384   17.4839    0.0000 [ 0.5652; 0.7638 ] 
#>   eta3 ~ eta1      0.4585      0.0669    6.8487    0.0000 [ 0.2751; 0.6214 ] 
#>   eta3 ~ eta2      0.3052      0.0723    4.2179    0.0000 [ 0.1207; 0.4949 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0376   17.6492    0.0000 [ 0.5678; 0.7621 ] 
#>   eta1 =~ y12      0.6493      0.0433   14.9829    0.0000 [ 0.5347; 0.7588 ] 
#>   eta1 =~ y13      0.7613      0.0289   26.3258    0.0000 [ 0.6830; 0.8326 ] 
#>   eta2 =~ y21      0.5165      0.0480   10.7562    0.0000 [ 0.3858; 0.6341 ] 
#>   eta2 =~ y22      0.7554      0.0349   21.6226    0.0000 [ 0.6625; 0.8432 ] 
#>   eta2 =~ y23      0.7997      0.0396   20.1954    0.0000 [ 0.7034; 0.9081 ] 
#>   eta3 =~ y31      0.8223      0.0337   24.4103    0.0000 [ 0.7366; 0.9108 ] 
#>   eta3 =~ y32      0.6581      0.0342   19.2570    0.0000 [ 0.5586; 0.7353 ] 
#>   eta3 =~ y33      0.7474      0.0318   23.5152    0.0000 [ 0.6820; 0.8464 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0211   18.7637    0.0000 [ 0.3435; 0.4525 ] 
#>   eta1 <~ y12      0.3873      0.0172   22.4625    0.0000 [ 0.3430; 0.4321 ] 
#>   eta1 <~ y13      0.4542      0.0206   22.1008    0.0000 [ 0.4003; 0.5065 ] 
#>   eta2 <~ y21      0.3058      0.0258   11.8303    0.0000 [ 0.2360; 0.3697 ] 
#>   eta2 <~ y22      0.4473      0.0188   23.7445    0.0000 [ 0.3980; 0.4954 ] 
#>   eta2 <~ y23      0.4735      0.0225   21.0054    0.0000 [ 0.4198; 0.5364 ] 
#>   eta3 <~ y31      0.4400      0.0177   24.8491    0.0000 [ 0.3928; 0.4844 ] 
#>   eta3 <~ y32      0.3521      0.0169   20.8814    0.0000 [ 0.3008; 0.3880 ] 
#>   eta3 <~ y33      0.3999      0.0161   24.8661    0.0000 [ 0.3654; 0.4486 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0384   17.4839    0.0000 [ 0.5652; 0.7638 ] 
#>   eta3 ~ eta1       0.6634      0.0381   17.4187    0.0000 [ 0.5549; 0.7518 ] 
#>   eta3 ~ eta2       0.3052      0.0723    4.2179    0.0000 [ 0.1207; 0.4949 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0487    4.2080    0.0000 [ 0.0792; 0.3310 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03839715 17.483937 1.899207e-68
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.06694847  6.848652 7.454907e-12
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.07234721  4.217870 2.466207e-05
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5652045          0.7637729          0.5890487          0.7399287
#> 2          0.2751429          0.6213626          0.3167171          0.5797883
#> 3          0.1207343          0.4948732          0.1656611          0.4499464
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5970721          0.7443269          0.6080260          0.7351151
#> 2          0.3423066          0.5879622          0.3431872          0.5747811
#> 3          0.1371195          0.4494637          0.1899802          0.4410489
```
