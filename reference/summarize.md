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
#>  Random seed                        = -667070716
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
#>   eta2 ~ eta1      0.6713      0.0394   17.0535    0.0000 [ 0.6161; 0.7453 ] 
#>   eta3 ~ eta1      0.4585      0.0756    6.0617    0.0000 [ 0.3476; 0.6314 ] 
#>   eta3 ~ eta2      0.3052      0.0865    3.5298    0.0004 [ 0.1432; 0.4342 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_percentile   
#>   Loading        Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 =~ y11      0.6631      0.0429   15.4393    0.0000 [ 0.5959; 0.7426 ] 
#>   eta1 =~ y12      0.6493      0.0320   20.3125    0.0000 [ 0.5949; 0.6996 ] 
#>   eta1 =~ y13      0.7613      0.0382   19.9056    0.0000 [ 0.6772; 0.8186 ] 
#>   eta2 =~ y21      0.5165      0.0599    8.6276    0.0000 [ 0.3905; 0.6079 ] 
#>   eta2 =~ y22      0.7554      0.0358   21.0821    0.0000 [ 0.6936; 0.8163 ] 
#>   eta2 =~ y23      0.7997      0.0391   20.4394    0.0000 [ 0.7462; 0.8907 ] 
#>   eta3 =~ y31      0.8223      0.0311   26.3976    0.0000 [ 0.7578; 0.8614 ] 
#>   eta3 =~ y32      0.6581      0.0470   13.9940    0.0000 [ 0.5865; 0.7397 ] 
#>   eta3 =~ y33      0.7474      0.0419   17.8180    0.0000 [ 0.6756; 0.8154 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_percentile   
#>   Weight         Estimate  Std. error   t-stat.   p-value         95%        
#>   eta1 <~ y11      0.3956      0.0219   18.0864    0.0000 [ 0.3611; 0.4463 ] 
#>   eta1 <~ y12      0.3873      0.0192   20.1922    0.0000 [ 0.3558; 0.4259 ] 
#>   eta1 <~ y13      0.4542      0.0207   21.9521    0.0000 [ 0.4127; 0.4837 ] 
#>   eta2 <~ y21      0.3058      0.0301   10.1588    0.0000 [ 0.2361; 0.3497 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.6091    0.0000 [ 0.4122; 0.5014 ] 
#>   eta2 <~ y23      0.4735      0.0230   20.6155    0.0000 [ 0.4382; 0.5041 ] 
#>   eta3 <~ y31      0.4400      0.0189   23.2470    0.0000 [ 0.4035; 0.4701 ] 
#>   eta3 <~ y32      0.3521      0.0192   18.3658    0.0000 [ 0.3240; 0.3855 ] 
#>   eta3 <~ y33      0.3999      0.0232   17.2672    0.0000 [ 0.3609; 0.4392 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_percentile   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta2 ~ eta1       0.6713      0.0394   17.0535    0.0000 [ 0.6161; 0.7453 ] 
#>   eta3 ~ eta1       0.6634      0.0378   17.5317    0.0000 [ 0.6175; 0.7546 ] 
#>   eta3 ~ eta2       0.3052      0.0865    3.5298    0.0004 [ 0.1432; 0.4342 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_percentile   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         95%        
#>   eta3 ~ eta1          0.2049      0.0584    3.5096    0.0004 [ 0.1089; 0.2894 ] 
#> ________________________________________________________________________________

# Extract e.g. the loadings
res_summarize$Estimates$Loading_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat       p_value
#> 1 eta1 =~ y11  Common factor 0.6630699 0.04294698 15.439267  8.910744e-54
#> 2 eta1 =~ y12  Common factor 0.6492779 0.03196444 20.312508  9.968225e-92
#> 3 eta1 =~ y13  Common factor 0.7613458 0.03824780 19.905611  3.638284e-88
#> 4 eta2 =~ y21  Common factor 0.5164548 0.05986044  8.627647  6.262675e-18
#> 5 eta2 =~ y22  Common factor 0.7553877 0.03583081 21.082073  1.161827e-98
#> 6 eta2 =~ y23  Common factor 0.7996637 0.03912359 20.439425  7.461249e-93
#> 7 eta3 =~ y31  Common factor 0.8222773 0.03114971 26.397590 1.460315e-153
#> 8 eta3 =~ y32  Common factor 0.6580689 0.04702515 13.993977  1.696529e-44
#> 9 eta3 =~ y33  Common factor 0.7474241 0.04194769 17.818006  5.123008e-71
#>   CI_percentile.95%L CI_percentile.95%U
#> 1          0.5959098          0.7425953
#> 2          0.5948843          0.6995798
#> 3          0.6771736          0.8186403
#> 4          0.3904720          0.6079157
#> 5          0.6935968          0.8163282
#> 6          0.7462030          0.8907450
#> 7          0.7577687          0.8613535
#> 8          0.5865170          0.7397089
#> 9          0.6755771          0.8153546

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
#>  Random seed                        = -667070716
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
#>   eta2 ~ eta1      0.6713      0.0394   17.0535    0.0000 [ 0.5595; 0.7631 ] 
#>   eta3 ~ eta1      0.4585      0.0756    6.0617    0.0000 [ 0.2505; 0.6416 ] 
#>   eta3 ~ eta2      0.3052      0.0865    3.5298    0.0004 [ 0.0891; 0.5362 ] 
#> 
#> Estimated loadings:
#> ===================
#>                                                              CI_standard_t   
#>   Loading        Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 =~ y11      0.6631      0.0429   15.4393    0.0000 [ 0.5443; 0.7664 ] 
#>   eta1 =~ y12      0.6493      0.0320   20.3125    0.0000 [ 0.5675; 0.7328 ] 
#>   eta1 =~ y13      0.7613      0.0382   19.9056    0.0000 [ 0.6696; 0.8674 ] 
#>   eta2 =~ y21      0.5165      0.0599    8.6276    0.0000 [ 0.3569; 0.6664 ] 
#>   eta2 =~ y22      0.7554      0.0358   21.0821    0.0000 [ 0.6578; 0.8431 ] 
#>   eta2 =~ y23      0.7997      0.0391   20.4394    0.0000 [ 0.6947; 0.8970 ] 
#>   eta3 =~ y31      0.8223      0.0311   26.3976    0.0000 [ 0.7495; 0.9106 ] 
#>   eta3 =~ y32      0.6581      0.0470   13.9940    0.0000 [ 0.5354; 0.7786 ] 
#>   eta3 =~ y33      0.7474      0.0419   17.8180    0.0000 [ 0.6408; 0.8577 ] 
#> 
#> Estimated weights:
#> ==================
#>                                                              CI_standard_t   
#>   Weight         Estimate  Std. error   t-stat.   p-value         99%        
#>   eta1 <~ y11      0.3956      0.0219   18.0864    0.0000 [ 0.3345; 0.4476 ] 
#>   eta1 <~ y12      0.3873      0.0192   20.1922    0.0000 [ 0.3382; 0.4374 ] 
#>   eta1 <~ y13      0.4542      0.0207   21.9521    0.0000 [ 0.4050; 0.5120 ] 
#>   eta2 <~ y21      0.3058      0.0301   10.1588    0.0000 [ 0.2289; 0.3846 ] 
#>   eta2 <~ y22      0.4473      0.0240   18.6091    0.0000 [ 0.3867; 0.5110 ] 
#>   eta2 <~ y23      0.4735      0.0230   20.6155    0.0000 [ 0.4168; 0.5356 ] 
#>   eta3 <~ y31      0.4400      0.0189   23.2470    0.0000 [ 0.3926; 0.4905 ] 
#>   eta3 <~ y32      0.3521      0.0192   18.3658    0.0000 [ 0.3004; 0.3995 ] 
#>   eta3 <~ y33      0.3999      0.0232   17.2672    0.0000 [ 0.3387; 0.4585 ] 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>                                                               CI_standard_t   
#>   Total effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta2 ~ eta1       0.6713      0.0394   17.0535    0.0000 [ 0.5595; 0.7631 ] 
#>   eta3 ~ eta1       0.6634      0.0378   17.5317    0.0000 [ 0.5556; 0.7513 ] 
#>   eta3 ~ eta2       0.3052      0.0865    3.5298    0.0004 [ 0.0891; 0.5362 ] 
#> 
#> Estimated indirect effects:
#> ===========================
#>                                                                  CI_standard_t   
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value         99%        
#>   eta3 ~ eta1          0.2049      0.0584    3.5096    0.0004 [ 0.0565; 0.3583 ] 
#> ________________________________________________________________________________

# Extract the loading including both confidence intervals
res_summarize$Estimates$Path_estimates
#>          Name Construct_type  Estimate    Std_err    t_stat      p_value
#> 1 eta2 ~ eta1  Common factor 0.6713334 0.03936628 17.053514 3.291361e-65
#> 2 eta3 ~ eta1  Common factor 0.4585068 0.07563986  6.061709 1.346830e-09
#> 3 eta3 ~ eta2  Common factor 0.3051511 0.08645110  3.529754 4.159458e-04
#>   CI_standard_t.99%L CI_standard_t.99%U CI_standard_t.95%L CI_standard_t.95%U
#> 1          0.5595076          0.7630878          0.5839536          0.7386418
#> 2          0.2504572          0.6416238          0.2974287          0.5946524
#> 3          0.0891435          0.5362198          0.1428286          0.4825346
#>   CI_percentile.99%L CI_percentile.99%U CI_percentile.95%L CI_percentile.95%U
#> 1          0.5908370          0.7586808          0.6161396          0.7453393
#> 2          0.3366046          0.6475232          0.3475736          0.6313516
#> 3          0.1102884          0.4404683          0.1432178          0.4342156
```
