# Assess model

**\[maturing\]**

## Usage

``` r
assess(
  .object              = NULL, 
  .quality_criterion   = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                           "hqc", "mallows_cp", "ave",
                           "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "dg", "dl", "dml", "df",
                           "effects", "f2", "fl_criterion", "chi_square", "chi_square_df",
                           "cfi", "cn", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "htmt2", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB"),
  .only_common_factors = TRUE, 
  ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .quality_criterion:

  Character string. A single character string or a vector of character
  strings naming the quality criterion to compute. See the Details
  section for a list of possible candidates. Defaults to "*all*" in
  which case all possible quality criteria are computed.

- .only_common_factors:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- ...:

  Further arguments passed to functions called by `assess()`. See
  [args_assess_dotdotdot](https://floschuberth.github.io/cSEM/reference/args_assess_dotdotdot.md)
  for a complete list of available arguments.

## Value

A named list of quality criteria. Note that if only a single quality
criteria is computed the return value is still a list!

## Details

Assess a model using common quality criteria. See the [Postestimation:
Assessing a
model](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
article on the [cSEM
website](https://floschuberth.github.io/cSEM/index.html) for details.

The function is essentially a wrapper around a number of internal
functions that perform an "assessment task" (called a **quality
criterion** in cSEM parlance) like computing reliability estimates, the
effect size (Cohen's f^2), the heterotrait-monotrait ratio of
correlations (HTMT) etc.

By default every possible quality criterion is calculated
(`.quality_criterion = "all"`). If only a subset of quality criteria are
needed a single character string or a vector of character strings naming
the criteria to be computed may be supplied to `assess()` via the
`.quality_criterion` argument. Currently, the following quality criteria
are implemented (in alphabetical order):

- Average variance extracted (AVE); "ave":

  An estimate of the amount of variation in the indicators that is due
  to the underlying latent variable. Practically, it is calculated as
  the ratio of the (indicator) true score variances (i.e., the sum of
  the squared loadings) relative to the sum of the total indicator
  variances. The AVE is inherently tied to the common factor model. It
  is therefore unclear how to meaningfully interpret AVE results for
  constructs modeled as composites. It is possible to report the AVE for
  constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning. Calculation is
  done by
  [`calculateAVE()`](https://floschuberth.github.io/cSEM/reference/calculateAVE.md).

- Congeneric reliability; "rho_C", "rho_C_mm", "rho_C_weighted",
  "rho_C_weighted_mm":

  An estimate of the reliability assuming a congeneric measurement model
  (i.e., loadings are allowed to differ) and a test score (proxy) based
  on unit weights. There are four different versions implemented. See
  the [Methods and
  Formulae](https://floschuberth.github.io/cSEM/articles/Using-assess.html#methods)
  section of the [Postestimation: Assessing a
  model](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
  article on the [cSEM
  website](https://floschuberth.github.io/cSEM/index.html) for details.
  Alternative but synonymous names for `"rho_C"` are: composite
  reliability, construct reliability, reliability coefficient,
  Jöreskog's rho, coefficient omega, or Dillon-Goldstein's rho. For
  `"rho_C_weighted"`: (Dijkstra-Henselers) rhoA. `rho_C_mm` and
  `rho_C_weighted_mm` have no corresponding names. The former uses unit
  weights scaled by (w'Sw)^(-1/2) and the latter weights scaled by
  (w'Sigma_hat w)^(-1/2) where Sigma_hat is the model-implied indicator
  correlation matrix. The Congeneric reliability is inherently tied to
  the common factor model. It is therefore unclear how to meaningfully
  interpret congeneric reliability estimates for constructs modeled as
  composites. It is possible to report the congeneric reliability for
  constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning. Calculation is
  done by
  [`calculateRhoC()`](https://floschuberth.github.io/cSEM/reference/reliability.md).

- Distance measures; "dg", "dl", "dml":

  Measures of the distance between the model-implied and the empirical
  indicator correlation matrix. Currently, the geodesic distance
  (`"dg"`), the squared Euclidean distance (`"dl"`) and the the maximum
  likelihood-based distance function are implemented (`"dml"`).
  Calculation is done by
  [`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  [`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  and
  [`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md).

- Degrees of freedom, "df":

  Returns the degrees of freedom. Calculation is done by
  [`calculateDf()`](https://floschuberth.github.io/cSEM/reference/calculateDf.md).

- Effects; "effects":

  Total and indirect effect estimates. Additionally, the variance
  accounted for (VAF) is computed. The VAF is defined as the ratio of a
  variables indirect effect to its total effect. Calculation is done by
  [`calculateEffects()`](https://floschuberth.github.io/cSEM/reference/calculateEffects.md).

- Effect size; "f2":

  An index of the effect size of an independent variable in a structural
  regression equation. This measure is commonly known as Cohen's f^2.
  The effect size of the k'th independent variable in this case is
  defined as the ratio (R2_included - R2_excluded)/(1 - R2_included),
  where R2_included and R2_excluded are the R squares of the original
  structural model regression equation (R2_included) and the alternative
  specification with the k'th variable dropped (R2_excluded).
  Calculation is done by
  [`calculatef2()`](https://floschuberth.github.io/cSEM/reference/calculatef2.md).

- Fit indices; "chi_square", "chi_square_df", "cfi", "cn", "gfi", "ifi",
  "nfi", "nnfi", "rmsea", "rms_theta", "srmr":

  Several absolute and incremental fit indices. Note that their
  suitability for models containing constructs modeled as composites is
  still an open research question. Also note that fit indices are not
  tests in a hypothesis testing sense and decisions based on common
  one-size-fits-all cut-offs proposed in the literature suffer from
  serious statistical drawbacks. Calculation is done by
  [`calculateChiSquare()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateChiSquareDf()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateCFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateGFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateIFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateNNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateRMSEA()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateRMSTheta()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  and
  [`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md).

- Fornell-Larcker criterion; "fl_criterion":

  A rule suggested by Fornell and Larcker (1981) to assess discriminant
  validity. The Fornell-Larcker criterion is a decision rule based on a
  comparison between the squared construct correlations and the average
  variance extracted. FL returns a matrix with the squared construct
  correlations on the off-diagonal and the AVEs on the main diagonal.
  Calculation is done by
  [`calculateFLCriterion()`](https://floschuberth.github.io/cSEM/reference/calculateFLCriterion.md).

- Goodness of Fit (GoF); "gof":

  The GoF is defined as the square root of the mean of the R squares of
  the structural model times the mean of the variances in the indicators
  that are explained by their related constructs (i.e., the average over
  all lambda^2_k). For the latter, only constructs modeled as common
  factors are considered as they explain their indicator variance in
  contrast to a composite where indicators actually build the construct.
  Note that, contrary to what the name suggests, the GoF is **not** a
  measure of model fit in a Chi-square fit test sense. Calculation is
  done by
  [`calculateGoF()`](https://floschuberth.github.io/cSEM/reference/calculateGoF.md).

- Heterotrait-monotrait ratio of correlations (HTMT); "htmt":

  An estimate of the correlation between latent variables assuming tau
  equivalent measurement models. The HTMT is used to assess convergent
  and/or discriminant validity of a construct. The HTMT is inherently
  tied to the common factor model. If the model contains less than two
  constructs modeled as common factors and
  `.only_common_factors = TRUE`, `NA` is returned. It is possible to
  report the HTMT for constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning. Calculation is
  done by
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md).

- HTMT2; "htmt2":

  An estimate of the correlation between latent variables assuming
  congeneric measurement models. The HTMT2 is used to assess convergent
  and/or discriminant validity of a construct. The HTMT is inherently
  tied to the common factor model. If the model contains less than two
  constructs modeled as common factors and
  `.only_common_factors = TRUE`, `NA` is returned. It is possible to
  report the HTMT for constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning. Calculation is
  done by
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md).

- Model selection criteria: "aic", "aicc", "aicu", "bic", "fpe", "gm",
  "hq", "hqc", "mallows_cp":

  Several model selection criteria as suggested by Sharma et al. (2019)
  in the context of PLS. See:
  [`calculateModelSelectionCriteria()`](https://floschuberth.github.io/cSEM/reference/calculateModelSelectionCriteria.md)
  for details.

- Reliability: "reliability":

  As described in the [Methods and
  Formulae](https://floschuberth.github.io/cSEM/articles/Using-assess.html#methods)
  section of the [Postestimation: Assessing a
  model](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
  article on the [cSEM
  website](https://floschuberth.github.io/cSEM/index.html) there are
  many different estimators for the (internal consistency) reliability.
  Choosing `.quality_criterion = "reliability"` computes the three most
  common measures, namely: "Cronbach's alpha" (identical to "rho_T"),
  "Jöreskog's rho" (identical to "rho_C_mm"), and "Dijkstra-Henseler's
  rho A" (identical to "rho_C_weighted_mm"). Reliability is inherently
  tied to the common factor model. It is therefore unclear how to
  meaningfully interpret reliability estimates for constructs modeled as
  composites. It is possible to report the three common reliability
  estimates for constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning.

- R square and R square adjusted; "r2", "r2_adj":

  The R square and the adjusted R square for each structural regression
  equation. Calculated when running
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- Tau-equivalent reliability; "rho_T":

  An estimate of the reliability assuming a tau-equivalent measurement
  model (i.e. a measurement model with equal loadings) and a test score
  (proxy) based on unit weights. Tau-equivalent reliability is the
  preferred name for reliability estimates that assume a tau-equivalent
  measurement model such as Cronbach's alpha. The tau-equivalent
  reliability (Cronbach's alpha) is inherently tied to the common factor
  model. It is therefore unclear how to meaningfully interpret
  tau-equivalent reliability estimates for constructs modeled as
  composites. It is possible to report tau-equivalent reliability
  estimates for constructs modeled as composites by setting
  `.only_common_factors = FALSE`, however, result should be interpreted
  with caution as they may not have a conceptual meaning. Calculation is
  done by
  [`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md).

- Variance inflation factors (VIF); "vif":

  An index for the amount of (multi-)collinearity between independent
  variables of a regression equation. Computed for each structural
  equation. Practically, VIF_k is defined as the ratio of 1 over (1 -
  R2_k) where R2_k is the R squared from a regression of the k'th
  independent variable on all remaining independent variables.
  Calculated when running
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- Variance inflation factors for PLS-PM mode B (VIF-ModeB); "vifmodeB":

  An index for the amount of (multi-)collinearity between independent
  variables (indicators) in mode B regression equations. Computed only
  if `.object` was obtained using `.weight_approach = "PLS-PM"` and at
  least one mode was mode B. Practically, VIF-ModeB_k is defined as the
  ratio of 1 over (1 - R2_k) where R2_k is the R squared from a
  regression of the k'th indicator of block j on all remaining
  indicators of the same block. Calculation is done by
  [`calculateVIFModeB()`](https://floschuberth.github.io/cSEM/reference/calculateVIFModeB.md).

For details on the most important quality criteria see the [Methods and
Formulae](https://floschuberth.github.io/cSEM/articles/Using-assess.html#methods)
section of the [Postestimation: Assessing a
model](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
article on the on the [cSEM
website](https://floschuberth.github.io/cSEM/index.html).

Some of the quality criteria are inherently tied to the classical common
factor model and therefore only meaningfully interpreted within a common
factor model (see the [Postestimation: Assessing a
model](https://floschuberth.github.io/cSEM/articles/Using-assess.html)
article for details). It is possible to force computation of all quality
criteria for constructs modeled as composites by setting
`.only_common_factors = FALSE`, however, we explicitly warn to interpret
quality criteria in analogy to the common factor model in this case, as
the interpretation often does not carry over to composite models.

### Resampling

To resample a given quality criterion supply the name of the function
that calculates the desired quality criterion to
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)'s
`.user_funs` argument. See
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
for details.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md),
[`exportToExcel()`](https://floschuberth.github.io/cSEM/reference/exportToExcel.md)

## Examples

``` r
# ===========================================================================
# Using the three common factors dataset
# ===========================================================================
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., modeled as latent variable
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res <- csem(threecommonfactors, model)
a   <- assess(res) # computes all quality criteria (.quality_criterion = "all")
a
#> ________________________________________________________________________________
#> 
#>  Construct        AVE           R2          R2_adj    
#>  eta1           0.4803          NA            NA      
#>  eta2           0.4923        0.4507        0.4496    
#>  eta3           0.5559        0.4912        0.4892    
#> 
#> -------------- Common (internal consistency) reliability estimates -------------
#> 
#>  Construct Cronbachs_alpha   Joereskogs_rho   Dijkstra-Henselers_rho_A 
#>  eta1        0.7318           0.7339                0.7388          
#>  eta2        0.7281           0.7380                0.7647          
#>  eta3        0.7860           0.7884                0.7964          
#> 
#> ----------- Alternative (internal consistency) reliability estimates -----------
#> 
#>  Construct       RhoC         RhoC_mm    RhoC_weighted
#>  eta1           0.7339        0.7341        0.7388    
#>  eta2           0.7380        0.7361        0.7647    
#>  eta3           0.7884        0.7875        0.7964    
#> 
#>  Construct  RhoC_weighted_mm     RhoT      RhoT_weighted
#>  eta1           0.7388        0.7318        0.7288    
#>  eta2           0.7647        0.7281        0.7095    
#>  eta3           0.7964        0.7860        0.7820    
#> 
#> --------------------------- Distance and fit measures --------------------------
#> 
#>  Geodesic distance             = 0.006013595
#>  Squared Euclidean distance    = 0.01121567
#>  ML distance                   = 0.03203348
#> 
#>  Chi_square       = 15.9847
#>  Chi_square_df    = 0.6660294
#>  CFI              = 1
#>  CN               = 1137.78
#>  GFI              = 0.9920803
#>  IFI              = 1.005614
#>  NFI              = 0.9889886
#>  NNFI             = 1
#>  RMSEA            = 0
#>  RMS_theta        = 0.1050618
#>  SRMR             = 0.01578725
#> 
#>  Degrees of freedom       = 24
#> 
#> --------------------------- Model selection criteria ---------------------------
#> 
#>  Construct        AIC          AICc          AICu     
#>  eta2          -296.5459     205.5025      -294.5419  
#>  eta3          -332.8544     169.2264      -329.8454  
#> 
#>  Construct        BIC           FPE           GM      
#>  eta2          -288.1166      0.5526       511.4292   
#>  eta3          -320.2106      0.5139       517.6438   
#> 
#>  Construct        HQ            HQc       Mallows_Cp  
#>  eta2          -293.2383     -293.1793      3.0000    
#>  eta3          -327.8930     -327.7823      5.0000    
#> 
#> ----------------------- Variance inflation factors (VIFs) ----------------------
#> 
#>   Dependent construct: 'eta3'
#> 
#>  Independent construct    VIF value 
#>  eta1                      1.8205   
#>  eta2                      1.8205   
#> 
#> -------------------------- Effect sizes (Cohen's f^2) --------------------------
#> 
#>   Dependent construct: 'eta2'
#> 
#>  Independent construct       f^2    
#>  eta1                      0.8205   
#> 
#>   Dependent construct: 'eta3'
#> 
#>  Independent construct       f^2    
#>  eta1                      0.2270   
#>  eta2                      0.1005   
#> 
#> ----------------------- Discriminant validity assessment -----------------------
#> 
#>  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
#> 
#>           eta1      eta2 eta3
#> eta1 1.0000000 0.0000000    0
#> eta2 0.6782752 1.0000000    0
#> eta3 0.6668841 0.6124418    1
#> 
#> 
#>  Advanced heterotrait-monotrait ratio of correlations matrix (HTMT2 matrix)
#> 
#>           eta1      eta2 eta3
#> eta1 1.0000000 0.0000000    0
#> eta2 0.6724003 1.0000000    0
#> eta3 0.6652760 0.5958725    1
#> 
#> 
#>  Fornell-Larcker matrix
#> 
#>           eta1      eta2      eta3
#> eta1 0.4802903 0.4506886 0.4400530
#> eta2 0.4506886 0.4922660 0.3757225
#> eta3 0.4400530 0.3757225 0.5559458
#> 
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>   Total effect    Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1       0.6713          NA        NA        NA
#>   eta3 ~ eta1       0.6634          NA        NA        NA
#>   eta3 ~ eta2       0.3052          NA        NA        NA
#> 
#> Estimated indirect effects:
#> ===========================
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value
#>   eta3 ~ eta1          0.2049          NA        NA        NA
#> ________________________________________________________________________________

## The return value is a named list. Type for example:
a$HTMT
#> $htmts
#>           eta1      eta2 eta3
#> eta1 1.0000000 0.0000000    0
#> eta2 0.6782752 1.0000000    0
#> eta3 0.6668841 0.6124418    1
#> 
#> $quantiles
#> NULL
#> 
#> $nr_admissibles
#> NULL
#> 

# You may also just compute a subset of the quality criteria
assess(res, .quality_criterion = c("ave", "rho_C", "htmt"))
#> ________________________________________________________________________________
#> 
#>  Construct        AVE     
#>  eta1           0.4803    
#>  eta2           0.4923    
#>  eta3           0.5559    
#> 
#> ----------- Alternative (internal consistency) reliability estimates -----------
#> 
#>  Construct       RhoC     
#>  eta1           0.7339    
#>  eta2           0.7380    
#>  eta3           0.7884    
#> 
#> ----------------------- Discriminant validity assessment -----------------------
#> 
#>  Heterotrait-monotrait ratio of correlations matrix (HTMT matrix)
#> 
#>           eta1      eta2 eta3
#> eta1 1.0000000 0.0000000    0
#> eta2 0.6782752 1.0000000    0
#> eta3 0.6668841 0.6124418    1
#> 
#> ________________________________________________________________________________

## Resampling ---------------------------------------------------------------
# To resample a given quality criterion use csem()'s .user_funs argument
# Note: The output of the quality criterion needs to be a vector or a matrix.
#       Matrices will be vectorized columnwise.
res <- csem(threecommonfactors, model, 
            .resample_method = "bootstrap", 
            .R               = 40,
            .user_funs       = cSEM:::calculateSRMR
)

## Look at the resamples
res$Estimates$Estimates_resample$Estimates1$User_fun$Resampled[1:4, ]
#> [1] 0.02459846 0.02738783 0.02846445 0.02611109

## Use infer() to compute e.g., the 95% percentile confidence interval
res_infer <- infer(res, .quantity = "CI_percentile")

## The results are saved under the name "User_fun"
res_infer$User_fun 
#> $CI_percentile
#>            [,1]
#> 95%L 0.01823702
#> 95%U 0.03545628
#> 

## Several quality criteria can be resampled simultaneously
res <- csem(threecommonfactors, model, 
            .resample_method = "bootstrap",
            .R               = 40,
            .user_funs       = list(
              "SRMR" = cSEM:::calculateSRMR,
              "RMS_theta" = cSEM:::calculateRMSTheta
            ),
            .tolerance = 1e-04
)
res$Estimates$Estimates_resample$Estimates1$SRMR$Resampled[1:4, ]
#> [1] 0.02650658 0.02545304 0.02006635 0.03189040
res$Estimates$Estimates_resample$Estimates1$RMS_theta$Resampled[1:4]
#> [1] 0.10888203 0.10234856 0.10636301 0.09283051
```
