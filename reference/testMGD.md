# Tests for multi-group comparisons

**\[stable\]**

## Usage

``` r
testMGD(
 .object                = NULL,
 .alpha                 = 0.05,
 .approach_p_adjust     = "none",
 .approach_mgd          = c("all", "Klesel", "Chin", "Sarstedt", 
                            "Keil", "Nitzl", "Henseler", "CI_para","CI_overlap"),
 .output_type           = c("complete", "structured"),
 .parameters_to_compare = NULL,
 .eval_plan             = c("sequential", "multicore", "multisession"),                           
 .handle_inadmissibles  = c("replace", "drop", "ignore"),
 .R_permutation         = 499,
 .R_bootstrap           = 499,
 .saturated             = FALSE,
 .seed                  = NULL,
 .type_ci               = "CI_percentile",
 .type_vcv              = c("indicator", "construct"),
 .verbose               = TRUE
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

- .approach_p_adjust:

  Character string or a vector of character strings. Approach used to
  adjust the p-value for multiple testing. See the `methods` argument of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for a
  list of choices and their description. Defaults to "*none*".

- .approach_mgd:

  Character string or a vector of character strings. Approach used for
  the multi-group comparison. One of: "*all*", "*Klesel*", "*Chin*",
  "*Sarstedt*", "*Keil*, "*Nitzl*", "*Henseler*", "*CI_para*", or
  "*CI_overlap*". Default to "*all*" in which case all approaches are
  computed (if possible).

- .output_type:

  Character string. The type of output to return. One of "*complete*" or
  "*structured*". See the Value section for details. Defaults to
  "*complete*".

- .parameters_to_compare:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) indicating
  which parameters (i.e, path (`~`), loadings (`=~`), weights (`<~`), or
  correlations (`~~`)) should be compared across groups. Defaults to
  `NULL` in which case all weights, loadings and path coefficients of
  the originally specified model are compared.

- .eval_plan:

  Character string. The evaluation plan to use. One of "*sequential*",
  "*multicore*", or "*multisession*". In the two latter cases all
  available cores will be used. Defaults to "*sequential*".

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*drop*", "*ignore*", or "*replace*". If "*drop*", all
  replications/resamples yielding an inadmissible result will be dropped
  (i.e. the number of results returned will potentially be less than
  `.R`). For "*ignore*" all results are returned even if all or some of
  the replications yielded inadmissible results (i.e. number of results
  returned is equal to `.R`). For "*replace*" resampling continues until
  there are exactly `.R` admissible solutions. Defaults to "*replace*"
  to accommodate all approaches.

- .R_permutation:

  Integer. The number of permutations. Defaults to `499`

- .R_bootstrap:

  Integer. The number of bootstrap runs. Ignored if `.object` contains
  resamples. Defaults to `499`

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .type_ci:

  Character string. Which confidence interval should be calculated? For
  possible choices, see the `.quantity` argument of the
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  function. Only used if `.approch_mgd` is one of "*CI_para*" or
  "*CI_overlap*". Ignored otherwise. Defaults to "*CI_percentile*".

- .type_vcv:

  Character string. Which model-implied correlation matrix should be
  calculated? One of "*indicator*" or "*construct*". Defaults to
  "*indicator*".

- .verbose:

  Logical. Should information (e.g., progress bar) be printed to the
  console? Defaults to `TRUE`.

## Value

If `.output_type = "complete"` a list of class `cSEMTestMGD`.
Technically, `cSEMTestMGD` is a named list containing the following list
elements:

- `$Information`:

  Additional information.

- `$Klesel`:

  A list with elements, `Test_statistic`, `P_value`, and `Decision`

- `$Chin`:

  A list with elements, `Test_statistic`, `P_value`, `Decision`, and
  `Decision_overall`

- `$Sarstedt`:

  A list with elements, `Test_statistic`, `P_value`, `Decision`, and
  `Decision_overall`

- `$Keil`:

  A list with elements, `Test_statistic`, `P_value`, `Decision`, and
  `Decision_overall`

- `$Nitzl`:

  A list with elements, `Test_statistic`, `P_value`, `Decision`, and
  `Decision_overall`

- `$Henseler`:

  A list with elements, `Test_statistic`, `P_value`, `Decision`, and
  `Decision_overall`

- `$CI_para`:

  A list with elements, `Decision`, and `Decision_overall`

- `$CI_overlap`:

  A list with elements, `Decision`, and `Decision_overall`

If `.output_type = "structured"` a tibble (data frame) with the
following columns is returned.

- `Test`:

  The name of the test.

- `Comparision`:

  The parameter that was compared across groups. If "overall" the
  overall fit of the model was compared.

- `alpha%`:

  The test decision for a given "alpha" level. If `TRUE` the null
  hypotheses was rejected; if FALSE it was not rejected.

- `p-value_correction`:

  The p-value correction.

- `CI_type`:

  Only for the "CI_para" and the "CI_overlap" test. Which confidence
  interval was used.

- `Distance_metric`:

  Only for Test = "Klesel". Which distance metric was used.

## Details

This function performs various tests proposed in the context of
multigroup analysis.

The following tests are implemented:

- `.approach_mgd = "Klesel"`: Approach suggested by Klesel et al. (2019)
  :

  The model-implied variance-covariance matrix (either indicator
  (`.type_vcv = "indicator"`) or construct (`.type_vcv = "construct"`))
  is compared across groups. If the model-implied indicator or construct
  correlation matrix based on a saturated structural model should be
  compared, set `.saturated = TRUE`. To measure the distance between the
  model-implied variance-covariance matrices, the geodesic distance (dG)
  and the squared Euclidean distance (dL) are used. If more than two
  groups are compared, the average distance over all groups is used.

- `.approach_mgd = "Sarstedt"`: Approach suggested by Sarstedt et
  al. (2011) :

  Groups are compared in terms of parameter differences across groups.
  Sarstedt et al. (2011) tests if parameter k is equal across all
  groups. If several parameters are tested simultaneously it is
  recommended to adjust the significance level or the p-values (in cSEM
  correction is done by p-value). By default no multiple testing
  correction is done, however, several common adjustments are available
  via `.approach_p_adjust`. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  details. Note: the test has some severe shortcomings. Use with
  caution.

- `.approach_mgd = "Chin"`: Approach suggested by Chin and
  Dibbern (2010) :

  Groups are compared in terms of parameter differences across groups.
  Chin and Dibbern (2010) tests if parameter k is equal between two
  groups. If more than two groups are tested for equality, parameter k
  is compared between all pairs of groups. In this case, it is
  recommended to adjust the significance level or the p-values (in cSEM
  correction is done by p-value) since this is essentially a multiple
  testing setup. If several parameters are tested simultaneously,
  correction is by group and number of parameters. By default no
  multiple testing correction is done, however, several common
  adjustments are available via `.approach_p_adjust`. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  details.

- `.approach_mgd = "Keil"`: Approach suggested by Keil et al. (2000) :

  Groups are compared in terms of parameter differences across groups.
  Keil et al. (2000) tests if parameter k is equal between two groups.
  It is assumed, that the standard errors of the coefficients are equal
  across groups. The calculation of the standard error of the parameter
  difference is adjusted as proposed by Henseler et al. (2009) . If more
  than two groups are tested for equality, parameter k is compared
  between all pairs of groups. In this case, it is recommended to adjust
  the significance level or the p-values (in cSEM correction is done by
  p-value) since this is essentially a multiple testing setup. If
  several parameters are tested simultaneously, correction is by group
  and number of parameters. By default no multiple testing correction is
  done, however, several common adjustments are available via
  `.approach_p_adjust`. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  details.

- `.approach_mgd = "Nitzl"`: Approach suggested by Nitzl (2010) :

  Groups are compared in terms of parameter differences across groups.
  Similarly to Keil et al. (2000) , a single parameter k is tested for
  equality between two groups. In contrast to Keil et al. (2000) , it is
  assumed, that the standard errors of the coefficients are unequal
  across groups (Sarstedt et al. 2011) . If more than two groups are
  tested for equality, parameter k is compared between all pairs of
  groups. In this case, it is recommended to adjust the significance
  level or the p-values (in cSEM correction is done by p-value) since
  this is essentially a multiple testing setup. If several parameters
  are tested simultaneously, correction is by group and number of
  parameters. By default no multiple testing correction is done,
  however, several common adjustments are available via
  `.approach_p_adjust`. See
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
  details.

- `.approach_mgd = "Henseler"`: Approach suggested by Henseler (2007) :

  Groups are compared in terms of parameter differences across groups.
  In doing so, the bootstrap estimates of one parameter are compared
  across groups. In the literature, this approach is also known as
  PLS-MGA. Originally, this test was proposed as an one-sided test. In
  this function we perform a left-sided and a right-sided test to
  investigate whether a parameter differs across two groups. In doing
  so, the significance level is divided by 2 and compared to p-value of
  the left and right-sided test. Moreover, `.approach_p_adjust` is
  ignored and no overall decision is returned. For a more detailed
  description, see also Henseler et al. (2009) .

- `.approach_mgd = "CI_param"`: Approach mentioned in Sarstedt et
  al. (2011) :

  This approach is based on the confidence intervals constructed around
  the parameter estimates of the two groups. If the parameter of one
  group falls within the confidence interval of the other group and/or
  vice versa, it can be concluded that there is no group difference.
  Since it is based on the confidence intervals `.approach_p_adjust` is
  ignored.

- `.approach_mgd = "CI_overlap"`: Approach mentioned in Sarstedt et
  al. (2011) :

  This approach is based on the confidence intervals (CIs) constructed
  around the parameter estimates of the two groups. If the two CIs
  overlap, it can be concluded that there is no group difference. Since
  it is based on the confidence intervals `.approach_p_adjust` is
  ignored.

Use `.approach_mgd` to choose the approach. By default all approaches
are computed (`.approach_mgd = "all"`).

For convenience, two types of output are available. See the "Value"
section below.

By default, approaches based on parameter differences across groups
compare all parameters (`.parameters_to_compare = NULL`). To compare
only a subset of parameters provide the parameters in [lavaan model
syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) just like the
model to estimate. Take the simple model:

    model_to_estimate <- "
    Structural model
    eta2 ~ eta1
    eta3 ~ eta1 + eta2

    # Each concept os measured by 3 indicators, i.e., modeled as latent variable
    eta1 =~ y11 + y12 + y13
    eta2 =~ y21 + y22 + y23
    eta3 =~ y31 + y32 + y33
    "

If only the path from eta1 to eta3 and the loadings of eta1 are to be
compared across groups, write:

    to_compare <- "
    Structural parameters to compare
    eta3 ~ eta1

    # Loadings to compare
    eta1 =~ y11 + y12 + y13
    "

Note that the "model" provided to `.parameters_to_compare` does not need
to be an estimable model!

Note also that compared to all other functions in cSEM using the
argument, `.handle_inadmissibles` defaults to `"replace"` to accommodate
the Sarstedt et al. (2011) approach.

Argument `.R_permuation` is ignored for the `"Nitzl"` and the `"Keil"`
approach. `.R_bootstrap` is ignored if `.object` already contains
resamples, i.e. has class `cSEMResults_resampled` and if only the
`"Klesel"` or the `"Chin"` approach are used.

The argument `.saturated` is used by `"Klesel"` only. If
`.saturated = TRUE` the original structural model is ignored and
replaced by a saturated model, i.e. a model in which all constructs are
allowed to correlate freely. This is useful to test differences in the
measurement models between groups in isolation.

## References

Chin WW, Dibbern J (2010). “An Introduction to a Permutation Based
Procedure for Multi-Group PLS Analysis: Results of Tests of Differences
on Simulated Data and a Cross Cultural Analysis of the Sourcing of
Information System Services Between Germany and the USA.” In *Handbook
of Partial Least Squares*, 171–193. Springer Berlin Heidelberg.
[doi:10.1007/978-3-540-32827-8_8](https://doi.org/10.1007/978-3-540-32827-8_8)
.  
  
Henseler J (2007). “A new and simple approach to multi-group analysis in
partial least squares path modeling.” In Martens H, Næ s T (eds.),
*Proceedings of PLS'07 - The 5th International Symposium on PLS and
Related Methods*, 104–107. PLS, Norway: Matforsk, As.  
  
Henseler J, Ringle CM, Sinkovics RR (2009). “The use of partial least
squares path modeling in international marketing.” *Advances in
International Marketing*, **20**, 277–320.
[doi:10.1108/S1474-7979(2009)0000020014](https://doi.org/10.1108/S1474-7979%282009%290000020014)
.  
  
Keil M, Tan BC, Wei K, Saarinen T, Tuunainen V, Wassenaar A (2000). “A
cross-cultural study on escalation of commitment behavior in software
projects.” *MIS Quarterly*, **24**(2), 299–325.  
  
Klesel M, Schuberth F, Henseler J, Niehaves B (2019). “A Test for
Multigroup Comparison Using Partial Least Squares Path Modeling.”
*Internet Research*, **29**(3), 464–477.
[doi:10.1108/intr-11-2017-0418](https://doi.org/10.1108/intr-11-2017-0418)
.  
  
Nitzl C (2010). “Eine anwenderorientierte Einfuehrung in die Partial
Least Square (PLS)-Methode.” In *Arbeitspapier*, number 21. Universitaet
Hamburg, Institut fuer Industrielles Management, Hamburg.  
  
Sarstedt M, Henseler J, Ringle CM (2011). “Multigroup Analysis in
Partial Least Squares (PLS) Path Modeling: Alternative Methods and
Empirical Results.” In *Advances in International Marketing*, 195–218.
Emerald Group Publishing Limited.
[doi:10.1108/s1474-7979(2011)0000022012](https://doi.org/10.1108/s1474-7979%282011%290000022012)
.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md),
[`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ===========================================================================
# Basic usage
# ===========================================================================
model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model

EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  <~ sat1  + sat2  + sat3  + sat4
VAL  <~ val1  + val2  + val3  + val4
"

## Create list of virtually identical data sets
dat <- list(satisfaction[-3,], satisfaction[-5, ], satisfaction[-10, ])
out <- csem(dat, model, .resample_method = "bootstrap", .R = 40) 

## Test 
testMGD(out, .R_permutation = 40,.verbose = FALSE)

# Notes: 
#  1. .R_permutation (and .R in the call to csem) is small to make examples run quicker; 
#     should be higher in real applications.
#  2. Test will not reject their respective H0s since the groups are virtually
#     identical.
#  3. Only exception is the approach suggested by Sarstedt et al. (2011), a
#     sign that the test is unreliable.
#  4. As opposed to other functions involving the argument, 
#     '.handle_inadmissibles' the default is "replace" as this is
#     required by Sarstedt et al. (2011)'s approach.

# ===========================================================================
# Extended usage
# ===========================================================================
### Test only a subset ------------------------------------------------------
# By default all parameters are compared. Select a subset by providing a 
# model in lavaan model syntax:

to_compare <- "
# Path coefficients
QUAL ~ EXPE

# Loadings
EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
"

## Test 
testMGD(out, .parameters_to_compare = to_compare, .R_permutation = 20, 
        .R_bootstrap = 20, .verbose = FALSE)

### Different p_adjustments --------------------------------------------------
# To adjust p-values to accommodate multiple testing use .approach_p_adjust. 
# The number of tests to use for adjusting depends on the approach chosen. For
# the Chin approach for example it is the number of parameters to test times the
# number of possible group comparisons. To compare the results for different
# adjustments, a vector of p-adjustments may be chosen.

## Test 
testMGD(out, .parameters_to_compare = to_compare, 
        .approach_p_adjust = c("none", "bonferroni"),
        .R_permutation = 20, .R_bootstrap = 20, .verbose = FALSE)
} # }
```
