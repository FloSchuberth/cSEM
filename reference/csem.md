# Composite-based SEM

**\[stable\]**

## Usage

``` r
csem(
.data                  = NULL,
.model                 = NULL,
.approach_2ndorder     = c("2stage", "mixed"),
.approach_cor_robust   = c("none", "mcd", "spearman"),
.approach_nl           = c("sequential", "replace"),
.approach_paths        = c("OLS", "2SLS"),
.approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", 
                           "MINVAR", "GENVAR","GSCA", "PCA",
                           "unit", "bartlett", "regression"),
.conv_criterion        = c("diff_absolute", "diff_squared", "diff_relative"),
.disattenuate          = TRUE,
.dominant_indicators   = NULL,
.estimate_structural   = TRUE,
.id                    = NULL,
.instruments           = NULL,
.iter_max              = 100,
.normality             = FALSE,
.PLS_approach_cf       = c("dist_squared_euclid", "dist_euclid_weighted", 
                           "fisher_transformed", "mean_arithmetic",
                           "mean_geometric", "mean_harmonic",
                           "geo_of_harmonic"),
.PLS_ignore_structural_model = FALSE,
.PLS_modes                   = NULL,
.PLS_weight_scheme_inner     = c("path", "centroid", "factorial"),
.reliabilities         = NULL,
.starting_values       = NULL,
.resample_method       = c("none", "bootstrap", "jackknife"),
.resample_method2      = c("none", "bootstrap", "jackknife"),
.R                     = 499,
.R2                    = 199,
.handle_inadmissibles  = c("drop", "ignore", "replace"),
.user_funs             = NULL,
.eval_plan             = c("sequential", "multicore", "multisession"),
.seed                  = NULL,
.sign_change_option    = c("none", "individual", "individual_reestimate", 
                           "construct_reestimate"),
.tolerance             = 1e-05
)
```

## Arguments

- .data:

  A `data.frame` or a `matrix` of standardized or unstandardized data
  (indicators/items/manifest variables). Additionally, a `list` of data
  sets (data frames or matrices) is accepted in which case estimation is
  repeated for each data set. Possible column types or classes of the
  data provided are: "`logical`", "`numeric`" ("`double`" or
  "`integer`"), "`factor`" ("`ordered`" and/or "`unordered`"),
  "`character`" (will be converted to factor), or a mix of several
  types.

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) or a
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  list.

- .approach_2ndorder:

  Character string. Approach used for models containing second-order
  constructs. One of: "*2stage*", or "*mixed*". Defaults to "*2stage*".

- .approach_cor_robust:

  Character string. Approach used to obtain a robust indicator
  correlation matrix. One of: "*none*" in which case the standard
  Bravais-Pearson correlation is used, "*spearman*" for the Spearman
  rank correlation, or "*mcd*" via
  [`MASS::cov.rob()`](https://rdrr.io/pkg/MASS/man/cov.rob.html) for a
  robust correlation matrix. Defaults to "*none*". Note that many
  postestimation procedures (such as
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
  or [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md)
  implicitly assume a continuous indicator correlation matrix (e.g.
  Bravais-Pearson correlation matrix). Only use if you know what you are
  doing.

- .approach_nl:

  Character string. Approach used to estimate nonlinear structural
  relationships. One of: "*sequential*" or "*replace*". Defaults to
  "*sequential*".

- .approach_paths:

  Character string. Approach used to estimate the structural
  coefficients. One of: "*OLS*" or "*2SLS*". If "*2SLS*", instruments
  need to be supplied to `.instruments`. Defaults to "*OLS*".

- .approach_weights:

  Character string. Approach used to obtain composite weights. One of:
  "*PLS-PM*", "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*",
  "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", or
  "*regression*". Defaults to "*PLS-PM*".

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .disattenuate:

  Logical. Should composite/proxy correlations be disattenuated to yield
  consistent loadings and path estimates if at least one of the
  construct is modeled as a common factor? Defaults to `TRUE`.

- .dominant_indicators:

  A character vector of `"construct_name" = "indicator_name"` pairs,
  where `"indicator_name"` is a character string giving the name of the
  dominant indicator and `"construct_name"` a character string of the
  corresponding construct name. Dominant indicators may be specified for
  a subset of the constructs. Default to `NULL`.

- .estimate_structural:

  Logical. Should the structural coefficients be estimated? Defaults to
  `TRUE`.

- .id:

  Character string or integer. A character string giving the name or an
  integer of the position of the column of `.data` whose levels are used
  to split `.data` into groups. Defaults to `NULL`.

- .instruments:

  A named list of vectors of instruments. The names of the list elements
  are the names of the dependent (LHS) constructs of the structural
  equation whose explanatory variables are endogenous. The vectors
  contain the names of the instruments corresponding to each equation.
  Note that exogenous variables of a given equation **must** be supplied
  as instruments for themselves. Defaults to `NULL`.

- .iter_max:

  Integer. The maximum number of iterations allowed. If `iter_max = 1`
  and `.approach_weights = "PLS-PM"` one-step weights are returned. If
  the algorithm exceeds the specified number, weights of iteration step
  `.iter_max - 1` will be returned with a warning. Defaults to `100`.

- .normality:

  Logical. Should joint normality of \\\[\eta\_{1:p}; \zeta;
  \epsilon\]\\ be assumed in the nonlinear model? See (Dijkstra and
  Schermelleh-Engel 2014) for details. Defaults to `FALSE`. Ignored if
  the model is not nonlinear.

- .PLS_approach_cf:

  Character string. Approach used to obtain the correction factors for
  PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
  "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*",
  "*mean_harmonic*", "*geo_of_harmonic*". Defaults to
  "*dist_squared_euclid*". Ignored if `.disattenuate = FALSE` or if
  `.approach_weights` is not PLS-PM.

- .PLS_ignore_structural_model:

  Logical. Should the structural model be ignored when calculating the
  inner weights of the PLS-PM algorithm? Defaults to `FALSE`. Ignored if
  `.approach_weights` is not PLS-PM.

- .PLS_modes:

  Either a named list specifying the mode that should be used for each
  construct in the form `"construct_name" = mode`, a single character
  string giving the mode that should be used for all constructs, or
  `NULL`. Possible choices for `mode` are: "*modeA*", "*modeB*",
  "*modeBNNLS*", "*unit*", "*PCA*", a single integer or a vector of
  fixed weights of the same length as there are indicators for the
  construct given by `"construct_name"`. If only a single number is
  provided this is identical to using unit weights, as weights are
  rescaled such that the related composite has unit variance. Defaults
  to `NULL`. If `NULL` the appropriate mode according to the type of
  construct used is chosen. Ignored if `.approach_weight` is not PLS-PM.

- .PLS_weight_scheme_inner:

  Character string. The inner weighting scheme used by PLS-PM. One of:
  "*centroid*", "*factorial*", or "*path*". Defaults to "*path*".
  Ignored if `.approach_weight` is not PLS-PM.

- .reliabilities:

  A character vector of `"name" = value` pairs, where `value` is a
  number between 0 and 1 and `"name"` a character string of the
  corresponding construct name, or `NULL`. Reliabilities may be given
  for a subset of the constructs. Defaults to `NULL` in which case
  reliabilities are estimated by `csem()`. Currently, only supported for
  `.approach_weights = "PLS-PM"`.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

- .resample_method:

  Character string. The resampling method to use. One of: "*none*",
  "*bootstrap*" or "*jackknife*". Defaults to "*none*".

- .resample_method2:

  Character string. The resampling method to use when resampling from a
  resample. One of: "*none*", "*bootstrap*" or "*jackknife*". For
  "*bootstrap*" the number of draws is provided via `.R2`. Currently,
  resampling from each resample is only required for the studentized
  confidence interval ("*CI_t_interval*") computed by the
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  function. Defaults to "*none*".

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .R2:

  Integer. The number of bootstrap replications to use when resampling
  from a resample. Defaults to `199`.

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*drop*", "*ignore*", or "*replace*". If "*drop*", all
  replications/resamples yielding an inadmissible result will be dropped
  (i.e. the number of results returned will potentially be less than
  `.R`). For "*ignore*" all results are returned even if all or some of
  the replications yielded inadmissible results (i.e. number of results
  returned is equal to `.R`). For "*replace*" resampling continues until
  there are exactly `.R` admissible solutions. Depending on the
  frequency of inadmissible solutions this may significantly increase
  computing time. Defaults to "*drop*".

- .user_funs:

  A function or a (named) list of functions to apply to every resample.
  The functions must take `.object` as its first argument (e.g.,
  `myFun <- function(.object, ...) {body-of-the-function}`). Function
  output should preferably be a (named) vector but matrices are also
  accepted. However, the output will be vectorized (columnwise) in this
  case. See the examples section for details.

- .eval_plan:

  Character string. The evaluation plan to use. One of "*sequential*",
  "*multicore*", or "*multisession*". In the two latter cases all
  available cores will be used. Defaults to "*sequential*".

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .sign_change_option:

  Character string. Which sign change option should be used to handle
  flipping signs when resampling? One of "*none*","*individual*",
  "*individual_reestimate*", "*construct_reestimate*". Defaults to
  "*none*".

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

## Value

An object of class `cSEMResults` with methods for all postestimation
generics. Technically, a call to `csem()` results in an object with at
least two class attributes. The first class attribute is always
`cSEMResults`. The second is one of `cSEMResults_default`,
`cSEMResults_multi`, or `cSEMResults_2ndorder` and depends on the
estimated model and/or the type of data provided to the `.model` and
`.data` arguments. The third class attribute `cSEMResults_resampled` is
only added if resampling was conducted. For a details see the
[cSEMResults
helpfile](https://floschuberth.github.io/cSEM/reference/csem_results.md)
.

## Details

Estimate linear, nonlinear, hierarchical and multigroup structural
equation models using a composite-based approach. In cSEM any method or
approach that involves linear compounds (scores/proxies/composites) of
observables (indicators/items/manifest variables) is defined as
composite-based. See the [Get
started](https://floschuberth.github.io/cSEM/articles/cSEM.html) section
of the [cSEM website](https://floschuberth.github.io/cSEM/index.html)
for a general introduction to composite-based SEM and cSEM.

`csem()` estimates linear, nonlinear, hierarchical or multigroup
structural equation models using a composite-based approach.

### Data and model:

The `.data` and `.model` arguments are required. `.data` must be given a
`matrix` or a `data.frame` with column names matching the indicator
names used in the model description. Alternatively, a `list` of data
sets (matrices or data frames) may be provided in which case estimation
is repeated for each data set. Possible column types/classes of the data
provided are: "`logical`", "`numeric`" ("`double`" or "`integer`"),
"`factor`" ("`ordered`" and/or "`unordered`"), "`character`", or a mix
of several types. Character columns will be treated as (unordered)
factors.

Depending on the type/class of the indicator data provided cSEM computes
the indicator correlation matrix in different ways. See
[`calculateIndicatorCor()`](https://floschuberth.github.io/cSEM/reference/calculateIndicatorCor.md)
for details.

In the current version `.data` must not contain missing values. Future
versions are likely to handle missing values as well.

To provide a model use the [lavaan model
syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html). Note,
however, that cSEM currently only supports the "standard" lavaan model
syntax (Types 1, 2, 3, and 7 as described on the help page). Therefore,
specifying e.g., a threshold or scaling factors is ignored.
Alternatively, a standardized (possibly incomplete)
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list
may be supplied. See
[`parseModel()`](https://floschuberth.github.io/cSEM/reference/parseModel.md)
for details.

### Weights and path coefficients:

By default weights are estimated using the partial least squares path
modeling algorithm (`"PLS-PM"`). A range of alternative weighting
algorithms may be supplied to `.approach_weights`. Currently, the
following approaches are implemented

1.  (Default) Partial least squares path modeling (`"PLS-PM"`). The
    algorithm can be customized. See
    [`calculateWeightsPLS()`](https://floschuberth.github.io/cSEM/reference/calculateWeightsPLS.md)
    for details.

2.  Generalized structured component analysis (`"GSCA"`) and generalized
    structured component analysis with uniqueness terms (GSCAm). The
    algorithms can be customized. See
    [`calculateWeightsGSCA()`](https://floschuberth.github.io/cSEM/reference/calculateWeightsGSCA.md)
    and
    [`calculateWeightsGSCAm()`](https://floschuberth.github.io/cSEM/reference/calculateWeightsGSCAm.md)
    for details. Note that GSCAm is called indirectly when the model
    contains constructs modeled as common factors only and
    `.disattenuate = TRUE`. See below.

3.  Generalized canonical correlation analysis (*GCCA*), including
    `"SUMCORR"`, `"MAXVAR"`, `"SSQCORR"`, `"MINVAR"`, `"GENVAR"`.

4.  Principal component analysis (`"PCA"`)

5.  Factor score regression using sum scores (`"unit"`), regression
    (`"regression"`) or Bartlett scores (`"bartlett"`)

It is possible to supply starting values for the weighting algorithm via
`.starting_values`. The argument accepts a named list of vectors where
the list names are the construct names whose indicator weights the user
wishes to set. The vectors must be named vectors of
`"indicator_name" = value` pairs, where `value` is the starting weight.
See the examples section below for details.

Composite-indicator and composite-composite correlations are properly
disattenuated by default to yield consistent loadings, construct
correlations, and path coefficients if any of the concepts are modeled
as a common factor.

For *PLS-PM* disattenuation is done using *PLSc* (Dijkstra and Henseler
2015) . For *GSCA* disattenuation is done implicitly by using *GSCAm*
(Hwang et al. 2017) . Weights obtained by *GCCA*, *unit*, *regression*,
*bartlett* or *PCA* are disattenuated using Croon's approach (Croon
2002) . Disattenuation my be suppressed by setting
`.disattenuate = FALSE`. Note, however, that quantities in this case are
inconsistent estimates for their construct level counterparts if any of
the constructs in the structural model are modeled as a common factor!

By default path coefficients are estimated using ordinary least squares
(`.approach_path = "OLS"`). For linear models, two-stage least squares
(`"2SLS"`) is available, however, *only if* *instruments are internal*,
i.e., part of the structural model. Future versions will add support for
external instruments if possible. Instruments must be supplied to
`.instruments` as a named list where the names of the list elements are
the names of the dependent constructs of the structural equations whose
explanatory variables are believed to be endogenous. The list consists
of vectors of names of instruments corresponding to each equation. Note
that exogenous variables of a given equation **must** be supplied as
instruments for themselves.

If reliabilities are known they can be supplied as `"name" = value`
pairs to `.reliabilities`, where `value` is a numeric value between 0
and 1. Currently, only supported for "PLS-PM".

### Nonlinear models:

If the model contains nonlinear terms `csem()` estimates a polynomial
structural equation model using a non-iterative method of moments
approach described in Dijkstra and Schermelleh-Engel (2014) . Nonlinear
terms include interactions and exponential terms. The latter is
described in model syntax as an "interaction with itself", e.g.,
`xi^3 = xi.xi.xi`. Currently only exponential terms up to a power of
three (e.g., three-way interactions or cubic terms) are allowed:

1.  \- Single, e.g., `eta1`

2.  \- Quadratic, e.g., `eta1.eta1`

3.  \- Cubic, e.g., `eta1.eta1.eta1`

4.  \- Two-way interaction, e.g., `eta1.eta2`

5.  \- Three-way interaction, e.g., `eta1.eta2.eta3`

6.  \- Quadratic and two-way interaction, e.g., `eta1.eta1.eta3`

The current version of the package allows two kinds of estimation:
estimation of the reduced form equation (`.approach_nl = "replace"`) and
sequential estimation (`.approach_nl = "sequential"`, the default). The
latter does not allow for multivariate normality of all exogenous
variables, i.e., the latent variables and the error terms.

Distributional assumptions are kept to a minimum (an i.i.d. sample from
a population with finite moments for the relevant order); for higher
order models, that go beyond interaction, we work in this version with
the assumption that as far as the relevant moments are concerned certain
combinations of measurement errors behave as if they were Gaussian. For
details see: Dijkstra and Schermelleh-Engel (2014) .

### Models containing second-order constructs

Second-order constructs are specified using the operators `=~` and `<~`.
These operators are usually used with indicators on their right-hand
side. For second-order constructs the right-hand side variables are
constructs instead. If c1, and c2 are constructs forming or measuring a
higher-order construct, a model would look like this:

    my_model <- "
    # Structural model
    SAT  ~ QUAL
    VAL  ~ SAT

    # Measurement/composite model
    QUAL =~ qual1 + qual2
    SAT  =~ sat1 + sat2

    c1 =~ x11 + x12
    c2 =~ x21 + x22

    # Second-order construct (in this case a second-order composite build by common
    # factors)
    VAL <~ c1 + c2
    "

Currently, two approaches are explicitly implemented:

- (Default) `"2stage"`. The (disjoint) two-stage approach as proposed by
  Agarwal and Karahanna (2000) . Note that by default a correction for
  attenuation is applied if common factors are involved in modeling
  second-order constructs. For instance, the three-stage approach
  proposed by Van Riel et al. (2017) is applied in case of a
  second-order construct specified as a composite of common factors. On
  the other hand, if no common factors are involved the two-stage
  approach is applied as proposed by Schuberth et al. (2020) .

- `"mixed"`. The mixed repeated indicators/two-stage approach as
  proposed by Ringle et al. (2012) .

The repeated indicators approach as proposed by Joereskog and Wold
(1982) and the extension proposed by Becker et al. (2012) are not
directly implemented as they simply require a respecification of the
model. In the above example the repeated indicators approach would
require to change the model and to append the repeated indicators to the
data supplied to `.data`. Note that the indicators need to be renamed in
this case as `csem()` does not allow for one indicator to be attached to
multiple constructs.

    my_model <- "
    # Structural model
    SAT  ~ QUAL
    VAL  ~ SAT

    VAL ~ c1 + c2

    # Measurement/composite model
    QUAL =~ qual1 + qual2
    SAT  =~ sat1 + sat2
    VAL  =~ x11_temp + x12_temp + x21_temp + x22_temp

    c1 =~ x11 + x12
    c2 =~ x21 + x22
    "

According to the extended approach indirect effects of `QUAL` on `VAL`
via `c1` and `c2` would have to be specified as well.

### Multigroup analysis

To perform a multigroup analysis provide either a list of data sets or
one data set containing a group-identifier-column whose column name must
be provided to `.id`. Values of this column are taken as levels of a
factor and are interpreted as group identifiers. `csem()` will split the
data by levels of that column and run the estimation for each level
separately. Note, the more levels the group-identifier-column has, the
more estimation runs are required. This can considerably slow down
estimation, especially if resampling is requested. For the latter it
will generally be faster to use `.eval_plan = "multisession"` or
`.eval_plan = "multicore"`.

### Inference:

Inference is done via resampling. See
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
and [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
for details.

## Postestimation

- [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md):

  Assess results using common quality criteria, e.g., reliability, fit
  measures, HTMT, R2 etc.

- [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md):

  Calculate common inferential quantities, e.g., standard errors,
  confidence intervals.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html):

  Creates a plot of the model. For the help file see
  [`plot.cSEMResults_default()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_default.md).

- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md):

  Predict endogenous indicator scores and compute common prediction
  metrics.

- [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md):

  Summarize the results. Mainly called for its side-effect the print
  method.

- [`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md):

  Verify/Check admissibility of the estimates.

Tests are performed using the test-family of functions. Currently the
following tests are implemented:

- [`testCVPAT()`](https://floschuberth.github.io/cSEM/reference/testCVPAT.md):

  Cross-validated predictive ability test proposed by Liengaard et al.
  (2021)

- [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md):

  Bootstrap-based test for overall model fit based on Beran and
  Srivastava (1985) .

- [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md):

  Permutation-based test for measurement invariance of composites
  proposed by Henseler et al. (2016) .

- [`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md):

  Several (mainly) permutation-based tests for multi-group comparisons.

- [`testHausman()`](https://floschuberth.github.io/cSEM/reference/testHausman.md):

  Regression-based Hausman test to test for endogeneity.

Other miscellaneous postestimation functions belong do the do-family of
functions. Currently three do functions are implemented:

- [`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md):

  Performs an importance-performance matrix analysis (IPMA).

- [`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md):

  Perform a nonlinear effects analysis as described in e.g., Spiller et
  al. (2013)

- [`doRedundancyAnalysis()`](https://floschuberth.github.io/cSEM/reference/doRedundancyAnalysis.md):

  Perform a redundancy analysis (RA) as proposed by Hair et al. (2016)
  with reference to Chin (1998)

## References

Agarwal R, Karahanna E (2000). “Time Flies When You're Having Fun:
Cognitive Absorption and Beliefs about Information Technology Usage.”
*MIS Quarterly*, **24**(4), 665.  
  
Becker J, Klein K, Wetzels M (2012). “Hierarchical Latent Variable
Models in PLS-SEM: Guidelines for Using Reflective-Formative Type
Models.” *Long Range Planning*, **45**(5-6), 359–394.
[doi:10.1016/j.lrp.2012.10.001](https://doi.org/10.1016/j.lrp.2012.10.001)
.  
  
Beran R, Srivastava MS (1985). “Bootstrap Tests and Confidence Regions
for Functions of a Covariance Matrix.” *The Annals of Statistics*,
**13**(1), 95–115.
[doi:10.1214/aos/1176346579](https://doi.org/10.1214/aos/1176346579) .  
  
Chin WW (1998). “Modern Methods for Business Research.” In Marcoulides
GA (ed.), chapter The Partial Least Squares Approach to Structural
Equation Modeling, 295–358. Mahwah, NJ: Lawrence Erlbaum.  
  
Croon MA (2002). “Using predicted latent scores in general latent
structure models.” In Marcoulides GA, Moustaki I (eds.), *Latent
Variable and Latent Structure Models*, chapter 10, 195–224. Lawrence
Erlbaum. ISBN 080584046X, Pagination: 288.  
  
Dijkstra TK, Henseler J (2015). “Consistent and Asymptotically Normal
PLS Estimators for Linear Structural Equations.” *Computational
Statistics & Data Analysis*, **81**, 10–23.  
  
Dijkstra TK, Schermelleh-Engel K (2014). “Consistent Partial Least
Squares For Nonlinear Structural Equation Models.” *Psychometrika*,
**79**(4), 585–604.  
  
Hair JF, Hult GTM, Ringle C, Sarstedt M (2016). *A Primer on Partial
Least Squares Structural Equation Modeling (PLS-SEM)*. Sage
publications.  
  
Henseler J, Ringle CM, Sarstedt M (2016). “Testing Measurement
Invariance of Composites Using Partial Least Squares.” *International
Marketing Review*, **33**(3), 405–431.
[doi:10.1108/imr-09-2014-0304](https://doi.org/10.1108/imr-09-2014-0304)
.  
  
Hwang H, Takane Y, Jung K (2017). “Generalized structured component
analysis with uniqueness terms for accommodating measurement error.”
*Frontiers in Psychology*, **8**(2137), 1–12.  
  
Joereskog KG, Wold HO (1982). *Systems under Indirect Observation:
Causality, Structure, Prediction - Part II*, volume 139. North
Holland.  
  
Liengaard BD, Sharma PN, Hult GTM, Jensen MB, Sarstedt M, Hair JF,
Ringle CM (2021). “Prediction: Coveted, Yet Forsaken? Introducing a
Cross-Validated Predictive Ability Test in Partial Least Squares Path
Modeling.” *Decision Sciences*, **52**(2), 362–392.  
  
Ringle CM, Sarstedt M, Straub D (2012). “A Critical Look at the Use of
PLS-SEM in MIS Quarterly.” *MIS Quarterly*, **36**(1), iii–xiv.  
  
Schuberth F, Rademaker ME, Henseler J (2020). “Estimating and assessing
second-order constructs using PLS-PM: the case of composites of
composites.” *Industrial Management & Data Systems*, **120**(12),
2211-2241.
[doi:10.1108/imds-12-2019-0642](https://doi.org/10.1108/imds-12-2019-0642)
.  
  
Spiller SA, Fitzsimons GJ, Lynch JG, Mcclelland GH (2013). “Spotlights,
Floodlights, and the Magic Number Zero: Simple Effects Tests in
Moderated Regression.” *Journal of Marketing Research*, **50**(2),
277–288. [doi:10.1509/jmr.12.0420](https://doi.org/10.1509/jmr.12.0420)
.  
  
Van Riel ACR, Henseler J, Kemeny I, Sasovova Z (2017). “Estimating
hierarchical constructs using Partial Least Squares: The case of second
order composites of factors.” *Industrial Management & Data Systems*,
**117**(3), 459–477.
[doi:10.1108/IMDS-07-2016-0286](https://doi.org/10.1108/IMDS-07-2016-0286)
.

## See also

[`args_default()`](https://floschuberth.github.io/cSEM/reference/args_default.md),
[cSEMArguments](https://floschuberth.github.io/cSEM/reference/csem_arguments.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`foreman()`](https://floschuberth.github.io/cSEM/reference/foreman.md),
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md),
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md),
[`plot.cSEMResults_default()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_default.md),
[`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md),
[`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md),
[`testCVPAT()`](https://floschuberth.github.io/cSEM/reference/testCVPAT.md),
[`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md),
[`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md),
[`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md),
[`testHausman()`](https://floschuberth.github.io/cSEM/reference/testHausman.md)

## Examples

``` r
# ===========================================================================
# Basic usage
# ===========================================================================
### Linear model ------------------------------------------------------------
# Most basic usage requires a dataset and a model. We use the 
#  `threecommonfactors` dataset. 

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
res <- csem(threecommonfactors, model)

## Postestimation
verify(res)
#> ________________________________________________________________________________
#> 
#> Verify admissibility:
#> 
#>     admissible
#> 
#> Details:
#> 
#>   Code   Status    Description
#>   1      ok        Convergence achieved                                   
#>   2      ok        All absolute standardized loading estimates <= 1       
#>   3      ok        Construct VCV is positive semi-definite                
#>   4      ok        All reliability estimates <= 1                         
#>   5      ok        Model-implied indicator VCV is positive semi-definite  
#> ________________________________________________________________________________
summarize(res)  
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
#>   Path           Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1      0.6713          NA        NA        NA
#>   eta3 ~ eta1      0.4585          NA        NA        NA
#>   eta3 ~ eta2      0.3052          NA        NA        NA
#> 
#> Estimated loadings:
#> ===================
#>   Loading        Estimate  Std. error   t-stat.   p-value
#>   eta1 =~ y11      0.6631          NA        NA        NA
#>   eta1 =~ y12      0.6493          NA        NA        NA
#>   eta1 =~ y13      0.7613          NA        NA        NA
#>   eta2 =~ y21      0.5165          NA        NA        NA
#>   eta2 =~ y22      0.7554          NA        NA        NA
#>   eta2 =~ y23      0.7997          NA        NA        NA
#>   eta3 =~ y31      0.8223          NA        NA        NA
#>   eta3 =~ y32      0.6581          NA        NA        NA
#>   eta3 =~ y33      0.7474          NA        NA        NA
#> 
#> Estimated weights:
#> ==================
#>   Weight         Estimate  Std. error   t-stat.   p-value
#>   eta1 <~ y11      0.3956          NA        NA        NA
#>   eta1 <~ y12      0.3873          NA        NA        NA
#>   eta1 <~ y13      0.4542          NA        NA        NA
#>   eta2 <~ y21      0.3058          NA        NA        NA
#>   eta2 <~ y22      0.4473          NA        NA        NA
#>   eta2 <~ y23      0.4735          NA        NA        NA
#>   eta3 <~ y31      0.4400          NA        NA        NA
#>   eta3 <~ y32      0.3521          NA        NA        NA
#>   eta3 <~ y33      0.3999          NA        NA        NA
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
assess(res)
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

# Notes: 
#   1. By default no inferential quantities (e.g. Std. errors, p-values, or
#      confidence intervals) are calculated. Use resampling to obtain
#      inferential quantities. See "Resampling" in the "Extended usage"
#      section below.
#   2. `summarize()` prints the full output by default. For a more condensed
#       output use:
print(summarize(res), .full_output = FALSE)
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
#>   Path           Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1      0.6713          NA        NA        NA
#>   eta3 ~ eta1      0.4585          NA        NA        NA
#>   eta3 ~ eta2      0.3052          NA        NA        NA
#> 
#> Estimated loadings:
#> ===================
#>   Loading        Estimate  Std. error   t-stat.   p-value
#>   eta1 =~ y11      0.6631          NA        NA        NA
#>   eta1 =~ y12      0.6493          NA        NA        NA
#>   eta1 =~ y13      0.7613          NA        NA        NA
#>   eta2 =~ y21      0.5165          NA        NA        NA
#>   eta2 =~ y22      0.7554          NA        NA        NA
#>   eta2 =~ y23      0.7997          NA        NA        NA
#>   eta3 =~ y31      0.8223          NA        NA        NA
#>   eta3 =~ y32      0.6581          NA        NA        NA
#>   eta3 =~ y33      0.7474          NA        NA        NA
#> 
#> Estimated weights:
#> ==================
#>   Weight         Estimate  Std. error   t-stat.   p-value
#>   eta1 <~ y11      0.3956          NA        NA        NA
#>   eta1 <~ y12      0.3873          NA        NA        NA
#>   eta1 <~ y13      0.4542          NA        NA        NA
#>   eta2 <~ y21      0.3058          NA        NA        NA
#>   eta2 <~ y22      0.4473          NA        NA        NA
#>   eta2 <~ y23      0.4735          NA        NA        NA
#>   eta3 <~ y31      0.4400          NA        NA        NA
#>   eta3 <~ y32      0.3521          NA        NA        NA
#>   eta3 <~ y33      0.3999          NA        NA        NA
#> ________________________________________________________________________________

## Dealing with endogeneity -------------------------------------------------

# See: ?testHausman()

### Models containing second constructs--------------------------------------
## Take a look at the dataset
#?dgp_2ndorder_cf_of_c

model <- "
# Path model / Regressions 
c4   ~ eta1
eta2 ~ eta1 + c4

# Reflective measurement model
c1   <~ y11 + y12 
c2   <~ y21 + y22 + y23 + y24
c3   <~ y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38
eta1 =~ y41 + y42 + y43
eta2 =~ y51 + y52 + y53

# Composite model (second order)
c4   =~ c1 + c2 + c3
"

res_2stage <- csem(dgp_2ndorder_cf_of_c, model, .approach_2ndorder = "2stage")
res_mixed  <- csem(dgp_2ndorder_cf_of_c, model, .approach_2ndorder = "mixed")

# The standard repeated indicators approach is done by 1.) respecifying the model
# and 2.) adding the repeated indicators to the data set

# 1.) Respecify the model
model_RI <- "
# Path model / Regressions 
c4   ~ eta1
eta2 ~ eta1 + c4
c4   ~ c1 + c2 + c3

# Reflective measurement model
c1   <~ y11 + y12 
c2   <~ y21 + y22 + y23 + y24
c3   <~ y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38
eta1 =~ y41 + y42 + y43
eta2 =~ y51 + y52 + y53

# c4 is a common factor measured by composites
c4 =~ y11_temp + y12_temp + y21_temp + y22_temp + y23_temp + y24_temp +
      y31_temp + y32_temp + y33_temp + y34_temp + y35_temp + y36_temp + 
      y37_temp + y38_temp
"

# 2.) Update data set
data_RI <- dgp_2ndorder_cf_of_c
coln <- c(colnames(data_RI), paste0(colnames(data_RI), "_temp"))
data_RI <- data_RI[, c(1:ncol(data_RI), 1:ncol(data_RI))]
colnames(data_RI) <- coln

# Estimate
res_RI <- csem(data_RI, model_RI)
summarize(res_RI)
#> ________________________________________________________________________________
#> ----------------------------------- Overview -----------------------------------
#> 
#>  General information:
#>  ------------------------
#>  Estimation status                  = Not ok!
#>  Number of observations             = 500
#>  Weight estimator                   = PLS-PM
#>  Inner weighting scheme             = "path"
#>  Type of indicator correlation      = Pearson
#>  Path model estimator               = OLS
#>  Second-order approach              = NA
#>  Type of path model                 = Linear
#>  Disattenuated                      = Yes (PLSc)
#> 
#>  Construct details:
#>  ------------------
#>  Name  Modeled as     Order         Mode      
#> 
#>  eta1  Common factor  First order   "modeA"   
#>  c1    Composite      First order   "modeB"   
#>  c2    Composite      First order   "modeB"   
#>  c3    Composite      First order   "modeB"   
#>  c4    Common factor  First order   "modeA"   
#>  eta2  Common factor  First order   "modeA"   
#> 
#> ----------------------------------- Estimates ----------------------------------
#> 
#> Estimated path coefficients:
#> ============================
#>   Path           Estimate  Std. error   t-stat.   p-value
#>   c4 ~ eta1        0.0029          NA        NA        NA
#>   c4 ~ c1          0.2333          NA        NA        NA
#>   c4 ~ c2          0.4381          NA        NA        NA
#>   c4 ~ c3          0.5448          NA        NA        NA
#>   eta2 ~ eta1      0.0345          NA        NA        NA
#>   eta2 ~ c4        0.5128          NA        NA        NA
#> 
#> Estimated loadings:
#> ===================
#>   Loading           Estimate  Std. error   t-stat.   p-value
#>   eta1 =~ y41         0.9416          NA        NA        NA
#>   eta1 =~ y42         0.7374          NA        NA        NA
#>   eta1 =~ y43         0.5285          NA        NA        NA
#>   c1 =~ y11           0.8999          NA        NA        NA
#>   c1 =~ y12           0.7204          NA        NA        NA
#>   c2 =~ y21           0.8041          NA        NA        NA
#>   c2 =~ y22           0.7090          NA        NA        NA
#>   c2 =~ y23           0.6563          NA        NA        NA
#>   c2 =~ y24           0.6739          NA        NA        NA
#>   c3 =~ y31           0.5380          NA        NA        NA
#>   c3 =~ y32           0.6091          NA        NA        NA
#>   c3 =~ y33           0.6759          NA        NA        NA
#>   c3 =~ y34           0.4268          NA        NA        NA
#>   c3 =~ y35           0.3482          NA        NA        NA
#>   c3 =~ y36           0.6089          NA        NA        NA
#>   c3 =~ y37           0.4549          NA        NA        NA
#>   c3 =~ y38           0.6092          NA        NA        NA
#>   c4 =~ y11_temp      0.6249          NA        NA        NA
#>   c4 =~ y12_temp      0.4750          NA        NA        NA
#>   c4 =~ y21_temp      0.6808          NA        NA        NA
#>   c4 =~ y22_temp      0.5885          NA        NA        NA
#>   c4 =~ y23_temp      0.5291          NA        NA        NA
#>   c4 =~ y24_temp      0.5771          NA        NA        NA
#>   c4 =~ y31_temp      0.4817          NA        NA        NA
#>   c4 =~ y32_temp      0.5415          NA        NA        NA
#>   c4 =~ y33_temp      0.5586          NA        NA        NA
#>   c4 =~ y34_temp      0.3624          NA        NA        NA
#>   c4 =~ y35_temp      0.3157          NA        NA        NA
#>   c4 =~ y36_temp      0.5154          NA        NA        NA
#>   c4 =~ y37_temp      0.4129          NA        NA        NA
#>   c4 =~ y38_temp      0.5111          NA        NA        NA
#>   eta2 =~ y51         0.7991          NA        NA        NA
#>   eta2 =~ y52         0.8477          NA        NA        NA
#>   eta2 =~ y53         0.7304          NA        NA        NA
#> 
#> Estimated weights:
#> ==================
#>   Weight            Estimate  Std. error   t-stat.   p-value
#>   eta1 <~ y41         0.5052          NA        NA        NA
#>   eta1 <~ y42         0.3957          NA        NA        NA
#>   eta1 <~ y43         0.2836          NA        NA        NA
#>   c1 <~ y11           0.7392          NA        NA        NA
#>   c1 <~ y12           0.4648          NA        NA        NA
#>   c2 <~ y21           0.4487          NA        NA        NA
#>   c2 <~ y22           0.3168          NA        NA        NA
#>   c2 <~ y23           0.2773          NA        NA        NA
#>   c2 <~ y24           0.3452          NA        NA        NA
#>   c3 <~ y31           0.2758          NA        NA        NA
#>   c3 <~ y32           0.2653          NA        NA        NA
#>   c3 <~ y33           0.2202          NA        NA        NA
#>   c3 <~ y34           0.1587          NA        NA        NA
#>   c3 <~ y35           0.1682          NA        NA        NA
#>   c3 <~ y36           0.2495          NA        NA        NA
#>   c3 <~ y37           0.2784          NA        NA        NA
#>   c3 <~ y38           0.2238          NA        NA        NA
#>   c4 <~ y11_temp      0.1510          NA        NA        NA
#>   c4 <~ y12_temp      0.1148          NA        NA        NA
#>   c4 <~ y21_temp      0.1645          NA        NA        NA
#>   c4 <~ y22_temp      0.1422          NA        NA        NA
#>   c4 <~ y23_temp      0.1279          NA        NA        NA
#>   c4 <~ y24_temp      0.1395          NA        NA        NA
#>   c4 <~ y31_temp      0.1164          NA        NA        NA
#>   c4 <~ y32_temp      0.1309          NA        NA        NA
#>   c4 <~ y33_temp      0.1350          NA        NA        NA
#>   c4 <~ y34_temp      0.0876          NA        NA        NA
#>   c4 <~ y35_temp      0.0763          NA        NA        NA
#>   c4 <~ y36_temp      0.1246          NA        NA        NA
#>   c4 <~ y37_temp      0.0998          NA        NA        NA
#>   c4 <~ y38_temp      0.1235          NA        NA        NA
#>   eta2 <~ y51         0.3873          NA        NA        NA
#>   eta2 <~ y52         0.4109          NA        NA        NA
#>   eta2 <~ y53         0.3540          NA        NA        NA
#> 
#> Estimated construct correlations:
#> =================================
#>   Correlation    Estimate  Std. error   t-stat.   p-value
#>   eta1 ~~ c1       0.2882          NA        NA        NA
#>   eta1 ~~ c2       0.2527          NA        NA        NA
#>   eta1 ~~ c3       0.2871          NA        NA        NA
#>   c1 ~~ c2         0.5772          NA        NA        NA
#>   c1 ~~ c3         0.6242          NA        NA        NA
#>   c2 ~~ c3         0.7480          NA        NA        NA
#> 
#> Estimated indicator correlations:
#> =================================
#>   Correlation    Estimate  Std. error   t-stat.   p-value
#>   y11 ~~ y12       0.3459          NA        NA        NA
#>   y21 ~~ y22       0.4341          NA        NA        NA
#>   y21 ~~ y23       0.3805          NA        NA        NA
#>   y21 ~~ y24       0.3255          NA        NA        NA
#>   y22 ~~ y23       0.3260          NA        NA        NA
#>   y22 ~~ y24       0.3102          NA        NA        NA
#>   y23 ~~ y24       0.3043          NA        NA        NA
#>   y31 ~~ y32       0.1558          NA        NA        NA
#>   y31 ~~ y33       0.2728          NA        NA        NA
#>   y31 ~~ y34      -0.1472          NA        NA        NA
#>   y31 ~~ y35       0.1617          NA        NA        NA
#>   y31 ~~ y36       0.3372          NA        NA        NA
#>   y31 ~~ y37       0.0961          NA        NA        NA
#>   y31 ~~ y38       0.2059          NA        NA        NA
#>   y32 ~~ y33       0.2355          NA        NA        NA
#>   y32 ~~ y34       0.4146          NA        NA        NA
#>   y32 ~~ y35       0.2228          NA        NA        NA
#>   y32 ~~ y36       0.2184          NA        NA        NA
#>   y32 ~~ y37       0.2684          NA        NA        NA
#>   y32 ~~ y38       0.0736          NA        NA        NA
#>   y33 ~~ y34       0.2908          NA        NA        NA
#>   y33 ~~ y35      -0.0586          NA        NA        NA
#>   y33 ~~ y36       0.3445          NA        NA        NA
#>   y33 ~~ y37       0.2018          NA        NA        NA
#>   y33 ~~ y38       0.6233          NA        NA        NA
#>   y34 ~~ y35       0.1723          NA        NA        NA
#>   y34 ~~ y36       0.1704          NA        NA        NA
#>   y34 ~~ y37       0.0729          NA        NA        NA
#>   y34 ~~ y38       0.1916          NA        NA        NA
#>   y35 ~~ y36       0.1125          NA        NA        NA
#>   y35 ~~ y37      -0.1425          NA        NA        NA
#>   y35 ~~ y38       0.3285          NA        NA        NA
#>   y36 ~~ y37       0.1102          NA        NA        NA
#>   y36 ~~ y38       0.2499          NA        NA        NA
#>   y37 ~~ y38       0.0858          NA        NA        NA
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>   Total effect    Estimate  Std. error   t-stat.   p-value
#>   c4 ~ eta1         0.0029          NA        NA        NA
#>   c4 ~ c1           0.2333          NA        NA        NA
#>   c4 ~ c2           0.4381          NA        NA        NA
#>   c4 ~ c3           0.5448          NA        NA        NA
#>   eta2 ~ eta1       0.0359          NA        NA        NA
#>   eta2 ~ c1         0.1197          NA        NA        NA
#>   eta2 ~ c2         0.2247          NA        NA        NA
#>   eta2 ~ c3         0.2794          NA        NA        NA
#>   eta2 ~ c4         0.5128          NA        NA        NA
#> 
#> Estimated indirect effects:
#> ===========================
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1          0.0015          NA        NA        NA
#>   eta2 ~ c1            0.1197          NA        NA        NA
#>   eta2 ~ c2            0.2247          NA        NA        NA
#>   eta2 ~ c3            0.2794          NA        NA        NA
#> ________________________________________________________________________________

### Multigroup analysis -----------------------------------------------------

# See: ?testMGD()

# ===========================================================================
# Extended usage
# ===========================================================================
# `csem()` provides defaults for all arguments except `.data` and `.model`.
#   Below some common options/tasks that users are likely to be interested in.
#   We use the threecommonfactors data set again:

model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

### PLS vs PLSc and disattenuation
# In the model all concepts are modeled as common factors. If 
#   .approach_weights = "PLS-PM", csem() uses PLSc to disattenuate composite-indicator 
#   and composite-composite correlations.
res_plsc <- csem(threecommonfactors, model, .approach_weights = "PLS-PM")
res$Information$Model$construct_type # all common factors
#>            eta1            eta2            eta3 
#> "Common factor" "Common factor" "Common factor" 

# To obtain "original" (inconsistent) PLS estimates use `.disattenuate = FALSE`
res_pls <- csem(threecommonfactors, model, 
                .approach_weights = "PLS-PM",
                .disattenuate = FALSE
                )

s_plsc <- summarize(res_plsc)
s_pls  <- summarize(res_pls)

# Compare
data.frame(
  "Path"      = s_plsc$Estimates$Path_estimates$Name,
  "Pop_value" = c(0.6, 0.4, 0.35), # see ?threecommonfactors
  "PLSc"      = s_plsc$Estimates$Path_estimates$Estimate,
  "PLS"       = s_pls$Estimates$Path_estimates$Estimate
  )
#>          Path Pop_value      PLSc       PLS
#> 1 eta2 ~ eta1      0.60 0.6713334 0.5046062
#> 2 eta3 ~ eta1      0.40 0.4585068 0.3588557
#> 3 eta3 ~ eta2      0.35 0.3051511 0.2972680

### Resampling --------------------------------------------------------------
if (FALSE) { # \dontrun{
## Basic resampling
res_boot <- csem(threecommonfactors, model, .resample_method = "bootstrap")
res_jack <- csem(threecommonfactors, model, .resample_method = "jackknife")

# See ?resamplecSEMResults for more examples

### Choosing a different weightning scheme ----------------------------------

res_gscam  <- csem(threecommonfactors, model, .approach_weights = "GSCA")
res_gsca   <- csem(threecommonfactors, model, 
                   .approach_weights = "GSCA",
                   .disattenuate = FALSE
)

s_gscam <- summarize(res_gscam)
s_gsca  <- summarize(res_gsca)

# Compare
data.frame(
  "Path"      = s_gscam$Estimates$Path_estimates$Name,
  "Pop_value" = c(0.6, 0.4, 0.35), # see ?threecommonfactors
  "GSCAm"      = s_gscam$Estimates$Path_estimates$Estimate,
  "GSCA"       = s_gsca$Estimates$Path_estimates$Estimate
)} # }
### Fine-tuning a weighting scheme ------------------------------------------
## Setting starting values

sv <- list("eta1" = c("y12" = 10, "y13" = 4, "y11" = 1))
res <- csem(threecommonfactors, model, .starting_values = sv)

## Choosing a different inner weighting scheme 
#?args_csem_dotdotdot

res <- csem(threecommonfactors, model, .PLS_weight_scheme_inner = "factorial",
            .PLS_ignore_structural_model = TRUE)


## Choosing different modes for PLS
# By default, concepts modeled as common factors uses PLS Mode A weights.
modes <- list("eta1" = "unit", "eta2" = "modeB", "eta3" = "unit")
res   <- csem(threecommonfactors, model, .PLS_modes = modes)
summarize(res) 
#> ________________________________________________________________________________
#> ----------------------------------- Overview -----------------------------------
#> 
#>  General information:
#>  ------------------------
#>  Estimation status                  = Not ok!
#>  Number of observations             = 500
#>  Weight estimator                   = PLS-PM
#>  Inner weighting scheme             = "path"
#>  Type of indicator correlation      = Pearson
#>  Path model estimator               = OLS
#>  Second-order approach              = NA
#>  Type of path model                 = Linear
#>  Disattenuated                      = Yes (PLSc)
#> 
#>  Construct details:
#>  ------------------
#>  Name  Modeled as     Order         Mode      
#> 
#>  eta1  Common factor  First order   "unit"    
#>  eta2  Common factor  First order   "modeB"   
#>  eta3  Common factor  First order   "unit"    
#> 
#> ----------------------------------- Estimates ----------------------------------
#> 
#> Estimated path coefficients:
#> ============================
#>   Path           Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1      0.5700          NA        NA        NA
#>   eta3 ~ eta1      0.2869          NA        NA        NA
#>   eta3 ~ eta2      0.3840          NA        NA        NA
#> 
#> Estimated loadings:
#> ===================
#>   Loading        Estimate  Std. error   t-stat.   p-value
#>   eta1 =~ y11      0.8017          NA        NA        NA
#>   eta1 =~ y12      0.7907          NA        NA        NA
#>   eta1 =~ y13      0.8278          NA        NA        NA
#>   eta2 =~ y21      0.5143          NA        NA        NA
#>   eta2 =~ y22      0.7570          NA        NA        NA
#>   eta2 =~ y23      0.7999          NA        NA        NA
#>   eta3 =~ y31      0.8603          NA        NA        NA
#>   eta3 =~ y32      0.8216          NA        NA        NA
#>   eta3 =~ y33      0.8286          NA        NA        NA
#> 
#> Estimated weights:
#> ==================
#>   Weight         Estimate  Std. error   t-stat.   p-value
#>   eta1 <~ y11      0.4132          NA        NA        NA
#>   eta1 <~ y12      0.4132          NA        NA        NA
#>   eta1 <~ y13      0.4132          NA        NA        NA
#>   eta2 <~ y21      0.1593          NA        NA        NA
#>   eta2 <~ y22      0.4538          NA        NA        NA
#>   eta2 <~ y23      0.5722          NA        NA        NA
#>   eta3 <~ y31      0.3983          NA        NA        NA
#>   eta3 <~ y32      0.3983          NA        NA        NA
#>   eta3 <~ y33      0.3983          NA        NA        NA
#> 
#> ------------------------------------ Effects -----------------------------------
#> 
#> Estimated total effects:
#> ========================
#>   Total effect    Estimate  Std. error   t-stat.   p-value
#>   eta2 ~ eta1       0.5700          NA        NA        NA
#>   eta3 ~ eta1       0.5057          NA        NA        NA
#>   eta3 ~ eta2       0.3840          NA        NA        NA
#> 
#> Estimated indirect effects:
#> ===========================
#>   Indirect effect    Estimate  Std. error   t-stat.   p-value
#>   eta3 ~ eta1          0.2189          NA        NA        NA
#> ________________________________________________________________________________
```
