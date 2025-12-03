# cSEMArguments

An alphabetical list of all arguments used by functions of the `cSEM`
package including their description and defaults. Mainly used for
internal purposes (parameter inheritance). To list all arguments and
their defaults, use
[`args_default()`](https://floschuberth.github.io/cSEM/reference/args_default.md).
To list all arguments and their possible choices, use
`args_default(.choices = TRUE)`.

## Arguments

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .absolute:

  Logical. Should the absolute HTMT values be returned? Defaults to
  `TRUE` .

- .approach_gcca:

  Character string. The Kettenring approach to use for GCCA. One of
  "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*".
  Defaults to "*SUMCORR*".

- .approach_2ndorder:

  Character string. Approach used for models containing second-order
  constructs. One of: "*2stage*", or "*mixed*". Defaults to "*2stage*".

- .approach_alpha_adjust:

  Character string. Approach used to adjust the significance level to
  accommodate multiple testing. One of "*none*" or "*bonferroni*".
  Defaults to "*none*".

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

- .approach_mgd:

  Character string or a vector of character strings. Approach used for
  the multi-group comparison. One of: "*all*", "*Klesel*", "*Chin*",
  "*Sarstedt*", "*Keil*, "*Nitzl*", "*Henseler*", "*CI_para*", or
  "*CI_overlap*". Default to "*all*" in which case all approaches are
  computed (if possible).

- .approach_nl:

  Character string. Approach used to estimate nonlinear structural
  relationships. One of: "*sequential*" or "*replace*". Defaults to
  "*sequential*".

- .approach_predict:

  Character string. Which approach should be used to predictions? One of
  "*earliest*" and "*direct*". If "*earliest*" predictions for
  indicators associated to endogenous constructs are performed using
  only indicators associated to exogenous constructs. If "*direct*",
  predictions for indicators associated to endogenous constructs are
  based on indicators associated to their direct antecedents. Defaults
  to "*earliest*".

- .approach_p_adjust:

  Character string or a vector of character strings. Approach used to
  adjust the p-value for multiple testing. See the `methods` argument of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for a
  list of choices and their description. Defaults to "*none*".

- .approach_paths:

  Character string. Approach used to estimate the structural
  coefficients. One of: "*OLS*" or "*2SLS*". If "*2SLS*", instruments
  need to be supplied to `.instruments`. Defaults to "*OLS*".

- .approach_score_benchmark:

  Character string. How should the aggregation of the estimates of the
  truncated normal distribution be done for the benchmark predictions?
  Ignored if not OrdPLS or OrdPLSc is used to obtain benchmark
  predictions. One of "*mean*", "*median*", "*mode*" or "*round*". If
  "*round*", the benchmark predictions are obtained using the
  traditional prediction algorithm for PLS-PM which are rounded for
  categorical indicators. If "*mean*", the mean of the estimated
  endogenous indicators is calculated. If "*median*", the mean of the
  estimated endogenous indicators is calculated. If "*mode*", the
  maximum empirical density on the intervals defined by the thresholds
  is used. If `.treat_as_continuous = TRUE` or if all indicators are on
  a continuous scale, `.approach_score_benchmark` is ignored. Defaults
  to "*round*".

- .approach_score_target:

  Character string. How should the aggregation of the estimates of the
  truncated normal distribution for the predictions using OrdPLS/OrdPLSc
  be done? One of "*mean*", "*median*" or "*mode*". If "*mean*", the
  mean of the estimated endogenous indicators is calculated. If
  "*median*", the mean of the estimated endogenous indicators is
  calculated. If "*mode*", the maximum empirical density on the
  intervals defined by the thresholds is used. Defaults to "*mean*".

- .approach_weights:

  Character string. Approach used to obtain composite weights. One of:
  "*PLS-PM*", "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*",
  "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", or
  "*regression*". Defaults to "*PLS-PM*".

- .args_used:

  A list of function argument names whose value was modified by the
  user.

- .attrbutes:

  Character string. Variables used as attributes in IPMA.

- .benchmark:

  Character string. The procedure to obtain benchmark predictions. One
  of "*lm*", "*unit*", "*PLS-PM*", "*GSCA*", "*PCA*", "*MAXVAR*", or
  "*NA*". Default to "*lm*".

- .bias_corrected:

  Logical. Should the standard and the tStat confidence interval be
  bias-corrected using the bootstrapped bias estimate? If `TRUE` the
  confidence interval for some estimated parameter `theta` is centered
  at `2*theta - theta*_hat`, where `theta*_hat` is the average over all
  `.R` bootstrap estimates of `theta`. Defaults to `TRUE`

- .by_equation:

  Should the criteria be computed for each structural model equation
  separately? Defaults to `TRUE`.

- .C:

  A (J x J) composite variance-covariance matrix.

- .check_errors:

  Logical. Should the model to parse be checked for correctness in a
  sense that all necessary components to estimate the model are given?
  Defaults to `TRUE`.

- .choices:

  Logical. Should candidate values for the arguments be returned?
  Defaults to `FALSE`.

- .ci:

  A vector of character strings naming the confidence interval to
  compute. For possible choices see
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).

- .ci_colnames:

  Internal argument used by several print helper functions.

- .closed_form_ci:

  Logical. Should a closed-form confidence interval be computed?
  Defaults to `FALSE`.

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .csem_resample:

  A list resulting from a call to
  [`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md).

- .cv_folds:

  Integer. The number of cross-validation folds to use. Setting
  `.cv_folds` to `N` (the number of observations) produces leave-one-out
  cross-validation samples. Defaults to `10`.

- .data:

  A `data.frame` or a `matrix` of standardized or unstandardized data
  (indicators/items/manifest variables). Possible column types or
  classes of the data provided are: "`logical`", "`numeric`" ("`double`"
  or "`integer`"), "`factor`" ("`ordered`" and/or "`unordered`"),
  "`character`" (converted to factor), or a mix of several types.

- .dependent:

  Character string. The name of the dependent variable.

- .disattenuate:

  Logical. Should composite/proxy correlations be disattenuated to yield
  consistent loadings and path estimates if at least one of the
  construct is modeled as a common factor? Defaults to `TRUE`.

- .dist:

  Character string. The distribution to use for the critical value. One
  of *"t"* for Student's t-distribution or *"z"* for the standard normal
  distribution. Defaults to *"z"*.

- .distance:

  Character string. A distance measure. One of: "*geodesic*" or
  "*squared_euclidean*". Defaults to "*geodesic*".

- .df:

  Character string. The method for obtaining the degrees of freedom.
  Choices are "*type1*" and "*type2*". Defaults to "*type1*" .

- .dominant_indicators:

  A character vector of `"construct_name" = "indicator_name"` pairs,
  where `"indicator_name"` is a character string giving the name of the
  dominant indicator and `"construct_name"` a character string of the
  corresponding construct name. Dominant indicators may be specified for
  a subset of the constructs. Default to `NULL`.

- .E:

  A (J x J) matrix of inner weights.

- .effect:

  Internal argument used by helper printEffects().

- .estimate_structural:

  Logical. Should the structural coefficients be estimated? Defaults to
  `TRUE`.

- .eval_plan:

  Character string. The evaluation plan to use. One of "*sequential*",
  "*multicore*", or "*multisession*". In the two latter cases all
  available cores will be used. Defaults to "*sequential*".

- .filename:

  Character string. The file name.

- .first_resample:

  A list containing the `.R` resamples based on the original data
  obtained by resamplecSEMResults().

- .fit_measures:

  Logical. (EXPERIMENTAL) Should additional fit measures be included?
  Defaults to `FALSE`.

- .force:

  Logical. Should .object be resampled even if it contains resamples
  already?. Defaults to `FALSE`.

- .full_output:

  Logical. Should the full output of summarize be printed. Defaults to
  `TRUE`.

- .graph_attrs:

  Character string. Additional attributes that should be passed to the
  DiagrammeR syntax, e.g., c("rankdir=LR", "ranksep=1.0"). Defaults to
  *c("rankdir=LR")*.

- .H:

  The (N x J) matrix of construct scores.

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

- .id:

  Character string or integer. A character string giving the name or an
  integer of the position of the column of `.data` whose levels are used
  to split `.data` into groups. Defaults to `NULL`.

- .inference:

  Logical. Should critical values be computed? Defaults to `FALSE`.

- .independent:

  Character string. The name of the independent variable.

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

- .level:

  Character. Used in `plot.cSEMIPMA` to indicate whether IPMA should be
  done for constructs or indicators.

- .matrix1:

  A `matrix` to compare.

- .matrix2:

  A `matrix` to compare.

- .matrices:

  A list of at least two matrices.

- .metrics:

  Character string or a vector of character strings. Which prediction
  metrics should be displayed? One of: "*MAE*", "*RMSE*", "*Q2*",
  "*MER*", "*MAPE*, "*MSE2*", "*U1*", "*U2*", "*UM*", "*UR*", or "*UD*".
  Default to c("*MAE*", "*RMSE*", "*Q2*").

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) or a
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  list.

- .moderator:

  Character string. The name of the moderator variable.

- .modes:

  A vector giving the mode for each construct in the form
  `"name" = "mode"`. Only used internally.

- .ms_criterion:

  Character string. Either a single character string or a vector of
  character strings naming the model selection criterion to compute.
  Defaults to `"all"`.

- .n:

  Integer. The number of observations of the original data.

- .n_steps:

  Integer. A value giving the number of steps (the spotlights, i.e.,
  values of .moderator in surface analysis or floodlight analysis)
  between the minimum and maximum value of the moderator. Defaults to
  `100`.

- .normality:

  Logical. Should joint normality of \\\[\eta\_{1:p}; \zeta;
  \epsilon\]\\ be assumed in the nonlinear model? See (Dijkstra and
  Schermelleh-Engel 2014) for details. Defaults to `FALSE`. Ignored if
  the model is not nonlinear.

- .nr_comparisons:

  Integer. The number of comparisons. Defaults to `NULL`.

- .null_model:

  Logical. Should the degrees of freedom for the null model be computed?
  Defaults to `FALSE`.

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .object1:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .object2:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .only_common_factors:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- .only_structural:

  Should the the log-likelihood be based on the structural model?
  Ignored if `.by_equation == TRUE`. Defaults to `TRUE`.

- .original_arguments:

  The list of arguments used within
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .output_type:

  Character string. The type of output to return. One of "*complete*" or
  "*structured*". See the Value section for details. Defaults to
  "*complete*".

- .P:

  A (J x J) construct variance-covariance matrix (possibly
  disattenuated).

- .parameters_to_compare:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) indicating
  which parameters (i.e, path (`~`), loadings (`=~`), weights (`<~`), or
  correlations (`~~`)) should be compared across groups. Defaults to
  `NULL` in which case all weights, loadings and path coefficients of
  the originally specified model are compared.

- .path:

  Character string. Path of the directory to save the file to. Defaults
  to `NULL`.

- .path_coefficients:

  List. A list that contains the resampled and the original path
  coefficient estimates. Typically a part of a `cSEMResults_resampled`
  object. Defaults to `NULL`.

- .PLS_approach_cf:

  Character string. Approach used to obtain the correction factors for
  PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
  "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*",
  "*mean_harmonic*", "*geo_of_harmonic*". Defaults to
  "*dist_squared_euclid*". Ignored if `.disattenuate = FALSE` or if
  `.approach_weights` is not PLS-PM.

- .plot_correlations:

  Character string. Specify which correlations should be plotted, i.e.,
  between the exogenous constructs (`exo`), between the exogenous
  constructs and the indicators (`all`), or not at all (`none`).
  Defaults to `exo`.

- .plot_labels:

  Logical. Whether to display edge labels. Defaults to TRUE.

- .plot_package:

  Character string. Indicates which packages should be used for
  plotting.

- .plot_significances:

  Logical. Should p-values in the form of stars be plotted? Defaults to
  `TRUE`.

- .plot_structural_model_only:

  Logical. Should only the structural model, i.e., the constructs and
  their relationships be plotted? Defaults to `FALSE`.

- .plot_type:

  Character string. Indicates the type of plot that is produced.

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

- .probs:

  A vector of probabilities.

- .postestimation_object:

  An object resulting from a call to one of cSEM's postestimation
  functions (e.g.
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)).

- .quality_criterion:

  Character string. A single character string or a vector of character
  strings naming the quality criterion to compute. See the Details
  section for a list of possible candidates. Defaults to "*all*" in
  which case all possible quality criteria are computed.

- .quantity:

  Character string. Which statistic should be returned? One of "*all*",
  "*mean*", "*sd*", "*bias*", "*CI_standard_z*", "*CI_standard_t*",
  "*CI_percentile*", "*CI_basic*", "*CI_bc*", "*CI_bca*",
  "*CI_t_interval*" Defaults to "*all*" in which case all quantities
  that do not require additional resampling are returned, i.e., all
  quantities but "*CI_bca*", "*CI_t_interval*".

- .Q:

  A vector of composite-construct correlations with element names equal
  to the names of the J construct names used in the measurement model.
  Note Q^2 is also called the reliability coefficient.

- .reliabilities:

  A character vector of `"name" = value` pairs, where `value` is a
  number between 0 and 1 and `"name"` a character string of the
  corresponding construct name, or `NULL`. Reliabilities may be given
  for a subset of the constructs. Defaults to `NULL` in which case
  reliabilities are estimated by
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
  Currently, only supported for `.approach_weights = "PLS-PM"`.

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

- \`.resample_object\`:

  An R object of class `cSEMResults_resampled` obtained from
  [`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
  or by setting `.resample_method = "bootstrap"` or `"jackknife"` when
  calling
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .resample_sarstedt:

  A matrix containing the parameter estimates that could potentially be
  compared and an id column indicating the group adherence of each row.

- .r:

  Integer. The number of repetitions to use. Defaults to `1`.

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .R2:

  Integer. The number of bootstrap replications to use when resampling
  from a resample. Defaults to `199`.

- .R_bootstrap:

  Integer. The number of bootstrap runs. Ignored if `.object` contains
  resamples. Defaults to `499`

- .R_permutation:

  Integer. The number of permutations. Defaults to `499`

- .S:

  The (K x K) empirical indicator correlation matrix.

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .second_resample:

  A list containing `.R2` resamples for each of the `.R` resamples of
  the first run.

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

- .sim_points:

  Integer. How many samples from the truncated normal distribution
  should be simulated to estimate the exogenous construct scores?
  Defaults to "*100*".

- .stage:

  Character string. The stage the model is needed for. One of "*first*"
  or "*second*". Defaults to "*first*".

- .standardized:

  Logical. Should standardized scores be returned? Defaults to `TRUE`.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

- .steps_mod:

  A numeric vector. Steps used for the moderator variable in calculating
  the simple effects of an independent variable on the dependent
  variable. Defaults to `NULL`.

- .terms:

  A vector of construct names to be classified.

- .test_data:

  A matrix of test data with the same column names as the training data.

- .testtype:

  Character string. One of "*twosided*" (H1: The models do not perform
  equally in predicting indicators belonging to endogenous constructs)"
  and *onesided*" (H1: Model 1 performs better in predicting indicators
  belonging

- .title:

  Character string. Title of an object. Defaults to *""*.

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

- .treat_as_continuous:

  Logical. Should the indicators for the benchmark predictions be
  treated as continuous? If `TRUE` all indicators are treated as
  continuous and PLS-PM/PLSc is applied. If `FALSE` OrdPLS/OrdPLSc is
  applied. Defaults to `TRUE`.

- .type_gfi:

  Character string. Which fitting function should the GFI be based on?
  One of *"ML"* for the maximum likelihood fitting function, *"GLS"* for
  the generalized least squares fitting function or *"ULS"* for the
  unweighted least squares fitting function (same as the squared
  Euclidean distance). Defaults to *"ML"*.

- .type_ci:

  Character string. Which confidence interval should be calculated? For
  possible choices, see the `.quantity` argument of the
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  function. Only used if `.approch_mgd` is one of "*CI_para*" or
  "*CI_overlap*". Ignored otherwise. Defaults to "*CI_percentile*".

- .type_htmt:

  Character string indicating the type of HTMT that should be
  calculated, i.e., the original HTMT ("*htmt*") or the HTMT2
  ("*htmt2*"). Defaults to "*htmt*"

- .type_vcv:

  Character string. Which model-implied correlation matrix should be
  calculated? One of "*indicator*" or "*construct*". Defaults to
  "*indicator*".

- .verbose:

  Logical. Should information (e.g., progress bar) be printed to the
  console? Defaults to `TRUE`.

- .user_funs:

  A function or a (named) list of functions to apply to every resample.
  The functions must take `.object` as its first argument (e.g.,
  `myFun <- function(.object, ...) {body-of-the-function}`). Function
  output should preferably be a (named) vector but matrices are also
  accepted. However, the output will be vectorized (columnwise) in this
  case. See the examples section for details.

- .value_independent:

  Integer. Only required for floodlight analysis; The value of the
  independent variable in case that it appears as a higher-order term.

- .values_moderator:

  A numeric vector. The values of the moderator in a the simple effects
  analysis. Typically these are difference from the mean (=0) measured
  in standard deviations. Defaults to `c(-2, -1, 0, 1, 2)`.

- .vcv_asymptotic:

  Logical. Should the asymptotic variance-covariance matrix be used,
  i.e., VCV(b0) - VCV(b1)= VCV(b1-b0), or should VCV(b1-b0) be computed
  directly? Defaults to `FALSE`.

- .vector1:

  A vector of numeric values.

- .vector2:

  A vector of numeric values.

- .W:

  A (J x K) matrix of weights.

- .what:

  Internal argument used by several print helper functions.

- .W_new:

  A (J x K) matrix of weights.

- .W_old:

  A (J x K) matrix of weights.

- .weighted:

  Logical. Should estimation be based on a score that uses the weights
  of the weight approach used to obtain `.object`?. Defaults to `FALSE`.

- .X:

  A matrix of processed data (scaled, cleaned and ordered).

- .X_cleaned:

  A data.frame of processed data (cleaned and ordered). Note:
  `X_cleaned` may not be scaled!
