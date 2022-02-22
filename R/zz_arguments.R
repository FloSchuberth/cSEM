#' cSEMArguments
#'
#' An alphabetical list of all arguments used by functions of the `cSEM` package
#' including their description and defaults.
#' Mainly used for internal purposes (parameter inheritance). To list all arguments
#' and their defaults, use [args_default()]. To list all arguments and
#' their possible choices, use `args_default(.choices = TRUE)`.
#'
#' @param .alpha An integer or a numeric vector of significance levels. 
#'   Defaults to `0.05`.
#' @param .absolute Logical. Should the absolute HTMT values be returned? 
#'   Defaults to `TRUE` .
#' @param .approach_gcca Character string. The Kettenring approach to use for GCCA. One of 
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*". Defaults to
#' "*SUMCORR*".
#' @param .approach_2ndorder Character string. Approach used for models containing
#'   second-order constructs. One of: "*2stage*", or "*mixed*". Defaults to "*2stage*".
#' @param .approach_alpha_adjust Character string. Approach used to adjust the 
#'   significance level to accommodate multiple testing. 
#'   One of "*none*" or "*bonferroni*". Defaults to "*none*". 
#' @param .approach_cor_robust Character string. Approach used to obtain a robust 
#'   indicator correlation matrix. One of: "*none*" in which case the standard 
#'   Bravais-Person correlation is used,
#'   "*spearman*" for the Spearman rank correlation, or
#'   "*mcd*" via \code{\link[MASS:cov.rob]{MASS::cov.rob()}} for a robust correlation matrix. 
#'   Defaults to "*none*". Note that many postestimation procedures (such as
#'   [testOMF()] or [fit()] implicitly assume a continuous  
#'   indicator correlation matrix (e.g. Bravais-Pearson correlation matrix).
#'   Only use if you know what you are doing.
#' @param .approach_mgd Character string or a vector of character strings. 
#'   Approach used for the multi-group comparison. One of: "*all*", "*Klesel*", "*Chin*", 
#'   "*Sarstedt*", "*Keil*, "*Nitzl*", "*Henseler*", "*CI_para*", or "*CI_overlap*". 
#'   Default to "*all*" in which case all approaches are computed (if possible).
#' @param .approach_nl Character string. Approach used to estimate nonlinear
#'   structural relationships. One of: "*sequential*" or "*replace*".
#'   Defaults to "*sequential*".
#' @param .approach_p_adjust Character string or a vector of character strings. 
#' Approach used to adjust the p-value for multiple testing. 
#' See the `methods` argument of \code{\link[stats:p.adjust]{stats::p.adjust()}} for a list of choices and
#' their description. Defaults to "*none*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. One of: "*OLS*" or "*2SLS*". If "*2SLS*", instruments
#'   need to be supplied to `.instruments`. Defaults to "*OLS*".
#' @param .approach_score_benchmark Character string. How should the aggregation
#'   of the estimates of the truncated normal distribution be done for the 
#'   benchmark predictions? Ignored if not OrdPLS or OrdPLSc is used to obtain benchmark predictions.
#'   One of "*mean*" or "*median*" or "*mode*" or "*round*". 
#'   If "*round*", the benchmark predictions are obtained using the traditional prediction
#'   algorithm for PLS-PM which are rounded for categorical indicators.
#'   If "*mean*", the mean of the estimated endogenous indicators is calculated. 
#'   If "*median*", the mean of the estimated endogenous indicators is calculated.
#'   If "*mode*", the maximum empirical density on the intervals defined by the thresholds
#'   is used.
#'   If `.treat_as_continuous = TRUE` or if all indicators are on a continuous scale,
#'   `.approach_score_benchmark` is ignored. Defaults to "*round*".
#' @param .approach_score_target Character string. How should the aggregation of the estimates of
#'   the truncated normal distribution for the predictions using OrdPLS/OrdPLSc be done?
#'   One of "*mean*" or "*median*" or "*mode*". 
#'   If "*mean*", the mean of the estimated endogenous indicators is calculated. 
#'   If "*median*", the mean of the estimated endogenous indicators is calculated.
#'   If "*mode*", the maximum empirical density on the intervals defined by the thresholds
#'   is used. Defaults to "*mean*".
#' @param .approach_weights Character string. Approach used to
#'   obtain composite weights. One of: "*PLS-PM*", "*SUMCORR*", "*MAXVAR*",
#'   "*SSQCORR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", 
#'   or "*regression*". Defaults to "*PLS-PM*".
#' @param .args_used A list of function argument names whose value was modified 
#'   by the user.
#'   
#' @param .attrbutes Character string. Variables used as attributes in IPMA.
#' @param .benchmark Character string. The procedure to obtain benchmark predictions.
#'   One of "*lm*", "*unit*", "*PLS-PM*", "*GSCA*", "*PCA*", or "*MAXVAR*".
#'   Default to "*lm*".
#' @param .bias_corrected Logical. Should the standard and the tStat
#'   confidence interval be bias-corrected using the bootstrapped bias estimate? 
#'   If `TRUE` the confidence interval for some estimated parameter `theta` 
#'   is centered at `2*theta - theta*_hat`,
#'   where `theta*_hat` is the average over all `.R` bootstrap estimates of `theta`.
#'   Defaults to `TRUE`
#' @param .by_equation Should the criteria be computed for each structural model
#'   equation separately? Defaults to `TRUE`.
#' @param .C A (J x J) composite variance-covariance matrix.
#' @param .check_errors Logical. Should the model to parse be checked for correctness 
#'   in a sense that all necessary components to estimate the model are given?
#'   Defaults to `TRUE`.
#' @param .choices Logical. Should candidate values for the arguments be returned?
#'   Defaults to `FALSE`.
#' @param .ci A vector of character strings naming the confidence interval to compute.
#'   For possible choices see [infer()].
#' @param .ci_colnames Internal argument used by several print helper functions.
#' @param .closed_form_ci Logical. Should a closed-form confidence interval be computed?
#'   Defaults to `FALSE`.
#' @param .conv_criterion Character string. The criterion to use for the convergence check.
#'   One of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*". Defaults
#'   to "*diff_absolute*".
#' @param .csem_model A (possibly incomplete) [cSEMModel]-list.
#' @param .csem_resample A list resulting from a call to [resamplecSEMResults()].
#' @param .cv_folds Integer. The number of cross-validation folds to use. Setting
#'   `.cv_folds` to `N` (the number of observations) produces 
#'   leave-one-out cross-validation samples. Defaults to `10`.
#' @param .data A `data.frame` or a `matrix` of standardized or unstandardized  
#'   data (indicators/items/manifest variables). Possible column types or classes 
#'   of the data provided are: "`logical`", "`numeric`" ("`double`" or "`integer`"), 
#'   "`factor`" ("`ordered`" and/or "`unordered`"), "`character`" (converted to factor),
#'   or a mix of several types.
#' @param .dependent Character string. The name of the dependent variable.
#' @param .disattenuate Logical. Should composite/proxy correlations 
#'   be disattenuated to yield consistent loadings and path estimates if at least
#'   one of the construct is modeled as a common factor? Defaults to `TRUE`.
#' @param .dist Character string. The distribution to use for the critical value.
#'  One of *"t"* for Student's t-distribution or *"z"* for the standard normal distribution.
#'  Defaults to *"z"*.
#' @param .distance Character string. A distance measure. One of: "*geodesic*"
#'   or "*squared_euclidean*". Defaults to "*geodesic*".
#' @param .df Character string. The method for obtaining the degrees of freedom.
#'   Choices are "*type1*" and "*type2*". Defaults to "*type1*" .
#' @param .dominant_indicators A character vector of `"construct_name" = "indicator_name"` pairs, 
#'   where `"indicator_name"` is a character string giving the name of the dominant indicator
#'   and `"construct_name"` a character string of the corresponding construct name.
#'   Dominant indicators may be specified for a subset of the constructs. 
#'   Default to `NULL`.
#' @param .E A (J x J) matrix of inner weights.
#' @param .effect Internal argument used by helper printEffects().
#' @param .estimate_structural Logical. Should the structural coefficients
#'   be estimated? Defaults to `TRUE`.
#' @param .eval_plan Character string. The evaluation plan to use. One of 
#'   "*sequential*", "*multiprocess*", or "*multisession*". In the two latter cases 
#'   all available cores will be used. Defaults to "*sequential*".
#' @param .filename Character string. The file name. Defaults to "results.xlsx".
#' @param .first_resample A list containing the `.R` resamples based on the original
#'   data obtained by resamplecSEMResults().
#' @param .fit_measures Logical. (EXPERIMENTAL) Should additional fit measures 
#'   be included? Defaults to `FALSE`. 
#' @param .force Logical. Should .object be resampled even if it contains resamples
#'   already?. Defaults to `FALSE`.
#' @param .full_output Logical. Should the full output of summarize be printed.
#'   Defaults to `TRUE`.
#' @param .H The (N x J) matrix of construct scores.
#' @param .handle_inadmissibles Character string. How should inadmissible results 
#'   be treated? One of "*drop*", "*ignore*", or "*replace*". If "*drop*", all
#'   replications/resamples yielding an inadmissible result will be dropped 
#'   (i.e. the number of results returned will potentially be less than `.R`). 
#'   For "*ignore*" all results are returned even if all or some of the replications
#'   yielded inadmissible results (i.e. number of results returned is equal to `.R`). 
#'   For "*replace*" resampling continues until there are exactly `.R` admissible solutions.
#'   Depending on the frequency of inadmissible solutions this may significantly increase
#'   computing time. Defaults to "*drop*".
#' @param .id Character string or integer. A character string giving the name or 
#'   an integer of the position of the column of `.data` whose levels are used
#'   to split `.data` into groups. Defaults to `NULL`.
#' @param .inference Logical. Should critical values be computed? Defaults to `FALSE`.
#' @param .independent Character string. The name of the independent variable.
#' @param .instruments A named list of vectors of instruments. The names
#'   of the list elements are the names of the dependent (LHS) constructs of the structural
#'   equation whose explanatory variables are endogenous. The vectors
#'   contain the names of the instruments corresponding to each equation. Note
#'   that exogenous variables of a given equation **must** be supplied as 
#'   instruments for themselves. Defaults to `NULL`.
#' @param .iter_max Integer. The maximum number of iterations allowed.
#'   If `iter_max = 1` and `.approach_weights = "PLS-PM"` one-step weights are returned. 
#'   If the algorithm exceeds the specified number, weights of iteration step 
#'   `.iter_max - 1`  will be returned with a warning. Defaults to `100`.
#' @param .level Character. Used in `plot.cSEMIPMA` to indicate whether IPMA should be done 
#' for constructs or indicators.    
#' @param .matrix1 A `matrix` to compare.
#' @param .matrix2 A `matrix` to compare.
#' @param .matrices A list of at least two matrices.
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] 
#'   or a [cSEMModel] list.
#' @param .moderator Character string. The name of the moderator variable.
#' @param .modes A vector giving the mode for each construct in the form `"name" = "mode"`. 
#'   Only used internally. 
#' @param .ms_criterion Character string. Either a single character string or a vector
#'   of character strings naming the model selection criterion to compute.
#'   Defaults to `"all"`.
#' @param .n Integer. The number of observations of the original data.
#' @param .n_steps Integer. A value giving the number of steps (the spotlights, i.e.,
#' values of .moderator in surface analysis or floodlight analysis) 
#' between the minimum and maximum value of the moderator. Defaults to `100`.
#' @param .normality Logical. Should joint normality of 
#' \eqn{[\eta_{1:p}; \zeta; \epsilon]}{[\eta_(1:p); \zeta; \epsilon]}
#'  be assumed in the nonlinear model? See \insertCite{Dijkstra2014}{cSEM} for details.
#'  Defaults to `FALSE`. Ignored if the model is not nonlinear.
#' @param .nr_comparisons Integer. The number of comparisons. Defaults to `NULL`.  
#' @param .null_model Logical. Should the degrees of freedom for the null model
#'   be computed? Defaults to `FALSE`.
#' @param .object An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param .only_common_factors Logical. Should only concepts modeled as common 
#'   factors be included when calculating one of the following quality critera: 
#'   AVE, the Fornell-Larcker criterion, HTMT, and all reliability estimates. 
#'   Defaults to `TRUE`.
#' @param .only_structural Should the the log-likelihood be based on the 
#'   structural model? Ignored if `.by_equation == TRUE`. Defaults to `TRUE`.
#' @param .original_arguments The list of arguments used within [csem()].
#' @param .output_type Character string. The type of output to return. One of
#'   "*complete*" or "*structured*". See the Value section for details. Defaults to
#'   "*complete*".
#' @param .P A (J x J) construct variance-covariance matrix (possibly disattenuated).
#' @param .parameters_to_compare A model in [lavaan model syntax][lavaan::model.syntax] indicating which 
#'   parameters (i.e, path (`~`), loadings (`=~`), weights (`<~`), or correlations (`~~`)) should be
#'   compared across groups. Defaults to `NULL` in which case all weights, loadings and 
#'   path coefficients of the originally specified model are compared.
#' @param .path Character string. Path of the directory to save the file to. Defaults
#'   to the current working directory.
#' @param .path_coefficients List. A list that contains the resampled and the original 
#' path coefficient estimates. Typically a part of a `cSEMResults_resampled` object.
#' Defaults to `NULL`. 
#' @param .PLS_approach_cf Character string. Approach used to obtain the correction
#'   factors for PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
#'   "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*", "*mean_harmonic*",
#'   "*geo_of_harmonic*". Defaults to "*dist_squared_euclid*". 
#'   Ignored if `.disattenuate = FALSE` or if `.approach_weights` is not PLS-PM.
#' @param .plot_package Character string. Indicates which packages should be used for plotting.
#' @param .plot_type Character string. Indicates the type of plot that is produced.
#' @param .PLS_ignore_structural_model Logical. Should the structural model be ignored
#'   when calculating the inner weights of the PLS-PM algorithm? Defaults to `FALSE`.
#'   Ignored if `.approach_weights` is not PLS-PM.
#' @param .PLS_modes Either a named list specifying the mode that should be used for
#'   each construct in the form `"construct_name" = mode`, a single character
#'   string giving the mode that should be used for all constructs, or `NULL`.
#'   Possible choices for `mode` are: "*modeA*", "*modeB*", "*modeBNNLS*", 
#'   "*unit*", "*PCA*", a single integer or 
#'   a vector of fixed weights of the same length as there are indicators for the
#'   construct given by `"construct_name"`. If only a single number is provided this is identical to
#'   using unit weights, as weights are rescaled such that the related composite 
#'   has unit variance.  Defaults to `NULL`.
#'   If `NULL` the appropriate mode according to the type
#'   of construct used is chosen. Ignored if `.approach_weight` is not PLS-PM.  
#' @param .PLS_weight_scheme_inner Character string. The inner weighting scheme
#'   used by PLS-PM. One of: "*centroid*", "*factorial*", or "*path*".
#'   Defaults to "*path*". Ignored if `.approach_weight` is not PLS-PM.
#' @param .probs A vector of probabilities.
#' @param .postestimation_object An object resulting from a call to one of cSEM's
#'   postestimation functions (e.g. [summarize()]).
#' @param .quality_criterion Character string. A single character string or a
#'   vector of character strings naming the quality criterion to compute. See 
#'   the Details section for a list of possible candidates. 
#'   Defaults to "*all*" in which case all possible quality criteria are computed.
#' @param .quantity Character string. Which statistic should be returned?
#'   One of "*all*", "*mean*", "*sd*", "*bias*", "*CI_standard_z*", "*CI_standard_t*",
#'   "*CI_percentile*", "*CI_basic*", "*CI_bc*", "*CI_bca*", "*CI_t_interval*"
#'   Defaults to "*all*" in which case all quantities that do not require
#'   additional resampling are returned, i.e., all quantities but "*CI_bca*", "*CI_t_interval*". 
#' @param .Q A vector of composite-construct correlations with element names equal to
#'   the names of the J construct names used in the measurement model. Note 
#'   Q^2 is also called the reliability coefficient.
#' @param .reliabilities A character vector of `"name" = value` pairs, 
#'   where `value` is a number between 0 and 1 and `"name"` a character string
#'   of the corresponding construct name, or `NULL`. Reliabilities
#'   may be given for a subset of the constructs. Defaults to `NULL` in which case
#'   reliabilities are estimated by `csem()`. Currently, only supported for
#'   `.approach_weights = "PLS-PM"`.
#' @param .resample_method Character string. The resampling method to use. One of: 
#'  "*none*", "*bootstrap*" or "*jackknife*". Defaults to "*none*".
#' @param .resample_method2 Character string. The resampling method to use when resampling
#'   from a resample. One of: "*none*", "*bootstrap*" or "*jackknife*". For 
#'   "*bootstrap*" the number of draws is provided via `.R2`. Currently, 
#'   resampling from each resample is only required for the studentized confidence
#'   intervall ("*CI_t_interval*") computed by the [infer()] function. Defaults to "*none*".  
#' @param `.resample_object` An R object of class `cSEMResults_resampled`
#'   obtained from [resamplecSEMResults()] or by setting `.resample_method = "bootstrap"`
#'   or `"jackknife"` when calling [csem()].
#' @param .resample_sarstedt A matrix containing the parameter estimates that 
#'   could potentially be compared and an id column indicating the group adherence
#'   of each row.
#' @param .r Integer. The number of repetitions to use. Defaults to `1`.
#' @param .R Integer. The number of bootstrap replications. Defaults to `499`.
#' @param .R2 Integer. The number of bootstrap replications to use when 
#'   resampling from a resample. Defaults to `199`.
#' @param .R_bootstrap Integer. The number of bootstrap runs. Ignored if `.object`
#'   contains resamples. Defaults to `499`
#' @param .R_permutation Integer. The number of permutations. Defaults to `499`
#' @param .S The (K x K) empirical indicator correlation matrix.
#' @param .saturated Logical. Should a saturated structural model be used? 
#'   Defaults to `FALSE`.
#' @param .second_resample A list containing `.R2` resamples for each of the `.R`
#'   resamples of the first run.
#' @param .seed Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
#'   case an arbitrary seed is chosen. Note that the scope of the seed is limited
#'   to the body of the function it is used in. Hence, the global seed will
#'   not be altered!
#' @param .sign_change_option Character string. Which sign change option should 
#' be used to handle flipping signs when resampling? One of "*none*","*individual*",
#' "*individual_reestimate*", "*construct_reestimate*". Defaults to "*none*".
#' @param .sim_points Integer. How many samples from the truncated normal distribution should
#'   be simulated to estimate the exogenous construct scores? Defaults to "*100*".
#' @param .stage Character string. The stage the model is needed for.
#'   One of "*first*" or "*second*". Defaults to "*first*".
#' @param .standardized Logical. Should standardized scores be returned? Defaults
#'   to `TRUE`.
#' @param .starting_values A named list of vectors where the
#'   list names are the construct names whose indicator weights the user
#'   wishes to set. The vectors must be named vectors of `"indicator_name" = value` 
#'   pairs, where `value` is the (scaled or unscaled) starting weight. Defaults to `NULL`.
#' @param .steps_mod A numeric vector. Steps used for the moderator variable in calculating 
#' the simple effects of an independent variable on the dependent variable. 
#' Defaults to `NULL`.
#' @param .terms A vector of construct names to be classified.
#' @param .test_data A matrix of test data with the same column names as the 
#'   training data.
#' @param .tolerance Double. The tolerance criterion for convergence. 
#'   Defaults to `1e-05`.
#' @param .treat_as_continuous Logical. Should the indicators for the benchmark predictions
#'   be treated as continuous? If `TRUE` all indicators are treated as continuous and PLS-PM/PLSc is applied. 
#'   If `FALSE` OrdPLS/OrdPLSc is applied. Defaults to `TRUE`.
#' @param .type_gfi Character string. Which fitting function should the GFI be based 
#'   on? One of *"ML"* for the maximum likelihood fitting function, *"GLS"* for 
#'   the generalized least squares fitting function or *"ULS"* for the
#'   unweighted least squares fitting function (same as the squared Euclidean distance). 
#'   Defaults to *"ML"*.
#' @param .type_ci Character string. Which confidence interval should be calculated? 
#'   For possible choices, see the `.quantity` argument of the [infer()] function.
#'   Only used if `.approch_mgd` is one of "*CI_para*" or "*CI_overlap*". Ignored otherwise.
#'   Defaults to "*CI_percentile*". 
#' @param .type_htmt Character string indicating the type of HTMT that should be 
#' calculated, i.e., the original HTMT ("*htmt*") or the HTMT2 ("*htmt2*").
#' Defaults to "*htmt*"
#' @param .type_vcv Character string. Which model-implied correlation 
#'  matrix should be calculated?
#'  One of "*indicator*" or "*construct*". Defaults to "*indicator*".   
#' @param .verbose Logical. Should information (e.g., progress bar) be printed 
#'   to the console? Defaults to `TRUE`.
#' @param .user_funs A function or a (named) list of functions to apply to every
#'   resample. The functions must take `.object` as its first argument (e.g., 
#'   `myFun <- function(.object, ...) {body-of-the-function}`).
#'   Function output should preferably be a (named)
#'   vector but matrices are also accepted. However, the output will be 
#'   vectorized (columnwise) in this case. See the examples section for details.
#' @param .value_independent Integer. Only required for floodlight analysis;
#' The value of the independent variable in case that it appears as a 
#' higher-order term.
#' @param .values_moderator A numeric vector. The values of the moderator in a 
#' the simple effects analysis. Typically these are difference from the mean (=0) 
#' measured in standard deviations. Defaults to `c(-2, -1, 0, 1, 2)`.
#' @param .vcv_asymptotic Logical. Should the asymptotic variance-covariance matrix be used, i.e., 
#' VCV(b0) - VCV(b1)= VCV(b1-b0), or should VCV(b1-b0) be computed directly? 
#'  Defaults to `FALSE`.
#' @param .vector1 A vector of numeric values.
#' @param .vector2 A vector of numeric values.
#' @param .W A (J x K) matrix of weights.
#' @param .what Internal argument used by several print helper functions.
#' @param .W_new A (J x K) matrix of weights.
#' @param .W_old A (J x K) matrix of weights.
#' @param .weighted Logical. Should estimation be based on a score that uses 
#'   the weights of the weight approach used to obtain `.object`?. Defaults to `FALSE`.
#' @param .X A matrix of processed data (scaled, cleaned and ordered).
#' @param .X_cleaned A data.frame of processed data (cleaned and ordered). Note: `X_cleaned`
#'   may not be scaled!
#'
#' @name csem_arguments
#' @aliases cSEMArguments
#' @keywords internal
NULL

#' Complete list of assess()'s ... arguments
#' 
#' A complete alphabetical list of all possible arguments accepted by `assess()`'s `...` 
#' (dotdotdot) argument.
#' 
#' Most arguments supplied to the `...` argument of `assess()` are only
#' accepted by a subset of the functions called by `assess()`. The following
#' list shows which argument is passed to which function:
#' \describe{
#' \item{.absolute}{Accepted by/Passed down to: [calculateHTMT()]}
#' \item{.alpha}{Accepted by/Passed down to: [calculateRhoT()], [calculateHTMT()], [calculateCN()]}
#' \item{.ci}{Accepted by/Passed down to: [calculateHTMT()]}
#' \item{.closed_form_ci}{Accepted by/Passed down to: [calculateRhoT()]}
#' \item{.handle_inadmissibles}{Accepted by/Passed down to: [calculateHTMT()]}
#' \item{.inference}{Accepted by/Passed down to: [calculateHTMT]}
#' \item{.null_model}{Accepted by/Passed down to: [calculateDf()]}
#' \item{.R}{Accepted by/Passed down to: [calculateHTMT()]}
#' \item{.saturated}{Accepted by/Passed down to: [calculateSRMR()], 
#'   [calculateDG()], [calculateDL()], [calculateDML()]and subsequently [fit()].}
#' \item{.seed}{Accepted by/Passed down to: [calculateHTMT()]}
#' \item{.type_gfi}{Accepted by/Passed down to: [calculateGFI()]}
#' \item{.type_vcv}{Accepted by/Passed down to: [calculateSRMR()], 
#'   [calculateDG()], [calculateDL()], [calculateDML()] and subsequently [fit()].}
#' }
#' @usage NULL
#' 
#' @inheritParams csem_arguments
#' @keywords internal

args_assess_dotdotdot <- function(
  .absolute            = TRUE,
  .alpha               = 0.05,
  .ci                  = c("CI_standard_z", "CI_standard_t", "CI_percentile", 
                           "CI_basic", "CI_bc", "CI_bca", "CI_t_interval"),
  .closed_form_ci      = FALSE,
  .handle_inadmissibles= c("drop", "ignore", "replace"),
  .inference           = FALSE,
  .null_model          = FALSE,
  .R                   = 499,
  .saturated           = FALSE,
  .seed                = NULL,
  .type_gfi            = c("ML", "GLS", "ULS"),
  .type_vcv            = "indicator"
) {NULL}
  
#' Show argument defaults or candidates
#'
#' Show all arguments used by package functions including default or candidate
#' values. For argument descriptions see: [csem_arguments].
#'
#' By default `args_default()`returns a list of default values by argument name.
#' If the list of accepted candidate values is required instead, use `.choices = TRUE`.
#'
#' @usage args_default(.choices = FALSE)
#' 
#' @inheritParams  csem_arguments
#' 
#' @return A named list of argument names and defaults or accepted candidates.
#'
#' @seealso [handleArgs()], [csem_arguments], [csem()], [foreman()]
#'
#' @export

args_default <- function(.choices = FALSE) {
  
  args <- list(
    .alpha                   = 0.05,
    .absolute                = TRUE,
    .approach_2ndorder       = c("2stage", "mixed"),
    .approach_alpha_adjust   = c("none", "bonferroni"),
    .approach_cor_robust     = c("none", "mcd", "spearman"),
    .approach_gcca           = c("SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR"),
    .approach_mgd            = c("all", "Klesel", "Chin", "Sarstedt", "Keil", "Nitzl", 
                                 "Henseler","CI_para","CI_overlap"),
    .approach_nl             = c("sequential", "replace"),
    .approach_p_adjust       = "none",
    .approach_paths          = c("OLS", "2SLS"),
    .approach_score_benchmark= c("mean", "median", "mode", "round"),
    .approach_score_target   = c("mean", "median", "mode"),
    .approach_weights        = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                                 "GSCA", "PCA", "unit", "bartlett", "regression"), 
    .arguments               = NULL,
    .attributes              = NULL,
    .benchmark               = c("lm", "unit", "PLS-PM", "GSCA", "PCA", "MAXVAR"),
    .bias_corrected          = TRUE,
    .by_equation             = TRUE,
    .C                       = NULL,
    .check_errors            = TRUE,
    .choices                 = FALSE,
    .ci                      = c("CI_standard_z", "CI_standard_t", "CI_percentile", 
                                 "CI_basic", "CI_bc", "CI_bca", "CI_t_interval"),
    .ci_colnames             = NULL,
    .closed_form_ci          = FALSE, 
    .conv_criterion          = c("diff_absolute", "diff_squared", "diff_relative"),
    .csem_model              = NULL,
    .csem_resample           = NULL,
    .cv_folds                = 10,
    .data                    = NULL,
    .dependent               = NULL,
    .disattenuate            = TRUE,
    .dist                    = c("z", "t"),
    .distance                = c("geodesic", "squared_euclidean"),
    .df                      = c("type1", "type2"),
    .dominant_indicators     = NULL,
    .E                       = NULL,
    .effect                  = NULL,
    .estimate_structural     = TRUE,
    .eval_plan               = c("sequential", "multiprocess","multisession"),
    .filename                = "results.xlsx",
    .fit_measures            = FALSE,
    .first_resample          = NULL,
    .force                   = FALSE,
    .full_output             = TRUE,
    .handle_inadmissibles    = c("drop", "ignore", "replace"),
    .H                       = NULL,
    .id                      = NULL,
    .inference               = FALSE,
    .independent             = NULL,
    .instruments             = NULL,
    .iter_max                = 100,
    .listMatrices            = NULL, 
    .level                   = c('construct',"indicator"),
    .matrix1                 = NULL,
    .matrix2                 = NULL,
    .matrices                = NULL,
    .model                   = NULL,
    .moderator               = NULL,
    .modes                   = NULL,
    .ms_criterion            = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                                 "hqc", "mallows_cp"),
    .n_steps                 = 100,
    .normality               = FALSE,
    .nr_comparisons          = NULL,
    .null_model              = FALSE,
    .object                  = NULL,
    .only_common_factors     = TRUE,
    .output_type             = c("complete", "structured"),
    .P                       = NULL,
    .parameters_to_compare   = NULL,
    .path                    = NULL,
    .path_coefficients       = NULL,
    .plot_package            = NULL,
    .plot_type               = NULL,
    .probs                   = NULL,
    .PLS_approach_cf         = c("dist_squared_euclid", "dist_euclid_weighted", 
                                 "fisher_transformed", "mean_arithmetic",
                                 "mean_geometric", "mean_harmonic",
                                 "geo_of_harmonic"),
    .PLS_ignore_structural_model = FALSE,
    .PLS_modes               = NULL,
    .PLS_weight_scheme_inner = c("path", "centroid", "factorial"),
    .postestimation_object   = NULL,
    .quality_criterion       = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
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
    
    .quantity                = c("all", "mean", "sd", "bias", "CI_standard_z", "CI_standard_t",
                                 "CI_percentile", "CI_basic", "CI_bc", "CI_bca", "CI_t_intervall"),
    .Q                       = NULL,
    .r                       = 1,
    .R                       = 499,
    .R2                      = 199,
    .R_bootstrap             = 499,
    .R_permutation           = 499,
    .reliabilities           = NULL,
    .resample_method         = c("none", "bootstrap", "jackknife"),
    .resample_method2        = c("none", "bootstrap", "jackknife"),
    .resample_object         = NULL,
    .resample_sarstedt       = NULL,
    .S                       = NULL,
    .saturated               = FALSE,
    .second_resample         = NULL,
    .seed                    = NULL,
    .sign_change_option      = c("none", "individual", "individual_reestimate",
                                 "construct_reestimate"),
    .sim_points              = 100,
    .stage                   = c("first", "second"),
    .standardized            = TRUE,
    .starting_values         = NULL,
    .steps_mod               = NULL,
    .terms                   = NULL,
    .test_data               = NULL,
    .tolerance               = 1e-05,
    .treat_as_continuous     = TRUE,
    .type_gfi                = c("ML", "GLS", "ULS"),
    .type_htmt               = c("htmt", "htmt2"),
    .type_vcv                = c("indicator", "construct"),
    .type_ci                 = c("CI_percentile","CI_standard_z","CI_standard_t",
                                 "CI_basic","CI_bc", "CI_bca"),
    .user_funs               = NULL,
    .value_independent       = 0,
    .values_moderator        = c(-2,-1,0,1,2),
    .vcv_asymptotic          = c(FALSE, TRUE),
    .verbose                 = TRUE,
    .W                       = NULL,
    .weighted                = FALSE,
    .what                    = NULL,
    .x                       = NULL,
    .X                       = NULL,
    .X_cleaned               = NULL,
    .y                       = NULL,
    .z                       = NULL
  )

  if(!.choices) {
    args <- lapply(args, function(x) eval(x)[1])
  }
    
  args_sorted <- args[sort(names(args))]
  
  return(args_sorted)
}

#' Internal: Handle arguments
#'
#' Internal helper function to handle arguments passed to any function within `cSEM`.
#'
#' @param .args_used A list of argument names and user picked values.
#'
#' @return The [args_default] list, with default values
#' changed to the values given by the user.
#' @keywords internal

handleArgs <- function(.args_used) {
  
  args_default  <- args_default()
  
  args_default_names  <- names(args_default)
  args_used_names <- names(.args_used)
  
  ## Are all argument names used valid?
  args_diff <- setdiff(args_used_names, args_default_names)
  
  if (length(args_diff) > 0) {
    stop("The following argument(s) are not known to cSEM: ",
         paste0("`", args_diff, '`', collapse = ", "), "\n",
         "Please type 'args_default()' to find all valid arguments and their defaults. Or check ?args_default.",
         call. = FALSE)
  }
  
  ## Are arguments valid
  # Note: for now, I will only check if string arguments are valid. In the
  # future we could refine and check if e.g. numeric values are within a reasonable
  # range or larger than zero if required etc.
  # choices_logical <- Filter(function(x) any(is.logical(x)), args_default(.choices = TRUE))
  # choices_numeric <- Filter(function(x) any(is.numeric(x)), args_default(.choices = TRUE))
  choices_character <- Filter(function(x) any(is.character(x)), args_default(.choices = TRUE))

  character_args <- intersect(names(choices_character), args_used_names)
  x <- Map(function(x, y) x %in% y, 
      x = .args_used[character_args], 
      y = choices_character[character_args]
      )

  lapply(seq_along(x), function(i) {
    if(isFALSE(x[[i]])) {
      n <- names(x[i])
      a <- args_default(.choices = TRUE)[[n]]
      stop(paste0("`", .args_used[n], "` is not a valid choice for ", "`", 
                  names(.args_used[n]),"`.\n"), 
           "Choices are: ", paste0("`", a[-length(a)],"`", collapse = ", "), 
           " or " , paste0("`", a[length(a)], "`"), call. = FALSE)
    }

  })
  ## Replace all arguments that were changed or explicitly given and keep
  #  the default values for the others
  
  for(i in args_used_names) {
    args_default[i] <- .args_used[i]
  }
  
  return(args_default)
}
