#' cSEMModel
#'
#' A standardized list containing model-related information. To convert a
#' a model written in \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#' to a `cSEMModel`-list use [parseModel()].
#'
#' An object of class `cSEMModel` is a standardized list containing the following
#' components (J constructs and K indicators)
#' \describe{
#'   \item{`$structural`}{A matrix mimicking the structural relationship between
#'      constructs. If constructs are only linearly related, `structural` is
#'      of dimension (J x J) with row- and column names equal to the construct
#'      names. If the structural model contains non-linear relationships
#'      `structural` is (J x (J + J\*)) where J\* is the number of
#'      nonlinear terms. Rows are ordered such that exogenous constructs are always
#'      first, followed by constructs that only depend on exogenous constructs and/or
#'      previously ordered constructs.}
#'   \item{`$measurement`}{A (J x K) matrix mimicking the measurement relationship
#'     between constructs and their related indicators. Rows are are in the same
#'     order as the `$structural` with row names equal to
#'     the construct names. The order of the columns is such that `$.measurement`
#'     forms a block diagonal matrix.}
#'   \item{`$error_cor`}{A (K x K) matrix mimicking the measurement error
#'     correlation relationship. The row and column order is identical to `$measurement`.}
#'   \item{`$construct_type`}{A named vector containing the names of each construct
#'    and their respective type (**"Common factor"** or **"Composite"**).}
#'   \item{`$model_type`}{The type of model (linear or nonlinear).}
  # \item{`$vars_endo`}{A vector of names of the endogenous constructs.}
  # \item{`$vars_exo`}{A vector of names of the exogenous constructs (includes
  #   possible interaction and exponential terms).}
  # \item{`$vars_explana`}{ A vector of names of the constructs that appear as
  #   explanatory variables in at least one structural equation (includes
  #   possible interaction and exponential terms).}
  # \item{`$explained_by_exo`}{A vector of names of the constructs that are
  #   solely explained by exogenous constructs.}
#' }
#' Note: it is possible to supply an incomplete `cSEMModel`-list
#' to all functions that require `.csem_model` as a mandatory argument. Currently,
#' only the structural and the measurement matrix are required.
#' However, specifying an incomplete cSEMModel list may lead to unexpected behavior 
#' and errors so do use this technique with caution.
#'
#' @seealso [parseModel]
#' @name csem_model
#' @aliases cSEMModel
#' @keywords internal
NULL

#' cSEMResults
#'
#' @return
#' An object of class `cSEMResults` with methods for all postestimation generics:
#' \itemize{
#'   \item{`assess.cSEMResults`   Asses results using common fit measures and indices.}
#'   \item{`summarize.cSEMResults`   Summarize the results.}
#'   \item{`test.cSEMResults`   Run tests.}
#'   \item{`verify.cSEMResults`   Verify admissibility (e.g., loadings larger than one).}
#' }
#' 
#' The structure and class of the resulting `cSEMResults` object is generally 
#' determined by the type of dataset(s) provided (`matrix`, `data.frame`, `list` 
#' of data)
#' 
#' \describe{
#' \item{Default}{If data is a single `matrix` or `data.frame` with no id-column, 
#' the result is a `list` with elements: 
#' \describe{
#'   \item{`$Estimates`}{A list containing a list of estimated quantities.}
#'   \item{`$Information`}{A list containing a list of additional information.}
#' }
#' The resulting object has classes `cSEMResults` and `cSEMResults_default`.
#' }
#' \item{Multi}{If the data provided contains an id-column to split the data by `G` group levels 
#' or if a list of `G` datasets is provided, the resulting object is a list of `G` 
#' lists, where `G` is equal to the number of groups or the number of datasets 
#' in the list of datasets provided. Each of the `G` list elements are itself 
#' a `cSEMResults_default` object. Hence its structure is identical to 
#' the structure described in `Default`.
#' 
#' The resulting object has classes `cSEMResults` and `cSEMResults_multi`. 
#' }
#' \item{2ndorder}{
#' A special output is generated if the model to estimate contains hierachical constructs
#' **and** the 2step (3step) approach is used to estimate the model. In this case
#' the resulting object is a list containing two elements `First_stage` and 
#' `Second_stage`.
#' 
#' Each list element is itself a `cSEMResults_default` object. Hence its structure is identical to 
#' the structure described in `Default`.
#' }
#' }
#' 
#' If `.resample_method = "bootstrap"` or `.resample_method = "jackknife"` resamples
#' are attached to each object. For objects of class `cSEMResults_default` the resamples are
#' attached to `.object$Estimates$Estimates_resamples`. For objects of type
#' `cSEMResults_multi` the same is done by group. All objects containing 
#' these elements gain the `cSEMResults_resampled` class.
#' 
#' As should be evident from the above, the core structure of the resulting object
#' is by and large similiar across all cases. Practically, most users wont have to worry
#' about the structure as postestimation functions have corresponding methods to
#' accomodate all cases.
#' 
#' @name csem_results
#' @aliases cSEMResults
#' @keywords internal
NULL

#' cSEMSummarize
#'
#' @return
#' An object of class `cSEMSummary`.
#' Technically `cSEMSummary` is a named list containing the following list elements:
#' \describe{
#'   \item{`...}{Not finished yet.}
#' }
#'
#' @name csem_summary
#' @aliases cSEMSummary
#' @keywords internal
NULL

#' cSEMTest
#'
#' @return
#' A standardized list of class `cSEMTest`. Technically `cSEMTest` is a named 
#' list containing the following list elements:
#' \describe{
#'   \item{`$Test_statistic`}{The value of test statistic(s).}
#'   \item{`$Critical_value`}{The critical value(s).}
#'   \item{`$Decision`}{The test decision. One of: **Reject** or **Do not reject**}
#'   \item{`$Number_admissibles`}{The number of admissible runs. See [verify()] for
#'     what constitues and inadmissible run.}
#' }
#'
#' @name csem_test
#' @aliases cSEMTest
#' @keywords internal
NULL