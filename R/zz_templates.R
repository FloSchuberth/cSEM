#' cSEMModel
#'
#' @details A standardized list containing model-related information. To convert a
#' a model written in [lavaan model syntax][lavaan::model.syntax] 
#' to a [cSEMModel] list use [parseModel()].
#'
#' @return An object of class `cSEMModel` is a standardized list containing the 
#' following components. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$structural`}{A matrix mimicking the structural relationship between
#'     constructs. If constructs are only linearly related, `structural` is
#'     of dimension (J x J) with row- and column names equal to the construct
#'     names. If the structural model contains nonlinear relationships
#'     `structural` is (J x (J + J\*)) where J\* is the number of
#'     nonlinear terms. Rows are ordered such that exogenous constructs are always
#'     first, followed by constructs that only depend on exogenous constructs and/or
#'     previously ordered constructs.}
#'   \item{`$measurement`}{A (J x K) matrix mimicking the measurement relationship
#'     between constructs and their related indicators. Rows are in the same
#'     order as the matrix `$structural` with row names equal to
#'     the construct names. The order of the columns is such that `$measurement`
#'     forms a block diagonal matrix.}
#'   \item{`$error_cor`}{A (K x K) matrix mimicking the measurement error
#'     correlation relationship. The row and column order is identical to 
#'     the column order of `$measurement`.}
#'   \item{`$construct_type`}{A named vector containing the names of each construct
#'     and their respective type ("Common factor" or "Composite").}
#'   \item{`$construct_order`}{A named vector containing the names of each construct
#'     and their respective order ("First order" or "Second order").}
#'   \item{`$model_type`}{The type of model ("Linear" or "Nonlinear").}
#'   \item{`$instruments`}{Only if instruments are supplied: a list of structural 
#'     equations relating endogenous RHS variables to instruments.}
#' }
#' It is possible to supply an incomplete list to `parseModel()`, resulting
#' in an incomplete `cSEMModel` list which can be passed
#' to all functions that require `.csem_model` as a mandatory argument. Currently,
#' only the structural and the measurement matrix are required.
#' However, specifying an incomplete `cSEMModel` list may lead to unexpected behavior 
#' and errors. Use with care.
#'
#' @seealso [parseModel]
#' @name csem_model
#' @aliases cSEMModel
#' @keywords internal
NULL

#' cSEMResults
#' 
#' A call to [csem()] results in an object with at least 
#' two class attributes. The first class attribute is always `cSEMResults` no matter
#' the type of data or model provided. 
#' The second is one of `cSEMResults_default`, `cSEMResults_multi`, or 
#' `cSEMResults_2ndorder` and depends on the estimated model and/or the type of 
#' data provided to the `.model` and `.data` arguments of [csem()].
#' The third class attribute `cSEMResults_resampled` is only added if resampling
#' was conducted.
#' 
#' Depending on the type of data and/or model provided three different output
#' types exists.
#' \describe{
#' \item{_default}{This will be the structure for the vaste majority of applications.
#'  If the data is a single `matrix` or `data.frame` with no id-column, 
#'  the result is a `list` with elements: 
#' \describe{
#'   \item{`$Estimates`}{A list containing a list of estimated quantities.}
#'   \item{`$Information`}{A list containing a list of additional information.}
#' }
#' The resulting object has classes `cSEMResults` and `cSEMResults_default`.
#' }
#' \item{_multi}{If the data provided is a single `matrix` or `data.frame` containing
#'  an id-column to split the data by `G` group levels 
#'  or if a list of `G` datasets is provided, the resulting object is a list of `G` 
#'  lists, where `G` is equal to the number of groups or the number of datasets 
#'  in the list of datasets provided. Each of the `G` list elements is itself 
#'  a `cSEMResults_default` object. Hence its structure is identical to 
#'  the structure described in `_default`.
#' 
#' The resulting object has classes `cSEMResults` and `cSEMResults_multi`. 
#' }
#' \item{_2ndorder}{
#'  A special output is generated if the model to estimate contains hierachical constructs
#'  **and** the "2stage" or "mixed" approach is used to estimate the model. In this case
#'  the resulting object is a list containing two elements `First_stage` and 
#' ` Second_stage`.
#' 
#' Each list element is itself a `cSEMResults_default` object. Hence its structure is identical to 
#' the structure described in `_default`.
#' }
#' }
#' 
#' If `.resample_method = "bootstrap"` or `.resample_method = "jackknife"`, resamples
#' are attached to each object. For objects of class `cSEMResults_default` the resamples are
#' attached to `.object$Estimates$Estimates_resample`. For objects of class
#' `cSEMResults_multi` the same is done by group. For objects of class
#' `cSEMResults_2ndorder` the resamples are attached to the 
#' `.object$Second_stage$Information$Resamples`. All objects containing 
#' these elements gain the `cSEMResults_resampled` class.
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