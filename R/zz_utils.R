#' cSEMArguments
#'
#' A list of all arguments including a description of their use and current defaults.
#' Mainly used for interal purposes (parameter inheritance).
#'
#' @param .data A data frame or a matrix containing the raw data.
#' @param .model A model in \href{http://lavaan.ugent.be/tutorial/syntax1.html}{lavaan model syntax}
#'   or a [cSEMModel] object.
#' @param .alpha A numeric vector of significance levels. Defaults to `c(0.1, 0.05, 0.01)`.
#' @param .approach_cf Character string. The approach to obtain the correction
#'   factor for PLSc. Possible choices are: "*dist_euclid*", "*dist_euclid_weighted*",
#'   "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*", "*mean_harmonic*",
#'   "*geo_of_harmonic*". Defaults to "*dist_euclid*". Ignored if `.disattenuate = FALSE`.
#' @param .approach_nl Character string. The approach used to handle non-linear
#'   structural relationships. Possible choices are: "*none*", or "*replace*".
#'   Defaults to "*none*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. Possible choices are: "*OLS*", "*2SLS*", or "*3SLS*".
#'   Defaults to "*OLS*".
#' @param .approach_weights Character string. The name of the approach used to
#'   obtain composite weights. Possible choices are: "*PLS*", "*SUMCOR*", "*MAXVAR*",
#'   "*SSQCOR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*fixed*", or "*unit*".
#'   Defaults to "*PLS*".
#' @param .C A (J x J) proxy variance-covariance matrix.
#' @param .criteria Character string. The criteria to use. One of "*SUMCOR*", "*MAXVAR*",
#'   "*SSQCOR*", "*MINVAR*" or "*GENVAR*".
#' @param .csem_model A (possibly incomplete) [cSEMModel] list.
#' @param .disattenuate Logical. Should proxy correlations be disattenuated
#'   if the construct is modeled as a common factor? Defaults to `TRUE`.
#' @param .dominant_indicators A character vector of `name = value` pairs, where `value` is 
#'   a character string giving the name of the dominant indicator and `name` 
#'   the corresponding construct name. Dominant indicators may be specified for 
#'   a subset of the constructs. 
#' @param .E A (J x J) matrix of inner weights.
#' @param .estimate_structural Logical. Should the structural (path) coefficients
#'   be estimated? Defaults to `TRUE`.
#' @param .group_var  Character string. The name of the column used to split the data into groups.
#' @param .H The (N x J) matrix of proxy values.
#' @param .ignore_structural_model Logical. Should the structural (path) model be ignored
#'   when calculating the inner weights of the PLS algorithm? Defaults to `FALSE`
#' @param .iter_max Integer. The maximum number of iterations of the PLS algorithm.
#'   If `iter_max = 1` one-step weights are returned. If the algorithm exceeds
#'   the specified number, weights of iteration step `.iter_max - 1`  will be returned
#'   with a warning. Defaults to `100`. The argument is ignored if
#'   `.approach_weight` is not PLS.
#' @param .mode A named vector giving the mode used to obtain new outer weights.
#' @param .normality Logical. Should joint normality be assumed if the model
#'   contains non-linear terms. For details see: \insertCite{Dijkstra2014;textual}{cSEM}. 
#'   Defaults to `TRUE`.
#' @param .permutations Integer. The number permutations used for step 2 and 3 of the test.
#' Defaults to `100`. Note however that the number should typically be a lot higher.
#' @param .PLS_weight_scheme_inner Character string. The inner weighting scheme
#'   used in PLS. Possible choices are: "*centroid*", "*factorial*", or "*path*".
#'   Defaults to "*centroid*". The argument is ignored if `.approach_weight` is not PLS.
#' @param .PLS_mode Either a named vector specifying the mode that should be used for
#'   each construct in the form `name = "mode"`, a single character
#'   string giving the mode that should be used for all constructs or `NULL`.
#'   Possible choices are: "*ModeA*" or "*ModeB*".
#'   If `NULL`, `cSEM` will choose the appropriate mode according to the type
#'   of construct used. The argument is ignored if `.approach_weight` is not PLS.
#' @param .Q A vector of reliability coefficients with element names equal to
#'   the names of the J construct names used in the measurement model.
#' @param .reliabilities A vector of `name = value` pairs of reliability coefficient. 
#'   Element names are equal to the names of the J construct names. Reliabilities
#'   may be given for a subset of the constructs. 
#' @param .S The (scaled) empirical (K x K) indicator covariance/correlation matrix,
#'   where K is the number of indicators.
#' @param .terms A vector of construct names to be classified.
#' @param .tolerance Double. The tolerance criterion for convergence of the PLS
#'   algorithm. Defaults to `1e-05`.
#'   The argument is ignored if `.approach_weight` is not PLS.
#' @param .verbose Logical. Should information be printed to the console ?
#' @param .W The (J x K) weight matrix, where J is the number of constructs and
#'   K the number of indicators.
#' @param ... Further arguments to be passed down to other functions.
#' See [args_default] for a complete list of possible arguments and their defaults.
#'
#'
#' @name csem_arguments
#' @aliases cSEMArguments csem_parameters
NULL

#' cSEMModel
#'
#' A standardized list containing model-related information. To convert a
#' a model written in \href{http://lavaan.ugent.be/tutorial/syntax1.html}{lavaan model syntax}
#' to a cSEMModel list use [parseModel].
#'
#' An object of class cSEMModel is a standardized list containing the following
#' components (assume in the following that there are J constructs and K indicators)
#' \describe{
#'   \item{`$structural`}{A matrix mimicking the structural relationship between
#'      constructs. If constructs are only linearly related, `structural` is
#'      of dimension (J x J) with row- and column names equal to the construct
#'      names. If the structural model contains non-linear relationships
#'      `structural` is (J x (J + J\*)) where J\* is the number of
#'      non-linear terms.}
#'   \item{`$structural_ordered`}{A matrix of the same dimension as `structural`
#'     with rows rearranged to satisfy the following criteria. 1. the construct
#'     in the first row depends on exogenous constructs only. 2. the constructs
#'     of the following rows depend only on exogenous constructs or those
#'     of previous rows. This is required to estimate non-linear structural equation
#'     model relationships using the replacement approach.
#'     }
#'   \item{`$measurement`}{A (J x K) matrix mimicking the measurement relationship
#'     between constructs and their related indicators with row names equal to
#'     the construct names and column names equal to the indicator names.}
#'   \item{`$error_cor`}{A (K x K) matrix mimicking the measurement error
#'     correlation relationship.}
#'   \item{`$construct_type`}{A data frame with two columns *"Name"* and *"Type"*
#'     containing the name of each construct and their respective type
#'     (**"Common factor"** or **"Composite"**).}
#'   \item{`$vars_endo`}{A vector of names of the endogenous constructs.}
#'   \item{`$vars_exo`}{A vector of names of the exogenous constructs (incudes
#'     possible interaction and exponential terms).}
#'   \item{`$vars_explana`}{ A vector of names of the constructs that appear as
#'     explanatory variables in at least one structural equation (incudes
#'     possible interaction and exponential terms).}
#'   \item{`$explained_by_exo`}{A vector of names of the constructs that are
#'     solely explained by exogenous constructs.}
#' }
#' Note: it is possible to supply an incomplete cSEMModel list
#' to all functions that require `.csem_model` as a mandatory argument. Currently,
#' only the structural and the measurement matrix are required.
#' However, specifying an incomplete cSEMModel list may lead to unexpected behaviour 
#' and errors so do use this technique with caution.
#'
#' @seealso [parseModel]
#' @name csem_model
#' @aliases cSEMModel
NULL

#' cSEMResults
#'
#' @return
#' An object of class `cSEMResults` for which the following methods exist:
#' \describe{
#'   \item{`print.cSEMResults`}{Prints a message to inform the user
#'   whether estimation has been successful and what functions may be used
#'   to examine the object.}
#'   \item{`summary.cSEMResults`}{Print and return a comprehensive summary of the results.}
#' }
#' Technically `cSEMResults` is a named list containing the following list elements:
#' \describe{
#'   \item{`$Estimates`}{A list containing the estimated quantities.}
#'   \item{`$Tests`}{A list containing test results. (not yet implemented)}
#'   \item{`$Meta_information`}{A list of additional information. (incomplete)}
#' }
#'
#' @name csem_results
#' @aliases cSEMResults
NULL
