#' cSEM: A package for composite-based structural equation modeling
#'
#' Estimate, analyse, test, and study linear, nonlinear, hierachical and
#' multigroup structural equation models using composite-based approaches and procedures including estimation 
#' techniques such as partial least squares path modelling (PLS) and its derivatives
#' (PLSc, ordPLSc, robustPLSc), generalized structured component analysis (GSCA), 
#' generalized structured component analyis with uniqueness terms (GSCAm), 
#' generalized canonical correlation analysis (GCCA) unit weights (sum score) and
#' fixed weights, as well as several tests and typical postestimation 
#' procedures (e.g., assess the model fit, compute direct, indirect and total effects).
#'
#' @importFrom magrittr %>%
#' @importFrom matrixStats rowProds
#' @importFrom rlang .data
#' @importFrom lifecycle deprecate_soft
#' @import crayon
#' @import cli
#' @import stats
#' @import utils
#' 
#' @keywords internal
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") globalVariables(c(".", "alpha", ""))