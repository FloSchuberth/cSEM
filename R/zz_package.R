#' cSEM: A package for composite-based structural equation modeling
#'
#' Estimate, analyze, test, and study linear, nonlinear, hierarchical and
#' multigroup structural equation models using composite-based approaches and procedures including estimation 
#' techniques such as partial least squares path modeling (PLS) and its derivatives
#' (PLSc, ordPLSc, robustPLSc), generalized structured component analysis (GSCA), 
#' generalized structured component analysis with uniqueness terms (GSCAm), 
#' generalized canonical correlation analysis (GCCA) unit weights (sum score) and
#' fixed weights, as well as several tests and typical postestimation 
#' procedures (e.g., assess the model fit, compute direct, indirect and total effects).
#'
#' @importFrom matrixStats rowProds
#' @importFrom rlang .data
#' @importFrom cli rule boxx ansi_align style_bold col_green col_red col_cyan
#'   col_magenta col_yellow
#' @importFrom lifecycle deprecated
#' @import stats
#' @import utils
#'
#' @keywords internal
"_PACKAGE"

## quiets concerns of R CMD check re: variables used in NSE contexts
if(getRversion() >= "2.15.1") globalVariables(c("alpha", "", "Construct_type", "pruned_tree"))