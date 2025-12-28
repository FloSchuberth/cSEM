#' Internal: Complete list of csem()'s ... arguments
#' 
#' A complete alphabetical list of all possible arguments accepted by `csem()`'s `...` 
#' (dotdotdot) argument.
#' 
#' Most arguments supplied to the `...` argument of `csem()` are only
#' accepted by a subset of the functions called by `csem()`. The following
#' list shows which argument is passed to which (internal) function:
#' \describe{
#' \item{.approach_cor_robust}{Accepted by/Passed down to: [calculateIndicatorCor()]}
#' \item{.conv_criterion}{Accepted by/Passed down to: [calculateWeightsPLS()],
#'   [calculateWeightsGSCA()], [calculateWeightsGSCAm()] and subsequently 
#'   [checkConvergence()].}
#' \item{.dominant_indicators}{Accepted by/Passed down to: [setDominantIndicator()]}
#' \item{.estimate_structural}{Accepted by/Passed down to: [foreman()]}
#' \item{.iter_max}{Accepted by/Passed down to: [calculateWeightsPLS()],
#'   [calculateWeightsGSCA()], [calculateWeightsGSCAm()]}
#' \item{.PLS_modes, .PLS_ignore_structural_model, .PLS_weight_scheme_inner, .PLS_approach_cf}{
#'   Accepted by/Passed down to: [calculateWeightsPLS()]}
#' \item{.tolerance}{Accepted by/Passed down to: [calculateWeightsPLS()],
#'   [calculateWeightsGSCA()], [calculateWeightsGSCAm()], [calculateWeightsUnit()]}
#' }
#' 
#' @usage NULL
#' 
#' @inheritParams csem_arguments
#' @keywords internal

args_csem_dotdotdot <- function(
  .approach_cor_robust     = c("none", "mcd", "spearman"),
  .conv_criterion          = c("diff_absolute", "diff_squared", "diff_relative"),
  .dominant_indicators     = NULL,
  .estimate_structural     = TRUE,
  .iter_max                = 100,
  .PLS_modes               = NULL,
  .PLS_ignore_structural_model = FALSE,
  .PLS_weight_scheme_inner     = c("path", "centroid", "factorial"),
  .PLS_approach_cf         = c("dist_squared_euclid", "dist_euclid_weighted", 
                               "fisher_transformed", "mean_arithmetic",
                               "mean_geometric", "mean_harmonic",
                               "geo_of_harmonic"),
  .tolerance               = 1e-05
) {NULL}