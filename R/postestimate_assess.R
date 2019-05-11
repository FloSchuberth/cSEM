#' Assess model
#'
#' Assess a model using common evaluation criteria and fit measures.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article on the
#' \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} for details.
#' 
#' The function is essentially a wraper around a number of internal functions
#' that perform an "assessment task" (called a **quality criterion** in \pkg{cSEM}
#' parlance) like computing the (congeneric) reliability,
#' the effect size, the heterotrait-monotrait ratio of correlations (HTMT) etc.
#' 
#' By default every possible quality criterion is calculated (`.what = "all"`). 
#' If only a subset of quality criteria needs to be computed a single character string
#' or a vector of character strings naming the quantity to compute may be 
#' supplied to `assess()` via the `.what` argument.
#' 
#' Some of the quality criteria are inherently tied to the classical common
#' factor model and therefore only meaningfully interpreted within a common
#' factor model (see the 
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article for details). 
#' It is possible to force computation of all quality criteria for constructs 
#' modeled as composites by setting `.only_common_factors = FALSE`, however, 
#' we explicitly warn to interpret quality criteria in this case with caution, 
#' as they may not even have a conceptual meaning. 
#'
#' @inheritParams csem_arguments
#' @param ... Further arguments passed to functions called by `assess()`.
#'   See [args_assess_dotdotdot] for a complete list of available arguments.
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @return A named list of quality criteria. Note that if only a single quality
#'   criteria is computed the return value is still a list!
#' @export

assess <- function(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  ...
  ){
  UseMethod("assess")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_default <- function(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  ...
  ){
  
  ## Get relevant objects
  con_types <-.object$Information$Model$construct_type
  names_cf  <- names(con_types[con_types == "Common factor"])
  P <- .object$Estimates$Construct_VCV
  
  # SRMR
  res_srmr  <- calculateSRMR(.object, ...)
  # dG
  res_dg    <- calculateDG(.object, ...)
  # dL
  res_dL    <- calculateDL(.object, ...)
  # dML
  res_dml   <- calculateDML(.object, ...)
  # GoF
  res_gof   <- calculateGoF(.object, .only_common_factors)
  
  # AVE
  res_ave   <- calculateAVE(.object, .only_common_factors)
  # RhoC
  res_rhoc  <- calculateRhoC(.object, .only_common_factors)
  # RhoT
  res_rhot  <- calculateRhoT(.object, .only_common_factors, ...)
  
  # Effect size
  res_esize <- calculateEffectSize(.object)
  # VIFModeB
  res_vifmodeb <- calculateVIFModeB(.object)
  
  # Redundancy analysis (RA)
  res_ra <- calculateRA(.object)
  
  # HTMT 
  res_htmt  <- calculateHTMT(.object, .only_common_factors)
  # Fornell-Larcker
  if(.only_common_factors) {
    P <- P[names_cf, names_cf]
  }

  FL_matrix <- cov2cor(P)^2
  diag(FL_matrix) <- res_ave 
  
  ## Output --------------------
  out <- list(
    "SRMR" = res_srmr,
    "dG"   = res_dg,
    "dL"   = res_dL,
    "dML"  = res_dml,
    "GoF"  = res_gof,
    "AVE"  = res_ave,
    "RA"   = res_ra,
    "R2"   = .object$Estimates$R2,
    "R2_adj" = .object$Estimates$R2adj,
    "VIF"  = .object$Estimates$VIF,
    "RhoC" = res_rhoc,
    "RhoT" = res_rhot,
    "HTMT" = res_htmt,
    "FL_matrix" = FL_matrix,
    "Effect size" = res_esize,
    "VIF_modeB" = res_vifmodeb
  )
  
  class(out) <- "cSEMAssess_default"
  return(out)
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_multi <- function(
  .object              = NULL,
  .only_common_factors = TRUE,
  ...
  ){
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, assess.cSEMResults_2ndorder, 
           .only_common_factors = TRUE, ...)
  } else {
    lapply(.object, assess.cSEMResults_default, 
           .only_common_factors = TRUE, ...)
  }
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_2ndorder <- function(.object, ...){
  
  stop2("Currently, second-order models are not supported by `assess()`.")
}