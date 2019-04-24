#' Assess model
#'
#' Assess model using common evaluation criteria and fit measures. 
#' 
#' For details see the vignette on assess.
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @return (TODO)
#' @export

assess <- function(.object, .only_common_factors = TRUE, ...) {
  UseMethod("assess")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_default <- function(.object, .only_common_factors = TRUE, ...){
  
  ## Get relevant objects
  con_types <-.object$Information$Model$construct_type
  names_cf  <- names(con_types[con_types == "Common factor"])
  P <- .object$Estimates$Construct_VCV
  
  # SRMR
  res_srmr  <- calculateSRMR(.object)
  # dG
  res_dg    <- calculateDG(.object)
  # dL
  res_dL    <- calculateDL(.object)
  # dML
  res_dml   <- calculateDML(.object)
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

assess.cSEMResults_multi <- function(.object, ...){
  
  paste("not yet implemented")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_2ndorder <- function(.object, ...){
  
  paste("not yet implemented")
}