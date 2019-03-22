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

assess <- function(.object, ...) {
  UseMethod("assess")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_default <- function(.object, ...){
  
  # SRMR
  res_srmr  <- calculateSRMR(.object)
  # dG
  res_dg    <- calculateDG(.object)
  # dL
  res_dL    <- calculateDL(.object)
  # dML
  res_dml   <- calculateDML(.object)
  # HTMT 
  res_htmt  <- calculateHTMT(.object, ...)
  # AVE
  res_ave   <- calculateAVE(.object, ...)
  # RhoC
  res_rhoc  <- calculateRhoC(.object, ...)
  # RhoT
  res_rhot  <- calculateRhoT(.object, ...)
  # Effect size
  res_esize <- calculateEffectSize(.object)
  
  out <- list(
    "AVE"  = res_ave,
    "RhoC" = res_rhoc,
    "RhoT" = res_rhot,
    "HTMT" = res_htmt,
    "SRMR" = res_srmr,
    "dG"   = res_dg,
    "dL"   = res_dL,
    "dML"  = res_dml,
    "Effect size" = res_esize
  )
  
  class(out) <- "cSEMAssess_default"
  return(out)
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_multi <- function(.object){
  
  paste("not yet implemented")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_2ndorder <- function(.object){
  
  paste("not yet implemented")
}