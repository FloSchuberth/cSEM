#' Verify admissibility
#'
#' Verify admissibility of the estimated quantities for a given model.
#' 
#' Results based on an estimated model exhibiting one of the 
#' following defects are deemed inadmissible: non-convergence, loadings and/or 
#' reliabilities larger than 1, a construct VCV and/or a
#' model-implied VCV matrix that is not positive (semi-)definite.
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @return A logical vector indicating which (if any) problem occurred. Status codes are:
#'  \itemize{
#' \item 1: The algorithm has not not converged.
#' \item 2: At least one absolute standardized loading estimate is larger than 1,
#    which implies either a negative variance of the measurement error or
#    a correlation larger than 1.
#' \item 3: The construct VCV is not positive semi-definite.
#' \item 4: The model-implied indicator VCV is not positive semi-definite.
#' \item 5: At least one construct reliability is larger than 1. 
#' }
#'
#' @export

verify <- function(.object){
  
  if(.object$Information$Model$model_type != "Linear"){
    stop("Currently, the status function only works for linear models.",
         call. = FALSE)}
  
  stat <- c("1" = FALSE, "2" = FALSE, "3" = FALSE, "4" = FALSE, "5" = FALSE)
  
  if(!(is.null(.object$Information$Weight_info$Convergence_status) || 
       .object$Information$Weight_info$Convergence_status)) {
    stat["1"] <- TRUE
  }
  
  if(max(abs(.object$Estimates$Cross_loadings)) > 1) {
    stat["2"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(.object$Estimates$Construct_VCV)) {
    stat["3"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(fit(.object))) {
    stat["4"] <- TRUE
  }
  
  if(max(.object$Estimates$Construct_reliabilities)>1) {
    stat["5"] <- TRUE
  }
  
  class(stat) <- "cSEMVerify"
  return(stat) 
}
