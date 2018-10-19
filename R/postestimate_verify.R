#' Verify admissibility
#'
#' Verify admissibility of the estimated quantities for a given model. 
#' 
#' Results based on an estimated model exhibiting one of the 
#' following defects are deemed inadmissible: non-convergence, loadings and/or 
#' reliabilities larger than 1, a construct VCV and/or a
#' model-implied VCV matrix that is not positive (semi-)definite.
#' 
#' For models containing second order constructs estimation is by default done
#' in a two/three stage procedure. In this case both the first and the second/third
#' stage model quantities are checked.
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @return A logical vector indicating which (if any) problem occurred. 
#'   A `TRUE` indicates that the problem has occured. For models containg second order
#'   constructs a list of two such vectors (one for the first and one for the second stage)
#'   is returned. Status codes are:
#'  \itemize{
#' \item 1: The algorithm has not converged.
#' \item 2: At least one absolute standardized loading estimate is larger than 1,
#    which implies either a negative variance of the measurement error or
#    a correlation larger than 1.
#' \item 3: The construct VCV is not positive semi-definite.
#' \item 4: At least one construct reliability is larger than 1. 
#' \item 5: The model-implied indicator VCV is not positive semi-definite. 
#'   This is only checked for linear models.
#' }
#'
#' @export

verify <- function(.object) {
  UseMethod("verify")
}

#' @describeIn verify (TODO)
#' @export
 
verify.cSEMResults_default <- function(.object){
  
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
  
  if(max(.object$Estimates$Construct_reliabilities)>1) {
    stat["4"] <- TRUE
  }
  
  if(.object$Information$Model$model_type == "Linear" && 
     !matrixcalc::is.positive.semi.definite(fit(.object))) {
    stat["5"] <- TRUE
  }
  
  class(stat) <- "cSEMVerify_default"
  return(stat) 
}

#' @describeIn verify (TODO)
#' @export

verify.cSEMResults_multi <- function(.object){
  
  lapply(.object, verify.cSEMResults_default)
}

#' @describeIn verify (TODO)
#' @export

verify.cSEMResults_2ndorder <- function(.object){
  
  stat1 <- stat2 <- c("1" = FALSE, "2" = FALSE, "3" = FALSE, "4" = FALSE, "5" = FALSE)
  
  ## First stage checks
  x1i <- .object$First_stage$Information
  x1e <- .object$First_stage$Estimates
  
  if(!(is.null(x1i$Weight_info$Convergence_status) || x1i$Weight_info$Convergence_status)) {
    stat1["1"] <- TRUE
  }
  
  if(max(abs(x1e$Cross_loadings)) > 1) {
    stat1["2"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(x1e$Construct_VCV)) {
    stat1["3"] <- TRUE
  }
  
  if(max(x1e$Construct_reliabilities)>1) {
    stat1["4"] <- TRUE
  }
  
  if(x1i$Model$model_type == "Linear" && 
     !matrixcalc::is.positive.semi.definite(fit.cSEMResults_default(.object$First_stage))) {
    stat1["5"] <- TRUE
  }
  
  ## Second stage checks
  x2i <- .object$Second_stage$Information
  x2e <- .object$Second_stage$Estimates  
  
  if(!(is.null(x2i$Weight_info$Convergence_status) || x2i$Weight_info$Convergence_status)) {
    stat2["1"] <- TRUE
  }
  
  if(max(abs(x2e$Cross_loadings)) > 1) {
    stat2["2"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(x2e$Construct_VCV)) {
    stat2["3"] <- TRUE
  }
  
  if(max(x1e$Construct_reliabilities)>1) {
    stat2["4"] <- TRUE
  }
  
  if(x1i$Model$model_type == "Linear" && 
     !matrixcalc::is.positive.semi.definite(fit.cSEMResults_default(.object$Second_stage))) {
    stat2["5"] <- TRUE
  }
  
  stat <- list("First_stage" = stat1, "Second_stage" = stat2)
  
  class(stat) <- "cSEMVerify_2ndorder"
  return(stat) 
}
