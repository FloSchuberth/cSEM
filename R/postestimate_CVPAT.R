#' Perform a Cross-Validated Predictive Ability Test (CVPAT)
#'
#'\lifecycle{maturing}
#'
#' Perform a Cross-Validated Predictive Ability Test (CVPAT) as described in 
#' \insertCite{Liengaard2020}{cSEM}.
#'
#' @return An object of class `cSEMCVPAT` with print and plot methods.
#'   Technically, `cSEMCVPAT` is a 
#'   named list containing the following list elements:
#'
#' \describe{
#'   \item{`$teststatistic`}{The test statistic for H0: The predictive performance of
#'                            both models is identical.}
#'   \item{`$pvalue`}{The p-value for H0: The predictive performance of
#'                            both models is identical.}
#'     }
#'   
#' @usage postestimate_CVPAT(
#' .object1  = NULL,
#' .object2  = NULL,
#' .seed     = NULL,
#' .cv_folds = 10)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults], [exportToExcel()]
#' 
#' @references
#'   \insertAllCited{}
#'
#' @example inst/examples/example_predict.R
#' 
#' @export

postestimate_CVPAT <- function(
  .object1  = NULL,
  .object2  = NULL,
  .seed     = NULL,
  .cv_folds = 10){
  
  ##Errors and warnings---------------------------------------------------------
  #Stop if one object is second order
  if(inherits(.object1, "cSEMResults_2ndorder")|inherits(.object2, "cSEMResults_2ndorder")) {
    stop2('Currently, `predict()` is not implemented for models containing higher-order constructs.')
  }
  
  #Stop if one object is not admissible
  if(sum(verify(.object1))!=0|sum(verify(.object2))!=0) {
    stop2('One csem object is not admissible.')
  }
  
  #Stop if one object does not have a structural model
  if(all(.object1$Information$Model$structural == 0)|all(.object2$Information$Model$structural == 0)) {
    stop2("`predict()` requires a structural model.")
  }
  
  #Stop if one object contains a non-linear model
  if(.object1$Information$Model$model_type != 'Linear'|.object2$Information$Model$model_type != 'Linear'){
    stop2('Currently, `predict()` works only for linear models.')
  }
  
  #Stop if both objects are based on different datasets
  if(!all(.object1$Information$Data == .object2$Information$Data)){
    stop2('The objects are not based on the same dataset.')
  }
  
  #Stop if no seed is provided
  if(is.null(.seed)){
    stop2('postestimate_CVPAT() requires a seed.')
  }
  
  #Perform out-of-sample predictions for both models
  predict1 <- predict(.object = .object1, .benchmark = "lm", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed)
  
  predict2 <- predict(.object = .object2, .benchmark = "lm", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed)
  
  L1 <- unlist(predict1$Residuals_target[[1]])
  L2 <- unlist(predict2$Residuals_target[[1]])
  
  D_bar <- mean(L2 - L1)
  
  test_stat <- D_bar/sqrt(var(L2-L1)/length(L1))
  
  p_value <- 2*(1-pt(abs(test_stat), length(L1)-1))
  
  out <- list(
    "test_statistic" = test_stat,
    "p-value"        = p_value
  )
  class(out) = "cSEMCVPAT"
  out
  
}