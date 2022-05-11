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
#' @usage testCVPAT(
#' .object1  = NULL,
#' .object2  = NULL,
#' .seed     = NULL,
#' .cv_folds = 10)
#'
#' @inheritParams csem_arguments
#' @param .object1 An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param .object2 An R object of class [cSEMResults] resulting from a call to [csem()].
#'
#' @seealso [csem], [cSEMResults], [exportToExcel()]
#' 
#' @references
#'   \insertAllCited{}
#'
#' @example inst/examples/example_predict.R
#' 
#' @export

testCVPAT <- function(
  .object1  = NULL,
  .object2  = NULL,
  .seed     = NULL,
  .cv_folds = 10){
  
  ##Errors and warnings---------------------------------------------------------
  #Stop if one object is not of class "cSEMResults"
  if(!inherits(.object1, "cSEMResults")|!inherits(.object2, "cSEMResults")) {
    stop2('The objects should be of class "cSEMResults".')
  }
  
  
  #Stop if one object is second order
  #TestCVPAT() uses predict() which is not implemented for models containing higher-order constructs
  if(inherits(.object1, "cSEMResults_2ndorder")|inherits(.object2, "cSEMResults_2ndorder")) {
    stop2('Currently, `testCVPAT()` is not implemented for models containing higher-order constructs.')
  }
  
  #Stop if one object is not admissible
  if(sum(verify(.object1))!=0|sum(verify(.object2))!=0) {
    stop2('One cSEM object is not admissible.')
  }
  
  #Stop if one object does not have a structural model
  if(all(.object1$Information$Model$structural == 0)|all(.object2$Information$Model$structural == 0)) {
    stop2("`testCVPAT()` requires a structural model.")
  }
  
  #Stop if one object contains a non-linear model
  if(.object1$Information$Model$model_type != 'Linear'|.object2$Information$Model$model_type != 'Linear'){
    stop2('Currently, `testCVPAT()` works only for linear models.')
  }
  
  #Stop if both objects are based on different datasets
  if(!all(.object1$Information$Data == .object2$Information$Data)){
    stop2('The objects are not based on the same dataset.')
  }
  
  #Stop if no seed is provided
  if(is.null(.seed)){
    warning2('`testCVPAT()` requires a seed. Since no seed has been provided a random seed is used.')
    .seed = sample(1:1000,1)
  }
  
  
  #Perform out-of-sample predictions for both models
  predict1 <- predict(.object = .object1, .benchmark = "lm", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed)
  
  predict2 <- predict(.object = .object2, .benchmark = "lm", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed)
  
  L1 <- c(as.matrix(predict1$Residuals_target[[1]]))
  L2 <- c(as.matrix(predict2$Residuals_target[[1]]))
  
  D <- L2 - L1
  
  D_bar <- mean(D)
  
  test_stat <- D_bar/sqrt(var(D)/length(L1))
  
  p_value <- 2*(1-pt(abs(test_stat), length(L1)-1))
  
  #missing: degrees of freedom, decision, information: seed, number cv_folds, alpha
  out <- list(
    "test_statistic" = test_stat,
    "p-value"        = p_value
  )
  class(out) = "cSEMCVPAT"
  out
  
}