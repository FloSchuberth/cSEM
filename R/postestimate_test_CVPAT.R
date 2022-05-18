#' Perform a Cross-Validated Predictive Ability Test (CVPAT)
#'
#'\lifecycle{maturing}
#'
#' Perform a Cross-Validated Predictive Ability Test (CVPAT) as described in 
#' \insertCite{Liengaard2020}{cSEM}. The predictive performance of two models 
#' based on the same dataset is compared. In doing so, the average difference in 
#' losses in predictions is compared for both models. 
#'
#' @return An object of class `cSEMCVPAT` with print and plot methods.
#'   Technically, `cSEMCVPAT` is a 
#'   named list containing the following list elements:
#'   
#'   \describe{
#'   \item{'$Information'}{Additional information.}
#'   }
#'   
#' @usage testCVPAT(
#' .object1  = NULL,
#' .object2  = NULL,
#' .approach_predict = c("earliest", "direct"),
#' .seed     = NULL,
#' .cv_folds = 10)
#'
#' @inheritParams csem_arguments
#' @param .object1 An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param .object2 An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param .approach_predict Character string. Which approach should be used to 
#'  predictions? One of "*earliest*" and "*direct*". If "*earliest*" predictions
#'  for indicators associated to endogenous constructs are performed using only
#'  indicators associated to exogenous constructs. If "*direct*", predictions for 
#'  indicators associated to endogenous constructs are based on indicators associated
#'  to their direct antecedents. Defaults to "*earliest*".
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
  .approach_predict = c("earliest", "direct"),
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
  if(!all(.object1$Information$Data[,colnames(.object1$Information$Data) %in% colnames(.object2$Information$Data)] == .object2$Information$Data[, colnames(.object2$Information$Data) %in% colnames(.object1$Information$Data)])){
      stop2('The objects are not based on the same dataset.')
  }

  
  #Stop if no seed is provided
  if(is.null(.seed)){
    warning2('`testCVPAT()` requires a seed. Since no seed has been provided a random seed is used.')
    .seed = sample(1:9999,1)
  }
  
  
  #Perform out-of-sample predictions for both models
  predict1 <- predict(.object = .object1, .benchmark = "NA", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed, 
                      .approach_predict = .approach_predict)
  
  predict2 <- predict(.object = .object2, .benchmark = "NA", 
                      .cv_folds = .cv_folds, .r = 1, .seed = .seed,
                      .approach_predict = .approach_predict)
  
  L1 <- rowMeans(predict1$Residuals_target[[1]]^2)
  L2 <- rowMeans(predict2$Residuals_target[[1]]^2)
  
  N <- length(L1)
  
  D <- L2 - L1
  
  D_bar <- mean(D)
  
  test_stat <- D_bar/sqrt(var(D)/N)
  
  p_value <- 2*(1-pt(abs(test_stat), N-1))
  
  #missing: degrees of freedom, decision, information: seed, number cv_folds, alpha
  out <- list(
    "Information"          = list(
      "Prediction 1"       = predict1,
      "Prediction 2"       = predict2,
      "Sample Size"        = N,
      "Seed"               = .seed,
      "Degrees_of_Freedom" = N-1
    ),
    "test_statistic" = test_stat,
    "p_value"        = p_value,
    "degrees of freedom" = N-1
  )
  class(out) = "cSEMTestCVPAT"
  out
  
}