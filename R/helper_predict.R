#' Internal: Calculate prediction metrics
#' 
#' Currently, the following prediction measures are available:
#' 
#' \itemize{
#' \item Mean absolute error
#' \item Mean absolute percentage error
#' \item Mean squared error
#' \item Root mean squared error
#' \item Theil's forecast accuracy
#' \item Theil's forecast quality
#' \item Bias proportion of MSE
#' \item Regression proportion of MSE
#' \item Disturbance proportion of MSE
#' }
#' 
#'
#' @return A vector of the prediction measures for the observed variables
#' belonging to endogenous constructs
#'   
#' @inheritParams csem_arguments
#'
#'
#' @keywords internal
#' 



calculateMAE <- function(resid){
  return(apply(na.omit(resid), 2, function(x) mean(abs(x))))
}

calculateRMSE <- function(resid){
  return(apply(na.omit(resid), 2, function(x) sqrt(mean(x^2))))
}

calculateMissclassification <- function(resid){
  return(apply(resid, 2, function(x) length(x[abs(x)>0.00001])/length(x)))
}

calculateMAPE <- function(resid, act){
  return(sapply(colnames(act), function(x) mean(abs(resid[,x]/act[,x]))))
}

calculateMSE2 <- function(pred, act, resid){
  return(sapply(colnames(act), function(x) (mean(pred[,x]) - mean(act[,x]))^2 +
           var(resid[,x])))
}

calculateU1 <- function(act, mse2){
  return(sapply(colnames(act), function(x) sqrt(mse2[x])/sqrt(mean(act[,x]^2))))
}

calculateU2 <- function(act, resid){
  return(sapply(colnames(act), function(x) sqrt(sum(resid[,x]^2))/sqrt(sum(act[,x]^2))))
}

calculateUM <- function(pred, act, mse2){
  return(sapply(colnames(act), function(x) (mean(pred[,x]) - mean(act[,x]))^2/mse2[x]))
}

calculateUR <- function(pred, act, mse2){
  return(sapply(colnames(act), function(x) suppressWarnings((sd(pred[,x]) - cor(pred[,x], act[,x])*
                                          sd(act[,x]))^2/mse2[x])))
}

calculateUD <- function(pred, act, mse2){
  return(sapply(colnames(act), function(x) suppressWarnings((1 - cor(pred[,x], act[,x])^2)*
                  var(act[,x])^2/mse2[x])))
}

calculateq2 <- function(res, MB){
  q2_predict <- c()
   for(i in colnames(res)) {
     q2_predict[i] <- 1- sum((na.omit(res[, i]) - mean(na.omit(res[, i])))^2) /
      sum((na.omit(MB[, i]) - mean(na.omit(MB[, i])))^2)
   }
  return(q2_predict)
}