#' Get construct scores
#'
#' Get the standardized or unstandardized construct scores.
#' 
#' @usage getConstructScores(
#'  .object        = NULL,
#'  .standardized  = TRUE
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @return A list of three with elements `Construct_scores`, `W_used`, 
#'   `Indicators_used`.
#' 
#' @seealso [csem()], [cSEMResults]
#' @export

getConstructScores <- function(.object = NULL, .standardized = TRUE){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, getConstructScores, .standardized = .standardized)
    return(out)
    
  } else if(inherits(.object, "cSEMResults_default")) {
    
    if(.standardized) {
      ## Get scaled indicator scores [E(x) = 0; Var(x) = 1]
      list("Construct_scores" = .object$Estimates$Construct_scores,
           "W_used"           = .object$Estimates$Weight_estimates,
           "Indicators_used"  = .object$Information$Data)
      
    } else {
      ## Get scaled indicator scores [E(x) = 0; Var(x) = 1]
      indicators_std <- .object$Information$Data
      
      ## Unscale
      indicators_unstd <- t(t(indicators_std) * attr(indicators_std, 'scaled:scale') + 
                                attr(indicators_std, 'scaled:center')) 
      
      W_std <- .object$Estimates$Weight_estimates
      
      W_unstd <- W_std
      
      for(i in colnames(W_std)){
        W_unstd[,i] <- W_std[,i]/sd(indicators_unstd[,i])  
      }
      
      # Scale weight so that their sum is equal to 1
      W_normed <- W_unstd
      W_normed <- solve(diag(rowSums(W_normed)))%*%W_normed 
      dimnames(W_normed) <- dimnames(W_unstd)
      
      ## Calculate Unstandarized construct scores
      scores_unstandardized <- indicators_unstd %*% t(W_normed)
      colnames(scores_unstandardized) <- rownames(.object$Estimates$Weight_estimates)
      
      list("Construct_scores" = scores_unstandardized,
           "W_used"           = W_normed,
           "Indicators_used"  = indicators_unstd)
    }

  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Currently, `getConstructScores()` is not implemented for models", 
          " containing second-order constructs.")
  } else {
    stop2("`.object` must be of class cSEMResults")
  }
}