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
#' @return A matrix of construct scores.
#' 
#' @seealso [csem()], [cSEMResults]
#' @export

getConstructScores <- function(.object, .standardized = TRUE){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, getConstructScores, .standardized = .standardized)
    
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    
    ## Get scaled indicator scores [E(x) = 0; Var(x) = 1]
    indicators_scaled <- .object$Information$Data
    
    ## Original data
    .object$Information$Arguments$.data
    
    # Unscale
    indicators_unscaled = t(t(indicators_scaled) * attr(indicators_scaled, 'scaled:scale') + 
      attr(indicators_scaled, 'scaled:center')) 
    
    #Unstandardized scores equal X_i%*%w_i/sum(w_i)
    
    w_transformed <- t(.object$Estimates$Weight_estimates) %*% solve(diag(rowSums(.object$Estimates$Weight_estimates)))
    
    scores_unstandardized <- indicators_unscaled%*%w_transformed
    
    colnames(scores_unstandardized) =  rownames(.object$Estimates$Weight_estimates)
    
    return(scores_unstandardized)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Currently not implemented for models containing second-order constructs")
  } else {
    stop2("Dont know how to handle objects of class: ", paste0(class(.object), collapse = ", "))
  }
}