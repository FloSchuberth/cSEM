#' Do a redundancy analysis
#'
#' Perform a redundancy analysis (RA) as proposed by \insertCite{Hair2016;textual}{cSEM}
#' with reference to \insertCite{Chin1998;textual}{cSEM}.
#'
#' RA is confined to PLS-PM, specifically PLS-PM with at least one construct
#' whose weights are obtained by mode B. In cSEM this is the case if the construct
#' is modeled as a composite or if argument `.PLS_modes` was explicitly set to 
#' mode B for at least one construct.
#' Hence RA is only conducted if `.approach_weights = "PLS-PM"` and if at least
#' one construct's mode is mode B.
#'  
#' The principal idea of RA is to take two different measures of the 
#' same construct and regress the scores obtained for each measure on each
#' other. If they are similar they are likely to measure the same "thing"
#' which is then taken as evidence that both measurement models actually
#' measure what they are supposed to measure (validity). 
#' 
#' There are several issues with the terminology  and the reasoning behind this logic.
#' RA is therefore only implemented since reviewers are likely to demand
#' its computation, however, its actual application for validity assessment 
#' is discouraged.
#' 
#' Currently, the function is not applicable to models containing second-order 
#' constructs.
#' 
#' @usage doRedundancyAnalysis(.object = NULL)
#'
#' @return A named numeric vector of correlations. If 
#'   the weighting approach used to obtain `.object` is not `"PLS-PM"` or 
#'   non of the PLS outer modes was mode B, the function silently returns `NA`.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @export

doRedundancyAnalysis <- function(.object = NULL) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doRedundancyAnalysis)
    
    class(out) <- c("cSEMRedundancyAnalysis", "cSEMRedundancyAnalysis_multi")
    return(out)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Currently, `doRedundancyAnalysis()` is not implemented for models", 
          "containing higher-order constructs.")
  } else {
    
    if(.object$Information$Arguments$.approach_weights == "PLS-PM") {
      
      modes <- .object$Information$Weight_info$Modes
      modesB <- modes[modes == "modeB"]
      
      if(length(modesB) > 0) {
        
        args <- .object$Information$Arguments
        
        # Functions resampleData() and testMICOM require the .id argument.
        # It is therefore present in the Arguments list although the 
        # data set in the Arguments list does not contain the id column anymore.
        # Therefore .id needs to be set to NULL
        args[[".id"]] <- NULL
        new_modes <- as.list(modes) 
        
        beta <- c()
        for(j in names(modesB)) {
          new_modes_j <- new_modes
          new_modes_j[j] <- "modeA"
          
          args[[".PLS_modes"]] <- new_modes_j
          
          res_reflective <- do.call(csem, args)
          
          Y <- res_reflective$Estimates$Construct_scores
          X <- .object$Estimates$Construct_scores
          
          beta[j] <- c(solve(t(X[,j]) %*% X[,j]) %*% t(X[,j]) %*% Y[,j])
        }
      } else {
        beta <- NA
      }
    } else {
      beta <- NA
    }
    
    class(beta) <- "cSEMRedundancyAnalysis"
    return(beta) 
  }
}