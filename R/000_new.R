### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


# Calculates the Heterotrait-Monotrait factor correlations, see Henseler et al. (2015)


HTMT = function(.object,.only_common_factors=TRUE){
  
  # Adapted from matrixpls
  

  if(.only_common_factors==TRUE){
  # Extract names of the common factors, the HTMT is only calculated for common factors
  cf_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Common factor"])
  
  # Indicators connected to a common factor
  cf_measurement=.object$Information$Model$measurement[cf_names,
                                        colSums(.object$Information$Model$measurement[cf_names,])!=0]
  } else{ # in case of composites
    cf_names=names(.object$Information$Model$construct_type)
    cf_measurement=.object$Information$Model$measurement[cf_names,
                                                         colSums(.object$Information$Model$measurement[cf_names,])!=0]
  }
  ind_names=colnames(cf_measurement)
  S_relevant=.object$Estimates$Indicator_VCV[ind_names,ind_names]
  
  1-diag(nrow(S_relevant))
  
  
  # calculate average correlation of the indicators of a block 
  average_correlation_per_block=cf_measurement%*%(S_relevant-diag(diag(S_relevant)))%*%t(cf_measurement)/
    cf_measurement%*%(1-diag(nrow(S_relevant))) %*%t(cf_measurement)
  
  # Choose constructs that are measured by at least 2 indicators
  i = which(rowSums(cf_measurement) > 1)
  relevant_average_block_correlations <- average_correlation_per_block[i,i]
  
  if(length(i)<2){
    if(.only_common_factors==TRUE){
    stop("The HTMT can only be calculated in case of two common factors with at least two indicators per common factor.")
    } else {
      stop("The HTMT can only be calculated in case of two constructs with at least two indicators per construct.")
    }
      
      }
  
  htmt <- relevant_average_block_correlations*lower.tri(relevant_average_block_correlations) /
    sqrt(diag(relevant_average_block_correlations) %o% diag(relevant_average_block_correlations))
  
  htmt
  
}

#' Effects
#'
#' Compute direct, indirect and total effects.
#'
#' @usage effects(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export
#'
effects.cSEMResults <- function(.object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  ## Endogenous (lhs) variables
  vars_endo <- .object$Information$Model$vars_endo
  
  ## Matrix of direct effects:
  direct <- .object$Estimates$Path_estimates
  
  ## Matrix of total total effects: B = direct
  # Note: eta = B x eta + zeta
  #       (I - B)*eta = zeta
  
  B_star <- diag(nrow(direct)) - direct
  total <- solve(B_star) - diag(nrow(direct))
  
  ## Matrix of indirect effects:
  indirect <- total - direct
  
  list(direct = direct[vars_endo, , drop = FALSE], 
       indirect = indirect[vars_endo, , drop = FALSE], 
       total = total[vars_endo, , drop = FALSE])
}
