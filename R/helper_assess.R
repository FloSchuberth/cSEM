#' HTMT
#'
#' Compute the heterotrait-monotrait ratio of correlations (TODO)
#'
#' @usage HTMT(.object)
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
HTMT <- function(.object,.only_common_factors=TRUE){
  
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
#' SRMR
#'
#' Compute the SRMR (TODO)
#'
#' @usage SRMR(.object)
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

SRMR <- function(.object = args_default()$.object) {
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object)
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 

#' RMS
#'
#' Compute the RMS (TODO)
#'
#' @usage RMS(.object)
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

RMS <- function() {
  
}