#' HTMT
#'
#' Compute the heterotrait-monotrait ratio of correlations \insertCite{Henseler2015}{cSEM}.
#'
#' @usage HTMT(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#' @references 
#'\insertAllCited{}
#'
#' @export
#'
HTMT <- function(
  .object              = args_default()$.object,
  .only_common_factors = args_default()$.only_common_factors
  ){

  ## Get relevant quantities
  m <- .object$Information$Model
 
  if(isTRUE(.only_common_factors)) {
    cf_names <-names(m$construct_type[m$construct_type == "Common factor"])
    
    ## Stop if there are no common factors
    if(length(cf_names) < 2) {
      stop("Computation of the HTMT requires at least two common factors, unless `.only_common_factors = FALSE`.",
           call. = FALSE)
    }
  } else {
    cf_names <- names(m$construct_type)
  }
  
  cf_measurement <- m$measurement[cf_names, colSums(m$measurement[cf_names, ]) != 0, drop = FALSE]
  
  ## HTMT can only be calculated for constructs with more than one indicator
  x <- rowSums(cf_measurement) > 1
  cf_measurement <- cf_measurement[x, colSums(cf_measurement[x, ]) != 0, drop = FALSE]
  
  ## At least two multi-indicator constructs required
  if(length(which(x)) < 2) {
    stop("Computation of the HTMT requires at least two multi indicator constructs.",
        call. = FALSE)
  }
  
  i_names <- colnames(cf_measurement)
  S       <- .object$Estimates$Indicator_VCV[i_names, i_names]
  
  ## Average correlation of the indicators of a block 
  avrg_cor <- cf_measurement %*% (S - diag(diag(S))) %*% t(cf_measurement) /
    cf_measurement %*% (1 - diag(nrow(S))) %*% t(cf_measurement)
  
  ## Compute HTMT
  out <- avrg_cor*lower.tri(avrg_cor) / sqrt(diag(avrg_cor) %o% diag(avrg_cor))
  
  # Return
  return(out)
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

SRMR <- function(.matrix1 = args_default()$.matrix1,
                 .matrix2 = args_default()$.matrix2) {
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  
  S         <- .matrix1
  Sigma_hat <- .matrix2

  
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