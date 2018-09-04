SRMR <- function(.object = args_default()$.object) {
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object)
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 

fitRMS <- function() {
  
}