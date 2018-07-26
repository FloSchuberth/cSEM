# Distance functions

### Squared euclidean distance

dL <- function(.object) {
  
  S     <- .object$Estimates$Indicator_VCV
  Sigma <- fitted(.object)
  0.5*sum((S - Sigma)[lower.tri(S, diag = FALSE)]^2)
}

### Geodesic distance

dG <- function(.object) {
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fitted(.object)
  Eigen     <- eigen(solve(S) %*% Sigma_hat)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   
  # not sure if logarithm naturalis is used or logarithm with base 10. 
  
  0.5 * sum(logEigenvaluessq)
}

SRMR <- function(.object) {
  
  # The SRMR as calculated by us is always based on the the difference between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fitted(.object)
  
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag=T)))
}


### dML: the fitting function used in FIML

dML <- function(.object){
  
  nobs <- .object$Information$Number_of_observations
  S         <- .object$Estimates$Indicator_VCV
  p         <- dim(S)[1]
  Sigma_hat <- fitted(.object)
  (nobs - 1)*(log(det(Sigma_hat)) + sum(diag(S %*% solve(Sigma_hat))) - log(det(S)) - p)
}