# Distance functions

### Squared euclidean distance

dL <- function(.A,.B) {
  
  # S     <- .object$Estimates$Indicator_VCV
  # Sigma <- fitted(.object)
  0.5*sum((.A - .B)[lower.tri(.A, diag = FALSE)]^2)
}

### Geodesic distance

dG <- function(.A,.B) {
  
  # S         <- .object$Estimates$Indicator_VCV
  # Sigma_hat <- fitted(.object)
  Eigen     <- eigen(solve(.A) %*% .B)
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