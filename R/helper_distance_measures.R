# Distance functions

### Squared euclidean distance

dL <- function(.matrix1=args_default()$.matrix1,.matrix2=args_default()$.matrix2) {
  
  if(!identical(dim(.matrix1),dim(.matrix2))){
    stop("The provided matrices are not identical.")
    }

  0.5*sum((.matrix1 - .matrix2)[lower.tri(.matrix1, diag = FALSE)]^2)
}

### Geodesic distance

dG <- function(.matrix1=args_default()$.matrix1,.matrix2=args_default()$.matrix2) {
  
  if(!identical(dim(.matrix1),dim(.matrix2))){
    stop("The provided matrices are not identical.")
  }
  
  Eigen     <- eigen(solve(.matrix1) %*% .matrix2)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   
  # not sure if logarithm naturalis is used or logarithm with base 10. 
  
  0.5 * sum(logEigenvaluessq)
}

SRMR <- function(.object=args_default()$.object) {
  
  # The SRMR as calculated by us is always based on the the difference between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fitted(.object)
  
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag=T)))
}


### dML: the fitting function used in FIML

dML <- function(.object=args_default()$.object){
  
  nobs <- .object$Information$Number_of_observations
  S         <- .object$Estimates$Indicator_VCV
  p         <- dim(S)[1]
  Sigma_hat <- fitted(.object)
  (nobs - 1)*(log(det(Sigma_hat)) + sum(diag(S %*% solve(Sigma_hat))) - log(det(S)) - p)
}