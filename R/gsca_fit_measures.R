#' Calculate the SRMR
#'
#' Calculates SRMR...
#'
#' Some more description...
#'
#' @usage calculateWeightsGSCA(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateSRMR <- function(.object = args_default()$.object) {
  
  Z <- .object$Information$Data
  S <- .object$Estimates$Indicator_VCV
  B <- t(.object$Estimates$Path_estimates)
  C <- t(.object$Estimates$Loading_estimates)
  W <- t(.object$Estimates$Weight_estimates)
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W) # number of indicators
  J = ncol(W) # number of constructs
  T = K + J
  
  A <- cbind(C, B)
  V <- cbind(diag(K), W)
  
  Q = t(t(V) - t(A) %*% t(W))
  Omega = solve(Q %*% t(Q)) %*% Q
  e = t(Q) %*% t(Z) # Residuals
  
  Xi = stats::corp(t(e))
  Sigma_hat = Omega %*% Xi %*% t(Omega)
  diagprod <- matrix(0, nrow = K, ncol = K) 
  
  for (k in 1:K) {
    for (q in 1:k) {
      diagprod[k,q] = S[k,k] * S[q,q]
    }
  }
  C_diff = (S - Sigma_hat)/diagprod
  SRMR = sqrt(2*sum(C_diff[lower.tri(C_diff, diag = T)]^2)/((K+1)*K))
  return(SRMR)
}

#' Calculate the GFI
#'
#' Calculates GFI
#'
#' Some more description...
#'
#' @usage calculateWeightsGSCA(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateGFI <- function(.object = args_default()$.object) {
  
  Z <- .object$Information$Data
  S <- .object$Estimates$Indicator_VCV
  B <- t(.object$Estimates$Path_estimates)
  C <- t(.object$Estimates$Loading_estimates)
  W <- t(.object$Estimates$Weight_estimates)
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W) # number of indicators
  J = ncol(W) # number of constructs
  T = K + J
  
  A <- cbind(C, B)
  V <- cbind(diag(K), W)
  
  Q = t(t(V) - t(A) %*% t(W))
  Omega = solve(Q %*% t(Q)) %*% Q
  e = t(Q) %*% t(Z) # Residuals 
  
  Xi = stats::cor(t(e)) 
  Sigma_hat = Omega %*% Xi %*% t(Omega)
  C_diff = S - Sigma_hat
  
  GFI = 1- matrixcalc::matrix.trace(C_diff %*% C_diff)/(matrixcalc::matrix.trace(S %*% S))
  return(GFI)
}