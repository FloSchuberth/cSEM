#' Calculate weights using GSCA
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsGSCA(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsGSCA <- function(
  .data                        = args_default()$.data,
  .model                       = args_default()$.model,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  
  ## Calculation of the actual weights, coefficients and loadings 
  
  Z <- X # Z is the data matrix in GSCA, data are already standardized
  W0 <- t(csem_model$measurement) # Matrix of the weighted relation model
  B0 <- t(csem_model$structural) # Matrix of the structural model
  C0 <- csem_model$measurement # Matrix of the measurement model if all indicators are reflective
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W0) # number of indicators
  J = ncol(W0) # number of constructs
  T = K + J
  
  A0 = cbind(C0, B0) # J rows, T columns
  
  ## get indices of those entries of A0 and W0 which are not zero 
  
  vecA <- matrix(A0, ncol = 1)
  vecW <- matrix(W0, ncol = 1)
  aindex <- which(vecA!=0,arr.ind = T)[,1]
  windex <- which(vecW!=0,arr.ind = T)[,1]
  
  ## set initial values of non-zero entries in A0 and W0 to random values
  vecA[aindex,] = runif(length(aindex), min = 0, max = 1)
  vecW[windex,] = runif(length(windex), min = 0, max = 1)
  # A = A0
  # W = W0
  A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
  W = matrix(vecW, ncol = J, nrow = K, byrow = FALSE)
  
  ## preparation of ALS-algorithm
  V = cbind(diag(K), W)
  Gamma = Z %*% W
  Psi = Z %*% V
  vecPsi = matrix(Psi, ncol = 1)
  
  iter = 0
  f0 = 100000
  imp = 100000
  
  ## actual ALS-algorithm begins right here
  while (iter < .iter_max & imp > .tolerance) {
    
    iter = iter + 1
    
    ## step 1: update A (loadings and path coefficients)
    
    Phi = kronecker(diag(T), Gamma)
    Phi = Phi[,aindex]
    a_hat = solve((t(Phi)%*%Phi))%*%t(Phi)%*%vecPsi
    vecA[aindex,1] = a_hat
    A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
    vecA <- matrix(A, ncol = 1)
    
    ## step 2: update W (weights)
    
    for (j in 1:J) {
      t = K + j
      w0 = W[,j, drop = FALSE]
      windex_j = which(w0!=0,arr.ind = T)[,1]
      m <- rep(0, ncol(A))
      m[t] = 1
      a = A[j,, drop = FALSE]
      beta = m - a
      H1 = diag(J)
      H1[j,j] = 0
      H2 = diag(T)
      H2[t,t] = 0
      Delta = W%*%H1%*%A - V%*%H2
      vecZDelta = matrix(Z%*%Delta, ncol = 1)
      XI = kronecker(t(beta), Z)
      XI = XI[,windex_j, drop = FALSE]
      theta_j_hat = solve(t(XI) %*% XI)%*%t(XI)%*%vecZDelta
      zw = Z[,windex_j, drop = FALSE] %*% theta_j_hat
      theta = sqrt(N)*theta_j_hat/norm(zw, type="2") 
      W[windex_j,j] = theta
      V[windex_j,t] = theta
      
    } 
    
    ## for-loop is over. Thus, the matrix A, W and V are updated and consequently
    ## all unknown parameters are estimated. Now, the value of the minimization
    ## criterion is calculated to check whether convergence is already valid
    
    ## updated matrices Gamma and Psi
    Gamma = Z %*% W
    Psi = Z %*% V
    dif = Psi - Gamma %*% A # This is the matrix of the optimization criterion 
    f = matrixcalc::matrix.trace(t(dif) %*% dif)
    imp = f0 - f # decrease of the minimization criterion
    f0 = f
    vecPsi = matrix(Psi, ncol = 1)
  }
  
  B <- A[,(K+1):T] # final matrix of path coefficients
  C <- A[,1:K] # final matrix of loadings
  
  # check for convergence
  if(imp <= .tolerance){
    Conv_status = TRUE  
  } else{ # the maximal number of iterations was reached without convergence 
    Conv_status = FALSE
  }
  
  ## Calculate proxies/scores
  H <- Z %*% W
  
  ## Calculate proxy covariance matrix
  proxyCV = t(W) %*% S %*% W
  
  ## Calculate global and local measures of fit
  FIT = 1 - matrixcalc::matrix.trace(t(dif) %*% dif)/matrixcalc::matrix.trace(t(Psi) %*% Psi)
  unknown = length(aindex) + length(windex) # total number of unknown parameters in the model
  AFIT = 1 - (1 - FIT)*(N*K)/(N*K-unknown)
  
  FIT_M = 1 - matrixcalc::matrix.trace(t(Z - Gamma %*% C) %*% (Z - Gamma %*% C))/(N*K) 
  FIT_S = 1 - matrixcalc::matrix.trace(t(Gamma - Gamma %*% B) %*% (Gamma - Gamma %*% B))/(N*J)
  
  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = t(B),
      "Loading_estimates"      = t(C),
      "Weight_estimates"       = t(W),
      "Construct_scores"       = H,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = proxyCV,
      "FIT"                    = FIT,
      "AFIT"                   = AFIT,
      "FIT_M"                  = FIT_M,
      "FIT_S"                  = FIT_S
    ),
    "Information" = list(
      "Data"          = Z,
      "Model"         = csem_model,
      "Arguments"     = as.list(match.call())[-1],
      "Weight_info"   = list(
        "Number_iterations"  = iter,
        "Convergence_status" = Conv_status
      )
    )
  )
  
  class(out) <- "cSEMResults"
  attr(out, "single") <- TRUE
  invisible(out)
  
} # END calculateWeightsGSCA


#' Calculate weights using GSCA and variance-covariance matrices
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsGSCAVCV(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsGSCAVCV <- function(
  .data                        = args_default()$.data,
  .model                       = args_default()$.model,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  
  ## Calculation of the actual weights, coefficients and loadings 
  
  Z <- X # Z is the data matrix in GSCA, data are already standardized
  W0 <- t(csem_model$measurement) # Matrix of the weighted relation model
  B0 <- t(csem_model$structural) # Matrix of the structural model
  C0 <- csem_model$measurement # Matrix of the measurement model if all indicators are reflective
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W0) # number of indicators
  J = ncol(W0) # number of constructs
  T = K + J
  
  ## in this approach, the variance-covariance matrix of the indicators is used instead
  ## of the data matrix Z. For this, a square root decomposition of M is needed
  ## which is obtained via a cholesky decomposition
  M <- (N-1) * S
  H <- chol(M)  ## upper triangular
  H <- t(H)  ## lower triangular: H %*% t(H) = t(Z) %*% Z = M
  
  A0 = cbind(C0, B0) # J rows, T columns
  
  ## get indices of those entries of A0 and W0 which are not zero 
  
  vecA <- matrix(A0, ncol = 1)
  vecW <- matrix(W0, ncol = 1)
  aindex <- which(vecA!=0,arr.ind = T)[,1]
  windex <- which(vecW!=0,arr.ind = T)[,1]
  
  ## set initial values of non-zero entries in A0 and W0 to random values
  vecA[aindex,] = runif(length(aindex), min = 0, max = 1)
  vecW[windex,] = runif(length(windex), min = 0, max = 1)
  # A = A0
  # W = W0
  A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
  W = matrix(vecW, ncol = J, nrow = K, byrow = FALSE)
  
  ## preparation of ALS-algorithm
  V = cbind(diag(K), W)
  Gamma = t(H) %*% W
  Psi = t(H) %*% V
  vecPsi = matrix(Psi, ncol = 1)
  
  iter = 0
  f0 = 100000
  imp = 100000
  
  ## actual ALS-algorithm begins right here
  while (iter < .iter_max & imp > .tolerance) {
    
    iter = iter + 1
    
    ## step 1: update A (loadings and path coefficients)
    
    Phi = kronecker(diag(T), Gamma)
    Phi = Phi[,aindex]
    a_hat = solve((t(Phi)%*%Phi))%*%t(Phi)%*%vecPsi
    vecA[aindex,1] = a_hat
    A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
    vecA <- matrix(A, ncol = 1)
    
    ## step 2: update W (weights)
    
    for (j in 1:J) {
      t = K + j
      w0 = W[,j, drop = FALSE]
      windex_j = which(w0!=0,arr.ind = T)[,1]
      m <- rep(0, ncol(A))
      m[t] = 1
      a = A[j,, drop = FALSE]
      beta = m - a
      H1 = diag(J)
      H1[j,j] = 0
      H2 = diag(T)
      H2[t,t] = 0
      Delta = W%*%H1%*%A - V%*%H2
      vecZDelta = matrix(t(H)%*%Delta, ncol = 1)
      XI = kronecker(t(beta), t(H))
      XI = XI[,windex_j, drop = FALSE]
      theta_j_hat = solve(t(XI) %*% XI)%*%t(XI)%*%vecZDelta
      zw = Z[,windex_j, drop = FALSE] %*% theta_j_hat
      theta = sqrt(N)*theta_j_hat/norm(zw, type="2") 
      W[windex_j,j] = theta
      V[windex_j,t] = theta
      
    } 
    
    ## for-loop is over. Thus, the matrix A, W and V are updated and consequently
    ## all unknown parameters are estimated. Now, the value of the minimization
    ## criterion is calculated to check whether convergence is already valid
    
    ## updated matrices Gamma and Psi
    Gamma = t(H) %*% W
    Psi = t(H) %*% V
    dif = Psi - Gamma %*% A # This is the matrix of the optimization criterion 
    f = matrixcalc::matrix.trace(t(dif) %*% dif)
    imp = f0 - f # decrease of the minimization criterion
    f0 = f
    vecPsi = matrix(Psi, ncol = 1)
  }
  
  B <- A[,(K+1):T] # final matrix of path coefficients
  C <- A[,1:K] # final matrix of loadings
  
  # check for convergence
  if(imp <= .tolerance){
    Conv_status = TRUE  
  } else{ # the maximal number of iterations was reached without convergence 
    Conv_status = FALSE
  }
  
  ## Calculate proxies/scores
  proxies <- Z %*% W
  
  ## Calculate proxy covariance matrix
  proxyCV = t(W) %*% S %*% W
  
  ## Calculate global and local measures of fit
  FIT = 1 - matrixcalc::matrix.trace(t(dif) %*% dif)/matrixcalc::matrix.trace(t(Psi) %*% Psi)
  unknown = length(aindex) + length(windex) # total number of unknown parameters in the model
  AFIT = 1 - (1 - FIT)*(N*K)/(N*K-unknown)
  
  FIT_M = 1 - matrixcalc::matrix.trace(t(t(H) - Gamma %*% C) %*% (t(H) - Gamma %*% C))/(N*K) 
  FIT_S = 1 - matrixcalc::matrix.trace(t(Gamma - Gamma %*% B) %*% (Gamma - Gamma %*% B))/(N*J)
  
  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = t(B),
      "Loading_estimates"      = t(C),
      "Weight_estimates"       = t(W),
      "Construct_scores"       = proxies,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = proxyCV,
      "FIT"                    = FIT,
      "AFIT"                   = AFIT,
      "FIT_M"                  = FIT_M,
      "FIT_S"                  = FIT_S
    ),
    "Information" = list(
      "Data"          = Z,
      "Model"         = csem_model,
      "Arguments"     = as.list(match.call())[-1],
      "Weight_info"   = list(
        "Number_iterations"  = iter,
        "Convergence_status" = Conv_status
      )
    )
  )
  
  class(out) <- "cSEMResults"
  attr(out, "single") <- TRUE
  invisible(out)
  
} # END calculateWeightsGSCAVCV