require(cSEM)

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
  
  Z <- X # Z is the data matrix in GSCA
  W0 <- t(csem_model$measurement) # Matrix of the weighted relation model
  B0 <- t(csem_model$structural) # Matrix of the structural model -> abklären, wie die aussieht!!
  C0 <- csem_model$measurement # Matrix of the measurement model; wenn nur reflektive Indikatoren, was ist sonst??
  
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
    # Frage: Konvergenz ggf. über checkConvergence berechnen? was soll später ausgegeben werden?
    
    ## updated matrix Gamma and Psi
    Gamma = Z %*% W
    Psi = Z %*% V
    dif = Psi - Gamma %*% A # This is the matrix of the optimization criterion 
    f = matrixcalc::matrix.trace(t(dif) %*% dif)
    imp = f0 - f
    f0 = f
    vecPsi = matrix(Psi, ncol = 1)
  }
  
  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = "here path",
      "Loading_estimates"      = "here loadings",
      "Weight_estimates"       = W$W,
      "Construct_scores"       = H,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = C,
      "Construct_VCV"          = P,
      "Correction_factors"     = correction_factors,
      "FIT"                     = "here FIT"
    ),
    "Information" = list(
      "Data"          = X,
      "Model"         = csem_model,
      "Arguments"     = as.list(match.call())[-1],
      "Weight_info"   = list(
        "Number_iterations"  = W$Iterations,
        "Convergence_status" = W$Conv_status
      )
    )
  )
  
  class(out) <- "cSEMResults"
  attr(out, "single") <- TRUE
  invisible(out)
  
} # END calculateWeightsGSCA