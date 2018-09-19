#' Calculate weights using GSCAm
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

calculateWeightsGSCAm <- function(
  .data                        = args_default()$.data,
  .model                       = args_default()$.model,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  
  # Use the usual GSCA procedure to obtain initial values for Gamma, C and B
  GSCA_results <-  calculateWeightsGSCA(
    .data                     = X,
    .model                    = csem_model,
    .iter_max                 = .iter_max,
    .tolerance                = .tolerance
  )
  
  W <- t(GSCA_results$Estimates$Weight_estimates) # Matrix of the weighted relation model
  B <- t(GSCA_results$Estimates$Path_estimates) # Matrix of the structural model
  C <- t(GSCA_results$Estimates$Loading_estimates) # Matrix of the measurement model if all indicators are reflective
  Gamma <- GSCA_results$Estimates$Construct_scores # Matrix of scores of latent variables
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W) # number of indicators
  J = ncol(W) # number of constructs
  T = K + J
  
  A = cbind(C, B) # J rows, T columns
  D <- diag(runif(K, min = 0, max = 1)) # matrix of unique loadings, random inital values
  U <- matrix(0, N, K) # matrix of unique variables
  
  ## get indices of those entries of A which are not zero 
  
  vecA <- matrix(A, ncol = 1)
  aindex <- which(vecA!=0,arr.ind = T)[,1]
  vecW <- matrix(W, ncol = 1)
  windex <- which(vecW!=0,arr.ind = T)[,1]
  
  ## Preparation of the calculation of the actual weights, coefficients and loadings 
  
  ## prepare data for GSCAm
  Z <- X # Z is the data matrix in GSCAm, data are already standardized
  # for GSCAm, normalized data are needed
  vZ = sqrt(diag(t(Z) %*% Z))
  factorZ = matrix(vZ, nrow=N, ncol=length(vZ), byrow=TRUE)
  Z = Z/factorZ # normalized matrix Z
  
  ## prepare latent variables (constructs) for GSCAm
  v = sqrt(diag(t(Gamma) %*% Gamma))
  factor = matrix(v, nrow=N, ncol=length(v), byrow=TRUE)
  Gamma = Gamma/factor # normalized matrix Gamma of latent variables
  
  # calculate initial values for the unique variables in U
  # analogous to step 3 of the modified ALS-algorithm
  Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, N-J, drop = FALSE]
  s <- svd(t(Gamma_orth) %*% Z %*% D)
  U_tilde = s$u %*% t(s$v)
  U = Gamma_orth %*% U_tilde # this is the initial matrix U
  
  iter = 0
  f0 = 100000
  imp = 100000
  
  ## actual modified ALS-algorithm of GSCAm begins right here

  while (iter < .iter_max & imp > .tolerance) {
    
    iter = iter + 1
    
    ## step 1: update Gamma (latent variable scores) and W (weights)
    
    H = A - cbind(matrix(0, nrow = J, ncol = K), diag(J))
    M = cbind(diag(K), matrix(0, nrow = K, ncol = J))
    W = M %*% t(H) %*% solve(H %*% t(H))
    Gamma = (Z - U %*% D) %*% W
    
    # normalize every latent variable and adjust the corresponding weights
    v = sqrt(diag(t(Gamma) %*% Gamma)) # length of v is J
    factor = matrix(v, nrow=N, ncol=length(v), byrow=TRUE)
    Gamma = Gamma/factor # normalized matrix of latent variables
    W = W/matrix(v, nrow=K, ncol=length(v), byrow=TRUE)
    
    ## step 2: update A (loadings and path coefficients)
    
    Sm = cbind(U %*% D, matrix(0, nrow = N, ncol = J))
    Psi = cbind(Z, Gamma)
    PsiStar = Psi - Sm
    vecPsi = matrix(Psi, ncol = 1)
    vecPsiStar = matrix(PsiStar, ncol = 1)
    
    Phi = kronecker(diag(T), Gamma)
    Phi = Phi[,aindex]
    a_hat = solve((t(Phi)%*%Phi))%*%t(Phi)%*%vecPsiStar
    vecA[aindex,1] = a_hat
    A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
    vecA <- matrix(A, ncol = 1)
    
    ## step 3: update U (unique parts) for fixed Gamma, A and D
    ## based on Trendafilov et al.'s procedure (2013)
    
    # Gamma_orth is an N by N-J orthonormal basis matrix of the null space of Gamma
    Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, N-J, drop = FALSE]
    s <- svd(t(Gamma_orth) %*% Z %*% D)
    U_tilde = s$u %*% t(s$v)
    U = Gamma_orth %*% U_tilde
  
    ## step 4: update D (unique loadings) for fixed Gamma, A and U 
    
    D = diag(diag(t(U) %*% Z))
    
    ## updated matrices Psi and Sm
    Psi = cbind(Z, Gamma)
    Sm = cbind(U %*% D, matrix(0, nrow = N, ncol = J))
    dif = Psi - Gamma %*% A - Sm # This is the matrix of the optimization criterion 
    f = matrixcalc::matrix.trace(t(dif) %*% dif)
    imp = f0 - f # decrease of the minimization criterion
    f0 = f
    PsiStar = Psi - Sm
    vecPsi = matrix(Psi, ncol = 1)
    vecPsiStar = matrix(PsiStar, ncol = 1)
    
  }
  
  B <- A[,(K+1):T] # final matrix of path coefficients
  C <- A[,1:K] # final matrix of loadings
  
  # check for convergence
  if(imp <= .tolerance){
    Conv_status = TRUE  
  } else{ # the maximal number of iterations was reached without convergence 
    Conv_status = FALSE
  }
  
  ## Calculate final proxies/scores (oder ist das einfach Gamma?)
  proxies <- (Z - U %*% D) %*% W
  
  ## Calculate proxy covariance matrix (S is the VCV-matrix of the indicators)
  proxyCV = t(W) %*% S %*% W
  
  ## Calculate global and local measures of fit
  FIT = 1 - matrixcalc::matrix.trace(t(dif) %*% dif)/matrixcalc::matrix.trace(t(Psi) %*% Psi)
  
  FIT_M = 1 - matrixcalc::matrix.trace(t(Z - Gamma %*% C - U %*% D) %*% (Z - Gamma %*% C - U %*% D))/K
  FIT_S = 1 - matrixcalc::matrix.trace(t(Gamma - Gamma %*% B) %*% (Gamma - Gamma %*% B))/J
  
  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = t(B),
      "Loading_estimates"      = t(C),
      "Weight_estimates"       = t(W),
      "Construct_scores"       = proxies,
      "Unique_variables"       = U,
      "Unique_loadings"        = D,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = proxyCV,
      "FIT"                    = FIT,
      "FIT_M"                  = FIT_M,
      "FIT_S"                  = FIT_S
    ),
    "Information" = list(
      "Data"          = X,
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
