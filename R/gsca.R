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
  
  if(.model$model_type != "linear")  {
    stop("Not yet implemented")
  }
  
  ### Make the function more autonomous ========================================
  ## Convert model to "cSEMModel" format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    csem_model <- parseModel(.model)
  } else {
    csem_model <- .model
  }
  
  ## Prepare, standardize, check, and clean data if not already in this format
  if(!(class(.data) == "cSEMData")) {
    if(is.matrix(.data) && isSymmetric.matrix(.data)) {
      S <- .data
    } else {
      #  If function is used externally, data is not automatically scaled!
      X <- processData(.data = .data, .model = csem_model)
      S <- stats::cov(X)
    }
  }
  
  ## End preparation of data
  
  ## Calculation of the actual weights, coefficients and loadings 
  
  Z <- X # Z is the data matrix in GSCA
  W <- t(csem_model$measurement) # Matrix of the weighted relation model
  B0 <- t(csem_model$structural) # Matrix of the structural model -> abklären, wie die aussieht!!
  C0 <- csem_model$measurement # Matrix of the measurement model; wenn nur reflektive Indikatoren, was ist sonst??
  
  N = nrow(Z) # number of observations per indicator
  K = nrow(W) # number of indicators
  J = ncol(W) # number of constructs
  T = K + J
  
  A0 = cbind(C0, B0) # J rows, T columns
  
  ## get indices of those entries of A0 and W0 which are not zero 
  
  vecA <- matrix(A0, ncol = 1)
  vecW <- matrix(W0, ncol = 1)
  aindex <- which(vecA!=0,arr.ind = T)[,1]
  windex <- which(vecW!=0,arr.ind = T)[,1]
  
  ## standardize data if not happened yet
  Z <- scale(Z) # Frage: wird bei scale durch die Wurzel der normalen oder der mod. Stichprobenvarianz geteilt?
  
  ## set initial values of non-zero entries in A0 and W0 to random values
  vecA[aindex,] = rnorm(length(aindex), mean = 0, sd = 1)
  vecW[windex,] = rnorm(length(windex), mean = 0, sd = 1)
  A = A0
  W = W0
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
  while (iter < .iter_max && imp > .tolerance) {
    
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
      w0 = W[,j]
      windex_j = which(w0!=0,arr.ind = T)[,1]
      m <- rep(0, ncol(A))
      m[t] = 1
      a = A[j,]
      beta = m - a
      H1 = diag(J)
      H1[j,j] = 0
      H2 = diag(T)
      H2[t,t] = 0
      Delta = W%*%H1%*%A - V%*%H2
      vecZDelta = matrix(Z%*%Delta, ncol = 1)
      XI = kronecker(t(beta), Z)
      XI = XI[,windex_j]
      theta_j_hat = solve(t(XI) %*% XI)%*%t(XI)%*%vecZDelta
      zw = Z[,windex_j]%*%theta_j_hat
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
    f = tr(t(dif) %*% dif)
    imp = f0 - f
    vecPsi = matrix(Psi, ncol = 1)
  }
  
  ### Output -------------------------------------------------------------------
  # was ist von Interesse, was soll ausgegeben oder ggf. noch berechnet werden?
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = if(.estimate_structural) {
        estim_results$Path_estimates
      } else {
        estim_results
      },
      "Loading_estimates"      = Lambda * csem_model$measurement,
      "Weight_estimates"       = W,
      "Inner_weight_estimates" = NULL,
      "Construct_scores"       = H,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = C,
      "Construct_VCV"          = P,
      "Cross_loadings"         = Lambda,
      "Construct_reliabilities"= Q^2,
      "Correction_factors"     = correction_factors,
      "R2"                     = if(.estimate_structural) {
        estim_results$R2
      } else {
        estim_results
      }
    ),
    "Information" = list(
      "Data"          = X,
      "Model"         = csem_model,
      "Arguments"     = as.list(match.call())[-1],
      "Weight_info"   = list(
        "Modes"              = W$Modes,
        "Number_iterations"  = W$Iterations,
        "Convergence_status" = W$Conv_status
      )
    )
  )
  
  class(out) <- "cSEMResults"
  invisible(out)
  
  
  
} # END calculateWeightsGSCA