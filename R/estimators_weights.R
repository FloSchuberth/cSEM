#' Calculate composite weights using PLS-PM
#'
#' Calculates composite weights using the PLS-PM algorithm. (TODO)
#'
#' More details here. (TODO)
#'
#' @usage calculateWeightsPLS(
#'   .S                           = args_default()$.S,
#'   .csem_model                  = args_default()$.csem_model,
#'   .conv_criterion              = args_default()$.conv_criterion,
#'   .iter_max                    = args_default()$.iter_max,
#'   .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model
#'   .PLS_modes                   = args_default()$.PLS_modes,
#'   .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
#'   .tolerance                   = args_default()$.tolerance
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A list with the elements
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{A (J x J) matrix of inner weights.}
#'   \item{`$Modes`}{A named vector of Modes used for the outer estimation.}
#'   \item{`$Conv_status`}{The convergence status. `TRUE` if the algorithm has converged 
#'     and `FASLE` otherwise. If one-step weights are used via `.iter_max = 1` 
#'     or a non-iterative procedure was used, the convergence status is set to `NULL`.}
#'   \item{`$Iterations`}{The number of iterations used.}
#' }
#' @export
#'

calculateWeightsPLS <- function(
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes                   = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
  .tolerance                   = args_default()$.tolerance
) {

  ### Preparation ==============================================================
  ## Get/set the modes for the outer estimation
  modes <- ifelse(.csem_model$construct_type == "Common factor", "ModeA", "ModeB")
  
  if(!is.null(.PLS_modes)) {
    # Error if other than "ModeA" or "ModeB"
    if(!all(.PLS_modes %in% c("ModeA", "ModeB"))) {
      stop(paste0("`", setdiff(.PLS_modes, c("ModeA", "ModeB")), "`", collapse = ", "),
           " in `.PLS_modes` is an unknown mode.", call. = FALSE)
    }
    # Error if construct names provided do not match the constructs of the model
    if(length(names(.PLS_modes)) != 0 &&
       length(setdiff(names(.PLS_modes), names(.csem_model$construct_type))) > 0) {
      stop(paste0("`", setdiff(names(.PLS_modes), names(csem_model$construct_type)), 
                  "`", collapse = ", ")," in `.PLS_modes` is an unknown construct name.", call. = FALSE)
    }
    # If only "ModeA" or "ModeB" is provided without set all of the modes to that mode.
    if(length(names(.PLS_modes)) == 0) {
      if(length(.PLS_modes) == 1) {
        modes <- rep(.PLS_modes, length(.csem_model$construct_type))
        names(modes) <- names(.csem_model$construct_type)
      } else {
        stop("Only a vector of `name = value` pairs or a single mode may be provided to `.PLS_modes`.", 
             call. = FALSE)
      }
    } else {
      # Replace modes if necessary and keep the others at their defaults
      modes[names(.PLS_modes)] <- .PLS_modes
    }
  }

  ### Calculation/Iteration ====================================================
  W <- .csem_model$measurement
  # Scale weights
  W <- scaleWeights(.S = .S, .W = W)

  W_iter       <- W
  tolerance    <- .tolerance
  iter_max     <- .iter_max
  iter_counter <- 0

  repeat {
    # Counter
    iter_counter <- iter_counter + 1

    # Inner estimation
    E <- calculateInnerWeightsPLS(
      .S                           = .S,
      .W                           = W_iter,
      .csem_model                  = .csem_model,
      .PLS_ignore_structural_model = .PLS_ignore_structural_model,
      .PLS_weight_scheme_inner     = .PLS_weight_scheme_inner
    )
    # Outer estimation

    W <- calculateOuterWeightsPLS(
      .S        = .S,
      .W        = W_iter,
      .E        = E,
      .modes    = modes
    )

    # Scale weights
    W <- scaleWeights(.S, W)

    # Check for convergence
    conv <- checkConvergence(W, W_iter, 
                             .conv_criterion = .conv_criterion, 
                             .tolerance = .tolerance)
    
    if(conv) {
      # Set convergence status to TRUE as algorithm has converged
      conv_status = TRUE
      break # return iterative PLS-PM weights
      
    } else if(iter_counter == iter_max & iter_max == 1) {
      # Set convergence status to NULL, NULL is used if no algorithm is used
      conv_status = NULL
      break # return one-step PLS-PM weights
      
    } else if(iter_counter == iter_max & iter_max > 1) {
      # Set convergence status to FALSE, as algorithm has not converged
      conv_status = FALSE
      break
      
    } else {
      W_iter <- W
    }
  }

  # Return
  l <- list("W" = W, "E" = E, "Modes" = modes, "Conv_status" = conv_status,
            "Iterations" = iter_counter)
  return(l)

  ### For maintenance: ### -----------------------------------------------------
  # W (J x K) := (Block-)diagonal matrix of weights with the same dimension as .csem_model$measurement
  # C (J x J) := Empirical composite correlation matrix
  # E (J x J) := Matrix of inner weights
  # D (J x J) := An "adjacency" matrix with elements equal to one, if proxy j and j' are adjacent
} # END calculateWeightsPLS

#' Calculate composite weights using one of Kettenrings's approaches
#'
#' Calculates weights according to on of the the five criteria
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*", and "*GENVAR*"
#' suggested by \insertCite{Kettenring1971;textual}{cSEM}.
#'
#' Some more description...
#'
#' @usage calculateWeightsKettenring(
#'   .S           = args_default()$.S, 
#'   .csem_model  = args_default()$.csem_model,   
#'   .approach    = args_default()$.approach
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#' 

calculateWeightsKettenring <- function(
  .S              = args_default()$.S,
  .csem_model     = args_default()$.csem_model,
  .approach       = args_default()$.approach
) {
  
  ### Preparation ==============================================================
  ## Check if all constructs are modeled as composites
  if("Common factor" %in% .csem_model$construct_type){
    stop("Currently Kettenring's approaches are only allowed for pure composite models.", 
         call. = FALSE)
  }

  ## Get relevant objects 
  W <- .csem_model$measurement
  
  construct_names <- rownames(W)
  indicator_names <- colnames(W)
  indicator_names_ls <- lapply(construct_names, function(x) {
    indicator_names[W[x, ] == 1]
  })
  
  ## Calculate S_{ii}^{-1/2}}
  sqrt_S_jj_list <- lapply(indicator_names_ls, function(x) {
    solve(expm::sqrtm(.S[x, x, drop = FALSE]))
  })
  
  ## Put in a diagonal matrix
  H <- as.matrix(Matrix::bdiag(sqrt_S_jj_list))
  dimnames(H) <- list(indicator_names, indicator_names)
  
  if(.approach %in% c("MAXVAR", "MINVAR")) {
    # Note: The MAXVAR implementation is based on a MATLAB code provided 
    #       by Theo K. Dijkstra
    Rs <- H %*% .S %*% H
    
    ## Calculate eigenvalues and eigenvectors of Rs and ...
    u <- switch (.approach,
          "MINVAR" = {
            # ... select the eigenvector corresponding to the smallest eigenvalue
            as.matrix(eigen(Rs)$vectors[, length(Rs[,1])], nocl = 1)
          },
          "MAXVAR" = {
            # ... select the eigenvector corresponding to the largest eigenvalue
            as.matrix(eigen(Rs)$vectors[, 1], nocl = 1)
          }
    )
    rownames(u) <- indicator_names

    ## Calculate weight
    W <- lapply(indicator_names_ls, function(x) {
      
      w <- as.vector((H[x, x] %*% u[x, ]) / sqrt(sum(c(u[x, ])^2)))
      
      return(w)
    })
  } else if(.approach %in% c("SUMCORR", "GENVAR", "SSQCORR")) {
    
    # Maximize the sum of the composite correlations sum(w_i%*%S_ii%*%t(w_i)), 
    #   see Asendorf (2015, Appendix B.3.2 Empirical, p. 295)
    
    ## Define function to be minimized:
    # R       := Empirical indicator correlation matrix (S)
    # RDsqrt  := Blockdiagonal matrix with S_{ii}^{-1/2} on its diagonal (H)
    # nameInd := List with indicator names per block (indicator_names_ls)
    # name    := Vector of composites names (construct_names)
    
    fn <- switch (.approach,
      "SUMCORR" = {
        function(wtilde, R, RDsqrt, nameInd, nameLV) {
          -(t(wtilde) %*% RDsqrt %*% R %*% RDsqrt %*% t(t(wtilde)))
        }
      },
      "GENVAR"  = {
        function(wtilde, R, RDsqrt, nameInd, nameLV) {
          det(-(t(wtilde) %*% RDsqrt %*% R %*% RDsqrt %*% t(t(wtilde))))
        }
      },
      "SSQCORR" = {
        function(wtilde, R, RDsqrt, nameInd, nameLV) {
          norm(-(t(wtilde) %*% RDsqrt %*% R %*% RDsqrt %*% t(t(wtilde))), 
               type = "F")
        }
      }
    )
    
    # Define constraint function: variance of each composite is one
    heq <- function(wtilde, R, RDsqrt, nameInd, nameLV) {
      names(wtilde) = unlist(nameInd)
      h <- rep(NA, length(nameLV))
      
      for (i in 1:length(h)) {
        h[i] <- t(wtilde[nameInd[[i]]]) %*% t(t(wtilde[nameInd[[i]]])) - 1
      }
      h
    }
    
    ## Optimization
    res <- alabama::auglag(
      par     = runif(length(indicator_names)),
      fn      = fn,
      R       = .S,
      RDsqrt  = H,
      nameInd = indicator_names_ls,
      nameLV  = construct_names,
      heq     = heq
    )
    
    # Get wtilde
    wtilde <- res$par
    names(wtilde) <- indicator_names
    
    # Calculate weights w_i = R_{ii}^{-1/2} wtilde_i
    W <- lapply(indicator_names_ls, function(x) {
      
      w <- as.vector(H[x, x] %*% wtilde[x])

      return(w)
    })
  }
  
  ## Format, name and scale
  W <- t(as.matrix(Matrix::bdiag(W)))
  dimnames(W) <- list(construct_names, indicator_names)
  W <- scaleWeights(.S, W)

  ## Prepare for output
  modes <- rep(.approach, length(construct_names))
  names(modes) <- construct_names
  
  ## Prepare output and return
  l <- list(
    "W"           = W, 
    "E"           = NULL, 
    "Modes"       = modes,
    "Conv_status" = if(.approach %in% c("MINVAR", "MAXVAR")) NULL else res$convergence == 0,  
    "Iterations"  = if(.approach %in% c("MINVAR", "MAXVAR")) 0 else res$counts[1]
    )
  
  return(l)
  
} # END calculateWeightsKettenring


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
  .S                           = args_default()$.S,
  .model                       = args_default()$.model,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  
  ## Calculation of the actual weights, coefficients and loadings 
  
  Z <- .data # Z is the data matrix in GSCA, data are already standardized
  W0 <- t(.model$measurement) # Matrix of the weighted relation model
  B0 <- t(.model$structural) # Matrix of the structural model
  C0 <- .model$measurement # Matrix of the measurement model if all indicators are reflective
  
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
  S <- .S
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
      "Model"         = .model,
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


#' Calculate composite weights using fixed weights
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsFixed(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsFixed <- function(
  .data         = NULL,
  .model        = NULL
) {

  stop("Not yet implemented")

} # END calculateWeightsFixed


#' Calculate composite weights using unit weights
#'
#' Calculate unit weights for all blocks, i.e., each indicator of a block is
#' equally weighted.
#'
#' @usage calculateWeightsUnit(
#'  .S                 = args_default()$.S,
#'  .csem_model        = args_default()$.csem_model
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'
calculateWeightsUnit = function(
  .S                 = args_default()$.S,
  .csem_model        = args_default()$.csem_model
){
  
  W <- .csem_model$measurement
  W <- scaleWeights(.S, W)
  
  modes        <- rep("unit", nrow(W))
  names(modes) <- rownames(W)
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = modes, "Conv_status" = NULL,
            "Iterations" = 0)
  return(l)  
}
