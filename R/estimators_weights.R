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
  ## Get/set the default modes for the outer estimation
  modes <- as.list(ifelse(.csem_model$construct_type == "Common factor", "modeA", "modeB"))
  
  if(!is.null(.PLS_modes)) {
    
    # Error if other than "modeA", "modeB", "unit", a number, or a vector
    # of numbers of the same length as there are indicators for block j
    modes_check <- sapply(.PLS_modes, function(x) all(x %in% c("modeA", "modeB", "unit") | is.numeric(x)))
    if(!all(modes_check)) {
      stop2("The following error occured in the `calculateWeightsPLS()` function:\n",
            paste0("`", .PLS_modes[!modes_check], "`", collapse = " and "),
            " in `.PLS_modes` is an unknown mode.")
    }
    
    # Error if construct names provided do not match the constructs of the model
    if(length(names(.PLS_modes)) != 0 &&
       length(setdiff(names(.PLS_modes), names(.csem_model$construct_type))) > 0) {
      stop2("The following error occured in the `calculateWeightsPLS()` function:\n",
        paste0("`", setdiff(names(.PLS_modes), names(.csem_model$construct_type)), 
               "`", collapse = ", ")," in `.PLS_modes` is an unknown construct name.")
    }
    
    # If only "modeA", "modeB" or "unit" is provided set all of the modes to that mode.
    if(length(names(.PLS_modes)) == 0) {
      if(length(.PLS_modes) == 1) {
        modes <- as.list(rep(.PLS_modes, length(.csem_model$construct_type)))
        names(modes) <- names(.csem_model$construct_type)
      } else {
        stop2("The following error occured in the `calculateWeightsPLS()` function:\n",
          "Only a list of `name = value` pairs or a single mode may be provided to `.PLS_modes`.")
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

  # Get the modes
  modes <- sapply(modes, function(x) 
    ifelse(is.numeric(x) & length(x) > 1, "fixed", 
           ifelse(is.numeric(x) & length(x) == 1, "unit", x)))
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
  # ## Check if all constructs are modeled as composites
  # if("Common factor" %in% .csem_model$construct_type){
  #   stop("Currently Kettenring's approaches are only allowed for pure composite models.", 
  #        call. = FALSE)
  # }

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
    capture.output(res <- alabama::auglag(
      par     = runif(length(indicator_names)),
      fn      = fn,
      R       = .S,
      RDsqrt  = H,
      nameInd = indicator_names_ls,
      nameLV  = construct_names,
      heq     = heq
    ))
    
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
  .X                           = args_default()$.X,
  .csem_model                  = args_default()$.csem_model,
  .disattenuate                = args_default()$.disattenuate,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  ## Notes
  # - If disattenuate = TRUE we have Z = Z - UD. Everywhere else in the 
  #   package we assume that S is the correlation matrix of the indicators
  #   where correlation can be polyserial, bravais-person, spearman etc. 
  #   The way GSCAm is implemented now requires the calculation of cor(Z) (e.g. 
  #   step 2) Currently i use R's cor() function. If we want to allow for different
  #   correlation types we need to use our calculateIndicatorCor() function instead
  # - Currently implementation assumes that B0 is a symmetric matrix. This 
  #   is not the case if we include non-linear model. This has to be checked
  #   once the algorithm is implemented for linear models.
  #
  ### For maintenance: ---------------------------------------------------------
  # Hwang and Takane (2004, 2010, 2014, 2017) use a completely different notation
  #   compared to the standard LISREL/Bollen-type notation.
  ## N              := Number of observations
  ## J              := Total number of indicators
  ## P              := Total number of constructs
  ## TT             := K + J (T in H&T)
  ## Z (N x J)      := Matrix of indicator values (=data)
  ## W (J x P)      := (Block-)diagonal matrix of weight relations (transposed
  ##                    compared to .csem_model$measurement!!)
  ## C (P x J)      := (Block-)diagonal matrix of measurement relations (Lambda in cSEM)
  ## B (P x P)      := Matrix of path coefficients (B' in cSEM). B is not necessarily
  ##                   upper-triangular (i.e may contain feedback loops)
  ## A (P x TT)     := "RHS-builder matrix". If postmultiplied on Gamma each column 
  ##                   of the resulting matrix is the RHS of an regression equation 
  ##                   If A is known each column represents the predicted values. 
  ##                   The columns of Gamma*A therefore correspond to X*beta_t in
  ##                   standard regression notation (where depending on the
  ##                   regression equation t, beta_t contains loadings or
  ##                   path coefficients)
  ## V (J x TT)     := "LHS-builder matrix". If postmultiplied on Z, the resulting
  ##                   (N x TT) matrix Psi contains the LHS of an regression equation
  ##                   in its columns. Depending on the regression equation this 
  ##                   is either an indicator or a proxy.
  ## U (N x J)      := Matrix if unique terms/variables
  ## D (J x J)      := Diagonal matrix of unique parameters/loadings
  ## S (N x TT)     := Matrix of J columns containing predicted unique terms 
  ##                   P zero columns.
  ## Gamma (N x P)  := Gamma = ZW - S; Matrix of proxies/scores. If 
  ##                   disattenuate = FALSE S is a zero matrix.
  ## Psi (N x TT)   := Psi = ZV = Z[I ; Gamma]. Matrix of indicator and construct 
  ##                   scores ("LHS" variables" in regression)
  ## 
  ## Define matrices
  
  Z  <- .X # Z is the data matrix in GSCA, data are already standardized
  W  <- t(.csem_model$measurement) # Matrix of the weighted relation model
  C  <- .csem_model$measurement # Matrix of the measurement model (non-zero only for common factors)
  C[which(.csem_model$construct_type == "Composite"), ]  <- 0
  B  <- t(.csem_model$structural) # Matrix of the structural model
  
  N  <- nrow(Z) # number of observations
  J  <- nrow(W) # number of indicators
  P  <- ncol(W) # number of constructs
  TT <- J + P
  
  Ij <- diag(J)
  Ip <- diag(P)
  A  <- cbind(C, B) # (P x TT)
  # V  <- cbind(Ij, W) # (J x TT)
  
  U  <- matrix(0, nrow = N, ncol = J) # matrix of unique terms
  D  <- diag(runif(J, min = 0, max = 1))
  
  Gamma <- (Z - U %*% D) %*% W
  # Calculate initial values for the unique variables in U
  # analogous to step 3 of the modified ALS-algorithm
  Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (P+1):N, drop = FALSE]
  s          <- svd(t(Gamma_orth) %*% Z %*% D)
  U_tilde    <- s$u %*% t(s$v)
  U          <- Gamma_orth %*% U_tilde # this is the initial matrix U 
  
  # Prepare repeat-loop
  iter_counter <- 0
  f0           <- 100000
  imp          <- 100000
  
  repeat {
    # Counter
    iter_counter <- iter_counter + 1
    
    ## Step 1: Update W and Gamma for given A, U and D ---------------------------
    H <- A - cbind(matrix(0, nrow = P, ncol = J), Ip) # (P x T)
    M <- cbind(Ij, matrix(0, nrow = J, ncol = P) )
    
    # Estimate new weights 
    W <- t(solve(H %*% t(H)) %*% H %*% t(M)) # (J x P)
    
    # Select only (block)-diagonal elements
    W <- W * t(.csem_model$measurement)
    
    # Build proxies/scores
    if(.disattenuate) {
      Z <- Z - U %*% D
    } 
    vcv_Z <- cor(Z)
    
    # Scale weights
    W <- t(scaleWeights(vcv_Z, t(W)))
    
    Gamma <- Z %*% W
    
    ## Step 2: Update A for given Gamma, U and D ---------------------------------
    
    # The matrix Theta in H&T (204, 2014, 2017) can get pretty large if N is large.
    # Its more efficient to update C und B in A seperately
    
    # Step 2a: Update B (the structural coefficients for a given W)
    vcv_gamma <- t(W) %*% vcv_Z %*% W
    vars_endo <- which(colSums(B) != 0)
    
    beta <- lapply(vars_endo, function(y) {
      x    <- which(B[, y, drop = FALSE] != 0)
      coef <- solve(vcv_gamma[x, x, drop = FALSE]) %*% vcv_gamma[x, y, drop = FALSE]
      
      # Since Var(dep_Var) = 1 we have R2 = Var(X coef) = t(coef) %*% X'X %*% coef
      # r2   <- t(coef) %*% .P[indep_var, indep_var, drop = FALSE] %*% coef
      # names(r2) <- x
    })
    B[B != 0] <- unlist(beta)
    
    # Step 2b: Update C (the loadings for a given W)
    vars_cf <- which(.csem_model$construct_type == "Common factor")
    if(length(vars_cf) > 0) {
      
      cov_gamma_indicators <- t(W) %*% vcv_Z
      Y <- which(colSums(C[vars_cf, ]) !=0)
      loadings <- lapply(Y, function(y) {
        x    <-  which(C[vars_cf, y] != 0)
        coef <- solve(vcv_gamma[x, x, drop = FALSE]) %*% cov_gamma_indicators[x, y, drop = FALSE]
      })
      
      # Transform
      tC <- t(C)
      tC[tC != 0] <- unlist(loadings)
      C <- t(tC)
    }
    
    # Build A
    A <- cbind(C, B)
    
    ### GSCAm steps
    if(.disattenuate) {
      ## Step 3: update U (unique parts) for given Gamma, A and D
      ## Based on Trendafilov et al.'s procedure (2013)
      # Gamma_orth is an N by N-P orthonormal basis matrix of the null space of Gamma
      Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (P+1):N, drop = FALSE]
      s          <- svd(t(Gamma_orth) %*% Z %*% D)
      U_tilde    <- s$u %*% t(s$v)
      U          <- Gamma_orth %*% U_tilde
      
      ## step 4: update D (unique loadings) and Z for fixed Gamma, A and U 
      D <- diag(diag(t(U) %*% Z))
      Z <- Z - U %*% D
    }
    
    ## Check convergence
    Psi <- cbind(Z, Gamma) # or Psi_star in GSCAm
    
    # Optimization criterion
    optim_crit <- Psi - Gamma %*% A  
    f   <- matrixcalc::matrix.trace(t(optim_crit) %*% optim_crit) # since trace(X'X) = SS(X)
    imp <- f0 - f # how much did the minimization criterion decrease
    
    if(imp < .tolerance) {
      # Set convergence status to TRUE as algorithm has converged
      conv_status <- TRUE
      break 
      
    } else if(iter_counter > .iter_max) {
      # Set convergence status to FALSE, as algorithm has not converged
      conv_status <- FALSE
      break
      
    } else {
      # Reset f0 and continue
      f0  <- f
    }
  }
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = "gsca", "Conv_status" = conv_status,
            "Iterations" = iter_counter)
  return(l)
  
} # END calculateWeightsGSCA

# calculateWeightsGSCA2 <- function(
#   .X                           = args_default()$.X,
#   .S                           = args_default()$.S,
#   .csem_model                  = args_default()$.csem_model,
#   .disattenuate                = args_default()$.disattenuate,
#   .conv_criterion              = args_default()$.conv_criterion,
#   .iter_max                    = args_default()$.iter_max,
#   .tolerance                   = args_default()$.tolerance
# ) {
#   ### Calculation (ALS algorithm) ==============================================
#   W0 <- Lambda0 <- .csem_model$measurement # Weight relation model
#   Lambda0[which(.csem_model$construct_type == "Composite"), ]  <- 0
#   B0 <- .csem_model$structural # Structural model
#   vars_endo    <- rownames(B0)[rowSums(B0) != 0]
#   vars_cf      <- names(which(.csem_model$construct_type == "Common factor"))
# 
#   J <- nrow(W0)
#   K <- ncol(W0)
#   JK <- J + K
#   
#   # Scale weights
#   W <- scaleWeights(.S = .S, .W = W0)
#   
#   W_iter       <- W
#   tolerance    <- .tolerance
#   iter_max     <- .iter_max
#   iter_counter <- 0
#   
#   repeat {
#     # Counter
#     iter_counter <- iter_counter + 1
#     
#     # Step 1a: Estimate B (the structural coefficients for a given W)
#     C <- W_iter %*% .S %*% t(W_iter)
#     B <- lapply(vars_endo, function(y) {
#       x <-  colnames(B0[y, B0[y, ] != 0, drop = FALSE])
#       coef <- solve(C[x, x, drop = FALSE]) %*% C[x, y, drop = FALSE]
#       
#       # Since Var(dep_Var) = 1 we have R2 = Var(X coef) = t(coef) %*% X'X %*% coef
#       # r2   <- t(coef) %*% .P[indep_var, indep_var, drop = FALSE] %*% coef
#       # names(r2) <- x
#     })
#     # Transform
#     tB <- t(B0)
#     tB[tB != 0] <- unlist(B)
#     
#     # Step 1b: Estimate Lambda (the loadings for a given W)
#     if(length(vars_cf) > 0) {
#       Lambda_tilde <- W_iter %*% .S
#       dep_vars <- names(which(colSums(Lambda0[vars_cf, ]) !=0))
#       Lambda <- lapply(dep_vars, function(y) {
#         x <-  which(Lambda0[vars_cf, y] != 0)
#         coef <- solve(C[x, x, drop = FALSE]) %*% Lambda_tilde[x, y, drop = FALSE]
#       })
#       # Transform
#       tLambda <- t(Lambda0)
#       tLambda[tLambda != 0] <- unlist(Lambda)
#     } else {
#       tLambda <- t(Lambda0)
#     }
# 
#     # Build matrices A and V used in GSCA
#     
#     # A <- rbind(tLambda, t(tB))
#     A <- cbind(t(tLambda), tB)
#     V <- rbind(diag(K), W_iter)
#     
#     # The following code is based on the ASGSCA package (licensed
#     # under GPL-3). Notation is adapted to be conform with the notation of the 
#     # cSEM package
#     
#     tr_w <- 0
#     W <- W_iter
#     for(j in 1:J){
#       t <- K + j
#       windex_j <- which(W[j, ] != 0)
#       m <- matrix(0, 1, JK)
#       m[t] <- 1
#       # a <- A[, j]
#       a <- A[j, ]
#       beta <- m - a
#       H1 <- diag(J)
#       H2 <- diag(JK)
#       H1[j,j] <- 0
#       H2[t,t] <- 0
#       Delta <- t(A) %*% H1 %*% W - H2 %*% V 
#       Sp <- .S[windex_j , windex_j]
#       if(length(windex_j) != 0) {        
#         
#         theta <- MASS::ginv(as.numeric(beta%*%t(beta))*.S[windex_j,windex_j]) %*%
#           t(beta %*% Delta %*% .S[,windex_j])
#         
#         # Update the weights based on the estimated parameters and standardize
#         W0[j, windex_j] <- theta
#         W <- scaleWeights(.S, W0)
#       }
#     }
#     
#     # Check for convergence
#     conv <- checkConvergence(W, W_iter, 
#                              .conv_criterion = .conv_criterion, 
#                              .tolerance = .tolerance)
#     
#     if(conv) {
#       # Set convergence status to TRUE as algorithm has converged
#       conv_status = TRUE
#       break # return iterative PLS-PM weights
#       
#     } else if(iter_counter == iter_max & iter_max == 1) {
#       # Set convergence status to NULL, NULL is used if no algorithm is used
#       conv_status = NULL
#       break # return one-step PLS-PM weights
#       
#     } else if(iter_counter == iter_max & iter_max > 1) {
#       # Set convergence status to FALSE, as algorithm has not converged
#       conv_status = FALSE
#       break
#       
#     } else {
#       W_iter <- W
#     }
#   }
# 
#   H <- A - cbind(matrix(0, nrow = J, ncol = K), diag(J)) # (P x T)
#   M <- cbind(diag(K), matrix(0, nrow = K, ncol = J) )
#   
#   Wnew <- M %*% t(H) %*% solve(H %*% t(H))
#   solve(H %*% t(H)) %*% H %*% t(M)
#   Zp <- .X[, 1:2]
#   wp <- Wnew[1:2, 1, drop = FALSE]
#   wp %*% (1/sqrt(t(wp) %*% (1/99*t(Zp) %*% Zp) %*% wp))
#   
#   var(Zp %*% wp)
#   cov(Zp)
#   Wnew2 <-  solve(H %*% t(H)) %*% H %*% t(.X %*% M)
#   all.equal(Wnew, t(Wnew2))
#   Wnew[1:5, ]
#   t(Wnew2)[1:5, ]
#   wp <- Wnew[, 1, drop = FALSE]
#   varproxy <- diag(t(Wnew) %*% .S %*% Wnew)
#   wp/sqrt(varproxy[1])
#   scaleWeights(t(Wnew), .S)
#   # GSCAm
#   if(.disattenuate) {
#     # W <- t(GSCA_results$Estimates$Weight_estimates) # Matrix of the weighted relation model
#     # B <- t(GSCA_results$Estimates$Path_estimates) # Matrix of the structural model
#     # C <- t(GSCA_results$Estimates$Loading_estimates) # Matrix of the measurement model if all indicators are reflective
#     # Gamma <- GSCA_results$Estimates$Construct_scores # Matrix of scores of latent variables
# 
#     Z <- .X # Z is the data matrix in GSCAm, data are already standardized
#     S <- .S # S is the empirical indicator covariance/correlation matrix
#     W <- t(W)
#     B <- tB
#     C <- t(W) %*% S * .csem_model$measurement
#     Gamma <- Z %*% W
# 
#     N <- nrow(Z) # number of observations per indicator
#     K <- nrow(W) # number of indicators
#     J <- ncol(W) # number of constructs
#     T <- K + J
# 
#     A <- cbind(C, B) # J rows, T columns
#     D <- diag(runif(K, min = 0, max = 1)) # matrix of unique loadings, random inital values
#     U <- matrix(0, N, K) # matrix of unique variables
# 
#     ## get indices of those entries of A which are not zero
# 
#     vecA <- matrix(A, ncol = 1)
#     aindex <- which(vecA!=0,arr.ind = T)[,1]
#     vecW <- matrix(W, ncol = 1)
#     windex <- which(vecW!=0,arr.ind = T)[,1]
# 
#     ## Preparation of the calculation of the actual weights, coefficients and loadings
# 
#     ## prepare data for GSCAm
#     # for GSCAm, normalized data are needed
#     vZ = sqrt(diag(t(Z) %*% Z))
#     factorZ = matrix(vZ, nrow=N, ncol=length(vZ), byrow=TRUE)
#     Z = Z/factorZ # normalized matrix Z
# 
#     ## prepare latent variables (constructs) for GSCAm
#     v = sqrt(diag(t(Gamma) %*% Gamma))
#     factor = matrix(v, nrow=N, ncol=length(v), byrow=TRUE)
#     Gamma = Gamma/factor # normalized matrix Gamma of latent variables
# 
#     # calculate initial values for the unique variables in U
#     # analogous to step 3 of the modified ALS-algorithm
#     Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, N-J, drop = FALSE]
#     s <- svd(t(Gamma_orth) %*% Z %*% D)
#     U_tilde = s$u %*% t(s$v)
#     U = Gamma_orth %*% U_tilde # this is the initial matrix U
# 
#     # iter = 0
#     iter_counter <- 1
#     f0 = 100000
#     imp = 100000
#     repeat {
#       # Counter
#       iter_counter <- iter_counter + 1
# 
#       # iter = iter + 1
# 
#       ## step 1: update Gamma (latent variable scores) and W (weights)
# 
#       H = A - cbind(matrix(0, nrow = J, ncol = K), diag(J))
#       M = cbind(diag(K), matrix(0, nrow = K, ncol = J))
#       W = M %*% t(H) %*% solve(H %*% t(H))
#       Gamma = (Z - U %*% D) %*% W
# 
#       # normalize every latent variable and adjust the corresponding weights
#       v = sqrt(diag(t(Gamma) %*% Gamma)) # length of v is J
#       factor = matrix(v, nrow=N, ncol=length(v), byrow=TRUE)
#       Gamma = Gamma/factor # normalized matrix of latent variables
#       W = W/matrix(v, nrow=K, ncol=length(v), byrow=TRUE)
# 
#       ## step 2: update A (loadings and path coefficients)
# 
#       Sm = cbind(U %*% D, matrix(0, nrow = N, ncol = J))
#       Psi = cbind(Z, Gamma)
#       PsiStar = Psi - Sm
#       vecPsi = matrix(Psi, ncol = 1)
#       vecPsiStar = matrix(PsiStar, ncol = 1)
# 
#       Phi = kronecker(diag(T), Gamma)
#       Phi = Phi[,aindex]
#       a_hat = solve((t(Phi)%*%Phi))%*%t(Phi)%*%vecPsiStar
#       vecA[aindex,1] = a_hat
#       A = matrix(vecA, ncol = T, nrow = J, byrow = FALSE)
#       vecA <- matrix(A, ncol = 1)
# 
#       ## step 3: update U (unique parts) for fixed Gamma, A and D
#       ## based on Trendafilov et al.'s procedure (2013)
# 
#       # Gamma_orth is an N by N-J orthonormal basis matrix of the null space of Gamma
#       Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, N-J, drop = FALSE]
#       s <- svd(t(Gamma_orth) %*% Z %*% D)
#       U_tilde = s$u %*% t(s$v)
#       U = Gamma_orth %*% U_tilde
# 
#       ## step 4: update D (unique loadings) for fixed Gamma, A and U
# 
#       D = diag(diag(t(U) %*% Z))
# 
#       ## updated matrices Psi and Sm
#       Psi = cbind(Z, Gamma)
#       Sm = cbind(U %*% D, matrix(0, nrow = N, ncol = J))
#       dif = Psi - Gamma %*% A - Sm # This is the matrix of the optimization criterion
#       f = matrixcalc::matrix.trace(t(dif) %*% dif)
#       imp = f0 - f # decrease of the minimization criterion
#       f0 = f
#       PsiStar = Psi - Sm
#       vecPsi = matrix(Psi, ncol = 1)
#       vecPsiStar = matrix(PsiStar, ncol = 1)
# 
# 
#       # Check for convergence
#       conv <- checkConvergence(W, W_iter,
#                                .conv_criterion = .conv_criterion,
#                                .tolerance = .tolerance)
# 
#       if(conv) {
#         # Set convergence status to TRUE as algorithm has converged
#         conv_status = TRUE
#         break # return iterative PLS-PM weights
# 
#       } else if(iter_counter == iter_max & iter_max == 1) {
#         # Set convergence status to NULL, NULL is used if no algorithm is used
#         conv_status = NULL
#         break # return one-step PLS-PM weights
# 
#       } else if(iter_counter == iter_max & iter_max > 1) {
#         # Set convergence status to FALSE, as algorithm has not converged
#         conv_status = FALSE
#         break
# 
#       } else {
#         W_iter <- W
#       }
#     }
#   }
#   
#   # Return
#   l <- list("W" = W, "E" = NULL, "Modes" = "gsca", "Conv_status" = conv_status,
#             "Iterations" = iter_counter)
#   return(l)
#   
#   ### For maintenance: ---------------------------------------------------------
#   # Hwang and Takane (2004, 2010, 2014, 2017) use a completely different notation
#   #   compared to the standard LISREL/Bollen-type notation. This implementation
#   #   uses the LISREL notation. This table translates some of the notation
#   ## N              := Number of observations
#   ## K              := Total number of indicators (J in H&T)
#   ## J              := Total number of constructs (P in H&T)
#   ## TT             := K + J (T in H&T)
#   ## X (N x K)      := Matrix of indicator values (=data) (Z in H&T)
#   ## W (J x K)      := (Block-)diagonal matrix of weight relations (W' in H&T)
#   ## Lambda (J x K) := (Block-)diagonal matrix of loadings (C in H&T)
#   ## B (J x J)      := Matrix of path coefficients (B' in H&T). B is not necessarily
#   ##                   lower-triangular (i.e may contain feedback loops)
#   ## H (N x J)      := Matrix of proxies/scores  
#   ## V (K x TT)     := Matrix of indicator and construct relations
#   ## Psi (N x TT)   := Matrix of all indicator values and construct scores
# 
# } # END calculateWeightsGSCA

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
