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
    
    # Error if other than "ModeA", "ModeB", "unit", a number, or a vector
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
        paste0("`", setdiff(names(.PLS_modes), names(csem_model$construct_type)), 
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
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .tolerance                   = args_default()$.tolerance
) {
  ### Calculation (ALS algorithm) ==============================================
  W0 <- Lambda0 <- .csem_model$measurement # Weight relation model
  Lambda0[which(.csem_model$construct_type == "Composite"), ]  <- 0
  B0 <- .csem_model$structural # Structural model
  vars_endo    <- rownames(B0)[rowSums(B0) != 0]
  vars_cf      <- names(which(.csem_model$construct_type == "Common factor"))

  J <- nrow(W0)
  K <- ncol(W0)
  JK <- J + K
  
  # Scale weights
  W <- scaleWeights(.S = .S, .W = W0)
  
  W_iter       <- W
  tolerance    <- .tolerance
  iter_max     <- .iter_max
  iter_counter <- 0
  
  repeat {
    # Counter
    iter_counter <- iter_counter + 1
    
    # Step 1a: Estimate B (the structural coefficients for a given W)
    C <- W_iter %*% .S %*% t(W_iter)
    B <- lapply(vars_endo, function(y) {
      x <-  colnames(B0[y, B0[y, ] != 0, drop = FALSE])
      coef <- solve(C[x, x, drop = FALSE]) %*% C[x, y, drop = FALSE]
      
      # Since Var(dep_Var) = 1 we have R2 = Var(X coef) = t(coef) %*% X'X %*% coef
      # r2   <- t(coef) %*% .P[indep_var, indep_var, drop = FALSE] %*% coef
      # names(r2) <- x
    })
    # Transform
    tB <- t(B0)
    tB[tB != 0] <- unlist(B)
    
    # Step 1b: Estimate Lambda (the loadings for a given W)
    if(length(vars_cf) > 0) {
      Lambda_tilde <- W_iter %*% .S
      dep_vars <- names(which(colSums(Lambda0[vars_cf, ]) !=0))
      Lambda <- lapply(dep_vars, function(y) {
        x <-  which(Lambda0[vars_cf, y] != 0)
        coef <- solve(C[x, x, drop = FALSE]) %*% Lambda_tilde[x, y, drop = FALSE]
      })
      # Transform
      tLambda <- t(Lambda0)
      tLambda[tLambda != 0] <- unlist(Lambda)
    } else {
      tLambda <- t(Lambda0)
    }

    # Build matrices A and V used in GSCA
    
    A <- rbind(tLambda, t(tB))
    V <- rbind(diag(K), W_iter)
    
    # The following code is based on the ASGSCA package (licensed
    # under GPL-3). Notation is adapted to be conform with the notation of the 
    # cSEM package
    
    tr_w <- 0
    for(j in 1:J){
      t <- K + j
      windex_j <- which(W_iter[j, ] != 0)
      m <- matrix(0, 1, JK)
      m[t] <- 1
      a <- A[, j]
      beta <- m - a
      H1 <- diag(J)
      H2 <- diag(JK)
      H1[j,j] <- 0
      H2[t,t] <- 0
      Delta <- A %*% H1 %*% W_iter - H2 %*% V 
      Sp <- .S[windex_j , windex_j]
      if(length(windex_j) != 0) {        
        
        theta <- MASS::ginv(as.numeric(beta%*%t(beta))*.S[windex_j,windex_j]) %*%
          t(beta %*% Delta %*% .S[,windex_j])
        
        # Update the weights based on the estimated parameters and standardize
        W0[j, windex_j] <- theta
        W <- scaleWeights(.S, W0)
      }
    }
    
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
  l <- list("W" = W, "E" = NULL, "Modes" = "gsca", "Conv_status" = conv_status,
            "Iterations" = iter_counter)
  return(l)
  
  ### For maintenance: ---------------------------------------------------------
  # Hwang and Takane (2004, 2010, 2014, 2017) use a completely different notation
  #   compared to the standard LISREL/Bollen-type notation. This implementation
  #   uses the LISREL notation. This table translates some of the notation
  ## N              := Number of observations
  ## K              := Total number of indicators (J in H&T)
  ## J              := Total number of constructs (P in H&T)
  ## TT             := K + J (T in H&T)
  ## X (N x K)      := Matrix of indicator values (=data) (Z in H&T)
  ## W (J x K)      := (Block-)diagonal matrix of weight relations (W' in H&T)
  ## Lambda (J x K) := (Block-)diagonal matrix of loadings (C in H&T)
  ## B (J x J)      := Matrix of path coefficients (B' in H&T). B is not necessarily
  ##                   lower-triangular (i.e may contain feedback loops)
  ## H (N x J)      := Matrix of proxies/scores  
  ## V (K x TT)     := Matrix of indicator and construct relations
  ## Psi (N x TT)   := Matrix of all indicator values and construct scores

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
