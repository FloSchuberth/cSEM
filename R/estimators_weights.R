#' Calculate composite weights using PLS
#'
#' Calculates composite weights using the PLS algorithm.
#'
#' More description here
#'
#' @usage calculateWeights(
#'   .data                        = NULL,
#'   .model                       = NULL,
#'   .PLS_modes                   = NULL,
#'   .tolerance                   = NULL,
#'   .iter_max                    = NULL,
#'   .PLS_ignore_structural_model = NULL
#'   .PLS_weight_scheme_inner     = NULL,
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{A (J x J) matrix of inner weights.}
#'   \item{`$Modes`}{A named vector of Modes used for the outer estimation.}
#'   \item{`$Conv_status`}{The convergence status. One of `TRUE`or `FASLE`. If 
#'     one-step weights are used via `.iter_max = 1` the convergence status
#'     is set to `NULL`.}
#'   \item{`$Iterations`}{The number of iterations used.}
#' }
#' @export
#'

calculateWeightsPLS <- function(
  .data                        = NULL,
  .model                       = NULL,
  .PLS_modes                   = NULL,
  .tolerance                   = NULL,
  .iter_max                    = NULL,
  .PLS_ignore_structural_model = NULL,
  .PLS_weight_scheme_inner     = NULL
) {


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

  ### Preparation ==============================================================
  ## Get/set the modes for the outer estimation

  if(is.null(.PLS_modes)) {
    
    modes <- ifelse(csem_model$construct_type == "Common factor", "ModeA", "ModeB")

  } else if(all(.PLS_modes %in% c("ModeA", "ModeB"))) {
    
    if(setequal(names(.PLS_modes), names(csem_model$construct_type))) {
      modes <- .PLS_modes
      modes <- modes[names(csem_model$construct_type)]
    } else if(length(.PLS_modes) == 1) {
      modes <- rep(.PLS_modes, length(csem_model$construct_type))
      names(modes) <- names(csem_model$construct_type)
    } else if(length(setdiff(names(.PLS_modes), names(csem_model$construct_type))) > 0) {
      stop(paste0("`", setdiff(names(.PLS_modes), names(csem_model$construct_type)), 
                  "`", collapse = ", ")," in `.PLS_modes` is an unknown construct name.", call. = FALSE)
    } else  {
      stop("Mode ", paste0("`", setdiff(names(csem_model$construct_type), names(.PLS_modes)), 
                  "`", collapse = ", ")," in `.PLS_modes` is missing.", call. = FALSE)
    }
  } else {
    stop(paste0("`", setdiff(.PLS_modes, c("ModeA", "ModeB")), "`", collapse = ", "),
         " in `.PLS_modes` is an unknown mode.", call. = FALSE)
  }

  ### Calculation/Iteration ====================================================
  W <- csem_model$measurement
  # Scale weights
  W <- scaleWeights(.S = S, .W = W)

  W_iter       <- W
  tolerance    <- .tolerance
  iter_max     <- .iter_max
  iter_counter <- 0

  repeat {
    # Counter
    iter_counter <- iter_counter + 1

    # Inner estimation
    E <- calculateInnerWeightsPLS(
      .S                           = S,
      .W                           = W_iter,
      .csem_model                  = csem_model,
      .PLS_ignore_structural_model = .PLS_ignore_structural_model,
      .PLS_weight_scheme_inner     = .PLS_weight_scheme_inner
    )
    # Outer estimation

    W <- calculateOuterWeightsPLS(
      .S        = S,
      .W        = W_iter,
      .E        = E,
      .modes    = modes
    )

    # Scale weights
    W <- scaleWeights(S, W)

    # Check for convergence
    if(max(abs(W_iter - W)) < .tolerance) {
      # Set convergence status to TRUE as algorithm has converged
      conv_status = TRUE
      break # return iterative PLS weights
    } else if(iter_counter == iter_max & iter_max == 1) {
      # Set convergence status to NULL, NULL is used if no algorithm is used
      conv_status = NULL
      break # return one-step PLS weights
      
    } else if(iter_counter == iter_max & iter_max > 1) {
      # Set convergence status to FALSE, as algorithm has not converged
      conv_status = FALSE
      warning("The PLS algorithm did not converge after ", iter_max, " steps. ",
              "Last weights are returned.", 
              call. = FALSE)
      
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

#' Calculate the inner weights for PLS
#'
#' Calculates inner weights
#'
#' More description here.
#'
#' @usage calculateInnerWeightsPLS(
#'   .S                           = NULL,
#'   .W                           = NULL,
#'   .csem_model                  = NULL,
#'   .PLS_ignore_structural_model = NULL,
#'   .PLS_weight_scheme_inner     = NULL
#' )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x J) matrix `E` of inner weights.
#'
calculateInnerWeightsPLS <- function(
  .S                           = NULL,
  .W                           = NULL,
  .csem_model                  = NULL,
  .PLS_ignore_structural_model = NULL,
  .PLS_weight_scheme_inner     = NULL
) {

  # Composite correlation matrix (C = V(H))
  C <- .W %*% .S %*% t(.W)

  # Due to floting point errors may not be symmetric anymore. In order
  # prevent that replace the lower triangular elements by the upper
  # triangular elements

  C[lower.tri(C)] <- t(C)[lower.tri(C)]

  # Path model relationship
  tmp <- rownames(.csem_model$measurement)
  E   <- .csem_model$structural[tmp, tmp]
  D   <- E + t(E)

  ## (Inner) weightning scheme:
  if(.PLS_weight_scheme_inner == "path" & .PLS_ignore_structural_model) {
    .PLS_ignore_structural_model <- FALSE
    warning("Structural model required for the path weighting scheme.\n",
            ".PLS_ignore_structural_model = TRUE was changed to FALSE.", 
            call. = FALSE)
  }

  if(.PLS_ignore_structural_model) {
    switch (.PLS_weight_scheme_inner,
            "centroid"  = {E <- sign(C) - diag(1, nrow = nrow(.W))},
            "factorial" = {E <- C - diag(1, nrow = nrow(.W))}
    )
  } else {
    switch (.PLS_weight_scheme_inner,
            "centroid"  = {E <- sign(C) * D},
            "factorial" = {E <- C * D},
            "path"      = {
              ## All construct that have an arrow pointing to construct j
              ## are called predecessors of j; all arrows going from j to other
              ## constructs are called followers of j

              ## Path weighting scheme:
              ## Take construct j:
              #  - For all predecessors of j: set the inner weight of
              #    predescessor i to the correlation of i with j.
              #  - For all followers of j: set the inner weight of follower i
              #    to the coefficient of a multiple regression of j on
              #    all followers i with i = 1,...,I.

              E_temp <- E
              ## Assign predescessor relation
              E[t(E_temp) == 1] <- C[t(E_temp) == 1]

              ## Assign follower relation
              for(j in tmp) {

                followers <- E_temp[j, ] == 1

                if(sum(followers) > 0) {

                  E[j, followers] <-  solve(C[followers, followers, drop = FALSE]) %*%
                    C[j, followers]
                }
              }
            }
    )
  }

  return(E)
} # END calculateInnerWeights

#' Calculate the outer weights for PLS
#'
#' Calculates outer weights
#'
#' More description here..
#'
#' @usage calculateOuterWeightsPLS(
#'    .S      = NULL,
#'    .W      = NULL,
#'    .E      = NULL,
#'    .modes  = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A (J x K) matrix of outer weights.
#'

calculateOuterWeightsPLS <- function(
  .S      = NULL,
  .W      = NULL,
  .E      = NULL,
  .modes  = NULL
) {
  # Covariance/Correlation matrix between each proxy and all indicators (not
  # only the related ones). Note: Cov(H, X) = WS, since H = XW'.
  W <- .W
  proxy_indicator_cor <- .E %*% W %*% .S

  for(i in 1:nrow(W)) {
    block      <- rownames(W[i, , drop = FALSE])
    indicators <- W[block, ] != 0

    if(.modes[block] == "ModeA") {
      ## Mode A - Regression of each indicator on its corresponding proxy
      # Select only
      W[block, indicators] <- proxy_indicator_cor[block, indicators]

    } else if(.modes[block] == "ModeB") {
      ## Mode B - Regression of each proxy on all its indicator
      # W_j = S_jj^{-1} %*% Cov(eta_j, X_j)
      W[block, indicators] <- solve(.S[indicators, indicators]) %*% proxy_indicator_cor[block, indicators]

    } # END ModeB
  }
  return(W)
} # END calculateOuterWeights

#' Calculate weights using one of Kettenrings's approaches
#'
#' Calculates weights according to on of the the five criteria
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*", and "*GENVAR*"
#' suggested by \insertCite{Kettenring1971;textual}{cSEM}.
#'
#' Some more description...
#'
#' @usage calculateWeightsKettenring(.data, .model, .approach)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsKettenring <- function(
  .data           = NULL,
  .model          = NULL,
  .approach       = NULL
) {
  
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
  
  ### Preparation ==============================================================
  ## Check if all constructs are modeled as composites
  if("Common factor" %in% csem_model$construct_type){
    stop("Currently Kettenring's approaches are only allowed for pure composite models.", 
         call. = FALSE)
  }

  ## Get relevant objects 
  W <- csem_model$measurement
  
  construct_names <- rownames(W)
  indicator_names <- colnames(W)
  indicator_names_ls <- lapply(construct_names, function(x) {
    indicator_names[W[x, ] == 1]
  })
  
  ## Calculate S_{ii}^{-1/2}}
  sqrt_S_jj_list <- lapply(indicator_names_ls, function(x) {
    solve(expm::sqrtm(S[x, x]))
  })
  
  ## Put in a diagonal matrix
  H <- as.matrix(Matrix::bdiag(sqrt_S_jj_list))
  dimnames(H) <- list(indicator_names, indicator_names)
  
  if(.approach %in% c("MAXVAR", "MINVAR")) {
    # Note: The MAXVAR implementation is based on a MATLAB code provided 
    #       by Theo K. Dijkstra
    Rs <- H %*% S %*% H
    
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
      
      # Correct if negative
      if(sum(w) < 0) w <- -w
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
      R       = S,
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
      
      # Correct if negative
      if(sum(w) < 0) w <- -w
      return(w)
    })
  }
  
  ## Format, name and scale
  W <- t(as.matrix(Matrix::bdiag(W)))
  dimnames(W) <- list(construct_names, indicator_names)
  W <- scaleWeights(S, W)

  ## Prepare for output
  modes <- rep(.approach, length(construct_names))
  names(modes) <- construct_names
  
  ## Prepare output and return
  l <- list(
    "W"           = W, 
    "E"           = NULL, 
    "Modes"       = modes,
    "Conv_status" = if(.approach %in% c("MINVAR", "MAXVAR")) NULL else res$convergence == 0,  
    "Iterations"  = if(.approach %in% c("MINVAR", "MAXVAR")) NULL else res$counts[1]
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
  .data         = NULL,
  .model        = NULL
) {

  stop("Not yet implemented")

} # END calculateWeightsGSCA

#' Calculate weights using fixed weights
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

#' Calculate unit weights for all blocks, i.e., each indicator of a block is equally weighted 
#' This approach might be useful for Kroon correction.
#'
#' Calculates weights...
#'
#' Some more description...
#'
#' @usage calculateWeightsUnit(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'
calculateWeightsUnit = function(
  .data         = NULL,
  .model        = NULL
){
  
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
  
  W=.model$measurement
  W=scaleWeights(S,W)
  
  modes=rep('unit',nrow(W))
  names(modes)=rownames(W)
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = modes, "Conv_status" = TRUE,
            "Iterations" = 0)
  return(l)  
}
