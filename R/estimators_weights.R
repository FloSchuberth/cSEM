#' Calculate composite weights using PLS
#'
#' Calculates composite weights using the PLS algorithm. (TODO)
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
#'   \item{`$Conv_status`}{The convergence status. One of `TRUE`or `FASLE`. If 
#'     one-step weights are used via `.iter_max = 1` or no iterative algorithm was used,
#'     the convergence status is set to `NULL`.}
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
       length(setdiff(names(.PLS_modes), names(csem_model$construct_type))) > 0) {
      stop(paste0("`", setdiff(names(.PLS_modes), names(csem_model$construct_type)), 
                  "`", collapse = ", ")," in `.PLS_modes` is an unknown construct name.", call. = FALSE)
    }
    # If only "ModeA" or "ModeB" is provided without set all of the modes to that mode.
    if(length(names(.PLS_modes)) == 0) {
      if(length(.PLS_modes) == 1) {
        modes <- rep(.PLS_modes, length(csem_model$construct_type))
        names(modes) <- names(csem_model$construct_type)
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
      break # return iterative PLS weights
      
    } else if(iter_counter == iter_max & iter_max == 1) {
      # Set convergence status to NULL, NULL is used if no algorithm is used
      conv_status = NULL
      break # return one-step PLS weights
      
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
  .data         = NULL,
  .model        = NULL
) {

  stop("Not yet implemented")

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
