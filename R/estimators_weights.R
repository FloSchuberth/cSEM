#' Calculate composite weights using PLS-PM
#'
#' Calculate composite weights using the partial least squares path modeling 
#' (PLS-PM) algorithm \insertCite{Wold1975}{cSEM}.
#'
#' @usage calculateWeightsPLS(
#'   .data                        = args_default()$.data,
#'   .S                           = args_default()$.S,
#'   .csem_model                  = args_default()$.csem_model,
#'   .conv_criterion              = args_default()$.conv_criterion,
#'   .iter_max                    = args_default()$.iter_max,
#'   .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
#'   .PLS_modes                   = args_default()$.PLS_modes,
#'   .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
#'   .starting_values             = args_default()$.starting_values,
#'   .tolerance                   = args_default()$.tolerance
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A named list. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{A (J x J) matrix of inner weights.}
#'   \item{`$Modes`}{A named vector of modes used for the outer estimation.}
#'   \item{`$Conv_status`}{The convergence status. `TRUE` if the algorithm has converged 
#'     and `FALSE` otherwise. If one-step weights are used via `.iter_max = 1` 
#'     or a non-iterative procedure was used, the convergence status is set to `NULL`.}
#'   \item{`$Iterations`}{The number of iterations required.}
#' }
#' 
#' @references
#'   \insertAllCited{}

calculateWeightsPLS <- function(
  .data                        = args_default()$.data,
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes                   = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
  .starting_values             = args_default()$.starting_values,
  .tolerance                   = args_default()$.tolerance
) {

  ### Preparation ==============================================================
  ## Get/set the default modes for the outer estimation
  modes <- as.list(ifelse(.csem_model$construct_type == "Common factor", "modeA", "modeB"))
  
  if(!is.null(.PLS_modes)) {
    
    # Error if other than "modeA", "modeB", "unit", "modeBNNLS", "PCA", a number, or a vector
    # of numbers of the same length as there are indicators for block j
    modes_check <- sapply(.PLS_modes, 
                          function(x) all(x %in% c("modeA", "modeB", "unit", 
                                                   "modeBNNLS", "PCA") | is.numeric(x)))
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
    
    # If only "modeA", "modeB", "modeBNNLS", "PCA", or "unit" is provided set all 
    # of the modes to that mode.
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
  
  # if starting values are provided
  if(!is.null(.starting_values)){
    W = setStartingValues(.W = W, .starting_values = .starting_values)
    }
  
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
      .data     = .data,
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

#' Calculate composite weights using GCCA
#'
#' Calculates composite weights according to one of the the five criteria
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*", and "*GENVAR*"
#' suggested by \insertCite{Kettenring1971;textual}{cSEM}.
#'
#' @usage calculateWeightsKettenring(
#'   .S              = args_default()$.S, 
#'   .csem_model     = args_default()$.csem_model,   
#'   .approach_gcca  = args_default()$.approach_gcca
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A named list. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{`NULL`}
#'   \item{`$Modes`}{The GCCA mode used for the estimation.}
#'   \item{`$Conv_status`}{The convergence status. `TRUE` if the algorithm has converged 
#'     and `FALSE` otherwise. For `.approach_gcca = "MINVAR"` or `.approach_gcca = "MAXVAR"`
#'     the convergence status is `NULL` since both are closed-form estimators.}
#'   \item{`$Iterations`}{The number of iterations required. 0 for 
#'     `.approach_gcca = "MINVAR"` or `.approach_gcca = "MAXVAR"`}
#' }
#' 
#' @references
#'   \insertAllCited{}

calculateWeightsKettenring <- function(
  .S              = args_default()$.S,
  .csem_model     = args_default()$.csem_model,
  .approach_gcca  = args_default()$.approach_gcca
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
  
  if(.approach_gcca %in% c("MAXVAR", "MINVAR")) {
    # Note: The MAXVAR implementation is based on a MATLAB code provided 
    #       by Theo K. Dijkstra
    Rs <- H %*% .S %*% H
    
    ## Calculate eigenvalues and eigenvectors of Rs and ...
    u <- switch (.approach_gcca,
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
  } else if(.approach_gcca %in% c("SUMCORR", "GENVAR", "SSQCORR")) {
    
    # Maximize the sum of the composite correlations sum(w_i%*%S_ii%*%t(w_i)), 
    #   see Asendorf (2015, Appendix B.3.2 Empirical, p. 295)
    
    ## Define function to be minimized:
    # R       := Empirical indicator correlation matrix (S)
    # RDsqrt  := Blockdiagonal matrix with S_{ii}^{-1/2} on its diagonal (H)
    # nameInd := List with indicator names per block (indicator_names_ls)
    # name    := Vector of composites names (construct_names)
    
    fn <- switch (.approach_gcca,
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
    utils::capture.output(res <- alabama::auglag(
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
  modes <- rep(.approach_gcca, length(construct_names))
  names(modes) <- construct_names
  
  ## Prepare output and return
  l <- list(
    "W"           = W, 
    "E"           = NULL, 
    "Modes"       = modes,
    "Conv_status" = if(.approach_gcca %in% c("MINVAR", "MAXVAR")) NULL else res$convergence == 0,  
    "Iterations"  = if(.approach_gcca %in% c("MINVAR", "MAXVAR")) 0 else res$counts[1]
    )
  
  return(l)
  
} # END calculateWeightsKettenring


#' Calculate composite weights using GSCA
#'
#' Calculate composite weights using generalized structure component analysis (GSCA). 
#' The first version of this approach was presented in \insertCite{Hwang2004;textual}{cSEM}. 
#' Since then, several advancements have been proposed. The latest version 
#' of GSCA can been found in \insertCite{Hwang2014;textual}{cSEM}. This is the version 
#' \pkg{cSEM}s implementation is based on.
#'
#' @usage calculateWeightsGSCA(
#'   .X                           = args_default()$.X,
#'   .S                           = args_default()$.S,
#'   .csem_model                  = args_default()$.csem_model,
#'   .conv_criterion              = args_default()$.conv_criterion,
#'   .GSCA_modes                  = args_default()$.GSCA_modes,
#'   .iter_max                    = args_default()$.iter_max,
#'   .starting_values             = args_default()$.starting_values,
#'   .tolerance                   = args_default()$.tolerance
#'    )
#'    
#' @inheritParams csem_arguments
#' @importFrom MASS ginv
#' @return A named list. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{`NULL`}
#'   \item{`$Modes`}{A named vector of Modes used for the outer estimation, for GSCA
#'     the mode is automatically set to "gsca".}
#'   \item{`$Conv_status`}{The convergence status. `TRUE` if the algorithm has converged 
#'     and `FALSE` otherwise.}
#'   \item{`$Iterations`}{The number of iterations required.}
#' }
#' 
#' @references
#'   \insertAllCited{}

calculateWeightsGSCA <- function(
  .X                           = args_default()$.X,
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .GSCA_modes                  = args_default()$.GSCA_modes,
  .iter_max                    = args_default()$.iter_max,
  .starting_values             = args_default()$.starting_values,
  .tolerance                   = args_default()$.tolerance
) {
  
  ### Calculation (ALS algorithm) ==============================================
  W0 <- Lambda0 <- .csem_model$measurement # Weight relation model
  B0 <- .csem_model$structural # Structural model
  vars_endo    <- rownames(B0)[rowSums(B0) != 0]
  # vars_cf      <- names(which(.csem_model$construct_type == "Common factor"))
  
  J <- nrow(W0)
  K <- ncol(W0)
  JK <- J + K


  # Composite Type ---------------------------------------------------------
  modes <- ifelse(.csem_model$construct_type == "Composite", yes = "CCMP", no = .csem_model$construct_type)
    if (!is.null(.GSCA_modes)) {
      stopifnot(
        "Invalid .GSCA_modes selected. Only NCMP or CCMP are supported" = all(sapply(
          .GSCA_modes,
          function(x) x %in% c('NCMP', 'CCMP')
        ))
      )

      if (length(names(.GSCA_modes)) != 0) {
        stopifnot(
          "Invalid construct name listed in .GSCA_modes" = length(setdiff(
            names(.GSCA_modes),
            names(.csem_model$construct_type)
          )) ==
            0
        )
        modes[names(.GSCA_modes)] <- .GSCA_modes
      } else if (length(names(.GSCA_modes)) == 0) {
        stopifnot("When passing a global setting to .GSCA_modes, only one choice may be passed" = length(.GSCA_modes) == 1)

        # If only an unnamed element is given (NCMP or CCMP), then set all composites to .GSCA_modes
        if (length(.GSCA_modes) == 1) {
          modes <- ifelse(modes == "CCMP", yes = .GSCA_modes, no = modes)
        }
      }
    }
  Lambda0[which(modes == "CCMP"), ]  <- 0
  vars_NCMP <- names(modes)[modes == "NCMP"]
  vars_cf_ncmp <- names(modes)[modes %in% c("Common factor", "NCMP")]

  
  
  # Starting Values --------------------------------------------------------
  if(!is.null(.starting_values)){;
    W0 = setStartingValues(.W = W0, .starting_values = .starting_values)
  }
  
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
      x <-  (colSums(B0[y, , drop = FALSE])  !=  0) |> 
        which() |> 
        names()
      # The old approach could create problems when there was only one exogenous and one endogenous variable
      # x <-  colnames(B0[y, B0[y, ] != 0, drop = FALSE])
      coef <- MASS::ginv(C[x, x, drop = FALSE]) %*% C[x, y, drop = FALSE]
      
      # Since Var(dep_Var) = 1 we have R2 = Var(X coef) = t(coef) %*% X'X %*% coef
      # r2   <- t(coef) %*% .P[indep_var, indep_var, drop = FALSE] %*% coef
      # names(r2) <- x
    })
    # Transform
    tB <- t(B0)
    tB[tB != 0] <- unlist(B)

    # Step 1b: Estimate Lambda (the loadings for a given W)
    if(any(modes %in% "NCMP")) {
      Lambda_tilde <- W_iter %*% .S
      dep_vars <- (colSums(Lambda0[vars_cf_ncmp, , drop = FALSE]) != 0) |> 
        which() |> 
        names()
      Lambda <- lapply(dep_vars, function(y) {
        x <- (rowSums(Lambda0[vars_cf_ncmp, y, drop = FALSE]) != 0) |> 
          which() |> 
          names()

        coef <- MASS::ginv(C[x, x, drop = FALSE]) %*% Lambda_tilde[x, y, drop = FALSE]
      })
      # Transform
      tLambda <- t(Lambda0)
      tLambda[tLambda != 0] <- unlist(Lambda)
    } else {
      tLambda <- t(Lambda0)
    } 
    
    # Build matrices A and V used in GSCA
    
    # A <- rbind(tLambda, t(tB))
    A <- cbind(t(tLambda), tB)
    V <- rbind(diag(K), W_iter)
    
    # The following code is based on the ASGSCA package (licensed
    # under GPL-3). Notation is adapted to be conform with the notation of the
    # cSEM package
    
    tr_w <- 0
    W <- W_iter
    for(j in 1:J){
      t <- K + j
      windex_j <- (colSums(W[j, , drop = FALSE]) != 0) |>
        which() |>
        names()
      m <- matrix(0, 1, JK)
      m[t] <- 1
      # a <- A[, j]
      a <- A[j, ]
      beta <- m - a
      H1 <- diag(J)
      H2 <- diag(JK)
      H1[j,j] <- 0
      H2[t,t] <- 0
      Delta <- t(A) %*% H1 %*% W - H2 %*% V
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

  if (all(modes %in% "CCMP")) {
    Lambda <- W %*% .S * .csem_model$measurement
  } else if (any(modes %in% "CCMP")) {
    CCMP_Lambda <- W %*% .S * .csem_model$measurement
    Lambda <- t(tLambda)
    Lambda[names(modes)[modes == "CCMP"], ] <- CCMP_Lambda[
      names(modes)[modes == "CCMP"],
    ]
  } else {
    Lambda <- t(tLambda)
  }
  
  # Return
  l <- list(
    "W" = W,
    "C" = Lambda, 
    "B" = t(tB), 
    "E" = NULL,
    "Construct_scores" = .X %*% t(W),
    "Unique_loading_estimates" = NULL,
    "Unique_scores"  = NULL,
    "Modes" = "gsca",
    "Conv_status" = conv_status,
    "Iterations" = iter_counter
  )
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

#' Calculate weights using GSCAm
#'
#' Calculate composite weights using generalized structured component analysis
#' with uniqueness terms (GSCAm) proposed by \insertCite{Hwang2017;textual}{cSEM}.
#' 
#' If there are only constructs modeled as common factors
#' calling [csem()] with `.appraoch_weights = "GSCA"` will automatically call 
#' [calculateWeightsGSCAm()] unless `.disattenuate = FALSE`.
#' 
#' **Note**: Here, we assume that there is only one unique loading per indicator. 
#' 
#' @usage calculateWeightsGSCAm(
#'   .X                           = args_default()$.X,
#'   .csem_model                  = args_default()$.csem_model,
#'   .conv_criterion              = args_default()$.conv_criterion,
#'   .iter_max                    = args_default()$.iter_max,
#'   .starting_values             = args_default()$.starting_values,
#'   .tolerance                   = args_default()$.tolerance
#'    )
#'    
#' @inheritParams csem_arguments
#' @importFrom MASS ginv
#' @return A list with the elements
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$C`}{The (J x K) matrix of estimated loadings.}
#'   \item{`$B`}{The (J x J) matrix of estimated path coefficients.}
#'   \item{`$E`}{`NULL`}
#'   \item{`$Modes`}{A named vector of Modes used for the outer estimation, for GSCA
#'     the mode is automatically set to 'gsca'.}
#'   \item{`$Conv_status`}{The convergence status. `TRUE` if the algorithm has converged 
#'     and `FALSE` otherwise.}
#'   \item{`$Iterations`}{The number of iterations required.}
#' }
#' 
#' @references
#'   \insertAllCited{}
#'
calculateWeightsGSCAm <- function(
  .X                           = args_default()$.X,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .starting_values             = args_default()$.starting_values,
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
  #   is not the case if we include nonlinear model. This has to be checked
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
  
  # Currently only pure common factor models are supported
  # if(any(.csem_model$construct_type == "Composite")) {
  #   stop2("The following error occured in the `calculateWeightsGSCAm()` function:\n",
  #         "GSCAm only applicable to pure common factor models.")
  # }
  
  Z  <- .X # Z is the data matrix in GSCA, data are already standardized
  W  <- t(.csem_model$measurement) # Matrix of the weighted relation model
  C  <- .csem_model$measurement # Matrix of the measurement model (non-zero only for common factors)
  # C[which(.csem_model$construct_type == "Composite"), ]  <- 0
  B  <- t(.csem_model$structural) # Matrix of the structural model
  indicator_type <- .csem_model$construct_type
  N  <- nrow(Z) # number of observations
  J  <- nrow(W) # number of indicators
  P  <- ncol(W) # number of constructs
  TT <- J + P
  
  Ij <- diag(J)
  Ip <- diag(P)
  A  <- cbind(C, B) # (P x TT)
  # V  <- cbind(Ij, W) # (J x TT)

  # Normalize Z
  normalization_factor <- sqrt(N - 1)
  Z <- Z / normalization_factor
  
  # if starting values are provided
  if(!is.null(.starting_values)){
    W <- setStartingValues(.W = t(W), .starting_values = .starting_values) |> 
      t()
  }

  # Normalized Gamma
  Gamma <- Z %*% W
  
  list_UD <- updateUD(
      D = diag(J),
      Eta_normed = Gamma,
      .indicator_type = indicator_type,
      n_constructs = P,
      n_indicators = J,
      n_case = N,
      Z_normed = Z
    )

  U <- list_UD[["U"]]
  D <- list_UD[["D"]]
  # Calculate initial values for the unique variables in U
  # analogous to step 3 of the modified ALS-algorithm
  # D          <- Ij
  
  # Implementation based on the updated MATLAB code provided by Heungsun in 
  # private communication to deal with big dataset (20.10.2021)
  # temp <- svd(t(Gamma)%*%Gamma)
  # gd2 <- diag(temp$d)
  # gv <- temp$v
  
  # GU <- Gamma%*%gv%*%solve(sqrt(gd2))
  
  # M3 <- Z%*%t(D) - GU%*%(t(GU)%*%Z)%*%t(D)
  
  # if(N>J){
  #   temp <- svd(t(M3)%*%M3)
  #   d2 <- diag(temp$d)
  #   v <- temp$v
    
  #   u <- M3%*%v%*%solve(sqrt(d2))
  #   U <- u[,1:J,drop=FALSE]%*%t(v)
  # } else {
  #   temp = svd(M3)
  #   u <- temp$u
  #   v <- temp$v
  #   U <- u[,1:J,drop = FALSE]%*%t(v)
  # }
  
  # Old GSCAm version which faces problem in case of large datasets
  # Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (P+1):N, drop = FALSE]
  # s          <- svd(D %*% t(Z) %*% Gamma_orth)
  # U_tilde    <- s$v %*% t(s$u)
  # U          <- Gamma_orth %*% U_tilde # this is the initial matrix U
  
  
  D <- diag(diag(t(U) %*% Z))
  est  <-  A[which(A != 0)]
  est0 <- est + 1
  iter_counter <- 0
  
  while (
    (!checkConvergence(
      .W_new = est,
      .W_old = est0,
      .tolerance = .tolerance,
      .conv_criterion = .conv_criterion
    )) &&
      (iter_counter < .iter_max)
  ) {
    iter_counter <- iter_counter + 1
    est0 <- est
    X <- Z - U %*% D
    WW <- crossprod(
      x = C,
      y = MASS::ginv(
        (tcrossprod(C) + diag(P) - (2 * B) + (tcrossprod(B)))
      )
    )
    # WW <- t(C) %*% solve(C %*% t(C) + Ip - 2 * B + B %*% t(B))
    for (p in 1:P) {
      # windex_p <- which(W[, p] != 0)
      windex_p <- (rowSums(W[, p, drop = FALSE]) != 0) |>
        which() |>
        names()
      w <- WW[windex_p, p]
      Xp <- X[, windex_p, drop = FALSE]
      # If construct p is a single-indicator construct, dividing w by its norm
      # sometimes doesnt yield 1 exactly due to the way norm() works. This causes
      # problems in verify() since we check if loadings are > 1, which they are

      # if the weight is 1.000000...0002.
      if (length(w) == 1) {
        w <- 1
      } else {
        w <- w / norm(Xp %*% w, type = "2")
      }
      W[windex_p, p] <- w
      Gamma[, p] <- Xp %*% w
    }

    # Step 2a: Update B (the structural coefficients for a given W)
    vcv_gamma <- crossprod(Gamma)
    # vcv_gamma <- t(Gamma) %*% Gamma
    vars_endo <- colnames(B)[colSums(B) != 0]

    beta <- lapply(vars_endo, function(y) {
      x <- (rowSums(B[, y, drop = FALSE]) != 0) |>
        which() |>
        names()
      coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*%
        vcv_gamma[x, y, drop = FALSE]
    })
    B[B != 0] <- unlist(beta)

    # Step 2b: Update C (the loadings for a given W)
    # vars_cf <- which(.csem_model$construct_type == "Common factor")
    # if (length(vars_cf) > 0) {
    cov_gamma_indicators <- crossprod(x = Gamma, y = X)
    # cov_gamma_indicators <- t(Gamma) %*% X
    # Y <- which(colSums(C[vars_cf, ]) != 0)
    Y <- (colSums(C) != 0) |>
      which() |>
      names()
    loadings <- lapply(Y, function(y) {
      x <- (rowSums(C[, y, drop = FALSE]) != 0) |>
        which() |>
        names()
      coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*%
        cov_gamma_indicators[x, y, drop = FALSE]
    })

    # Transform
    tC <- t(C)
    tC[tC != 0] <- unlist(loadings)
    C <- t(tC)
    # }

    list_UD <- updateUD(
      D = D,
      Eta_normed = Gamma,
      .indicator_type = indicator_type,
      n_constructs = P,
      n_indicators = J,
      n_case = N,
      Z_normed = Z
    )

    U <- list_UD[["U"]]
    D <- list_UD[["D"]]

    # Implementation based on the updated MATLAB code provided by Heungsun in
    # private communication to deal with big dataset (20.10.2021)
    # temp <- svd(t(Gamma) %*% Gamma)
    # gd2 <- diag(temp$d)
    # gv <- temp$v

    # GU <- Gamma %*% gv %*% solve(sqrt(gd2))
    # M3 <- Z %*% D - GU %*% (t(GU) %*% Z) %*% t(D)

    # if (N > J) {
    #   temp <- svd(t(M3) %*% M3)
    #   d2 <- diag(temp$d)
    #   v <- temp$v

    #   u <- M3 %*% v %*% solve(sqrt(d2))
    #   U <- u[, 1:J, drop = FALSE] %*% t(v)
    # } else {
    #   temp = svd(M3)
    #   u <- temp$u
    #   v <- temp$v
    #   U <- u[, 1:J, drop = FALSE] %*% t(v)
    # }

    # Old GSCAm implementation
    #
    # Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (P+1):N, drop = FALSE]
    # s          <- svd(D %*% t(Z) %*% Gamma_orth)
    # U_tilde    <- s$v %*% t(s$u)
    # U          <- Gamma_orth %*% U_tilde # this is the initial matrix U

    # D <- diag(diag(t(U) %*% Z))

    # Build A
    A <- cbind(C, B)
    est <- A[which(A != 0)]
  }

  
# Output Formatting ------------------------------------------------------
  # isTRUE(identical(diag(D^2), diag(D)^2))
  D_diag <- diag(D)
  names(D_diag) <- colnames(Z)

  # Retrieve standardized unique_scores
  Unique_scores <- U * normalization_factor
  colnames(Unique_scores) <- colnames(Z)

  # Return
  l <- list(
    "W" = t(W),
    "C" = C,
    "B" = t(B), 
    # Recall that .X is the (unmodified) standardized data passed to this optimization function
    "Construct_scores" = (.X - (Unique_scores %*% D)) %*% W,
    "Unique_loading_estimates" = D_diag,
    "Unique_scores" = Unique_scores,
    "E" = NULL,
    "Modes" = "gsca (gsca_m)",
    "Conv_status" = ifelse(iter_counter > .iter_max, FALSE, TRUE),
    "Iterations" = iter_counter
  )
  return(l)
  
} # END calculateWeightsGSCAm

#' Calculate weights using Integrated Generalised Structured Component Analysis (IGSCA)
#'
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#'
#' In the example section, the specified model is based on the tutorial I-GSCA model associated with GSCA Pro \insertCite{hwangetal2023StructuralEquationModelingAMultidisciplinaryJournal}{cSEM}.
#'
#' **Note**: Here, we assume that there is only one unique loading per indicator.
#'
#' @inheritParams csem
#' @inheritParams csem_arguments
#'
#' @author Michael S. Truong
#' @return List of 4 matrices that make up a fitted I-GSCA Model:
#' * (1) Weights
#' * (2) Loadings
#' * (3) Squared Unique Loadings D^2
#' * (4) Path Coefficients.
#' * (5) Unique Component of Indicators (for common factors) DU
#'
#' @importFrom MASS ginv
#'
#' @references
#'   \insertAllCited{}
#' @examples
#' \dontrun{
#' # Specify the model according to GSCA Pro's example
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 +
#'                       Behavior7 + Behavior8 +  Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2
#'
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 +
#'                     Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 +
#'                 Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 +
#'                 Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 +
#'                  Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 +
#'                      Conscientious6 + Conscientious7 + Conscientious8 +
#'                      Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 +
#'                         Openness7 + Openness8 + Openness9 + Openness10
#'
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion +
#'                      Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#'
#' data(LeDang2022)
#'
#' csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "GSCA",
#' .dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion =
#' "sum_diff_absolute")
#' }
calculateWeightsIGSCA <- function(
  .data,
  .S = args_default()$.S,
  .csem_model = args_default()$.csem_model,
  .conv_criterion = args_default()$.conv_criterion,
  .GSCA_modes = args_default()$.GSCA_modes,
  .iter_max = args_default()$.iter_max,
  .tolerance = args_default()$.tolerance,
  .starting_values = args_default()$.starting_values
) {
  ## Initialize Computational Variables -----------------------------------------------------
  csemify <- parseModel(.model = .csem_model)
  #  Initial Indicators Z0
  Z0 <- .data[, csemify$indicators, drop = FALSE]

  # Igsca assumes \eta \times B, so the rows should be from
  # and the columns should be to. This is in contrast to
  # the rest of cSEM
  B0 <- t(csemify$structural)

  C0 <- csemify$measurement
  W0 <- t(csemify$measurement)

  # whether a construct is a latent or composite variable (con_type)
  con_type <- csemify$construct_type

  # Whether an indicator corresponds to a latent or composite variable (indicator_type)
  # Default value must be "Composite" because of how the measurement matrix
  #  returned by parseModel uses 0 to denote both composite variables and the
  # absence of any corresponding construct variable.
  indicator_type <- vector(
    mode = "character",
    length = ncol(csemify$measurement)
  )
  names(indicator_type) <- csemify$indicators

  composite_indicators <- C0[
    names(csemify$construct_type[csemify$construct_type == "Composite"]),
    ,
    drop = FALSE
  ] |>
    colSums() |>
    sapply(\(x) x == 1)

  indicator_type[composite_indicators] <- "Composite"

  common_factor_indicators <- C0[
    names(csemify$construct_type[csemify$construct_type == "Common factor"]),
    ,
    drop = FALSE
  ] |>
    colSums() |>
    sapply(\(x) x == 1)

  indicator_type[common_factor_indicators] <- "Common factor"
  n_case <- nrow(Z0)
  n_indicators <- ncol(Z0)
  n_constructs <- ncol(W0)
  n_total_var <- n_indicators + n_constructs

  # Normalize data matrix to make normalized Z
  normalization_factor <- sqrt(n_case - 1)
  Z <- Z0 / normalization_factor
  w_index <- which(c(W0) == 1)
  b_index <- which(c(B0) == 1)

  ### Initial Values ------------------------------------
  GSCA_starting_values <- calculateWeightsGSCA(
    .X = .data,
    .S = .S,
    .csem_model = .csem_model,
    .conv_criterion = .conv_criterion,
    .GSCA_modes = .GSCA_modes,
    .iter_max = 10,
    .tolerance = .tolerance,
    .starting_values = .starting_values
  )

  W <- t(GSCA_starting_values$W)
  B <- t(GSCA_starting_values$B)
  C <- GSCA_starting_values$C
  V <- cbind(diag(n_indicators), W)

  # Create initial values for U and D, using normalized data (Z) and construct scores (Eta) --------------

  # Create Initial Normalized Construct Scores Matrix
  Eta <- Z %*% W

  # Initalize unique loadings

  list_UD <- updateUD(
    D = diag(n_indicators),
    Eta_normed = Eta,
    .indicator_type = indicator_type,
    n_constructs = n_constructs,
    n_indicators = n_indicators,
    n_case = n_case,
    Z_normed = Z
  )

  U <- list_UD[["U"]]
  D <- list_UD[["D"]]

  # if starting values are provided
  if (!is.null(.starting_values)) {
    W <- setStartingValues(.W = t(W), .starting_values = .starting_values) |>
      t()
  }

  # Define and Check Modes: 'Common factor', 'NCMP', and 'CCMP'
  modes <- ifelse(
    .csem_model$construct_type == "Composite",
    yes = "CCMP",
    no = .csem_model$construct_type
  )
  if (!is.null(.GSCA_modes)) {
    stopifnot(
      "Invalid .GSCA_modes selected. Only NCMP or CCMP are supported" = all(sapply(
        .GSCA_modes,
        function(x) x %in% c('NCMP', 'CCMP')
      ))
    )

    if (length(names(.GSCA_modes)) != 0) {
      stopifnot(
        "Invalid construct name listed in .GSCA_modes" = length(setdiff(
          names(.GSCA_modes),
          names(.csem_model$construct_type)
        )) ==
          0
      )
      modes[names(.GSCA_modes)] <- .GSCA_modes
    } else if (length(names(.GSCA_modes)) == 0) {
      stopifnot(
        "When passing a global setting to .GSCA_modes, only one choice may be passed" = length(
          .GSCA_modes
        ) ==
          1
      )

      # If only an unnamed element is given (NCMP or CCMP), then set all composites to .GSCA_modes
      if (length(.GSCA_modes) == 1) {
        modes <- ifelse(modes == "CCMP", yes = .GSCA_modes, no = modes)
      }
    }
  }
  # Set Loadings of canonical composites to 0 just in case it wasn't handled properly in the initial values
  C[which(modes == "CCMP"), ] <- 0
  # Create c_index which should only exist for nomological composites and common factors
  C0[which(modes == "CCMP"), ] <- 0
  c_index <- which(c(C0) == 1)

  ## Alternating Least Squares Algorithm -------------------------------------

  ### Optimization Preparation -----------------------------------------------
  # Set the initial estimates based on either the structural model or the loadings
  # if there's no structural model
  if (length(b_index) > 0) {
    est <- B[b_index]
  } else {
    # Because canonical composites do not estimate loadings, here the weights are generally used to determine convergenece, instead of loadings
    est <- W[w_index]
    # est <- C[c_index]
  }
  est0 <- est + 1
  it <- 0

  while (
    (!checkConvergence(
      .W_new = est,
      .W_old = est0,
      .tolerance = .tolerance,
      .conv_criterion = .conv_criterion
    )) &&
      (it <= .iter_max)
  ) {
    # Update Counter Variables
    it <- it + 1
    est0 <- est

    ### Update Weights --------------------------------------------------

    # X is the indicators Z with measurement error removed UD.
    X <- Z - (U %*% D)

    # WW is important for updating the Theta for common factors
    WW <- crossprod(
      x = C,
      y = MASS::ginv(
        (tcrossprod(C) + diag(n_constructs) - (2 * B) + (tcrossprod(B)))
      )
    )

    # WW <- t(C) %*% solve((C %*% t(C) + diag(n_constructs) - (2 * B) + (B %*% t(B))))

    #### for-loop update per construct ---------------------------------------
    A <- cbind(C, B)

    # After each cycle, the Gamma, W and V matrices are updated
    for (eta_idx in seq_len(n_constructs)) {
      tot <- n_indicators + eta_idx
      windex_eta_idx <- (W0[, eta_idx, drop = FALSE] == 1)
      X_eta_idx <- X[, windex_eta_idx, drop = FALSE]

      if (con_type[eta_idx] == "Composite") {
        Theta <-
          updateCompositeTheta(
            W = W,
            A = A,
            X = X,
            V = V,
            eta_idx = eta_idx,
            tot = tot,
            n_constructs = n_constructs,
            n_total_var = n_total_var,
            windex_eta_idx = windex_eta_idx,
            .S = .S
          )
      } else if (con_type[eta_idx] == "Common factor") {
        Theta <- WW[windex_eta_idx, eta_idx, drop = FALSE]
      } else {
        stop("con_type should only either be `Composite` or `Common factor`")
      }

      normed_weights <- Theta / norm(X_eta_idx %*% Theta, "2")
      W[windex_eta_idx, eta_idx] <- normed_weights
      V[windex_eta_idx, tot] <- normed_weights
    }
    Eta <- X %*% W # Trying to save on compute time

    ### Update Loadings, Path Coefficients and Uniqueness Terms ----------
    list_CB <- updateCB(
      X = X,
      Eta = Eta,
      Lambda = C,
      B = B,
      n_indicators = n_indicators,
      .indicator_type = indicator_type,
      n_constructs = n_constructs,
      n_case = n_case,
      lambda_index = c_index,
      b_index = b_index,
      modes = modes
    )

    list2env(list_CB[c("C", "B")], envir = environment())

    ## Uniqueness Scores and Loadings Update ------------------------------------
    list_UD <- updateUD(
      D = D,
      Eta_normed = Eta,
      .indicator_type = indicator_type,
      n_constructs = n_constructs,
      n_indicators = n_indicators,
      n_case = n_case,
      Z_normed = Z
    )

    U <- list_UD[["U"]]
    D <- list_UD[["D"]]

    if (length(b_index) > 0) {
      est <- B[b_index]
    } else {
      est <- W[w_index]
      # Because canonical composites do not estimate loadings, here the weights are generally used to determine convergenece, instead of loadings
      # est <- C[c_index]
    }
  }

  ## Output Formatting -------------------------------------------------------

  # Compute loadings for Canonical Composites
  if (any(modes %in% "CCMP")) {
    CCMP_C <- t(W) %*% .S * .csem_model$measurement
    C[names(modes)[modes == "CCMP"], ] <- CCMP_C[
      names(modes)[modes == "CCMP"],
      ,
      drop = FALSE
    ]
  }

  D_diag <- diag(D)
  names(D_diag) <- colnames(Z)

  # Get the standardized unique scores back from normalized U and zero out the unique scores for composite indicators
  U[, indicator_type == "Composite"] <- 0
  Unique_scores <- U * normalization_factor
  colnames(Unique_scores) <- colnames(Z)

  return(
    list(
      "W" = t(W), # W is J X P, so W^T is P \times J. As shown in the examples ?csem, res$Estimates$Weight_estimates should be P \times J
      "C" = C, # C is P \times J. As shown in the examples ?csem, res$Estimates$Loading_estimates should be P \times J
      "B" = t(B), # B is From \times To; so t(B) is To \times From. As shown in the examples ?csem, res$Estimates$Path_estimates should be To \times From
      # Recall that Z0 is the original standardized data
      "Construct_scores" = (Z0 - (Unique_scores %*% D)) %*% W,
      "Unique_loading_estimates" = D_diag,
      "Unique_scores" = Unique_scores,
      "Modes" = "gsca (igsca)",
      "Conv_status" = ifelse(it > .iter_max, FALSE, TRUE),
      "Iterations" = it
    )
  )
}


#' Calculate composite weights using unit weights
#'
#' Calculate unit weights for all blocks, i.e., each indicator of a block is
#' equally weighted.
#'
#' @usage calculateWeightsUnit(
#'  .S                 = args_default()$.S,
#'  .csem_model        = args_default()$.csem_model,
#'  .starting_values   = args_default()$.starting_values
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A named list. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{`NULL`}
#'   \item{`$Modes`}{The mode used. Always "unit".}
#'   \item{`$Conv_status`}{`NULL` as there are no iterations}
#'   \item{`$Iterations`}{0 as there are no iterations}
#' }
#'
calculateWeightsUnit = function(
  .S                 = args_default()$.S,
  .csem_model        = args_default()$.csem_model,
  .starting_values   = args_default()$.starting_values
){
  
  W <- .csem_model$measurement
  
  # if starting values are provided
  if(!is.null(.starting_values)){
    W = setStartingValues(.W = W, .starting_values = .starting_values)
  }
  
  W <- scaleWeights(.S, W)
  
  modes        <- rep("unit", nrow(W))
  names(modes) <- rownames(W)
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = modes, "Conv_status" = NULL,
            "Iterations" = 0)
  return(l)  
}

#' Calculate composite weights using principal component analysis (PCA)
#'
#' Calculate weights for each block by extracting the first principal component
#' of the indicator correlation matrix S_jj for each blocks, i.e., weights
#' are the simply the first eigenvector of S_jj.
#'
#' @usage calculateWeightsPCA(
#'  .S                 = args_default()$.S,
#'  .csem_model        = args_default()$.csem_model
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A named list. J stands for the number of constructs and K for the number
#' of indicators.
#' \describe{
#'   \item{`$W`}{A (J x K) matrix of estimated weights.}
#'   \item{`$E`}{`NULL`}
#'   \item{`$Modes`}{The mode used. Always "PCA".}
#'   \item{`$Conv_status`}{`NULL` as there are no iterations}
#'   \item{`$Iterations`}{0 as there are no iterations}
#' }
#'

calculateWeightsPCA = function(
  .S                 = args_default()$.S,
  .csem_model        = args_default()$.csem_model
){
  
  W <- .csem_model$measurement
  
  for(j in 1:nrow(W)) {
    block      <- rownames(W[j, , drop = FALSE])
    indicators <- W[block, ] != 0
    
    temp <- psych::principal(r = .S[indicators, indicators], nfactors = 1)
    # Note that temp$weights is the same as:
    #  1. Calculating the eigenvectors of S[indicators, indicators] and
    #     then taking the first eigenvector. This is the unscaled weight w
    #  2. Obtain the scaled weight by w_s = w / sqrt(w' S[i, i] w)
    # 
    #  R Code:
    #    tt <- eigen(.S[indicators, indicators])
    #    w  <- tt$vectors[, 1]
    #    w  <- w / sqrt(w %*% .S[indicators, indicators] %*% w)
    
    W[block, indicators] <- c(temp$weights)
  }
  
  modes        <- rep("PCA", nrow(W))
  names(modes) <- rownames(W)
  
  # Return
  l <- list("W" = W, "E" = NULL, "Modes" = modes, "Conv_status" = NULL,
            "Iterations" = 0)
  return(l)  
}
