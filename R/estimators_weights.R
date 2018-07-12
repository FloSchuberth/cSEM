#' Calculate composite weights
#'
#' Calculates weights
#'
#' More description here
#'
#' @usage calculateWeights(
#'   .data                     = NULL,
#'   .model                    = NULL,
#'   .PLS_mode                 = NULL,
#'   .tolerance                = NULL,
#'   .iter_max                 = NULL,
#'   .PLS_weight_scheme_inner  = NULL,
#'   .ignore_structural_model  = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of estimated weights.
#'
#' @export
#'

calculateWeightsPLS <- function(
  .data                     = NULL,
  .model                    = NULL,
  .PLS_mode                 = NULL,
  .tolerance                = NULL,
  .iter_max                 = NULL,
  .PLS_weight_scheme_inner  = NULL,
  .ignore_structural_model  = NULL
) {


  ### Make the function more autonomous ========================================
  ## Convert model to "cSEMModel" format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    csem_model <- parseModel(.model)
  } else {
    csem_model <- .model
  }

  # ## Prepare, standardize, check, and clean data if not already in this format
  # if(!(class(.data) == "cSEMData")) {
  #   if(is.matrix(.data) && isSymmetric.matrix(.data)) {
  #     S <- .data
  #   } else {
  #     X <- processData(.data = .data, .model = csem_model) 
  #     S <- stats::cor(X)
  #   }
  # }

  ### Preparation ==============================================================
  ## Get/set the modes for the outer estimation

  if(is.null(.PLS_mode)) {
    modes <- apply(csem_model$construct_type, 1, function(x) {

      if(x["Type"] == "Common factor") {
        "ModeA"
      } else {
        "ModeB"
      }
    })
    names(modes) <- csem_model$construct_type$Name

  } else if(all(.PLS_mode %in% c("ModeA", "ModeB"))) {
    if(setequal(names(.PLS_mode), csem_model$construct_type$Name)) {
      modes <- .PLS_mode
    } else if (length(.PLS_mode) == 1) {
      modes <- sapply(csem_model$construct_type$Name, function(x) .PLS_mode)
    } else {
      stop("At least one mode has not been specified.")
    }
  } else {
    stop(paste(setdiff(.PLS_mode, c("ModeA", "ModeB")), collapse = ", "),
         " is an unknown Mode.")
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
      .S                       = S,
      .W                       = W_iter,
      .csem_model              = csem_model,
      .PLS_weight_scheme_inner = .PLS_weight_scheme_inner,
      .ignore_structural_model = .ignore_structural_model
    )
    # Outer estimation

    W <- calculateOuterWeightsPLS(
      .S        = S,
      .W        = W_iter,
      .E        = E,
      .PLS_mode = modes
    )

    # Scale weights
    W <- scaleWeights(S, W)

    # Check for convergence
    if(max(abs(W_iter - W)) < .tolerance) {
      break # return iterative PLS weights
    } else if(iter_counter == iter_max & iter_max == 1) {
      break # return one-step PLS weights
    } else if(iter_counter == iter_max & iter_max > 1) {
      warning("Iteration did not converge after ", iter_max, " steps. ",
              "Last weights are returned.")
    } else {
      W_iter <- W
    }
  }

  # Return
  l <- list("W" = W, "E" = E, "Modes" = modes)
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
#'   .S                        = NULL,
#'   .W                        = NULL,
#'   .csem_model               = NULL,
#'   .PLS_weight_scheme_inner  = NULL,
#'   .ignore_structural_model  = NULL
#' )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x J) matrix of inner weights.
#'
calculateInnerWeightsPLS <- function(
  .S                        = NULL,
  .W                        = NULL,
  .csem_model               = NULL,
  .PLS_weight_scheme_inner  = NULL,
  .ignore_structural_model  = NULL
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
  if(.PLS_weight_scheme_inner == "path" & .ignore_structural_model) {
    .ignore_structural_model <- FALSE
    warning("Structural model required for the path weighting scheme.\n",
            ".ignore_structural_model = TRUE was changed to FALSE.")
  }

  if(.ignore_structural_model) {
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
#'    .mode   = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A (J x K) matrix of outer weights.
#'

calculateOuterWeightsPLS <- function(
  .S              = NULL,
  .W              = NULL,
  .E              = NULL,
  .PLS_mode       = NULL
) {
  # Covariance/Correlation matrix between each proxy and all indicators (not
  # only the related ones). Note: Cov(H, X) = WS, since H = XW'.
  W <- .W
  proxy_indicator_cor <- .E %*% W %*% .S

  for(i in 1:nrow(W)) {
    block      <- intersect(rownames(W[i, , drop = FALSE]), names(.PLS_mode))
    indicators <- W[block, ] != 0

    if(.PLS_mode[block] == "ModeA") {
      ## Mode A - Regression of each indicator on its corresponding proxy
      # Select only
      W[block, indicators] <- proxy_indicator_cor[block, indicators]

    } else if(.PLS_mode[block] == "ModeB") {
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
#' "*SUMCOR*", "*MAXVAR*", "*SSQCOR*", "*MINVAR*", and "*GENVAR*"
#' suggested by \insertCite{Kettenring1971;textual}{cSEM}.
#'
#' Some more description...
#'
#' @usage calculateWeightsKettenring(.data, .model, .criteria)
#'
#' @inheritParams csem_arguments
#'
#' @inherit calculateWeightsPLS return
#'

calculateWeightsKettenring <- function(
  .data           = NULL,
  .model          = NULL,
  .criteria       = NULL
) {

  stop("Not yet implemented")

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
