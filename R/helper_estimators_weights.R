#' Internal: Calculate the inner weights for PLS-PM
#'
#' PLS-PM forms "inner" composites as a weighted sum of its *I* related composites.
#' These inner weights are obtained using one of the following schemes \insertCite{Lohmoeller1989}{cSEM}:
#' \describe{
#'   \item{`centroid`}{According to the centroid weighting scheme each inner weight used
#'     to form composite *j* is either 1 if the correlation between composite *j* and 
#'     its via the structural model related composite *i = 1, ..., I* is positive 
#'     and -1 if it is negative.}
#'   \item{`factorial`}{According to the factorial weighting scheme each inner weight used
#'     to form inner composite *j* is equal to the correlation between composite *j* 
#'     and its via the structural model related composite *i = 1, ..., I*.}
#'   \item{`path`}{Lets call all construct that have an arrow pointing to construct *j*
#'     **predecessors of j** and all arrows going from j to other constructs **followers of j**.
#'     According the path weighting scheme, inner weights are computed as follows.
#'     Take construct *j*: 
#'     \itemize{
#'       \item For all predecessors of *j* set the inner weight of predecessor 
#'             *i* to the correlation of *i* with *j*.
#'       \item For all followers of *j* set the inner weight of follower *i* to 
#'             the coefficient of a multiple regression of *j* on all 
#'             followers *i* with *i = 1,...,I*.
#'    }}
#' }
#' Except for the path weighting scheme relatedness can come in two flavors.
#' If `.PLS_ignore_structural_model = TRUE` all constructs are considered related.
#' If `.PLS_ignore_structural_model = FALSE` (the default) only adjacent constructs
#' are considered. If `.PLS_ignore_structural_model = TRUE` and `.PLS_weight_scheme_inner = "path"`
#' a warning is issued and `.PLS_ignore_structural_model` is changed to `FALSE`.
#' 
#' @usage calculateInnerWeightsPLS(
#'   .S                           = args_default()$.S,
#'   .W                           = args_default()$.W,
#'   .csem_model                  = args_default()$.csem_model,
#'   .PLS_ignore_structural_model = args_default()$.PLS_ignore_structrual_model,
#'   .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner
#' )
#' @inheritParams csem_arguments
#'
#' @return The (J x J) matrix `E` of inner weights.
#' @keywords internal

calculateInnerWeightsPLS <- function(
  .S                           = args_default()$.S,
  .W                           = args_default()$.W,
  .csem_model                  = args_default()$.csem_model,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structrual_model,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner
) {
  
  # Composite correlation matrix (C = V(H))
  C <- .W %*% .S %*% t(.W)
  
  # Due to floting point errors may not be symmetric anymore. In order
  # prevent that replace the lower triangular elements by the upper
  # triangular elements
  
  C[lower.tri(C)] <- t(C)[lower.tri(C)]
  
  # Structural model relationship; if only correlations are specified
  # use these
  tmp <- rownames(.csem_model$structural)
  if(sum(rowSums(.csem_model$structural)) == 0) {
    if(.PLS_weight_scheme_inner == "path") {
      stop2("Structural model is required for the PLS path weighting scheme.\n",
            "Please change the inner weighting scheme or supply a path model.")
    }
    D <- .csem_model$cor_specified[tmp, tmp]
  } else {
    E <- .csem_model$structural[tmp, tmp]
    D <- E + t(E)
  }
  
  # Note: June 2019
  if(any(D == 2)) { # non recursive model
    # Set elements back to 1 
    D[D == 2] <- 1 
  }
  
  ## (Inner) weighting scheme:
  if(.PLS_weight_scheme_inner == "path" & .PLS_ignore_structural_model) {
    .PLS_ignore_structural_model <- FALSE
    warning("Structural model is required for the path weighting scheme.\n",
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
              #    predecessor i to the correlation of i with j.
              #  - For all followers of j: set the inner weight of follower i
              #    to the coefficient of a multiple regression of j on
              #    all followers i with i = 1,...,I.
              
              E_temp <- E
              ## Assign predecessor relation
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

#' Internal: Calculate the outer weights for PLS-PM
#'
#' Calculates outer weights in PLS-PM. Currently, the originally suggested mode A
#' and mode B are suggested. Additionally, non-negative least squares (modeBNNLS) and 
#' weights of principal component analysis (PCA) are implemented.  
#'
#' @usage calculateOuterWeightsPLS(
#'    .data   = args_default()$.data,  
#'    .S      = args_default()$.S,
#'    .W      = args_default()$.W,
#'    .E      = args_default()$.E,
#'    .modes  = args_default()$.modes
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @return A (J x K) matrix of outer weights.
#' @keywords internal

calculateOuterWeightsPLS <- function(
  .data   = args_default()$.data,
  .S      = args_default()$.S,
  .W      = args_default()$.W,
  .E      = args_default()$.E,
  .modes  = args_default()$.modes
) {
  # Covariance/Correlation matrix between each proxy and all indicators (not
  # only the related ones). Note: Cov(H, X) = WS, since H = XW'.
  W <- .W
  proxy_indicator_cor <- .E %*% W %*% .S
  
  # Scale the inner proxy. Inner weights are usually scaled such that the inner
  # proxies are standardized (mean = 0, var = 1)
  inner_proxy <- scale(.data %*% t(W) %*% t(.E))
  colnames(inner_proxy) = rownames(W)
  
  # Compute outer weights by block/ construct
  for(i in 1:nrow(W)) {
    block      <- rownames(W[i, , drop = FALSE])
    indicators <- W[block, ] != 0
    
    if(is.numeric(.modes[[block]]) & length(.modes[[block]]) > 1) {
      if(length(.modes[[block]]) == sum(indicators)) {
        ## Fixed weights - Each weight of "block" is fixed to a user-given value
        W[block, indicators] <- .modes[[block]]
      } else {
        stop2("Construct ", paste0("`", block, "` has ", sum(indicators), 
                                   " indicators but only ",
                                  length(.modes[[block]]), " fixed weights are provided.")) 
      }
    } else if(.modes[block] == "modeA") {
      ## Mode A - Regression of each indicator on its corresponding proxy
      W[block, indicators] <- proxy_indicator_cor[block, indicators]
      
    } else if(.modes[block] == "modeB") {
      ## Mode B - Regression of each proxy on all its indicator
      # W_j = S_jj^{-1} %*% Cov(eta_j, X_j)
      W[block, indicators] <- solve(.S[indicators, indicators]) %*% proxy_indicator_cor[block, indicators]
      
    } else if(.modes[block] == "modeBNNLS"){
      ## Mode BNNLS - Regression of each proxy on its indicators using non-negative LS
      # Note: .data is standardized, i.e., mean 0 and unit variance, inner proxy is also 
      #       standardized (standardization of the inner proxy 
      #       apparently has no effect though.)
      if (!requireNamespace("nnls", quietly = TRUE)) {
        stop2(
          "Package `nnls` needed for \"modeBNNLS\" to work. Use `install.packages(\"nnls\")` and rerun.")
      }
      temp <- nnls::nnls(A = .data[, indicators, drop = FALSE], b = inner_proxy[,block])
      W[block, indicators] <- temp$x
    
    } else if(.modes[block] == "PCA"){
      ## PCA - Weights to create the first principal component are used  (= the first eigenvector of
      ##       of S_jj).
      temp <- psych::principal(r = .S[indicators, indicators], nfactors = 1)
      W[block, indicators] <- c(temp$weights)
      
    } 
    # Set weights of single-indicator constructs to 1 (in order to avoid floating
    # point imprecision)
    if(sum(indicators) == 1){
      W[block, indicators] <- 1 
    } 
    
    # If .modes[block] == "unit" or a single value has been given, nothing needs
    # to happen since W[block, indicators] would be set to 1 (which it already is). 
  }
  return(W)
} # END calculateOuterWeights

#' Internal: Check convergence
#'
#' Check convergence of an algorithm using one of the following criteria:
#' \describe{
#'   \item{`diff_absolute`}{Checks if the largest elementwise absolute difference
#'                          between two matrices `.W_new` and `W.old` is 
#'                          smaller than a given tolerance.}
#'   \item{`diff_squared`}{Checks if the largest elementwise squared difference
#'                         between two matrices `.W_new` and `W.old` is 
#'                         smaller than a given tolerance.}
#'   \item{`diff_relative`}{Checks if the largest elementwise absolute rate of change
#'                          (new - old / new) for two matrices `.W_new` 
#'                          and `W.old` is smaller than a given tolerance.}
#'   \item{`sum_diff_absolute`}{Checks if the sum of the element-wise absolute
#'                              difference between two matrices `.W_new` and `W.old` is smaller than a
#'                              given tolerance}
#'   \item{`mean_diff_absolute`}{Checks if the mean of the element-wise absolute
#'                              difference between two matrices `.W_new` and `W.old` is smaller than a
#'                              given tolerance
#'   }
#' }
#'
#' @usage checkConvergence(
#'   .W_new          = args_default()$.W_new,
#'   .W_old          = args_default()$.W_old,
#'   .conv_criterion = args_default()$.conv_criterion,
#'   .tolerance      = args_default()$.tolerance
#'   )
#'
#' @inheritParams csem_arguments
#' 
#' @return `TRUE` if converged; `FALSE` otherwise.
#' @keywords internal

checkConvergence <- function(
  .W_new          = args_default()$.W_new,
  .W_old          = args_default()$.W_old,
  .conv_criterion = args_default()$.conv_criterion,
  .tolerance      = args_default()$.tolerance
  ){
  ## Check if correct value is provided:
  match.arg(.conv_criterion, args_default(.choices = TRUE)$.conv_criterion)
  
  switch (.conv_criterion,
    "diff_absolute" = {
      max(abs(.W_old - .W_new)) < .tolerance
    },
    "diff_squared"  = {
      max((.W_old - .W_new)^2) < .tolerance
    },
    "diff_relative" = {
      max(abs((.W_old[.W_new != 0] - .W_new[.W_new != 0]) /
                .W_new[.W_new != 0])) < .tolerance
    }, 
    "sum_diff_absolute" = {
      (sum(abs(.W_old - .W_new))) < .tolerance
    },
    "mean_diff_absolute" = {
      (mean(abs(.W_old - .W_new))) < .tolerance
    }
  )
}

#' Internal: Scale weights
#'
#' Scale weights such that the formed composite has unit variance.
#'
#' @usage scaleWeights(
#'   .S = args_default()$.S, 
#'   .W = args_default()$.W
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of scaled weights.
#' @keywords internal

scaleWeights <- function(
  .S = args_default()$.S, 
  .W = args_default()$.W
  ) {
  
  ## Calculate the variance of the proxies:
  var_proxies <- diag(.W %*% .S %*% t(.W))
  
  # Using the solve function is suboptimal as if one of the proxies' variances
  # is closed to 0, matrix inversion might not work. 
  # W_scaled <- solve(diag(sqrt(var_proxies), 
                         # nrow = length(var_proxies),
                         # ncol = length(var_proxies)
                         # )) %*% .W
  
  ## Scale the weights to ensure that the proxies have a variance of one  
  ### For GSCA_M and IGSCA models it multiplies the weights for indicators of common factors by a construct-specific constant
  W_scaled <- diag(1/sqrt(var_proxies)) %*%.W
  
  ## Assign rownames and colnames to the scaled weights and return
  rownames(W_scaled) <- rownames(.W)
  colnames(W_scaled) <- colnames(.W)
  
  return(W_scaled)
}

#' Internal: Set starting values
#'
#' Set the starting values.
#'
#' @usage setStartingValues(
#'   .W               = args_default()$.W,
#'   .starting_values = args_default()$.starting_values
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of starting values.
#' @keywords internal

setStartingValues = function(.W = args_default()$.W,
                             .starting_values = args_default()$.starting_values){

  if(!is.list(.starting_values)){
    stop2(
      "The following error occured in the `setStartingValues()` function:\n",
      "Starting values must be as a list."
      )
  }
  
  tmp <- setdiff(names(.starting_values), rownames(.W))
  
  if(length(tmp) != 0) {
    stop2(
      "The following error occured in the `setStartingValues()` function:\n",
      "Construct name(s): ", paste0("`", tmp, "`", collapse = ", "), 
      " provided to `.starting_values`", 
      ifelse(length(tmp) == 1, " is", " are"), " unknown.")
  }
  
  # Replace the original ones by the starting value
  for(i in names(.starting_values)) {
    ## Error if starting values for construct i have not been names
    if(is.null(names(.starting_values[[i]]))) {
      stop2(
        "The following error occured in the `setStartingValues()` function:\n",
        "Starting weights must be named."
        )
    }
    # tmp <- setdiff(names(.starting_values[[i]]), colnames(W[i,,drop=FALSE]))
    tmp <- setdiff(names(.starting_values[[i]]), colnames(.W[i,.W[i,]!=0,drop = FALSE]))
    
    
    if(length(tmp) != 0) {
      stop2(
        "The following error occured in the `setStartingValues()` function:\n",
        "Indicator name(s): ", paste0("`", tmp, "`", collapse = ", "), 
        " provided to `.starting_values`", 
        ifelse(length(tmp) == 1, " is", " are"), " unknown.")
    }
    
    .W[i,names(.starting_values[[i]])] = .starting_values[[i]]
    
  }
  
  return(.W)
}

#' Update Theta for Composite Variables in IGSCA
#' 
#' It is unintuitive that X is used here, seeing as how X = Z-UD; and we use X to update composite variables. 
#' However, from a non-computational point of view, it shouldn't matter because for the composite indicators, X_comp = Z_comp
#'
#' @param W Weights matrix
#' @param A Stacked matrix of loadings and path coefficients \eqn{\left[\Lambda \mid B \right]}
#' @param V Stacked matrix of identity matrix and weights \eqn{\left[I \mid W \right]}
#' @param X The matrix X is equal to \eqn{Z - UD}
#' @param windex_eta_idx Index of weights related to the indicators for the construct of interest
#' @param n_total_var Number of indicators and constructs
#' @param tot Index dependent on which construct variable we are examining
#' @param n_constructs Number of constructs
#' @param eta_idx Index of which construct we are examining
#' @inheritParams csem_arguments
#' @importFrom MASS ginv
#' @return Theta: A matrix that will later be used to update the weights for the composite variable.
#'
#'
updateCompositeTheta <-
  function(
    W,
    A,
    V,
    X,
    windex_eta_idx,
    n_total_var,
    tot,
    n_constructs,
    eta_idx,
    .S = args_default()$.S
  ) {
    # The following code is based on the ASGSCA package (licensed
    # under GPL-3). Notation is adapted to be conform with the notation of the
    # cSEM package
    e <- matrix(0, nrow = 1, ncol = n_total_var)
    e[tot] <- 1
    H1 <- diag(n_total_var)
    H2 <- diag(n_constructs)
    H1[tot, tot] <- 0
    H2[eta_idx, eta_idx] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)

    beta <- e - A[eta_idx, , drop = FALSE]

    Theta <- tcrossprod(
      x = MASS::ginv(
        as.numeric(beta %*% t(beta)) *
          .S[windex_eta_idx, windex_eta_idx, drop = FALSE]
      ),
      y = beta %*% t(Delta) %*% .S[, windex_eta_idx, drop = FALSE]
    )

    # Theta <- MASS::ginv(
    #   as.numeric(beta %*% t(beta)) * .S[windex_eta_idx, windex_eta_idx, drop = FALSE]
    # ) %*%
    #   t(beta %*% t(Delta) %*% .S[, windex_eta_idx, drop = FALSE])

    return(Theta)

    # Kronecker Method
    # vecZDelta <- c(X %*% Delta) 
    # XI <- kronecker(t(beta), X)
    # XI <- XI[, windex_eta_idx]
    # XI <- kroneckerC(t(beta), X, which(windex_eta_idx))
    # Theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta
  }

#' Update Loadings and Path-Coefficients for IGSCA
#'
#' @param X Indicators with measurement error removed
#' @param Eta Construct Scores
#' @param Lambda Loadings matrix
#' @param B Path coefficients matrix
#' @param n_indicators Number of indicators
#' @param n_constructs Number of oncstructs
#' @param lambda_index Index of loadings
#' @param b_index Index of Path Coefficients
#' @param n_case Number of Cases
#' @param .indicator_type Vector of whether each indicator corresponds to a common factor or composite
#' @param modes Named vector of whether the construct is a Common factor, nomological composite or canonical composite.
#' @importFrom MASS ginv
#' @return List of matrices:
#'
#' * (1) Estimated Loadings matrix (C)
#' * (2) Estimated Path Coefficients matrix (B)
#'
updateCB <-
  function(
    X,
    Eta,
    Lambda,
    B,
    .indicator_type,
    n_indicators,
    lambda_index,
    n_constructs,
    b_index,
    n_case,
    modes
  ) {
    # Loading Update ----------------------------------------------------------

    # Kronecker bypass
    # browser()
    vars_cf_ncmp <- names(modes)[modes %in% c("Common factor", "NCMP")]
    # cov_eta_indicators <- t(Eta) %*% X 
    cov_eta_indicators <- crossprod(Eta, X) 
    # cor_eta <- t(Eta) %*% Eta 
    cor_eta <- crossprod(Eta) 
    

    dep_vars <- (colSums(Lambda[vars_cf_ncmp, , drop = FALSE]) != 0) |> 
        which() |> 
        names()
    # This approach assumes that every factor/NCMP loads onto one indicator: no cross-loadings
    loadings <- lapply(dep_vars, function(y) {
      x <- (rowSums(Lambda[vars_cf_ncmp, y, drop = FALSE]) != 0) |> 
          which() |> 
          names()
      coef <- MASS::ginv(cor_eta[x, x, drop = FALSE]) %*% cov_eta_indicators[x, y, drop = FALSE]
    })
    # A future approach should consider avoiding c_index and using explicit names, for safety.
    Lambda[lambda_index] <- unlist(loadings, use.names =  FALSE)

    # Kronecker Approach and Assumes All Composites are Nomological
    # t1 <- c(X)
    # M1 <- kroneckerC(diag(n_indicators), Eta, c_index)
    # C[c_index] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)

    # Path Coefficients Update ------------------------------------------------
    vars_endo <- colnames(B)[colSums(B) != 0]
    beta <- lapply(vars_endo, function(y) {
      x <- (rowSums(B[, y, drop = FALSE]) != 0) |> 
        which() |> 
        names()
      coef <- MASS::ginv(cor_eta[x, x, drop = FALSE]) %*%
        cor_eta[x, y, drop = FALSE]
    })
    B[b_index] <- unlist(beta, use.names = FALSE)

    return(
      list(
        "C" = Lambda,
        "B" = B
      )
    )
  }


#' Update unique scores and unique loadings
#' 
#' Intended to be used within the alternating least squares algorithm for either GSCA_M or IGSCA. Assumes that the construct scores and data are normalized.
#'
#' @param D Unique loadings
#' @param Eta_normed Normalized data
#' @param Z_normed Normalized data
#' @param n_constructs Number of constructs
#' @param n_case Number of cases
#' @param n_indicators Number of indicators
#' @param .indicator_type Vector of whether each indicator corresponds to a common factor or composite
#' @returns List of 2 elements, normalized unique scores (`U`) and normalized unique loadings (`D`)
#' 
#'
updateUD <- function(D, Eta_normed, .indicator_type, n_constructs, n_case, n_indicators, Z_normed) {

  qr_eta <- qr(Eta_normed)
  # QtZ_null <- qr.qty(qr_eta, Z_normed)[(n_constructs + 1):n_case, , drop = FALSE]
  svd_mx <- svd(tcrossprod(x = D, y = qr.qty(qr_eta, Z_normed)[(n_constructs + 1):n_case, , drop = FALSE]))
  #  svd_mx <- svd(D %*% t(QtZ_null))
  Utilde <-   # (N-P) × J
  U <- qr.qy(qr_eta, rbind(matrix(0, n_constructs, n_indicators),  tcrossprod(x= svd_mx$v, y = svd_mx$u)))
  # U <- qr.qy(qr_eta, rbind(matrix(0, n_constructs, n_indicators), svd_mx$v %*% t(svd_mx$u)))

  # Old method based on Hwang et al. (2017) — O(N^2) memory and computation
  # Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
  #   (n_constructs + 1):n_case,
  #   drop = FALSE
  # ]
  # svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
  # Utilde <- svd_mx$v %*% t(svd_mx$u)
  # U <- Eta_Q2 %*% Utilde
  
  # U[, .indicator_type == "Composite"] <- 0

  # Update Unique Loadings

  # D <- diag(diag(t(U) %*% Z_normed))
  D <- diag(diag(crossprod(U, Z_normed)))
  D[.indicator_type == "Composite", .indicator_type == "Composite"] <- 0

  # Return output
  return(list("U" = U, "D" = D))
} 

#' Block Diagonalize Estimated Parameter Matrices to Facilitate Computation of FIT Statistics
#'
#' Block diagonalizes the estimated paramater matrices as shown on Equations
#' 3.28-3.29 on page 111 of \insertCite{Hwang2014;textual}{cSEM}. Should only be used on multi-group models
#'
#' @inheritParams tidy.cSEMResults
#'
#' @return cSEMResults in single-group data structure with block diagonalized parameter estimates
#' @importFrom Matrix bdiag
bdiagonalizeMultiGroupIgscaEstimates <- function(x) {
  if (!identical(names(x), c("Estimates", "Information"))) {
    # Multi-Group Code

    ## Extract Matrices ------------------------------------------------------
    # Extract Estimated Matrices
    # TODO: Test whether the following code works

    estimates_to_be_extracted <- list(
      "Path_estimates",
      "Loading_estimates",
      "Weight_estimates",
      "Construct_scores",
      "Unique_scores"
    )

    names(estimates_to_be_extracted) <- unlist(estimates_to_be_extracted)

    extraction <- lapply(
      X = estimates_to_be_extracted,
      function(matrix_name, multigroup_output) {
        extraction <- lapply(
          multigroup_output,
          function(onegroup_output, matrix_name) {
            return(onegroup_output[["Estimates"]][[matrix_name]])
          },
          matrix_name = matrix_name
        )
        return(extraction)
      },
      multigroup_output = x
    )

    # Extract Data
    extraction$Data <- lapply(x, function(onegroup_output) {
      return(onegroup_output[["Information"]][["Data"]])
    })

    # Remove Null Matrices
    extracts_to_remove <- lapply(extraction, \(x) {
      lapply(x, is.null) |> unlist() |> all()
    })

    extraction <- extraction[which(
      !unlist(lapply(extraction, \(x) lapply(x, is.null) |> unlist() |> all()))
    )]

    extraction <- mapply(
      function(extract, extract_name) {
        # We don't keep it as a sparse matrix because that might break
        # functionality with other functions unless much more of Matrix is
        # imported

        bdiaged <- Matrix::bdiag(extract)
        colnames(bdiaged) <- rep(
          colnames(extract[[1]]),
          times = length(extract)
        )
        if (
          !(extract_name %in% c("Data", "Construct_scores", "Unique_scores"))
        ) {
          rownames(bdiaged) <- rep(
            rownames(extract[[1]]),
            times = length(extract)
          )
        }

        return(as.matrix(bdiaged))
      },
      extract = extraction,
      extract_name = names(extraction),
      SIMPLIFY = FALSE
    )

    ## Create Surrogate Output -----------------------------------------------
    surrogate_out <- list()
    # Insert the diagonalized matrices into the surrogate
    ## It's not simple to take x[[1]] as the surrogate structure because it has
    ## many estimatates (such as reliabilities) that are specific to the group
    ## model

    for (extract_name in names(extraction)[which(
      names(extraction) != "Data"
    )]) {
      surrogate_out[["Estimates"]][[extract_name]] <- extraction[[extract_name]]
    }

    surrogate_out[["Information"]][["Data"]] <- extraction[["Data"]]

    return(surrogate_out)
  } else if (identical(names(x), c("Estimates", "Information"))) {
    # Single-Group Code
    stop("This function is only meant for multi-group models.")
  }
}