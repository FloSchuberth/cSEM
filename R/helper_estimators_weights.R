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

#' Internal: Check for complete starting values
#'
#' Check whether `.starting_values` completely specifies the starting weights,
#' i.e., whether every construct that has at least one free (nonzero-pattern)
#' weight is present in `.starting_values`, and, for each such construct, the
#' named starting-value vector covers exactly (`setequal()`) that construct's
#' pattern indicators. Partial specification (e.g., missing constructs or
#' missing indicators within a specified construct) returns `FALSE`.
#'
#' @usage hasCompleteStartingValues(
#'   .W               = args_default()$.W,
#'   .starting_values = args_default()$.starting_values
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A single logical.
#' @keywords internal

hasCompleteStartingValues <- function(.W = args_default()$.W,
                                       .starting_values = args_default()$.starting_values){

  if(!is.list(.starting_values)) {
    return(FALSE)
  }

  # Constructs that have at least one free (nonzero-pattern) weight
  constructs_with_free_weights <- rownames(.W)[rowSums(.W != 0) > 0]

  if(length(constructs_with_free_weights) == 0) {
    return(FALSE)
  }

  if(!all(constructs_with_free_weights %in% names(.starting_values))) {
    return(FALSE)
  }

  all(vapply(constructs_with_free_weights, function(i) {
    pattern_indicators <- colnames(.W)[.W[i, ] != 0]
    given_indicators   <- names(.starting_values[[i]])

    !is.null(given_indicators) && setequal(pattern_indicators, given_indicators)
  }, logical(1)))
}

#' Internal: Is `.starting_values` a multi-start specification?
#'
#' A multi-start specification is a non-empty list whose every element is
#' itself a list (i.e., a list of `.starting_values` sets), as opposed to a
#' single set (a named list of named numeric vectors).
#'
#' @inheritParams csem_arguments
#' @return A single logical.
#' @keywords internal

isMultiStartStartingValues <- function(.starting_values = args_default()$.starting_values) {
  is.list(.starting_values) &&
    length(.starting_values) > 0 &&
    all(vapply(.starting_values, is.list, logical(1)))
}

#' Internal: Validate a random multi-start specification
#'
#' Checks that `.starting_values` is a named numeric vector of length three with
#' names `n`, `min` and `max`, a positive integer `n`, finite `min`/`max`, and
#' `min < max`. Returns the parsed values.
#'
#' @inheritParams csem_arguments
#' @return A list with elements `n` (integer), `min` and `max` (numeric).
#' @keywords internal

checkRandomStartingSpec <- function(.starting_values = args_default()$.starting_values) {

  x  <- .starting_values
  nm <- names(x)

  if(length(x) != 3 || is.null(nm) || !setequal(nm, c("n", "min", "max"))) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "A numeric `.starting_values` must be a named vector of length 3 with names ",
      "`n`, `min` and `max`, e.g. `c(n = 10, min = -1, max = 1)`.")
  }

  n  <- x[["n"]]
  mn <- x[["min"]]
  mx <- x[["max"]]

  if(!all(is.finite(c(n, mn, mx)))) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "`n`, `min` and `max` in `.starting_values` must all be finite.")
  }

  if(n < 1 || n != round(n)) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "`n` in `.starting_values` must be a positive integer.")
  }

  if(mn >= mx) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "`min` must be strictly smaller than `max` in `.starting_values`.")
  }

  list(n = as.integer(n), min = mn, max = mx)
}

#' Internal: Generate random starting-value sets
#'
#' For every construct in `.measurement` that has at least one free
#' (nonzero-pattern) weight, sample `runif(., .min, .max)` starting weights for
#' its pattern indicators. Repeated `.n` times to produce a recursive list of
#' complete `.starting_values` sets suitable for GSCA multi-start.
#'
#' @param .measurement A (construct x indicator) 0/1 pattern matrix
#'   (`parseModel()$measurement`).
#' @param .n Integer. Number of random sets to generate.
#' @param .min Numeric. Lower bound passed to [stats::runif()].
#' @param .max Numeric. Upper bound passed to [stats::runif()].
#'
#' @return A list of length `.n`; each element is a named list of named numeric
#'   vectors.
#' @keywords internal

generateRandomStartingValues <- function(.measurement, .n, .min, .max) {

  free_constructs <- rownames(.measurement)[rowSums(.measurement != 0) > 0]

  lapply(seq_len(.n), function(i) {
    sv <- lapply(free_constructs, function(cn) {
      inds <- colnames(.measurement)[.measurement[cn, ] != 0]
      stats::setNames(stats::runif(length(inds), min = .min, max = .max), inds)
    })
    names(sv) <- free_constructs
    sv
  })
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
#' Ensures that the unique scores have a mean of 0.
#' 
#' Following page 879 of \insertCite{trendafilovExploratoryFactorAnalysis2011;textual}{cSEM}, as the ratio of the sample size to the number of indicators decreases, the unique scores and unique loadings approach 0. This means that both GSCAm and IGSCA will reduce to GSCA at small enough sample sizes.
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

  # Alternative version of Hwang et al. (2017) for Obtaining Unique Scores--------------------------------------------------
  # Author: Claude Fable 5 (Anthropic AI assistant) and revised by Michael S. Truong
  #
  # This alternative version forces U'1 = 0, so that the unique scores have a mean of 0,
  # which makes the Construct scores have a mean of 0 and is needed or else Path estimates
  # may be biased, especially as theta/n grows. By doing a QR decomposition of [1 | Eta_normed], it's ensured that 
  # the n_constructs + 2 columns of the Q matrix are orthogonal to the ones column. Orthogonality to the ones column also
  # means that the mean of each column is equal to 0.
  # 
  # When n_case is small relative to n_constructs and n_indicators, U may be under-identified and parameter bias may be much more apparent.
  # Solutions to this problem are discussed in Hwang et al. (2017), citing Unkel and Trendafilo (2013) and Trendafilo and Unkel (2011). 
  # It should be noted that likelihood-based covariance structured analysis has similar problems, but faces convergence problems earlier than IGSCA.
  #
  qr_eta <- qr(cbind(1, Eta_normed))
  # The first `k` columns of the Q matrix of [1 | Eta_normed] is equal to the rank of [1 | Eta_normed]. 
  # The full Q matrix is N x N. In the case that some of the constructs are collinear with each other, then indexing naively 
  # Will cause an error, particular when n_case is less than (rank_eta + 1). So the below approach is taken instead
  # Additionally, recall that the rank of [1 | Eta_normed] must always be either equal to or less than the number of rows or columns [1 | Eta_normed],
  #  therefore m will never be equal to a negative number
  rank_eta <- qr_eta$rank
  m        <- n_case - rank_eta
  if (m > 0) {
    #   rank_eta + seq_len(m) should equal to `n_constructs + 2: n_case`, unless there is collinearity in [1 | Eta_normed]
    svd_mx <- svd(tcrossprod(x = D, y = qr.qty(qr_eta, Z_normed)[rank_eta + seq_len(m), , drop = FALSE]))
    # rank_eta should equal to 1 + n_constructs (number of columns of [1 | Eta_normed]), unless there is collinearity
    # matrix(0, n_constructs, n_indicators) is zero-padding so that we take the matrix product of the parts of Q that we want.
    U <- qr.qy(qr_eta, rbind(matrix(0, rank_eta, n_indicators), tcrossprod(x = svd_mx$v, y = svd_mx$u)))
  } else if (m == 0) {
    # See trendafilovExploratoryFactorAnalysis2011 regarding how the unique scores becomes 0
    # Because the unique scores are 0, the unique loadings will also become zero --- meaning that GSCA-M converges to GSCA
    U <- matrix(0, n_case, n_indicators)
  } else if (m <= 0) {
    stop2("Unique score computation failed, sample size is likely too small or there is an error in the algorithm. Please report to developers as you should not be able to see this error.")
  }

  # R Optimized version of Hwang et al. (2017) for Obtaining Unique Scores----------------------------------------------------
  # qr_eta <- qr(Eta_normed)
  # svd_mx <- svd(tcrossprod(x = D, y = qr.qty(qr_eta, Z_normed)[(n_constructs + 1):n_case, , drop = FALSE]))
  ## Optimized Version of
  ## QtZ_null <- qr.qty(qr_eta, Z_normed)[(n_constructs + 1):n_case, , drop = FALSE]
  ## svd_mx <- svd(D %*% t(QtZ_null))
  # U <- qr.qy(qr_eta, rbind(matrix(0, n_constructs, n_indicators),  tcrossprod(x= svd_mx$v, y = svd_mx$u)))
  ## Optimized Version of
  ## Utilde <- svd_mx$v %*% t(svd_mx$u)  # (N-P) × J
  ## U <- qr.qy(qr_eta, rbind(matrix(0, n_constructs, n_indicators), svd_mx$v %*% t(svd_mx$u)))

  # Old method based on Hwang et al. (2017) for Obtaining Unique Scores — O(N^2) memory and computation-----------------------
  # Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
  #   (n_constructs + 1):n_case,
  #   drop = FALSE
  # ]
  # svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
  # Utilde <- svd_mx$v %*% t(svd_mx$u)
  # U <- Eta_Q2 %*% Utilde
  
  # U[, .indicator_type == "Composite"] <- 0

  # Update Unique Loadings------------------------------------------------------------------------------------------------------

  # D <- diag(diag(t(U) %*% Z_normed))
  D <- diag(diag(crossprod(U, Z_normed)))
  D[.indicator_type == "Composite", .indicator_type == "Composite"] <- 0

  # Return output
  return(list("U" = U, "D" = D))
} 

#' Constructs the different matrices of the GSCA-type Least Squares Objective function
#'
#'
#' This function should be called through `calculateFIT`, `calculateFIT_s`, etc... and is not intended to be user-facing.
#' 
#' @inheritParams csem_arguments
#'
#' @return List of Psi, Z, Eta, A, Lambda, B, S and UD matrices 
constructGSCAObjectiveParts <- function(.object) {
  
  # For multigroup models, block diagonalize each part (Z, Gamma, C, ...) so
  # that the single-group algebra below yields one overall FIT for the model.
  HasStructural <- TRUE
  isMULTIGROUP <- FALSE
  if (inherits(.object, "cSEMResults_multi")) {
    if (inherits(.object, "cSEMResults_2ndorder")) {
      stop2(
        "The following error occured in the calculateFIT() function:\n",
        "FIT statistics are not supported for second-order multigroup models."
      )
    }
    if (.object[[1]]$Information$Arguments$.approach_weights != "GSCA") {
      return(NA)
    }
    if (all(.object[[1]]$Information$Model$structural == 0)) {
      HasStructural <- FALSE
    } 
    .object <- bdiagGSCA(.object)
    isMULTIGROUP <- TRUE
  } else if (.object$Information$Arguments$.approach_weights != "GSCA") {
    return(NA)
  } else if (all(.object$Information$Model$structural == 0)) {
    HasStructural <- FALSE
  } 

  # As shown in Equation 4 and 6 the GSCA_m publication (Hwang et al., 2017)
  Eta <- .object$Estimates$Construct_scores
  Z <- .object$Information$Data
  Psi <- cbind(Z, Eta)
  Lambda <- .object$Estimates$Loading_estimates
  
  if (isTRUE(HasStructural)) {
    # Eta[, 1, drop = FALSE] %*% Lambda[1, , drop = FALSE]
    # Eta[, 1, drop = FALSE] %*% Lambda[1, 1, drop = FALSE] 

    # Eta[,1,drop = FALSE] %*% t(B)[1,2, drop = FALSE]
    # Eta[,3,drop = FALSE] %*% t(B)[3,4, drop = FALSE]

    # If there's a structural model
    # Path_estimates is row = to; column = from
    # t(Path_estimates) is row = from; column = to
    # Eta[,1, drop = FALSE] * t(B)[1,2]
    B <- .object$Estimates$Path_estimates
    A <- cbind(
      Lambda,
      t(B)
    )
  }
  else {
    # The rows of the Loading estimates correspond to the construct names
    B <- matrix(
        data = 0,
        nrow = nrow(Lambda),
        ncol = nrow(Lambda),
        dimnames = list(rownames(Lambda), rownames(Lambda))
      )
    A <- cbind(
      Lambda,
      B
    )
  }
  
  if (!all(is.na(.object$Estimates$Unique_scores))) {
    if (isTRUE(isMULTIGROUP)) {
      D <- .object$Estimates$Unique_loading_estimates
    } else {
      D <- diag(.object$Estimates$Unique_loading_estimates)
      rownames(D) <- names(.object$Estimates$Unique_loading_estimates)
      colnames(D) <- names(.object$Estimates$Unique_loading_estimates)
    }

    UD <- .object$Estimates$Unique_scores %*% D
    S <- cbind(
      UD,
      matrix(
        data = 0,
        nrow = nrow(Eta),
        ncol = ncol(Eta),
        dimnames = list(rownames(Eta), colnames(Eta))
      )
    )
  } else if (all(is.na(.object$Estimates$Unique_scores))) {
    # Unique_scores should be NULL when GSCA and not GSCA_m/I-GSCA is run
    UD <- matrix(
        data = 0,
        nrow = nrow(Z),
        ncol = ncol(Z),
        dimnames = list(rownames(Z), colnames(Z))
      )
    S <- cbind(
      UD,
      matrix(
        data = 0,
        nrow = nrow(Eta),
        ncol = ncol(Eta),
        dimnames = list(rownames(Eta), colnames(Eta))
      )
    )
  }

  return(list(Psi = Psi, Z = Z, Eta = Eta, Lambda = Lambda, B = B, A = A, S = S, UD = UD))

}

#' Calculate the matrix of indicator and construct errors for GSCA-type models
#'
#' Computes the residual matrix of the GSCA-type least squares objective
#' function, i.e. `Psi - (Eta %*% A) - S`, based on the output of
#' `constructGSCAObjectiveParts`. The columns corresponding to the indicators
#' contain the indicator errors (`Z - Eta %*% Lambda - UD`) and the columns
#' corresponding to the constructs contain the construct errors
#' (`Eta - Eta %*% t(B)`). For multigroup models the rows of all groups are
#' stacked, with the parameter matrices block diagonalized by `bdiagGSCA`.
#'
#' This function is not intended to be user-facing.
#'
#' @inheritParams csem_arguments
#'
#' @return A numeric matrix with one row per observation and one column per
#'   indicator and construct, or `NA` if `.object` was not estimated with
#'   `.approach_weights = "GSCA"`.
#' @keywords internal
calculateGSCAErrors <- function(.object) {

  parts <- constructGSCAObjectiveParts(.object)

  # constructGSCAObjectiveParts() returns NA for non-GSCA objects
  if (!is.list(parts)) {
    return(NA)
  }

  parts$Psi - (parts$Eta %*% parts$A) - parts$S
}

#' Calculate the degrees of freedom for GSCA-type models
#'
#' Computes the number of model parameters (`npar`), the total number of
#' degrees of freedom (`d_0`) and the residual degrees of freedom
#' (`d_1 = d_0 - npar`) of a GSCA-type model. Following
#' \insertCite{Hwang2014;textual}{cSEM} and the `gesca.mg` function of the
#' gesca package, the number of path, loading, weight and unique loading
#' estimates is multiplied by the number of groups for multigroup models.
#'
#' This function is not intended to be user-facing.
#'
#' @inheritParams csem_arguments
#'
#' @return A named numeric vector of length 3 with elements `npar`, `d_0` and
#'   `d_1`.
#' @keywords internal
calculateGSCADegreesOfFreedom <- function(.object) {

  if (inherits(.object, "cSEMResults_multi")) {
    model <- .object[[1]]$Information$Model
    n_total <- sum(vapply(.object, function(group) {
      nrow(group$Information$Data)
    }, numeric(1)))
    n_groups <- length(.object)
  } else {
    model <- .object$Information$Model
    n_total <- nrow(.object$Information$Data)
    n_groups <- 1
  }

  d_0 <- n_total * length(model$indicators)

  # Hwang, De Sarbo and Takane (2014); combined with the gesca.mg function from
  # the gesca package (June 18/2026) make it clear that the number of
  # parameters is multiplied by the number of groups. This makes sense, or else
  # there would be no penalty for fitting a multigroup model over a
  # single-group one.
  nPath <- sum(model$structural)
  nLoadings <- sum(model$measurement)
  nWeights <- nLoadings
  nUniqueLoadings <- sum(model$measurement[model$construct_type == "Common factor", ])

  npar <- (nPath + nLoadings + nWeights + nUniqueLoadings) * n_groups
  d_1 <- d_0 - npar

  c(npar = npar, d_0 = d_0, d_1 = d_1)
}

#' Block Diagonalize GSCA Parameter Estimates and Scores
#'
#' Block diagonalizes the estimated paramater matrices as shown on Equations
#' 3.28-3.29 on page 111 of \insertCite{Hwang2014;textual}{cSEM}. This is used to compute a single FIT-type statistic given a multigroup GSCA-type model.
#'
#' This function should be called through `constructGSCAObjectiveParts` and is not intended to be user-facing.
#' 
#' @inheritParams csem_arguments
#'
#' @return Single group cSEMResults of Block diagonalized path estimates, loadingestimates, weight estaimes, construct scores and unique scores for GSCA type models.
#' @importFrom Matrix bdiag
bdiagGSCA <- function(.object) {
  
  matrices_to_be_extracted <- list(
    Path_estimates = "Path_estimates",
    Loading_estimates = "Loading_estimates",
    Weight_estimates = "Weight_estimates",
    Construct_scores = "Construct_scores",
    Unique_scores = "Unique_scores",
    Unique_loading_estimates = "Unique_loading_estimates",
    Data = "Data"
  )

  bdiaged_list <- lapply(matrices_to_be_extracted, function(estName) {
    if (estName != "Data") {
      extraction <- lapply(.object, function(group) {
        group[["Estimates"]][[estName]]
      })

      if (estName == "Unique_loading_estimates") {
        extraction <- lapply(extraction, function(x) {
          mx <- diag(x)
          rownames(mx) <- names(x)
          colnames(mx) <- names(x)
          return(mx)
        })
      }
    } else if (estName == "Data") {
      extraction <- lapply(.object, function(group) {
        group[["Information"]][[estName]]
      })
    }

    if (all(unlist(lapply(extraction, is.null)))) {
      return(base::as.matrix(NA))
    }

    extraction_rownames <- lapply(extraction, rownames) |>
      unlist()
    extraction_colnames <- lapply(extraction, colnames) |>
      unlist()
    bdiaged <- Matrix::bdiag(extraction)
    rownames(bdiaged) <- extraction_rownames
    colnames(bdiaged) <- extraction_colnames

    return(base::as.matrix(bdiaged))
  })

  ## Create Surrogate Output -----------------------------------------------
  out <- list()

  out[["Estimates"]] <- bdiaged_list[names(bdiaged_list) != "Data"]

  out[["Information"]] <- bdiaged_list["Data"]

  return(out)
}