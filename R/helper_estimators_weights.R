#' Internal: Calculate the inner weights for PLS-PM
#'
#' PLS-PM forms "inner" composites as a weighted sum of its *I* related composites.
#' These inner weights are obtained using one of the following schemes:
#' \describe{
#'   \item{`centroid`}{According to the centroid scheme each inner weight used
#'     to form composite *j* is either 1 if the correlation between composite *j* and 
#'     its related composite *i = 1, ..., I* is positive 
#'     and -1 if it is negative.}
#'   \item{`factorial`}{According to the factorial scheme each inner weight used
#'     to form composite *j* is equal to the correlation between composite *j* 
#'     and its related composite *i = 1, ..., I*.}
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
#' Calculates outer weights using Mode A, Mode B, Unit or fixed weights.
#'
#' @usage calculateOuterWeightsPLS(
#'    .data   = args_default()$.data  
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
  
  # scale the inner proxy, inner weights are usually scaled that the inner
  # proxies are standardized
  inner_proxy <- scale(.data %*% t(W) %*% t(.E))
  colnames(inner_proxy) = rownames(W)
  
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
      # Select only
      W[block, indicators] <- proxy_indicator_cor[block, indicators]
      
    } else if(.modes[block] == "modeB") {
      ## Mode B - Regression of each proxy on all its indicator
      # W_j = S_jj^{-1} %*% Cov(eta_j, X_j)
      W[block, indicators] <- solve(.S[indicators, indicators]) %*% proxy_indicator_cor[block, indicators]
      
    } else if(.modes[block] == "modeBNNLS"){
      # .data is standardized, i.e., mean 0 and unit variance, inner proxy is also 
      #  standardized (standardization of the inner proxy has no effect)
      temp <- nnls::nnls(A = .data[,indicators,drop=FALSE], b = inner_proxy[,block])
      W[block, indicators] <- temp$x
    } else if(.modes[block] == "PCA"){
      temp <- psych::principal(r = .S[indicators, indicators], nfactors = 1)
      
      W[block, indicators] <- c(temp$weights)
      
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
    }
  )
}

#' Internal: Scale weights
#'
#' Scale weights such that the composite variance is one.
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
  
  ## Scale the weights to ensure that the proxies have a variance of one
  W_scaled <- solve(diag(sqrt(var_proxies))) %*% .W
  
  ## Assign rownames and colnames to the scaled weights and return
  rownames(W_scaled) <- rownames(.W)
  colnames(W_scaled) <- colnames(.W)
  
  return(W_scaled)
}
