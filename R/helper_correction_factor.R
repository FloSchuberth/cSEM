#' Calculate the correction factor
#'
#' Calculates the correction factor of PLSc.
#'
#' Currently seven methods are available:
#'
#' \itemize{
#' \item "dist_euclid" (default)
#' \item "dist_euclid_weighted"
#' \item "fisher_transformed"
#' \item "mean_geometric"
#' \item "mean_harmonic"
#' \item "mean_arithmetic"
#' \item "geo_of_harmonic" (not yet implemented)
#' }
#'
#' @inheritParams csem_arguments
#'
#' @return A numeric vector of correction factors with element names equal
#'   to the names of the J constructs used in the measurement model
#'

calculateCorrectionFactors <- function(
  .S            = NULL,
  .W            = NULL,
  .csem_model   = NULL,
  .approach_cf  = NULL
  ) {

  approach_cf <- .approach_cf

  ### Compute correction factors  ----------------------------------------------
  correction_factors <- vector(mode = "double", length = nrow(.W))
  for(j in 1:nrow(.W)) {

    ## Extract vector of weights of block j
    w_j <- .W[j, ] %>%
      .[. != 0] %>%
      as.matrix(.)

    ## Check if single indicator block or composite; If yes, set cf to 1
    if(nrow(w_j) == 1 | .csem_model$construct_type[j, ]["Type"] == "Composite") {
      correction_factors[j] <- 1

    } else {
      ## Extract relevant objects
      E_jj <- .csem_model$error_cor[rownames(w_j), rownames(w_j)]
      S_jj <- .S[rownames(w_j), rownames(w_j)]
      W_jj <- w_j %*% t(w_j)

      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of S_jj (result is vectorized)
      S_vect <- replace(S_jj, which(E_jj == 1), 0) %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[. != 0]

      ## Set indicator pairs whise measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of W_jj (vectorized)
      W_vect <- replace(W_jj, which(E_jj == 1), 0)  %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[. != 0]

      if(length(S_vect) == 0) {
        stop("At least one pair of indicators with uncorrelated measurement
             errors required.\n Please revise your model.",
             call. = FALSE)
      }

      ## Do the actual computation ---------------------------------------------
      switch (approach_cf,
              "dist_euclid"          = {
                cf <- sum(W_vect * S_vect) / sum(W_vect^2)
              },
              "dist_euclid_weighted" = {
                weights <- 1 / (1 - S_vect^2)
                cf      <- sum(W_vect * S_vect * weights) / sum(W_vect^2 * weights)
              },
              "fisher_transformed"   = {
                # Function to be minimized
                temp_fun <- function(.c, .W_vect, .S_vect){
                  sum((0.5*log((1 + .S_vect) / (1 - .S_vect)) -
                         0.5*log((1 + .c*.W_vect) / (1 - .c*W_vect)))^2)
                }

                # Optimaziation
                temp_optim <- optim(fn = temp_fun, par = 0.5, method = "BFGS",
                      .W_vect = W_vect, .S_vect = S_vect)
                cf <- temp_optim$par
              },
              "mean_geometric"       = {
                cf <- prod(S_vect/W_vect)^(1/length(S_vect))
              },
              "mean_arithmetic"      = {
                cf <- mean(S_vect/W_vect)
              },
              "mean_harmonic"        = {
                cf <- 1/mean(1/(S_vect/W_vect))
              },
              "geo_of_harmonic"      = {stop("not implemented yet")}
      )

      ## Compute absolute value and take the sqrt since cf = c^2
      correction_factors[j] <- sqrt(abs(cf))
    }
  }
  names(correction_factors) <- row.names(.W)
  return(correction_factors)
  ### For maintenance: ### -------------------------------------------------------------------------
  # w_j (K_j x 1)      := Column vector of indicator weights of block j
  # W_vect (1 x 2*K_j) := Vector of off-diagonal elements of w_jw'_j
  # S_vect (1 x 2*K_j) := Vector of off-diagonal elements of S_jj
}
