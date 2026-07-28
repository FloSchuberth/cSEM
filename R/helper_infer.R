#' Internal: Helper for infer()
#' 
#' Collection of various functions that compute an inferential quantity. 
#' 
#' Implementation and terminology of the confidence intervals is based on 
#' \insertCite{Hesterberg2015;textual}{cSEM} and 
#' \insertCite{Davison1997;textual}{cSEM}.
#' 
#' @inheritParams csem_arguments
#'  
#' @references
#'   \insertAllCited{} 
#'   
#' @name inference_helper
#' @rdname inference_helper
#' @keywords internal
MeanResample <- function(.first_resample) {
  
  lapply(.first_resample, function(x) {
    out        <- colMeans(x$Resampled)
    names(out) <- names(x$Original)
    out
  })
}
#' @rdname inference_helper
SdResample <- function(.first_resample, .resample_method, .n) {
  
  lapply(.first_resample, function(x) {
    out        <- matrixStats::colSds(x$Resampled)
    names(out) <- names(x$Original)
    if(.resample_method == "jackknife") {
      out <- out * (.n - 1)/sqrt(.n)
    }
    out
  })
}
#' @rdname inference_helper
BiasResample <- function(.first_resample, .resample_method, .n) {
  
  lapply(.first_resample, function(x) {
    out        <- colMeans(x$Resampled) - x$Original
    names(out) <- names(x$Original)
    if(.resample_method == "jackknife") {
      out <- out * (.n - 1)
    }
    out
  })
}
#' @rdname inference_helper
# Computes the *Standard CI with bootstrap SE's*.
# Critical quantiles can be based on both the `t`- or the 
# standard normal distribution (`z`). The former may perform better in
# small samples but there is no clear consenus on what the degrees of freedom
# should be. We use N - 1 ("type1").
StandardCIResample <- function(
  .first_resample, 
  .bias_corrected,
  .dist = c("z", "t"), 
  .df = c("type1", "type2"),
  .resample_method, 
  .n,
  .probs
) {
  # Standard CI with bootstrap SE's
  # Critical quantiles can be based on both the t- or the 
  # standard normal distribution (z). The former may perform better in
  # small samples but there is no clear consenus on what the degrees of freedom
  # should be.
  # CI: [theta_hat + c(alpha/2)*boot_sd ; theta_hat + c(1 - alpha/2)*boot_sd]
  #   = [theta_hat - c(1 - alpha/2)*boot_sd ; theta_hat + c(1 - alpha/2)*boot_sd]
  # if c() is based on a symmetric distribution.
  
  ## Compute standard deviation
  boot_sd <- SdResample(.first_resample, .resample_method = .resample_method, .n = .n)
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  ## Compute intervals
  # Notation:
  # w = .first_resample := the list containing the estimated values based on .R resamples
  #                        and the original values
  # y      := the estimated standard errors (based on .R resamples)
  # z      := a vector of probabilities
  out <- mapply(function(w, y) {
    lapply(.probs, function(z) {
      
      theta_star <- if(.bias_corrected) {
        2*w$Original - colMeans(w$Resampled) 
        # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
      } else {
        w$Original
      }
      
      if(.dist == "t") {
        df <- switch(.df,
                     "type1" = .n - 1,
                     "type2" = nrow(w) - 1
        )
        theta_star + qt(z, df = df) * y
      } else {
        theta_star + qnorm(z) * y
      }
    })
  }, w = .first_resample, y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      # rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      rownames(x) <- unlist(lapply(100*cl, function(x) 
        c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
      x
    })
  return(out)
}
#' @rdname inference_helper
#The function takes the distribution F* (the CDF) of the resamples as an estimator for
#the true distribution F of the statistic/estimator of interest. 
#Quantiles of the estimated distribution are then used as lower and upper bound.
PercentilCIResample <- function(.first_resample, .probs) {
  # Percentile CI 
  # Take the bootstrap distribution F* (the CDF) as an estimator for
  # the true distribution F. Use the quantiles of the estimated distribution
  # to estimate the CI for theta.
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  lapply(.first_resample, function(x) {
    out <- t(matrixStats::colQuantiles(x$Resampled, probs = .probs, drop = FALSE))
    colnames(out) <- names(x$Original)
    # rownames(out) <- paste0(c("L_", "U_"), rep(cl, each = 2))
    rownames(out) <- unlist(lapply(100*cl, function(x) 
      c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
    out
  })
}
#' @rdname inference_helper
BasicCIResample <- function(.first_resample, .bias_corrected, .probs) {
  # Basic CI 
  # Estimate the distribution of delta_hat = theta_hat - theta by the bootstrap 
  # distribution delta_hat_star = theta_hat_star - theta_hat.
  # Since P(a <= theta_hat - theta <= b) = P(theta_hat - b <= theta <= theta_hat - a)
  # (notice how the limits a and b switch places!!) 
  # Define q_p := p%-quantile of the bootstrap distribution of delta_hat_star
  # and
  #  b = q_(1 - alpha/2) 
  #  a = q_(alpha/2) 
  # then the CI is:
  # CI: [theta_hat - q_(1 - alpha/2); theta_hat - q_(alpha/2)]
  # Note that q_p is just the alpha percentile quantile shifted by theta_hat:
  # q_p = Q_p - theta_hat, where Q_P = F*^-1(p) := the p% quantile of the
  # distribution of theta_hat_star.
  # Therefore:
  # CI: [2*theta_hat - Q_(1 - alpha/2); 2*theta_hat - Q_(alpha/2)]
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  lapply(.first_resample, function(x) {
    out <- t(matrixStats::colQuantiles(x$Resampled, probs = .probs, drop = FALSE))
    
    theta_star <- if(.bias_corrected) {
      2*x$Original - colMeans(x$Resampled) 
      # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
    } else {
      x$Original
    }
    
    out <- t(2*theta_star - t(out))
    colnames(out) <- names(x$Original)
    out <- out[1:nrow(out) + rep(c(1, -1), times = nrow(out)/2), , drop = FALSE]
    # rownames(out) <- paste0(c("L_", "U_"), rep(cl, each = 2))
    rownames(out) <- unlist(lapply(100*cl, function(x) 
      c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
    out
  })
}
#' @rdname inference_helper
# The function computes a boostrap t-statisic (since it is roughly pivotal) 
# and constructs the CI based on bootstraped t-values and bootstraped/jackknife
# SE's
TStatCIResample <- function(
  .first_resample, 
  .second_resample, 
  .bias_corrected,
  .resample_method, 
  .resample_method2, 
  .n, 
  .probs
) {
  # Bootstraped t statistic
  # Bootstrap the t-statisic (since it is roughly pivotal) and compute
  # CI based on bootstraped t-values and bootstraped/jackknife SE
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  boot_sd_star <- lapply(.second_resample, SdResample, .resample_method = .resample_method2, .n = .n) %>% 
    purrr::transpose(.) %>% 
    lapply(function(y) do.call(rbind, y))
  
  boot_sd <- SdResample(.first_resample, .resample_method = .resample_method, .n = .n)
  
  boot_t_stat <- mapply(function(x, y) t(t(x$Resampled) - x$Original) / y,
                        x = .first_resample,
                        y = boot_sd_star)
  
  i <- 1
  y <- boot_sd[[1]]
  z <- .probs[1]
  out <- mapply(function(i, y) {
    lapply(.probs, function(z) {
      qt_boot <- t(matrixStats::colQuantiles(boot_t_stat[[i]], probs = z, drop = FALSE))
      
      theta_star <- if(.bias_corrected) {
        2*.first_resample[[i]]$Original - colMeans(.first_resample[[i]]$Resampled) # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
      } else {
        .first_resample[[i]]$Original
      }
      ## Confidence interval
      theta_star - qt_boot * y
    })
  }, i = seq_along(.first_resample), y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      x <- x[1:nrow(x) + rep(c(1, -1), times = nrow(x)/2), , drop = FALSE]
      # rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      rownames(x) <- unlist(lapply(100*cl, function(x) 
        c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
      x
    })
  
  names(out) <- names(.first_resample)
  return(out)
}
#' @rdname inference_helper
BcCIResample <- function(.first_resample, .probs) {
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  out <- lapply(.first_resample, function(x) {
    p0 <- colMeans(t(t(x$Resampled) <= x$Original))
    z0 <- qnorm(p0)
    
    bc <- lapply(1:length(z0), function(i) {
      bc_quantiles <- c()
      for(j in seq_along(.probs)) {
        q  <- pnorm(2*z0[i] + qnorm(.probs[j]))
        bc_quantiles <- c(bc_quantiles, quantile(x$Resampled[, i], probs = q))
      }
      bc_quantiles
    })
    names(bc) <- names(x$Original)
    bc
  }) %>% 
    lapply(function(x) t(do.call(rbind, x))) %>% 
    lapply(function(x) {
      # rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      rownames(x) <- unlist(lapply(100*cl, function(x) 
        c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
      x
    })
  
  return(out)
} 
#' @rdname inference_helper
BcaCIResample <- function(.object, .first_resample, .probs) {
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  ## influence values (estimated via jackknife)
  # Delete resamples and cSEMResults_resampled class tag 
  # to be able to run resampling again. 
  # if(any(class(.object) == "cSEMResults_2ndorder")) {
  #   .object$Second_stage$Information$ResamplesEstimates$Estimates_resample <- NULL
  #   .object$Second_stage$Estimates$Information_resample <- NULL
  # } else {
  #   .object$Estimates$Estimates_resample <- NULL
  #   .object$Estimates$Information_resample <- NULL
  # }
  class(.object) <- setdiff(class(.object), "cSEMResults_resampled")
  
  jack <- resamplecSEMResults(
    .object = .object,
    .resample_method = "jackknife"
  )
  
  jack_estimates <- if(any(class(jack) == "cSEMResults_2ndorder")) {
    jack$Second_stage$Information$Resamples$Estimates$Estimates1
  } else {
    jack$Estimates$Estimates_resample$Estimates1
  }
  aFun <- function(x) {
    1/6 * (sum(x^3) / sum(x^2)^1.5)
  }
  
  a <- lapply(jack_estimates, function(x) t(x$Original - t(x$Resampled))) %>% 
    lapply(function(x) apply(x, 2, aFun)) 
  
  p0 <- lapply(.first_resample, function(x) colMeans(t(t(x$Resampled) <= x$Original)))
  z0 <- lapply(p0, qnorm)
  
  out <- lapply(seq_along(.first_resample), function(i) {
    bca <- lapply(1:length(z0[[i]]), function(j) {
      q <- pnorm(z0[[i]][j] + (z0[[i]][j] + qnorm(.probs))/ (1 - a[[i]][j]*(z0[[i]][j] + qnorm(.probs)))) 
      quantile(.first_resample[[i]]$Resampled[, j], probs = q)
    })
    names(bca) <- names(.first_resample[[i]]$Original)
    bca
  })%>% 
    lapply(function(x) t(do.call(rbind, x))) %>% 
    lapply(function(x) {
      # rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      rownames(x) <- unlist(lapply(100*cl, function(x) 
        c(sprintf("%.6g%%L", x), sprintf("%.6g%%U", x))))
      x
    })
  
  names(out) <- names(.first_resample)
  out
}

# =============================================================================
# Asymptotic (analytic) inference helpers
# =============================================================================

#' Distribution-free asymptotic covariance matrix of correlations
#'
#' Large-sample variance-covariance matrix of the unique sample correlations in
#' `.S`, using the estimator of Steiger & Hakstian (1982). 
#'
#' The correlations are taken in `.S[lower.tri(.S)]` (column-major) order - the
#' same order used by the HTMT gradient - so that a gradient vector and this
#' matrix align by position for a delta-method confidence interval, i.e.
#' `Var(f) = t(g) %*% calculateCorVCV(.S, .data) %*% g`.
#'
#' @param .data A `(n x P)` matrix or data.frame of indicator data whose columns
#'   include (and are matched by name to) `rownames(.S)`.
#'
#' @return The `(L x L)` asymptotic covariance matrix of the correlations,
#'   `L = P(P - 1)/2`, in `lower.tri(.S)` order.
#'
#' @references
#'   Steiger, J. H., & Hakstian, A. R. (1982). The asymptotic distribution of
#'   elements of a correlation matrix: Theory and application. British Journal of
#'   Mathematical and Statistical Psychology, 35(2), 208-215.
#'
#' @keywords internal
calculateCorVCV <- function(.data) {
  n <- nrow(.data)
  p <- ncol(.data)
  size_n <- (p * p - p) / 2
  
  # Pre-calculate means and centered data
  data_means <- colMeans(.data)
  data_sd <- apply(.data, 2, sd)
  
  # Create indices for upper triangle
  indices <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  
  # Initialize result matrix
  vc_r <- matrix(0, nrow = size_n, ncol = size_n)
  
  # Pre-calculate all possible products of centered variables
  products <- array(0, dim = c(n, p, p))
  for(i in 1:p) {
    for(j in i:p) {
      products[, i, j] <- .data[, i] * .data[, j]
      products[, j, i] <- products[, i, j]
    }
  }
  
  # Function to get position in upper triangular matrix
  get_pos <- function(i, j) {
    if(i > j) {
      temp <- i
      i <- j
      j <- temp
    }
    return(p*(i-1) - i*(i-1)/2 + j - i)
  }
  
  # Calculate covariances using vectorized operations
  for(idx1 in 1:nrow(indices)) {
    x <- indices[idx1, 1]
    y <- indices[idx1, 2]
    
    for(idx2 in idx1:nrow(indices)) {
      z <- indices[idx2, 1]
      t <- indices[idx2, 2]
      
      sd <- sd(.data[,x]) * sd(.data[,y]) * sd(.data[,z]) * sd(.data[,t])
      
      
      # Calculate fourth-order moments using pre-computed products
      mu_xyzt <- mean(products[, x, y] * products[, z, t])  / sd
      
      mu_xxzt <- mean(products[, x, x] * products[, z, t]) / sd
      mu_yyzt <- mean(products[, y, y] * products[, z, t]) / sd
      mu_xyzz <- mean(products[, x, y] * products[, z, z]) / sd
      mu_xytt <- mean(products[, x, y] * products[, t, t]) / sd
      
      mu_xxtt <- mean(products[, x, x] * products[, t, t]) / sd
      mu_xxzz <- mean(products[, x, x] * products[, z, z]) / sd
      mu_yytt <- mean(products[, y, y] * products[, t, t]) / sd
      mu_yyzz <- mean(products[, y, y] * products[, z, z]) / sd
      
      # Get correlation coefficients
      rxy <- cor(.data[, x], .data[, y])
      rzt <- cor(.data[, z], .data[, t])
      
      # Calculate covariance
      cov_val <- (mu_xyzt - 0.5 * rxy * (mu_xxzt + mu_yyzt) - 
                    0.5 * rzt * (mu_xyzz + mu_xytt) + 
                    0.25 * rxy * rzt * (mu_xxzz + mu_xxtt + mu_yyzz + mu_yytt))
      
      pos1 <- get_pos(x, y)
      pos2 <- get_pos(z, t)
      
      vc_r[pos1, pos2] <- cov_val
      vc_r[pos2, pos1] <- cov_val  # Matrix is symmetric
    }
  }
  
  return(vc_r/nrow(.data))
}


#' Delta-method standard errors of the HTMT
#'
#' Combines the per-construct-pair HTMT gradients with the asymptotic covariance
#' matrix of the indicator correlations ([calculateCorVCV()]) into the
#' delta-method standard error of each HTMT, `sqrt(t(g) %*% Sigma %*% g)`. The
#' covariance matrix is formed once and reused for every pair. Gradient vectors
#' and covariance matrix share the `lower.tri(.S)` ordering, so they align by
#' position (both are built from the same `.S`).
#'
#' No lower bound is imposed on the variance: since `Sigma` is positive
#' semi-definite, `t(g) %*% Sigma %*% g` is non-negative by construction, so a
#' negative value would signal a genuine problem and is deliberately allowed to
#' surface as `NaN` (with a warning) rather than being silently clamped to zero.
#'
#' @param .gradients A named list of gradient vectors, one per construct pair,
#'   each aligned with `.S[lower.tri(.S)]`.
#' @param .S A `(P x P)` indicator correlation matrix.
#' @param .data A `(n x P)` matrix or data.frame of indicator data.
#'
#' @return A named numeric vector of standard errors, one per construct pair.
#'
#' @keywords internal
calculateHTMTasymptoticSE <- function(.gradients, .data) {
  Sigma <- calculateCorVCV(.data = .data)          # once, reused for all pairs
  vapply(
    .gradients,
    function(g) sqrt(drop(t(g) %*% Sigma %*% g)),           # sqrt( t(g) %*% Sigma %*% g )
    numeric(1)
  )
}
