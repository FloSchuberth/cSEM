# =============================================================================
# Common setup for the H-O composite + common factor PoC.
#
# Defines:
#   * dgp_syntax(): cSEM.DGP / lavaan model syntax for the population DGP
#   * true_params(): closed-form population parameters that the Stan model
#                    should recover (sign-normalized weights, β, λ, ψ, Σ_x)
#   * simulate_one(): draws one synthetic sample of size N from the DGP
#
# We work with a 3+3 indicator design (K=3 composite indicators, M=3
# reflective indicators), with the population standardized so all latent
# variables have unit variance.
# =============================================================================

suppressPackageStartupMessages({
  library(cSEM.DGP)
  library(cSEM)
  library(cmdstanr)
})

# -----------------------------------------------------------------------------
# Population DGP -- canonical composite c affecting common factor η.
# Indicator block sizes: K = 3, M = 3.
# -----------------------------------------------------------------------------
dgp_syntax <- function() {
"
# Structural model
eta ~ 0.5 * c

# Composite measurement model (unscaled weights in lavaan composite syntax)
c   <~ 0.6 * x1 + 0.7 * x2 + 0.8 * x3

# Reflective measurement model
eta =~ 0.8 * y1 + 0.7 * y2 + 0.9 * y3

# Within-block indicator correlations (composite indicators are arbitrarily
# correlated in a composite; we specify a moderately-correlated block)
x1 ~~ 0.4 * x2
x1 ~~ 0.5 * x3
x2 ~~ 0.6 * x3
"
}

# -----------------------------------------------------------------------------
# Closed-form population parameters that we will compare posterior draws to.
# cSEM.DGP normalizes data and indicators to unit variance, and normalizes
# the composite weights so Var(c) = 1.
# -----------------------------------------------------------------------------
true_params <- function() {
  # Indicator correlation matrix among composite indicators
  Sigma_x <- matrix(c(1.0, 0.4, 0.5,
                      0.4, 1.0, 0.6,
                      0.5, 0.6, 1.0), 3, 3)
  # Unscaled weights from the lavaan syntax
  w_unscaled <- c(0.6, 0.7, 0.8)
  # Scale so that w' Sigma_x w = 1 (cSEM.DGP convention)
  s <- as.numeric(sqrt(t(w_unscaled) %*% Sigma_x %*% w_unscaled))
  w <- w_unscaled / s

  beta     <- 0.5
  sigma_eta <- sqrt(1 - beta^2)  # Var(eta) = 1 in population
  lambda   <- c(0.8, 0.7, 0.9)
  sigma_y  <- sqrt(1 - lambda^2)

  list(w = w, beta = beta, sigma_eta = sigma_eta,
       lambda = lambda, sigma_y = sigma_y, Sigma_x = Sigma_x)
}

# -----------------------------------------------------------------------------
# Draw one synthetic sample using cSEM.DGP::generateData()
# -----------------------------------------------------------------------------
simulate_one <- function(N = 300, seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  dat <- cSEM.DGP::generateData(
    .model       = dgp_syntax(),
    .N           = N,
    .return_type = "matrix"
  )
  # cSEM.DGP returns indicators in alphabetical order; coerce columns
  dat <- as.matrix(dat)
  list(
    N = N, K = 3, M = 3,
    x = dat[, c("x1", "x2", "x3"), drop = FALSE],
    y = dat[, c("y1", "y2", "y3"), drop = FALSE]
  )
}

# -----------------------------------------------------------------------------
# Sign-fix posterior draws of w (and induced flip of beta) so they are
# comparable across replications. We follow the convention that w[1] > 0.
# -----------------------------------------------------------------------------
sign_fix_draws <- function(draws_df, K = 3) {
  flip <- draws_df[[paste0("w[1]")]] < 0
  for (k in seq_len(K)) {
    col <- paste0("w[", k, "]")
    draws_df[[col]][flip] <- -draws_df[[col]][flip]
  }
  # beta also flips because c -> -c sends β -> -β
  draws_df[["beta"]][flip] <- -draws_df[["beta"]][flip]
  draws_df
}

# -----------------------------------------------------------------------------
# Build Stan initial values from cSEM PLS estimates on the same data.
#
# Rationale: random inits sometimes leave chains stuck in a degenerate
# corner (λ_1 ≈ 0, β ≈ 0) from which they cannot escape. PLS gives a fast,
# consistent estimate of (w, β, λ) that is well inside the regular region.
# Using it as the starting point is the standard practical workflow for
# composite SEMs and does NOT compromise inference (the prior + likelihood
# decide the posterior; inits are just where the sampler starts).
# -----------------------------------------------------------------------------
init_from_pls <- function(sim, jitter_sd = 0.05, n_chains = 4) {
  dat_df <- as.data.frame(cbind(sim$x, sim$y))
  pls    <- cSEM::csem(dat_df, dgp_syntax())
  s      <- cSEM::summarize(pls)
  w_hat  <- s$Estimates$Weight_estimates$Estimate[s$Estimates$Weight_estimates$Construct_type == "Composite"]
  if (length(w_hat) != sim$K) w_hat <- s$Estimates$Weight_estimates$Estimate[1:sim$K]
  # Loadings: rows whose Name begins with "eta =~"
  lam_idx <- grep("^eta =~", s$Estimates$Loading_estimates$Name)
  lam_hat <- s$Estimates$Loading_estimates$Estimate[lam_idx]
  if (length(lam_hat) != sim$M) lam_hat <- rep(0.7, sim$M)
  beta_hat <- s$Estimates$Path_estimates$Estimate[1]
  # Ensure sign conventions: w[1] > 0 and λ[1] > 0; flip companion signs too
  if (w_hat[1] < 0)   { w_hat <- -w_hat; beta_hat <- -beta_hat }
  if (lam_hat[1] < 0) { lam_hat <- -lam_hat; beta_hat <- -beta_hat }
  Sigma_x_hat <- cor(sim$x)
  L_corr_hat  <- t(chol(Sigma_x_hat))
  sigma_y_hat <- sqrt(pmax(1 - lam_hat^2, 0.05))
  beta_hat    <- pmin(pmax(beta_hat, -0.95), 0.95)

  lapply(seq_len(n_chains), function(.c) {
    j <- rnorm(1, 0, jitter_sd)  # consume RNG so chains differ
    list(
      w_raw_head  = max(0.05, w_hat[1] + abs(rnorm(1, 0, jitter_sd))),
      w_raw_tail  = w_hat[-1] + rnorm(sim$K - 1, 0, jitter_sd),
      L_corr      = L_corr_hat,
      lambda_head = max(0.05, lam_hat[1] + abs(rnorm(1, 0, jitter_sd))),
      lambda_tail = lam_hat[-1] + rnorm(sim$M - 1, 0, jitter_sd),
      sigma_y     = sigma_y_hat + abs(rnorm(sim$M, 0, jitter_sd)),
      beta        = beta_hat + rnorm(1, 0, jitter_sd)
    )
  })
}
