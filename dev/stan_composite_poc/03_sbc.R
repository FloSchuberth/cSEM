# =============================================================================
# Simulation-Based Calibration (SBC) -- minimal version.
#
# SBC tests posterior calibration relative to the *prior+model*, not relative
# to a particular DGP. The recipe:
#   for each of R replications:
#     θ ~ prior; generate data y ~ p(·|θ); fit; compute rank of θ among
#     L thinned posterior draws.
#   ranks should be uniform on {0, 1, ..., L} if posterior is correctly
#   computed.
#
# CAVEATS for this model:
#   * The prior on the composite indicator correlation matrix Σ_x (LKJ(2))
#     and on w_raw (N(0,1)) jointly induce a non-trivial prior on the
#     weights w. We sample from this prior here directly to make ranks
#     interpretable.
#   * For SBC to be meaningful the model must be identified up to the sign
#     conventions we have already imposed (w[1] > 0, λ[1] > 0, |β| < 1).
#   * This is a coarse SBC. A definitive analysis would use many more
#     replications (e.g., 500+) and pay careful attention to the chain
#     thinning required to obtain i.i.d. posterior draws.
# =============================================================================
source("dev/stan_composite_poc/00_setup.R")

R_sbc   <- as.integer(Sys.getenv("R_SBC", 50))
L_draws <- as.integer(Sys.getenv("L_DRAWS", 50))
N_each  <- as.integer(Sys.getenv("N_SBC", 300))

mod <- cmdstan_model("dev/stan_composite_poc/composite_ho.stan")

# ---- Helper: simulate data from prior draw (uses cSEM.DGP's syntax route) ---
# To avoid coding the H-O DGP manually, we draw plausible parameter values
# from prior-like distributions and render them via cSEM.DGP. We constrain
# everything so that all relevant moments stay in the feasible region.
sim_from_prior <- function(N) {
  ok <- FALSE
  while (!ok) {
    w_unscaled <- abs(rnorm(3, 0.7, 0.15))           # positive weights
    r_pairs    <- runif(3, 0.2, 0.7)                 # indicator correlations
    Sigma_x <- matrix(c(1, r_pairs[1], r_pairs[2],
                        r_pairs[1], 1, r_pairs[3],
                        r_pairs[2], r_pairs[3], 1), 3, 3)
    beta_true   <- runif(1, -0.85, 0.85)
    lambda_true <- runif(3, 0.55, 0.95)              # all positive, identifiable
    if (any(lambda_true >= 0.99)) next
    if (min(eigen(Sigma_x, only.values = TRUE)$values) <= 0.05) next
    ok <- TRUE
  }
  syntax <- sprintf(
    "eta ~ %.4f * c\nc <~ %.4f * x1 + %.4f * x2 + %.4f * x3\n
     eta =~ %.4f * y1 + %.4f * y2 + %.4f * y3\n
     x1 ~~ %.4f * x2\n x1 ~~ %.4f * x3\n x2 ~~ %.4f * x3\n",
    beta_true,
    w_unscaled[1], w_unscaled[2], w_unscaled[3],
    lambda_true[1], lambda_true[2], lambda_true[3],
    r_pairs[1], r_pairs[2], r_pairs[3]
  )
  dat <- try(cSEM.DGP::generateData(.model = syntax, .N = N, .return_type = "matrix"),
             silent = TRUE)
  if (inherits(dat, "try-error")) return(NULL)
  # scale weights to unit-variance composite
  s <- as.numeric(sqrt(t(w_unscaled) %*% Sigma_x %*% w_unscaled))
  w_true <- w_unscaled / s
  list(
    sim = list(N = N, K = 3, M = 3,
               x = as.matrix(dat[, c("x1","x2","x3"), drop = FALSE]),
               y = as.matrix(dat[, c("y1","y2","y3"), drop = FALSE])),
    truth = c(w1 = w_true[1], w2 = w_true[2], w3 = w_true[3],
              beta = beta_true,
              lambda1 = lambda_true[1], lambda2 = lambda_true[2],
              lambda3 = lambda_true[3])
  )
}

monitored <- c("w[1]", "w[2]", "w[3]", "beta",
               "lambda[1]", "lambda[2]", "lambda[3]")

ranks <- matrix(NA_integer_, nrow = R_sbc, ncol = length(monitored))
colnames(ranks) <- monitored

t0 <- Sys.time()
for (r in seq_len(R_sbc)) {
  d <- sim_from_prior(N_each)
  if (is.null(d)) next
  inits <- try(init_from_pls(d$sim, n_chains = 2), silent = TRUE)
  if (inherits(inits, "try-error")) next

  fit <- try(mod$sample(
    data = list(N = d$sim$N, K = d$sim$K, M = d$sim$M, x = d$sim$x, y = d$sim$y),
    chains = 2, parallel_chains = 2,
    iter_warmup = 800, iter_sampling = ceiling(L_draws * 5),
    thin = 5, refresh = 0, seed = 7000 + r, adapt_delta = 0.95,
    init = inits, show_messages = FALSE, show_exceptions = FALSE
  ), silent = TRUE)
  if (inherits(fit, "try-error")) next

  draws <- posterior::as_draws_df(fit$draws(variables = monitored))
  draws <- sign_fix_draws(draws, K = 3)
  if (nrow(draws) < L_draws) next
  draws_sub <- draws[seq_len(L_draws), monitored, drop = FALSE]
  truth_named <- c(d$truth["w1"], d$truth["w2"], d$truth["w3"],
                   d$truth["beta"], d$truth["lambda1"], d$truth["lambda2"], d$truth["lambda3"])
  for (j in seq_along(monitored)) {
    ranks[r, j] <- sum(draws_sub[[monitored[j]]] < truth_named[j])
  }
  cat(sprintf("[sbc %3d/%d] ranks: %s   t=%.0fs\n", r, R_sbc,
              paste(ranks[r, ], collapse = ","),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

ok_rows <- complete.cases(ranks)
ranks_ok <- ranks[ok_rows, , drop = FALSE]
cat(sprintf("\n%d / %d valid replications.\n", nrow(ranks_ok), R_sbc))

# Crude check: chi-squared test for uniformity in B equal bins
B <- 5
chisq_p <- numeric(length(monitored)); names(chisq_p) <- monitored
for (j in seq_along(monitored)) {
  br <- cut(ranks_ok[, j], breaks = seq(-0.5, L_draws + 0.5, length.out = B + 1),
            include.lowest = TRUE)
  tab <- table(br)
  chisq_p[j] <- suppressWarnings(chisq.test(tab)$p.value)
}
cat("\nChi-squared p-values for uniform ranks (per parameter):\n")
print(round(chisq_p, 3))

saveRDS(list(ranks = ranks_ok, R_sbc = R_sbc, L = L_draws, N = N_each, chisq_p = chisq_p),
        "dev/stan_composite_poc/sbc_result.rds")
cat("\nSaved dev/stan_composite_poc/sbc_result.rds\n")
