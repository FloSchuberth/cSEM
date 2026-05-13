# =============================================================================
# Multi-replicate parameter recovery study.
#
# For R independent samples from the population DGP, fit the Stan model and
# record:
#   * posterior mean of each parameter of interest
#   * 90% credible interval bounds
#
# Aggregate over replications to estimate:
#   * MC bias        (mean of posterior means − truth)
#   * MC RMSE        (sqrt mean squared error of posterior means)
#   * 90% CI coverage (fraction of CIs containing truth)
#
# 90% nominal coverage should be ~0.9 if the posterior is well calibrated
# at the chosen N.  This is a frequentist calibration check, NOT a Bayesian
# coherence check (that would be SBC, see 03_sbc.R).
# =============================================================================
source("dev/stan_composite_poc/00_setup.R")

R_reps <- as.integer(Sys.getenv("REPS", 30))
N_each <- as.integer(Sys.getenv("N", 500))
cat(sprintf("Running %d replications with N=%d each\n", R_reps, N_each))

mod <- cmdstan_model("dev/stan_composite_poc/composite_ho.stan")

pars_of_interest <- c("w[1]", "w[2]", "w[3]",
                      "beta",
                      "lambda[1]", "lambda[2]", "lambda[3]",
                      "sigma_y[1]", "sigma_y[2]", "sigma_y[3]",
                      "R2_eta")
tp <- true_params()
truth <- c(tp$w,
           tp$beta,
           tp$lambda,
           tp$sigma_y,
           tp$beta^2)                  # R^2 since Var(η) = 1
names(truth) <- pars_of_interest

results <- vector("list", R_reps)
diag_log <- data.frame(rep = integer(), max_rhat = double(),
                       min_ess = double(), divergences = integer())

t_start <- Sys.time()
for (r in seq_len(R_reps)) {
  set.seed(1000 + r)
  sim <- simulate_one(N = N_each, seed = 1000 + r)
  inits <- try(init_from_pls(sim, n_chains = 2), silent = TRUE)
  if (inherits(inits, "try-error")) {
    cat(sprintf("[rep %d] PLS init failed, skipping\n", r))
    next
  }

  fit <- tryCatch(
    mod$sample(
      data            = list(N = sim$N, K = sim$K, M = sim$M, x = sim$x, y = sim$y),
      chains          = 2,
      parallel_chains = 2,
      iter_warmup     = 800,
      iter_sampling   = 800,
      refresh         = 0,
      seed            = 2000 + r,
      adapt_delta     = 0.95,
      init            = inits,
      show_messages   = FALSE,
      show_exceptions = FALSE
    ),
    error = function(e) { cat(sprintf("[rep %d] FAIL: %s\n", r, conditionMessage(e))); NULL }
  )
  if (is.null(fit)) next

  draws <- posterior::as_draws_df(fit$draws(variables = pars_of_interest))
  draws <- sign_fix_draws(draws, K = sim$K)

  summ <- fit$summary(variables = pars_of_interest)
  divs <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  diag_log <- rbind(diag_log, data.frame(
    rep = r,
    max_rhat = max(summ$rhat, na.rm = TRUE),
    min_ess  = min(summ$ess_bulk, na.rm = TRUE),
    divergences = divs
  ))

  results[[r]] <- list(
    rep = r,
    post_mean = sapply(pars_of_interest, function(p) mean(draws[[p]])),
    ci_lo     = sapply(pars_of_interest, function(p) quantile(draws[[p]], 0.05)),
    ci_hi     = sapply(pars_of_interest, function(p) quantile(draws[[p]], 0.95))
  )
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  rate <- elapsed / r
  remaining <- rate * (R_reps - r)
  cat(sprintf("[rep %3d/%d] max R-hat=%.3f  min ESS=%6.0f  divs=%d  elapsed=%5.0fs  ETA=%5.0fs\n",
              r, R_reps, max(summ$rhat, na.rm = TRUE),
              min(summ$ess_bulk, na.rm = TRUE), divs,
              elapsed, remaining))
}

results <- results[!sapply(results, is.null)]
cat(sprintf("\nCompleted %d / %d replications.\n\n", length(results), R_reps))

# ---- Aggregate ---------------------------------------------------------------
pm <- t(sapply(results, function(r) r$post_mean))
lo <- t(sapply(results, function(r) r$ci_lo))
hi <- t(sapply(results, function(r) r$ci_hi))
colnames(pm) <- colnames(lo) <- colnames(hi) <- pars_of_interest

bias  <- colMeans(pm) - truth
rmse  <- sqrt(colMeans(sweep(pm, 2, truth, "-")^2))
cover <- sapply(pars_of_interest, function(p) mean(lo[, p] <= truth[p] & truth[p] <= hi[, p]))
mc_se <- apply(pm, 2, sd) / sqrt(nrow(pm))

agg <- data.frame(
  parameter   = pars_of_interest,
  truth       = round(truth, 3),
  mean_postmean = round(colMeans(pm), 3),
  bias        = round(bias, 4),
  mc_se_bias  = round(mc_se, 4),
  rmse        = round(rmse, 4),
  coverage_90 = round(cover, 3)
)
cat("\n--- Aggregate recovery (90% nominal CI) ---\n")
print(agg, row.names = FALSE)

cat("\n--- Diagnostics across replications ---\n")
cat(sprintf("Max R-hat:  median=%.3f, max=%.3f\n",
            median(diag_log$max_rhat), max(diag_log$max_rhat)))
cat(sprintf("Min ESS:    median=%.0f, min=%.0f\n",
            median(diag_log$min_ess), min(diag_log$min_ess)))
cat(sprintf("Divergences: total=%d, max per rep=%d\n",
            sum(diag_log$divergences), max(diag_log$divergences)))

saveRDS(list(
  results = results, agg = agg, diag_log = diag_log,
  truth = truth, R_reps = R_reps, N_each = N_each
), "dev/stan_composite_poc/sim_study_result.rds")
cat("\nSaved dev/stan_composite_poc/sim_study_result.rds\n")
