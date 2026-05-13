# =============================================================================
# Single-replication diagnostic: draw one dataset, fit Stan, examine
# convergence and parameter recovery in detail.
# =============================================================================
source("dev/stan_composite_poc/00_setup.R")

set.seed(1)
N <- 500
sim <- simulate_one(N = N, seed = 1)
tp  <- true_params()

# ---- Sanity: cSEM PLS / composite estimator on the same data ------------
dat_df <- as.data.frame(cbind(sim$x, sim$y))
csem_fit <- cSEM::csem(dat_df, dgp_syntax())
cat("\n--- cSEM PLS estimates (sanity check) ---\n")
sm <- cSEM::summarize(csem_fit)
print(sm$Estimates$Weight_estimates[, c("Name", "Estimate")])
print(sm$Estimates$Path_estimates[, c("Name", "Estimate")])
print(sm$Estimates$Loading_estimates[, c("Name", "Estimate")])

# ---- Compile and fit Stan -----------------------------------------------
mod <- cmdstan_model("dev/stan_composite_poc/composite_ho.stan")

stan_data <- list(
  N = sim$N, K = sim$K, M = sim$M,
  x = sim$x, y = sim$y
)

fit <- mod$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 500,
  seed            = 12345,
  adapt_delta     = 0.95,
  init            = init_from_pls(sim, n_chains = 4)
)

cat("\n--- Stan diagnostics summary ---\n")
fit$cmdstan_diagnose()

cat("\n--- Posterior summary (key parameters) ---\n")
pars <- c("w[1]", "w[2]", "w[3]",
          "beta", "sigma_eta",
          "lambda[1]", "lambda[2]", "lambda[3]",
          "sigma_y[1]", "sigma_y[2]", "sigma_y[3]",
          "R2_eta")
print(fit$summary(variables = pars))

cat("\n--- Truth vs posterior mean ---\n")
post <- posterior::as_draws_df(fit$draws(variables = pars))
post <- sign_fix_draws(post, K = 3)

truth <- c(tp$w, tp$beta, tp$sigma_eta, tp$lambda, tp$sigma_y,
           tp$beta^2)  # R^2 = beta^2 / Var(eta) = beta^2 / 1 = beta^2
names(truth) <- pars
ests  <- sapply(pars, function(p) mean(post[[p]]))
qlo   <- sapply(pars, function(p) quantile(post[[p]], 0.05))
qhi   <- sapply(pars, function(p) quantile(post[[p]], 0.95))
tab <- data.frame(parameter = pars, truth = round(truth, 3),
                  post_mean = round(ests, 3),
                  ci_lo = round(qlo, 3), ci_hi = round(qhi, 3),
                  bias = round(ests - truth, 3),
                  in_ci = (truth >= qlo & truth <= qhi))
print(tab, row.names = FALSE)

saveRDS(list(fit_summary = fit$summary(),
             draws_subset = post,
             truth = tp,
             recovery_table = tab),
        file = "dev/stan_composite_poc/single_run_result.rds")
cat("\nSaved dev/stan_composite_poc/single_run_result.rds\n")
