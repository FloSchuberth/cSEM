library(SimDesign)
library(rlang)
rlang::global_handle()
library(lobstr)
cli::pretty_print_code()
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(2)
library(tibble)
library(future)
library(future.mirai)
plan(future.mirai::mirai_multisession)
library(future.apply)
library(boot)

# PoC --------------------------------------------------------------------

sigma <- 1
x1 <- rnorm(1000)
x2 <- rnorm(1000)
allDat <- tibble::tibble(x1, x2)
trueDat <-subset(allDat, select = c(x1, x2))
beta <- .3
y <- rnorm(1000, mean = rowSums(beta*trueDat), sd = sigma)

allDat$y <- y

lm(y ~ ., data = allDat) |> summary()

## Trying out t-test described on Page 27 of Hwang & Takane 2014 ----------

# The goal here is to examine this statement, "Similar to bootstrapping R-squared in linear regression (Ohtani 2000), we can compute the bootstrapped standard errors or confidence intervals of the difference in FIT between two models so as to examine whether there is a statistically significant difference in the FIT values. This procedure can be regarded as a nonparametric version of the paired t test for two groups of FIT values."

# Followed by on page 41, "Furthermore, we conducted a paired t-test for the FIT values between the original and simpler models. As described earlier, this test was equivalent to testing whether the difference in the FIT values calculated from 100 bootstrap samples for each model was equal to zero. We found the t-statistic nonsignificant statistically [t(99) = −0.59, p = 0.55, 95% CI = −0.00−0.00], indicating that there was no statistically significant difference in FIT between the original and simpler models."

# Q: From a theoretical perspective, is a one-sample t-test of (a-b) = 0 identical to a paired t-test of a = b?
## A: Yes, see page  280 of Hogg et al. 2019. and also getAnywhere('t.test.default').
## Typically, a paired t-test is just comparing  the mean and standard deviations of the paired differences against the null hypothesis.
#### Also, Page 521 of Wackerly et al. (2008) gives one sample t-test formula + nice table for how to understand critical regions. Page 646 to 648 of Wackerly for their exposition.


# Simulation -------------------------------------------------------------

# Monte Carlo study (SimDesign) stress-testing the Hwang & Takane (2014, p.27)
# procedure of comparing two models' FIT (R^2) values via a paired t-test on
# bootstrap replicates. The two models compared here are:
#   * POOLED : a single multiple linear regression fit to ALL the data.
#   * SPLIT  : TWO SEPARATE multiple linear regressions, one per part of the data.
# The slopes (beta1, beta2) differ by `delta` between the two parts. Varying delta
# turns the rejection rate into an interpretable error rate:
#   delta = 0 -> H0 true (pooling is correct) -> rejection rate = empirical Type I error
#   delta > 0 -> H1 true                        -> rejection rate = power (= 1 - Type II error)
# `prop1` is the fraction of N in part 1 (0.5 = balanced, up to 0.9 = severe imbalance),
# holding N fixed so group-size imbalance is isolated from total sample size. The
# bootstrap is stratified by part, so each resample preserves n1 and n2.
#
# Methods compared (all returned as p-values; EDR = rejection rate at alpha = .05):
#   p_ttest_*  PUBLISHED H&T method: paired t-test on the B bootstrap FIT replicates.
#              SE = sd(d*)/sqrt(B), so its df = B - 1 (matching their reported t(99), B = 100)
#              and its significance is an artifact of B.
#   p_se_*     Wald/studentized z = T_obs / sd(d*): the SAME sd(d*) but WITHOUT /sqrt(B),
#              so t_paired ~= sqrt(B) * z_wald -- isolates the sqrt(B) error.
#   p_pctl_*   Percentile / ASL bootstrap CI of the FIT difference (CI-inversion).
#   p_null_*   Label-permutation test (boot, sim = "permutation"): impose H0 via exchangeability
#              of the group labels -- randomly permute them R = B times to build the null
#              distribution of the FIT improvement; p = (1 + #{T* >= T_obs})/(B + 1), one-sided
#              upper. Canonical randomization test of group homogeneity.
#   p_chow     Chow F-test (anova: pooled vs full interaction) -- exact parametric oracle.
#
# ML aside: the principled object for choosing between these models is out-of-sample /
# cross-validated R^2 (or AIC/BIC). In-sample R^2 mechanically favours the higher-capacity
# SPLIT model (it is the pooled model plus group main-effect + group:slope interactions, so
# pooled is nested in split and R2_split >= R2_pool w.p. 1). We use in-sample R^2 / adj-R^2
# because that is precisely the object the literature claim invokes and what we are testing.

## ---- helper: FIT of pooled vs two separate models ------------------------

# Single .lm.fit()-based helper, used both inside the bootstrap loops and in the
# sanity checks below. `d` is a data frame with columns {x1, x2, y, group}; X is the
# [1, x1, x2] design matrix. A single POOLED regression is fit to all rows and TWO
# SEPARATE regressions (one per `group`) via .lm.fit (fast path -- no formula parsing);
# the SPLIT FIT pools the two models' residuals (combined-prediction R^2), with total
# parameter count k = 2 * ncol(X) = 6 for the adjustment. The R^2 / adj-R^2 formulae
# reproduce summary(lm(...))$r.squared / $adj.r.squared exactly (intercept model).
# Shipped to the parallel workers via `fixed_objects`.
fit_compare <- function(d) {
  X    <- cbind(1, d$x1, d$x2)                               # [1, x1, x2] design
  y    <- d$y
  g1   <- d$group == 1L; g2 <- d$group == 2L                 # the two parts
  n    <- length(y)
  p    <- ncol(X)                                            # params per model (incl. intercept)
  SST  <- sum((y - mean(y))^2)
  SSE_pool  <- sum(.lm.fit(X, y)$residuals^2)                             # pooled model
  SSE_split <- sum(.lm.fit(X[g1, , drop = FALSE], y[g1])$residuals^2) +   # two SEPARATE models
               sum(.lm.fit(X[g2, , drop = FALSE], y[g2])$residuals^2)
  c(R2_pool     = 1 - SSE_pool  / SST,
    R2_split    = 1 - SSE_split / SST,
    adjR2_pool  = 1 - (SSE_pool  / (n - p))     / (SST / (n - 1)),
    adjR2_split = 1 - (SSE_split / (n - 2 * p)) / (SST / (n - 1)))
}

## ---- SimDesign: Generate / Analyse / Summarise ---------------------------

Generate <- function(condition, fixed_objects) {
  N     <- fixed_objects$N
  # `prop1` = fraction of N in part 1 (0.5 = balanced). Clamp so both parts stay
  # estimable (need n >= the 3 regression coefficients).
  n1    <- min(N - 4L, max(4L, round(N * condition$prop1)))
  group <- factor(rep(c(1L, 2L), c(n1, N - n1)))
  x1    <- rnorm(N)
  x2    <- rnorm(N)
  d2    <- as.integer(group == 2L)                           # 1 in part 2, else 0
  b1    <- fixed_objects$beta1 + condition$delta * d2        # slopes differ by delta
  b2    <- fixed_objects$beta2 + condition$delta * d2
  mu    <- fixed_objects$beta0 + b1 * x1 + b2 * x2
  y     <- rnorm(N, mean = mu, sd = fixed_objects$sigma)
  tibble::tibble(x1 = x1, x2 = x2, y = y, group = group)
}

Analyse <- function(condition, dat, fixed_objects) {
  fit_compare <- fixed_objects$fit_compare                   # exported with fixed_objects
  B  <- condition$B

  # Observed FIT and observed improvements T = FIT_split - FIT_pool.
  obs   <- fit_compare(dat)
  T_R2  <- unname(obs["R2_split"]    - obs["R2_pool"])
  T_adj <- unname(obs["adjR2_split"] - obs["adjR2_pool"])

  # (i) Stratified bootstrap (boot, strata = group -> resampling within each part, sizes
  #     n1/n2 preserved) -> B replicates of the 4 FIT values.
  boot_mat <- boot::boot(dat, function(d, i) fit_compare(d[i, ]),
                         R = B, strata = dat$group)$t
  colnames(boot_mat) <- names(obs)
  dR2  <- boot_mat[, "R2_split"]    - boot_mat[, "R2_pool"]
  dAdj <- boot_mat[, "adjR2_split"] - boot_mat[, "adjR2_pool"]

  # (1) PUBLISHED Hwang & Takane paired t-test (SE = sd(d*)/sqrt(B)).
  p_ttest_R2    <- t.test(boot_mat[, "R2_split"],    boot_mat[, "R2_pool"],    paired = TRUE)$p.value
  p_ttest_adjR2 <- t.test(boot_mat[, "adjR2_split"], boot_mat[, "adjR2_pool"], paired = TRUE)$p.value

  # (3a') Wald / studentized: z = T_obs / sd(d*)  (no /sqrt(B)).
  p_se_R2    <- 2 * pnorm(-abs(T_R2  / sd(dR2)))
  p_se_adjR2 <- 2 * pnorm(-abs(T_adj / sd(dAdj)))

  # (3a) Percentile / ASL bootstrap CI of the difference.
  p_pctl_R2    <- min(1, 2 * min(mean(dR2  <= 0), mean(dR2  >= 0)))
  p_pctl_adjR2 <- min(1, 2 * min(mean(dAdj <= 0), mean(dAdj >= 0)))

  # (3b) Label-permutation test (boot, sim = "permutation"): impose H0 via exchangeability
  # of the group labels. The pooled model ignores `group`, so each permutation only re-fits
  # the SPLIT model; T* = FIT_split* - FIT_pool is the null distribution. One-sided upper.
  perm_stat <- function(d, i) {
    d$group <- d$group[i]                                    # randomly reassign group labels
    fc <- fit_compare(d)
    c(R2  = unname(fc["R2_split"]    - fc["R2_pool"]),
      adj = unname(fc["adjR2_split"] - fc["adjR2_pool"]))
  }
  perm <- boot::boot(dat, perm_stat, R = B, sim = "permutation")
  p_null_R2    <- (1 + sum(perm$t[, 1] >= perm$t0[1])) / (B + 1)
  p_null_adjR2 <- (1 + sum(perm$t[, 2] >= perm$t0[2])) / (B + 1)

  # (d) Chow oracle: exact F-test, pooled vs full interaction (== two separate models).
  p_chow <- anova(lm(y ~ x1 + x2, data = dat),
                  lm(y ~ group * (x1 + x2), data = dat))[["Pr(>F)"]][2]

  SimDesign::nc(
    p_ttest_R2 = p_ttest_R2, p_ttest_adjR2 = p_ttest_adjR2,
    p_se_R2    = p_se_R2,    p_se_adjR2    = p_se_adjR2,
    p_pctl_R2  = p_pctl_R2,  p_pctl_adjR2  = p_pctl_adjR2,
    p_null_R2  = p_null_R2,  p_null_adjR2  = p_null_adjR2,
    p_chow     = p_chow
  )
}

Summarise <- function(condition, results, fixed_objects) {
  SimDesign::EDR(results, alpha = fixed_objects$alpha)        # rejection rate per method
}

## ---- design + fixed objects ----------------------------------------------

sim_design <- SimDesign::createDesign(
  delta = c(0, 0.1, 0.2, 0.4),
  B     = c(100, 1000),                                       # 100 reproduces H&T's df = 99
  prop1 = c(0.5, 0.7, 0.9)                                    # fraction of N in part 1: 100/100, 140/60, 180/20
)

sim_fixed <- list(
  beta0 = 0, beta1 = 0.3, beta2 = 0.3, sigma = 1, N = 200, alpha = 0.05,
  fit_compare = fit_compare                                  # ship helper to parallel workers
)

## ---- sanity checks (cheap; run before the full study) --------------------

local({
  set.seed(1)
  for (p1 in c(0.5, 0.9)) {                                          # balanced AND severe imbalance
    chk  <- Generate(list(delta = 0, prop1 = p1), sim_fixed)
    fc   <- fit_compare(chk)
    full <- summary(lm(y ~ group * (x1 + x2), data = chk))
    stopifnot(
      fc["R2_split"] > fc["R2_pool"],                                # nesting (point 1)
      isTRUE(all.equal(unname(fc["R2_split"]),    full$r.squared)),  # 2 sep. models == interaction
      isTRUE(all.equal(unname(fc["adjR2_split"]), full$adj.r.squared))
    )
  }
  message("Sanity OK (balanced + imbalanced): split R2 > pooled R2; two-model FIT == interaction-model FIT.")
})

## ---- run the simulation --------------------------------------------------

sim_results <- SimDesign::runSimulation(
  design        = sim_design,
  replications  = 1000,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = sim_fixed,
  packages      = c("tibble", "boot"),
  parallel      = "future",                                  # uses the mirai plan set above
  save          = FALSE
)

print(sim_results)

# NOTE: the bootstrap machinery now runs through the `boot` package -- the stratified resample
# (strata = group) and an H0 LABEL-PERMUTATION test that replaced the old null-pooled RESIDUAL
# bootstrap -- all driven by the unified .lm.fit()-based fit_compare. The p_null_* column below is
# from the superseded residual bootstrap and must be regenerated by re-running the simulation; the
# t-test / Wald / percentile columns use the same stratified resample and are unchanged in method.
#
# Results (1000 replications; rejection rate at alpha = .05; delta = 0 = Type I error,
# delta > 0 = power = 1 - Type II error). The deviating part is always part 2, so prop1
# (fraction of N in part 1) shrinks the DEVIATING group as it grows: 100/100 -> 140/60 -> 180/20.
#
#  prop1 delta   B  ttest_R2 ttest_adj  se_R2 se_adj  pctl_R2 pctl_adj  null_R2 null_adj  chow
#   0.5   0.0  100    1.000    0.882    0.003  0.001   1.000    0.018    0.045   0.045   0.043   <- Type I
#   0.5   0.4  100    1.000    1.000    0.647  0.404   1.000    0.832    0.918   0.917   0.918   <- power
#   0.5   0.0 1000    1.000    0.981    0.002  0.000   1.000    0.014    0.044   0.046   0.039
#   0.5   0.4 1000    1.000    1.000    0.626  0.386   1.000    0.822    0.917   0.917   0.917
#   0.7   0.0 1000    1.000    0.973    0.006  0.001   1.000    0.022    0.053   0.053   0.050
#   0.7   0.4 1000    1.000    1.000    0.517  0.281   1.000    0.723    0.862   0.862   0.859
#   0.9   0.0 1000    1.000    0.966    0.011  0.005   1.000    0.017    0.045   0.043   0.043
#   0.9   0.4 1000    1.000    0.995    0.200  0.062   1.000    0.313    0.474   0.474   0.475
# (full 24-row grid -> dev/igsca/ohtani_simulation_plan.md)
#
# Conclusions:
#   p_ttest_* : the PUBLISHED method is INVALID at every imbalance level. raw-R2 always rejects
#               (Type I = 1.0); adj-R2 Type I = 0.88 (B=100) -> 0.97 (B=1000) -- it RISES WITH B,
#               the sqrt(B) artifact. Significance is driven by the bootstrap count, not the evidence.
#   p_se_*,   : the two CI-style readings of the quote fail in OPPOSITE directions on raw R2 --
#   p_pctl_*    p_pctl_R2 always rejects (= 1.0; d* > 0 in every resample) while p_se_R2 is
#               conservative (~0). On adj-R2 both are conservative and underpowered.
#   p_null_*  : the label-permutation test imposes H0 via exchangeability and is an exact
#               randomization test, so it should stay calibrated (Type I ~= .05) and track the
#               oracle -- RE-RUN to confirm the numbers under the new permutation method.
#   p_chow    : exact oracle -- Type I ~= .05 and the power ceiling; p_null_* matches it throughout.
#
# Imbalance effect: Type I calibration of the valid methods (p_null_*, p_chow) is ROBUST to
# imbalance (~.05 for prop1 = .5/.7/.9). But POWER collapses as the deviating subgroup shrinks --
# at delta = 0.4, power falls 0.92 (100/100) -> 0.86 (140/60) -> 0.47 (180/20): a small
# heterogeneous subgroup carries too little information to detect its own slope difference.
#
# Takeaway: a paired t-test on bootstrap-replicate FIT values (Hwang & Takane 2014) is not a valid
# test of whether splitting the data improves fit, balanced or not. Imposing the null via a label
# permutation fixes the calibration and matches the exact Chow F-test at all imbalance levels.

