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
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

# Simulation -------------------------------------------------------------

# Authored by Claude 4.8. Revised and verified by me.

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

fit_compare <- function(d) {
  X    <- cbind(1, d$x1, d$x2)                               # [1, x1, x2] design matrix
  y    <- d$y
  g1   <- d$group == 1 
  g2 <- d$group == 2
  n    <- length(y)
  p    <- ncol(X)                                            # params per model (incl. intercept)
  SST  <- sum((y - mean(y))^2)
  SSE_pool  <- sum(.lm.fit(X, y)$residuals^2)                             # pooled model
  SSE_split <- sum(.lm.fit(X[g1, , drop = FALSE], y[g1])$residuals^2) +   # two SEPARATE models
               sum(.lm.fit(X[g2, , drop = FALSE], y[g2])$residuals^2)
  c(
    R2_pool = 1 - SSE_pool / SST,
    R2_split = 1 - SSE_split / SST,
    adjR2_pool = 1 - (SSE_pool / (n - p)) / (SST / (n - 1)),
    adjR2_split = 1 - (SSE_split / (n - 2 * p)) / (SST / (n - 1))
  )
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

  # R^2 Statistics
  obs   <- fit_compare(dat) 
  T_R2  <- unname(obs["R2_split"]    - obs["R2_pool"]) 
  T_adj <- unname(obs["adjR2_split"] - obs["adjR2_pool"])

  # (i) Stratified bootstrap (boot, strata = group -> resampling within each part, sizes so ratio of n1/n2 preserved) -> B replicates of the 4 FIT values.
  boot_mat <- boot::boot(dat, function(d, i) fit_compare(d[i, ]),
                         R = B, strata = dat$group)$t
  colnames(boot_mat) <- names(obs)
  dR2  <- boot_mat[, "R2_split"]    - boot_mat[, "R2_pool"]
  dAdj <- boot_mat[, "adjR2_split"] - boot_mat[, "adjR2_pool"]

  # (1) Published Hwang & Takane (2014) paired t-test (SE = sd(d*)/sqrt(B)).
  p_ttest_R2    <- t.test(boot_mat[, "R2_split"],    boot_mat[, "R2_pool"],    paired = TRUE)$p.value
  p_ttest_adjR2 <- t.test(boot_mat[, "adjR2_split"], boot_mat[, "adjR2_pool"], paired = TRUE)$p.value

  # (3a') Wald / studentized: z = T_obs / sd(d*)  (no /sqrt(B)).
  # p_se_R2    <- 2 * pnorm(-abs(T_R2  / sd(dR2)))
  # p_se_adjR2 <- 2 * pnorm(-abs(T_adj / sd(dAdj)))

  # (3a) Percentile / ASL bootstrap CI of the difference.
  # p_pctl_R2    <- min(1, 2 * min(mean(dR2  <= 0), mean(dR2  >= 0)))
  # p_pctl_adjR2 <- min(1, 2 * min(mean(dAdj <= 0), mean(dAdj >= 0)))

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
  # p_chow <- anova(lm(y ~ x1 + x2, data = dat),
  #                 lm(y ~ group * (x1 + x2), data = dat))[["Pr(>F)"]][2]

  SimDesign::nc(
    p_ttest_R2 = p_ttest_R2,
    p_ttest_adjR2 = p_ttest_adjR2,
    # p_se_R2 = p_se_R2,
    # p_se_adjR2 = p_se_adjR2,
    # p_pctl_R2 = p_pctl_R2,
    # p_pctl_adjR2 = p_pctl_adjR2,
    p_null_R2 = p_null_R2,
    p_null_adjR2 = p_null_adjR2 #,
    # p_chow = p_chow
  )
}

Summarise <- function(condition, results, fixed_objects) {
  SimDesign::EDR(results, alpha = fixed_objects$alpha)        # rejection rate per method
}

## ---- run the simulation --------------------------------------------------

sim_results <- SimDesign::runSimulation(
  design        = sim_design,
  replications  = 1000,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = sim_fixed,
  packages      = c("tibble", "boot"),
  parallel = "mirai",
  ncores = 20,
  save          = FALSE
)

print(sim_results)

## ---- plot Type I / power curves ------------------------------------------

plot_df <- sim_results |>
  as.data.frame() |>
  dplyr::select(delta, B, prop1, starts_with("p_")) |>
  tidyr::pivot_longer(
    starts_with("p_"),
    names_to = "stat",
    values_to = "reject"
  ) |>
  tidyr::extract(stat, c("method", "metric"), regex = "^p_(.*)_(R2|adjR2)$") |>
  dplyr::mutate(
    method = dplyr::recode(
      method,
      ttest = "H&T",
      # ttest = "H&T paired t (/\u221AB)",
      # se = "Wald / SE",
      null = "Permutation"
    ),
    metric = dplyr::recode(metric, R2 = "R\u00B2", adjR2 = "adj-R\u00B2"),
    # factor versions for facet/linetype labels; keep delta numeric for the x-axis
    Bf = factor(
      B,
      levels = sort(unique(B)),
      labels = paste0("B = ", sort(unique(B)))
    ),
    prop1f = factor(
      prop1,
      levels = sort(unique(prop1)),
      labels = paste0("prop1 = ", sort(unique(prop1)))
    )
  )

# Facets: metric (rows) separates the scale-crushing raw R2 from adj-R2; prop1
# (cols) walks through group-size imbalance. linetype = B exposes the sqrt(B)
# Type-I climb of the published t-test directly within each panel.
p_curves <- ggplot(plot_df, aes(delta, reject, colour = method, linetype = Bf)) +
  geom_hline(yintercept = sim_fixed$alpha,            # alpha reference (Type I target)
             linetype = "dotted", colour = "black") +
  geom_hline(yintercept = .8, linetype = "dotted", colour = "black") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_grid(prop1f ~ metric) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = sort(unique(plot_df$delta))) +
  labs(
    x = expression(delta ~ "(slope difference; " * delta == 0 * " = Type I)"),
    y = "Rejection rate",
    colour = "Method", linetype = "Bootstraps",
    title = "Type I error (\u03B4 = 0) and power (\u03B4 > 0) by method",
    caption = "Dotted line = nominal \u03B1"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

print(p_curves)

# Save for the Typst write-up (dev/igsca/ohtani.typ embeds this in the Results section).

ggsave(
  here::here("dev/igsca/ohtani_curves.png"),
  p_curves,
  width = 5,
  height = 6,
  dpi = 300
)


## ---- APA results table -> Typst partial (tinytable, auto-generated) ------
# Mirrors facet_grid(prop1f ~ metric): metric = column spanners, prop1 = row
# groups, method x B = rows, delta = columns. EVERYTHING below is derived from
# plot_df / sim_results -- no grid levels are hard-coded -- so re-sourcing this
# script refreshes BOTH the figure and the table from the same data. tinytable
# emits Typst natively (and does the number formatting); ohtani.typ #includes
# the partial it writes.
library(tinytable)
library(stringr)

# Dimensions read straight off the plotted data (nothing about the grid hard-coded).
deltas    <- sort(unique(plot_df$delta))
metrics   <- unique(plot_df$metric)    # e.g. "R²", "adj-R²"
methods   <- unique(plot_df$method)    # first-appearance order == sim_results column order
props     <- sort(unique(plot_df$prop1))
rows_per_group <- length(methods) * length(unique(plot_df$B))
N <- sim_fixed$N

# Value-column order: metric-major, delta-minor (expand_grid varies its last arg fastest).
val_order <- tidyr::expand_grid(metric = metrics, delta = deltas) |>
  dplyr::mutate(key = paste0(metric, delta)) |>
  dplyr::pull(key)

# One row per (prop1, method, B); value columns crossed metric x delta, left numeric
# so tinytable formats them.
wide <- plot_df |>
  dplyr::mutate(key = paste0(metric, delta), method = factor(method, levels = methods)) |>
  dplyr::select(prop1, method, B, key, reject) |>
  tidyr::pivot_wider(names_from = key, values_from = reject) |>
  dplyr::arrange(prop1, method, B)

body <- wide |>
  dplyr::select(Method = method, B, dplyr::all_of(val_order)) |>
  dplyr::mutate(Method = as.character(Method), B = as.character(B))

# Display headers: two stub columns, then the bare delta values repeated under each
# metric (what delta indexes is explained in the note, keeping the columns narrow).
delta_labels <- as.character(deltas)
colnames(body) <- c("Method", "$B$", rep(delta_labels, times = length(metrics)))

# Column spanners: one per metric, over its block of delta columns (after the 2 stub cols).
spanners <- lapply(seq_along(metrics), function(k) {
  first <- 3 + (k - 1) * length(deltas)
  first:(first + length(deltas) - 1)
})
names(spanners) <- metrics

# Row groups: one per prop1, labelled with realised part sizes (matching Generate's clamp).
n1 <- pmin(N - 4, pmax(4, round(N * props)))
group_rows <- as.list(seq(from = 1, by = rows_per_group, length.out = length(props)))
names(group_rows) <- paste0("$n_1\\/n_2 = ", n1, "\\/", N - n1, "$")

cap <- paste0(
  "Empirical rejection rates over ", as.integer(sim_results$REPLICATIONS[1]),
  " Monte Carlo replications at $alpha = .05$, by model-comparison method and number of ",
  "bootstrap resamples $B$. The $delta = 0$ columns give the Type I error rate; ",
  "$delta > 0$ columns give power."
)
note <- paste0(
  "_Note._ Columns within each metric block index the slope difference $delta$. ",
  "Cells are the proportion of replications rejecting $H_0$ at $alpha = .05$; ",
  "$delta = 0$ is the Type I error rate and $delta > 0$ is power. Row groups give the ",
  "part-1/part-2 sample sizes; the deviating subgroup is part 2."
)

tbl <- tt(body, caption = cap, notes = note) |>
  format_tt(j = 3:ncol(body), digits = 3, num_fmt = "decimal", num_zero = TRUE) |>
  group_tt(j = spanners) |>
  group_tt(i = group_rows) |>
  style_tt(j = 2:ncol(body), align = "c")

typ_out <- save_tt(tbl, output = "typst")
# Share the existing #figure(table(...)) "Table" counter, and add a label for @tbl-results.
typ_out <- str_replace(typ_out, fixed('kind: "tinytable"'), "kind: table")
typ_out <- paste0(typ_out, "\n<tbl-results>\n")
writeLines(typ_out, here::here("dev/igsca/ohtani_results_table.typ"))