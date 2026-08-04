# Hard gate: the partition-family extree driver must reproduce ctree() ------
#
# After the 2026-07-14 refactor the COIN family rides partykit::ctree()
# directly (nothing left to prove there); what stays custom is grow_extree()
# -- the extree_fit() driver + generic argmax splitfun + candidate_
# partitions() used by the partition testers. This script drives that engine
# with cheap Gaussian fixtures (identity influence, no IGSCA fitting) whose
# tests are computed by the coin package -- an implementation INDEPENDENT of
# partykit's internal libcoin -- and asserts EXACT agreement with
# partykit::ctree() on split variables, cutpoints and partitions.
#
# Run: Rscript R/tree/igscaTree/ctree_equivalence_proof.R
# (exits non-zero on any failure). coin and MASS are PROOF-ONLY dependencies
# now; the production code needs neither. select_gauss below is also the
# template for any future conditional-inference selector that cannot ride
# igsca_ctree() (e.g. COIN selection mixed with a partition splitter).

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressMessages({library(partykit); library(coin); library(MASS)})
source(here::here("R/tree/igscaTree/igsca_tree.R"))

same_partition <- function(a, b) {           # from R/tree/HowTree/lavaan.R
  tab <- table(a, b)
  all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
}

# ---- Gaussian fixtures (proof-local plain functions) ----------------------

#' Strasser-Weber (1999) standardized quadratic statistic for one binary
#' partition of the influence h -- the ctree split criterion, reimplemented
#' independently of libcoin.
sw_quadstat <- function(h, goes_left) {
  h <- as.matrix(h); n <- nrow(h); nl <- sum(goes_left)
  if (nl == 0L || nl == n) return(NA_real_)
  hbar <- colMeans(h)
  Tt <- colSums(h[goes_left, , drop = FALSE])
  mu <- nl * hbar
  Vh <- crossprod(sweep(h, 2, hbar)) / n
  Sig <- nl * (n - nl) / (n - 1) * Vh
  d <- Tt - mu
  as.numeric(t(d) %*% MASS::ginv(Sig) %*% d)
}

#' Identity-influence selector (extree selectfun contract; exact formal
#' names required by extree_fit). h = y itself -- matches ctree's numeric
#' default -- tested with coin::independence_test, which computes exactly
#' what ctree's internal libcoin computes. model$estfun here is node-local
#' (rows = subset), the custom-path convention.
select_gauss <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
  mf <- model.frame(data)
  crit <- matrix(NA_real_, 2L, ncol(mf),
                 dimnames = list(c("statistic", "p.value"), names(mf)))
  dd <- data.frame(h = model$estfun[, 1], .z = NA)
  for (j in whichvar) {
    dd$.z <- mf[[j]][subset]
    if (length(unique(dd$.z)) < 2L) {
      next
    }
    it <- coin::independence_test(
      h ~ .z,
      data = dd,
      teststat = "quadratic",
      distribution = "asymptotic"
    )
    crit["statistic", j] <- log(max(
      as.numeric(coin::statistic(it)),
      .Machine$double.eps
    ))
    crit["p.value", j] <- log1p(-min(as.numeric(coin::pvalue(it)), 1 - 1e-12))
  }
  list(criteria = crit)
}

#' Matched splitter kernel: SW quadratic statistic on the node influence.
split_gauss <- function(model, mf, subset, goes_left, ctrl)
  sw_quadstat(model$estfun, goes_left)

## Our engine, matched knob-for-knob to ctree.
grow_ours <- function(dat, yvar, zvars, alpha, minbucket) {
  ctrl <- igsca_tree_control(
    alpha = alpha,
    bonferroni = TRUE,
    minbucket = minbucket,
    minsplit = 2L * minbucket,
    maxdepth = 10L,
    max_cuts = .Machine$integer.max
  )
  fml <- stats::as.formula(paste(yvar, "~", paste(zvars, collapse = " + ")))
  d <- extree_data(fml, data = dat, yx = "none", nmax = c(yx = Inf, z = Inf))
  mf <- model.frame(d)
  ctrl$collector <- new_collector()
  ## Gaussian trafo: influence = the response itself, refit-free.
  trafo <- function(
    subset,
    weights,
    info = NULL,
    estfun = TRUE,
    object = TRUE,
    ...
  ) {
    y <- mf[[yvar]][subset]
    list(
      coefficients = NA_real_,
      objfun = sum((y - mean(y))^2),
      object = NULL,
      estfun = matrix(y, ncol = 1),
      converged = TRUE,
      nobs = length(subset)
    )
  }
  nodes <- grow_extree(d, trafo, select_gauss, split_gauss, ctrl)
  fitted <- data.frame(
    "(fitted)" = fitted_node(nodes, mf),
    "(weights)" = rep.int(1L, nrow(mf)),
    check.names = FALSE
  )
  ret <- party(nodes, data = mf, fitted = fitted)
  class(ret) <- c("constparty", class(ret))
  ret
}

check_same <- function(ours, ct, label) {
  po <- predict(ours, type = "node")
  pc <- predict(ct, type = "node")
  vo <- sort(used_split_vars(ours))
  vc <- sort(setdiff(
    names(ct$data)[unique(na.omit(unlist(nodeapply(
      ct,
      ids = nodeids(ct),
      function(nd) {
        s <- split_node(nd)
        if (is.null(s)) NA_integer_ else varid_split(s)
      }
    ))))],
    character(0)
  ))
  ok <- same_partition(po, pc) && identical(vo, vc)
  cat(sprintf(
    "%-28s vars[%s|%s] partition_match=%s\n",
    label,
    paste(vo, collapse = ","),
    paste(vc, collapse = ","),
    ok
  ))
  stopifnot(ok)
}

alpha <- 0.05; minbucket <- 7L
for (seed in c(11, 22, 33)) {
  set.seed(seed)
  n <- 300
  ## design 1: numeric-only covariates, one informative
  d1 <- data.frame(z1 = rnorm(n), z2 = rnorm(n), z3 = runif(n))
  d1$y <- ifelse(d1$z2 > 0.3, 1.2, 0) + rnorm(n)
  ct1 <- ctree(
    y ~ z1 + z2 + z3,
    data = d1,
    control = ctree_control(
      teststat = "quadratic",
      testtype = "Bonferroni",
      splitstat = "quadratic",
      alpha = alpha,
      minbucket = minbucket,
      minsplit = 2L * minbucket,
      maxsurrogate = 0L
    )
  )
  check_same(
    grow_ours(d1, "y", c("z1", "z2", "z3"), alpha, minbucket),
    ct1,
    paste0("seed ", seed, " numeric")
  )
  ## design 2: mixed numeric + binary + 4-level factor
  d2 <- data.frame(
    z1 = rnorm(n),
    z2 = factor(sample(1:2, n, TRUE)),
    z3 = factor(sample(letters[1:4], n, TRUE))
  )
  d2$y <- ifelse(d2$z2 == "2", 0.9, 0) + rnorm(n)
  ct2 <- ctree(
    y ~ z1 + z2 + z3,
    data = d2,
    control = ctree_control(
      teststat = "quadratic",
      testtype = "Bonferroni",
      splitstat = "quadratic",
      alpha = alpha,
      minbucket = minbucket,
      minsplit = 2L * minbucket,
      maxsurrogate = 0L
    )
  )
  check_same(
    grow_ours(d2, "y", c("z1", "z2", "z3"), alpha, minbucket),
    ct2,
    paste0("seed ", seed, " mixed")
  )
}

## Null-uniformity smoke (forced stump, all-noise covariates of unequal
## cardinality): root_pick should be ~uniform across the 3 covariates.
set.seed(99)
B <- 400
picks <- character(B)
for (b in seq_len(B)) {
  n <- 120
  d0 <- data.frame(
    z1 = rnorm(n),
    z2 = factor(sample(1:2, n, TRUE)),
    z3 = factor(sample(1:8, n, TRUE))
  )
  d0$y <- rnorm(n)
  tr <- grow_ours(d0, "y", c("z1", "z2", "z3"), alpha = 1 - 1e-9, minbucket = 7L)
  # root pick = first split variable (forced by alpha ~ 1); NA if degenerate
  us <- used_split_vars(tr)
  picks[b] <- if (length(us)) us[1] else NA_character_
}
tab <- table(factor(picks, levels = c("z1", "z2", "z3")))
gof <- suppressWarnings(chisq.test(tab))
cat("null selection frequencies:\n")
print(prop.table(tab))
cat(sprintf("chi-square GOF p = %.4f (must exceed 0.001)\n", gof$p.value))
stopifnot(gof$p.value > 0.001)

cat("ALL EQUIVALENCE CHECKS PASS\n")
