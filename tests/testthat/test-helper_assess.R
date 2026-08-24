# Tests for the delta-method (asymptotic) confidence interval of the HTMT:
#   calculateHTMT(.inference = "asymptotic")
#     -> calculateHTMTcore()            (HTMT + analytic gradients)
#     -> calculateHTMTasymptoticSE()    (delta-method SE)
#     -> calculateCorVCV()              (Steiger & Hakstian 1982 acov of correlations)
#
# The pipeline is validated on three levels:
#   1. the algebra of each building block (numerical gradients, normal-theory
#      closed form for the correlation acov),
#   2. deterministic end-to-end invariances (CI assembly, affine transformations,
#      indicator ordering, submodel consistency),
#   3. statistical validity (Monte-Carlo coverage; skipped on CRAN).

model_3cf <- "
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
eta2 ~ eta1
eta3 ~ eta1 + eta2
"

res_3cf <- csem(threecommonfactors, model_3cf)

blocks_3cf <- list(
  eta1 = c("y11", "y12", "y13"),
  eta2 = c("y21", "y22", "y23"),
  eta3 = c("y31", "y32", "y33")
)
pairs_3cf <- utils::combn(names(blocks_3cf), 2, simplify = FALSE)

# Independent reference implementation: HTMT of one construct pair from a full
# indicator correlation matrix (arithmetic-mean version of Henseler et al. 2015).
htmt_ref <- function(S, A, B) {
  SA <- S[A, A]
  SB <- S[B, B]
  mean(S[A, B]) / sqrt(mean(SA[lower.tri(SA)]) * mean(SB[lower.tri(SB)]))
}

ht_asy <- calculateHTMT(res_3cf,
                        .type_htmt          = "htmt",
                        .absolute           = FALSE,
                        .alpha              = 0.05,
                        .inference = "asymptotic")

test_that("analytic HTMT gradients equal numerical gradients", {
  core <- cSEM:::calculateHTMTcore(res_3cf,
                                   .type_htmt           = "htmt",
                                   .absolute            = FALSE,
                                   .asymptotic          = TRUE,
                                   .only_common_factors = TRUE)
  grads <- attr(core, "gradients")
  expect_named(grads, c("eta1__eta2", "eta1__eta3", "eta2__eta3"))
  
  S    <- res_3cf$Estimates$Indicator_VCV
  ind  <- rownames(S)
  lt   <- lower.tri(S)
  rows <- ind[row(S)[lt]]
  cols <- ind[col(S)[lt]]
  h    <- 1e-6
  
  for (p in seq_along(pairs_3cf)) {
    A <- blocks_3cf[[pairs_3cf[[p]][1]]]
    B <- blocks_3cf[[pairs_3cf[[p]][2]]]
    g_num <- vapply(seq_along(rows), function(k) {
      Sp <- S; Sm <- S
      Sp[rows[k], cols[k]] <- Sp[cols[k], rows[k]] <- S[rows[k], cols[k]] + h
      Sm[rows[k], cols[k]] <- Sm[cols[k], rows[k]] <- S[rows[k], cols[k]] - h
      (htmt_ref(Sp, A, B) - htmt_ref(Sm, A, B)) / (2 * h)
    }, numeric(1))
    expect_equal(unname(grads[[p]]), g_num, tolerance = 1e-6)
  }
})

test_that("calculateCorVCV matches the normal-theory acov of correlations", {
  # Under multivariate normality the Steiger-Hakstian (1982) distribution-free
  # estimator must converge to the Pearson-Filon normal-theory covariance
  # (Steiger 1980, Eq. 2). Anchors the estimator to the literature rather than
  # to a re-implementation of itself.
  nt_cov <- function(R, j, k, h, m) {
    0.5 * R[j, k] * R[h, m] * (R[j, h]^2 + R[j, m]^2 + R[k, h]^2 + R[k, m]^2) +
      R[j, h] * R[k, m] + R[j, m] * R[k, h] -
      R[j, k] * (R[j, h] * R[j, m] + R[k, h] * R[k, m]) -
      R[h, m] * (R[j, h] * R[k, h] + R[j, m] * R[k, m])
  }
  Rpop <- matrix(c(1, .5, .3,
                   .5, 1, .4,
                   .3, .4, 1), 3, 3)
  set.seed(123)
  n <- 100000
  X <- MASS::mvrnorm(n, mu = rep(0,3), Sigma = Rpop, empirical = TRUE)
  #X <- matrix(rnorm(n * 3), n, 3) %*% chol(Rpop)
  
  # calculateCorVCV() assumes standardized input (csem standardizes internally)
  Sig_hat <- cSEM:::calculateCorVCV(scale(X)) * n
  
  # correlations in the pipeline's ordering: (1,2), (1,3), (2,3)
  idx    <- which(lower.tri(diag(3)), arr.ind = TRUE)
  Sig_nt <- outer(1:3, 1:3, Vectorize(function(a, b)
    nt_cov(Rpop, idx[a, 1], idx[a, 2], idx[b, 1], idx[b, 2])))
  
  expect_true(max(abs(Sig_hat - Sig_nt)) < 0.005)
  
  # basic properties on real (standardized) model data
  Sig <- cSEM:::calculateCorVCV(res_3cf$Information$Data)
  expect_true(isSymmetric(Sig))
  expect_true(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values) >
                -1e-12)
})

test_that("CI equals estimate +/- z * sqrt(t(g) Sigma g) and layout is correct", {
  S   <- res_3cf$Estimates$Indicator_VCV
  Sig <- cSEM:::calculateCorVCV(res_3cf$Information$Data)
  core <- cSEM:::calculateHTMTcore(res_3cf,
                                   .type_htmt           = "htmt",
                                   .absolute            = FALSE,
                                   .asymptotic          = TRUE,
                                   .only_common_factors = TRUE)
  grads <- attr(core, "gradients")
  
  est <- vapply(pairs_3cf, function(p)
    htmt_ref(S, blocks_3cf[[p[1]]], blocks_3cf[[p[2]]]), numeric(1))
  se  <- unname(vapply(grads, function(g) sqrt(drop(t(g) %*% Sig %*% g)),
                       numeric(1)))
  z   <- qnorm(1 - 0.05)
  
  # quantiles is 2 x 9 (column-major over the 3x3 HTMT matrix); the pairs
  # (eta1,eta2), (eta1,eta3), (eta2,eta3) sit in cells (2,1), (3,1), (3,2),
  # i.e. columns 2, 3 and 6
  qcols <- c(2, 3, 6)
  expect_equal(rownames(ht_asy$quantiles), c("90%L", "90%U"))
  expect_equal(unname(ht_asy$quantiles[1, qcols]), est - z * se,
               tolerance = 1e-10)
  expect_equal(unname(ht_asy$quantiles[2, qcols]), est + z * se,
               tolerance = 1e-10)
  
  # point estimates in the lower triangle, upper CI bound (HTMT > 0 here) in
  # the upper triangle, unit diagonal
  expect_equal(unname(ht_asy$htmts[lower.tri(ht_asy$htmts)]), est,
               tolerance = 1e-10)
  expect_equal(unname(t(ht_asy$htmts)[lower.tri(ht_asy$htmts)]), est + z * se,
               tolerance = 1e-10)
  expect_equal(unname(diag(ht_asy$htmts)), rep(1, 3))
  expect_null(ht_asy$nr_admissibles)
  expect_null(attr(ht_asy$htmts, "gradients"))  # consumed, not leaked
})

test_that("the CI is invariant to affine transformations of the indicators", {
  res_t <- csem(threecommonfactors * 3 + 5, model_3cf)
  ht_t  <- calculateHTMT(res_t,
                         .type_htmt          = "htmt",
                         .absolute           = FALSE,
                         .alpha              = 0.05,
                         .inference = "asymptotic")
  expect_equal(ht_t$quantiles, ht_asy$quantiles, tolerance = 1e-10)
})

test_that("the CI is invariant to indicator and block ordering", {
  # Reordering blocks/indicators permutes the correlation vector, so any
  # misalignment between the gradient ordering and the ordering of the
  # correlation acov matrix would change t(g) Sigma g.
  perm_cols  <- c("y31", "y32", "y33", "y11", "y12", "y13", "y21", "y22", "y23")
  model_perm <- "
  eta3 =~ y31 + y32 + y33
  eta1 =~ y11 + y12 + y13
  eta2 =~ y21 + y22 + y23
  eta2 ~ eta1
  eta3 ~ eta1 + eta2
  "
  res_p <- csem(threecommonfactors[, perm_cols], model_perm)
  ht_p  <- calculateHTMT(res_p,
                         .type_htmt          = "htmt",
                         .absolute           = FALSE,
                         .alpha              = 0.05,
                         .inference = "asymptotic")
  cn <- colnames(ht_asy$htmts)
  for (p in pairs_3cf) {
    cell   <- which(rownames(ht_asy$htmts) == p[2]) +
      (which(cn == p[1]) - 1) * nrow(ht_asy$htmts)
    cell_p <- which(rownames(ht_p$htmts) == p[2]) +
      (which(colnames(ht_p$htmts) == p[1]) - 1) * nrow(ht_p$htmts)
    expect_equal(unname(ht_p$quantiles[, cell_p]),
                 unname(ht_asy$quantiles[, cell]), tolerance = 1e-10)
  }
})

test_that("the CI of a pair is unaffected by indicators outside the pair", {
  # The gradient is zero for correlations not involving the pair's blocks, so
  # the delta-method variance must reduce exactly to the two-construct submodel.
  model_12 <- "
  eta1 =~ y11 + y12 + y13
  eta2 =~ y21 + y22 + y23
  eta2 ~ eta1
  "
  res_12 <- csem(threecommonfactors[, 1:6], model_12)
  ht_12  <- calculateHTMT(res_12,
                          .type_htmt          = "htmt",
                          .absolute           = FALSE,
                          .alpha              = 0.05,
                          .inference = "asymptotic")
  expect_equal(unname(ht_12$quantiles[, 2]), unname(ht_asy$quantiles[, 2]),
               tolerance = 1e-10)
})

test_that("guards of the asymptotic path fire", {
  # asymptotic CIs are not available for the HTMT2
  expect_warning(
    calculateHTMT(res_3cf, .type_htmt = "htmt2", .absolute = FALSE,
                  .inference = "asymptotic"),
    "not available for the HTMT2"
  )
  # a supplied .ci is ignored with a warning
  expect_warning(
    calculateHTMT(res_3cf, .type_htmt = "htmt", .absolute = FALSE,
                  .ci = "CI_percentile",
                  .inference = "asymptotic"),
    "ignored"
  )
  # .absolute = TRUE triggers the recommendation warning
  expect_warning(
    calculateHTMT(res_3cf, .type_htmt = "htmt", .absolute = TRUE,
                  .inference = "asymptotic"),
    "absolute"
  )
  # only a single .alpha is accepted
  expect_error(
    calculateHTMT(res_3cf, .type_htmt = "htmt", .absolute = FALSE,
                  .alpha = c(0.05, 0.10),
                  .inference = "asymptotic" ),
    "single numeric probability"
  )
})

test_that("the delta-method CI attains nominal coverage (Monte Carlo)", {
  skip_on_cran()
  
  # Tau-equivalent two-factor population (all loadings 0.7, factor correlation
  # 0.5): the population HTMT equals the factor correlation, so the true value
  # is known exactly. 400 replications, nominal two-sided level 90%; the
  # acceptance band is the nominal level +/- ~3 binomial standard errors.
  lam <- 0.7
  rho <- 0.5
  L   <- rbind(c(lam, lam, lam, 0, 0, 0),
               c(0, 0, 0, lam, lam, lam))
  Phi <- matrix(c(1, rho, rho, 1), 2, 2)
  Sig_pop <- t(L) %*% Phi %*% L
  diag(Sig_pop) <- 1
  ch <- chol(Sig_pop)
  
  model_2cf <- "
  eta1 =~ y1 + y2 + y3
  eta2 =~ y4 + y5 + y6
  eta2 ~ eta1
  "
  
  set.seed(2026)
  R_mc <- 400
  n_mc <- 300
  cover <- logical(R_mc)
  ests  <- ses <- numeric(R_mc)
  for (b in seq_len(R_mc)) {
    d <- matrix(rnorm(n_mc * 6), n_mc, 6) %*% ch
    colnames(d) <- paste0("y", 1:6)
    r  <- csem(as.data.frame(d), model_2cf)
    ht <- calculateHTMT(r, .type_htmt = "htmt", .absolute = FALSE,
                        .alpha = 0.05, .inference = "asymptotic")
    lo <- ht$quantiles[1, 2]
    hi <- ht$quantiles[2, 2]
    ests[b]  <- ht$htmts[2, 1]
    ses[b]   <- (hi - lo) / (2 * qnorm(0.95))
    cover[b] <- (lo <= rho) && (rho <= hi)
  }
  
  # coverage close to the nominal 90%
  expect_true(mean(cover) > 0.85 && mean(cover) < 0.95)
  # the HTMT estimates center on the true value
  expect_equal(mean(ests), rho, tolerance = 0.02)
  # the delta-method SE tracks the Monte-Carlo sd of the estimates
  expect_equal(mean(ses), sd(ests), tolerance = 0.1 * sd(ests))
})
