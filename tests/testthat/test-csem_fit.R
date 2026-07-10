# Model syntax taken from tests/testthat/test-assess.R (BergamiBagozzi2000)
model_GSCA <- "
# Measurement models
OrgPres <~ cei1 + cei2 + cei3
OrgIden <~ ma1 + ma2 + ma3
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_GSCAM <- "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3
OrgIden =~ ma1 + ma2 + ma3
AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_IGSCA <- "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3
OrgIden <~ ma1 + ma2 + ma3
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

## Single group fits ----------------------------------------------------------
res_gsca_ncmp <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCA,
  .approach_weights = "GSCA",
  .tolerance = 1e-5,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

res_gscam_ncmp <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCAM,
  .approach_weights = "GSCA",
  .tolerance = 1e-5,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

res_igsca_ncmp <- csem(
  .data = BergamiBagozzi2000,
  .model = model_IGSCA,
  .approach_weights = "GSCA",
  .tolerance = 1e-5,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

# Plain GSCA with default modes (CCMP for all composites). CCMP is the case
# where the loading-substitution logic matters.
res_gsca_ccmp <- csem(
  .data = BergamiBagozzi2000,
  .model = model_GSCA,
  .approach_weights = "GSCA",
  .tolerance = 1e-5,
  .conv_criterion = "sum_diff_absolute"
)

## Multigroup fit -------------------------------------------------------------
res_igsca_mg <- csem(
  .data = BergamiBagozzi2000,
  .model = model_IGSCA,
  .approach_weights = "GSCA",
  .id = "gender",
  .tolerance = 1e-5,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

res_single <- list(
  "GSCA_NCMP"  = res_gsca_ncmp,
  "GSCAM_NCMP" = res_gscam_ncmp,
  "IGSCA_NCMP" = res_igsca_ncmp,
  "GSCA_CCMP"  = res_gsca_ccmp
)

## Helpers ---------------------------------------------------------------------

# Indicators belonging to GSCA-M-treated (common factor) blocks; under
# .approach_weights = "GSCA" these are exactly the indicators with unique
# loadings (z2 in Cho et al. 2022), all others are composite indicators (z1).
z2_indicators <- function(res) {
  ct <- res$Information$Model$construct_type
  M  <- res$Information$Model$measurement
  colSums(M[names(ct)[ct == "Common factor"], colnames(M), drop = FALSE]) == 1
}

# Literal, independent implementation of the reproduced covariance matrix of
# Cho et al. (2022), Appendix 2, Eq. (A-11), used as an oracle for fit():
#   Sigma = G (I_K - T)^{-1} E(ee') ((I_K - T)^{-1})' G' + E(ss')
# with T = [0 C; 0 B], G = [I_J, 0], E(ss') = D^2 and the model-implied
# block-diagonal E(ee'): E(zeta zeta') = blockdiag(Phi_exo, diag(1 - R^2)),
# E(eps1 eps1') = within-composite-block residual of S, E(eps2 eps2') = 0.
sigma_cho_A11 <- function(res) {
  mod   <- res$Information$Model
  S     <- res$Estimates$Indicator_VCV
  C     <- t(res$Estimates$Loading_estimates)             # J x P
  Bfull <- res$Estimates$Path_estimates                   # P x P, to x from
  Pmat  <- res$Estimates$Construct_VCV
  J     <- nrow(C); P <- ncol(C); K <- J + P
  cons  <- colnames(C); inds <- rownames(C)

  d <- res$Estimates$Unique_loading_estimates
  if(is.null(d)) d <- setNames(rep(0, J), inds)
  d <- d[inds]; d[is.na(d)] <- 0

  # E(zeta zeta')
  Ezz <- matrix(0, P, P, dimnames = list(cons, cons))
  Ezz[mod$cons_exo, mod$cons_exo] <- Pmat[mod$cons_exo, mod$cons_exo]
  diag(Ezz)[match(mod$cons_endo, cons)] <-
    1 - rowSums(Bfull * Pmat)[mod$cons_endo]

  # E(eps eps'): block-diagonal for composite blocks, zero elsewhere
  IBi   <- solve(diag(P) - Bfull)
  V_eta <- IBi %*% Ezz %*% t(IBi)
  Eee   <- matrix(0, J, J, dimnames = list(inds, inds))
  resid <- S - C %*% V_eta %*% t(C)
  for(p in names(mod$construct_type)[mod$construct_type == "Composite"]) {
    idx <- inds[mod$measurement[p, inds] == 1]
    Eee[idx, idx] <- resid[idx, idx]
  }

  Tmat <- matrix(0, K, K)
  Tmat[1:J, (J + 1):K]       <- C
  Tmat[(J + 1):K, (J + 1):K] <- Bfull
  EeeK <- matrix(0, K, K)
  EeeK[1:J, 1:J]             <- Eee
  EeeK[(J + 1):K, (J + 1):K] <- Ezz
  G    <- cbind(diag(J), matrix(0, J, P))
  ITi  <- solve(diag(K) - Tmat)
  out  <- G %*% ITi %*% EeeK %*% t(ITi) %*% t(G) + diag(d^2)
  dimnames(out) <- dimnames(S)
  out
}

test_that("fit(.type_vcv = 'indicator') returns a proper matrix for GSCA-type objects", {
  for(i in names(res_single)) {
    res <- res_single[[i]]
    S   <- res$Estimates$Indicator_VCV

    Sigma <- fit(res, .saturated = FALSE, .type_vcv = "indicator")

    expect_true(is.matrix(Sigma), info = i)
    expect_identical(dim(Sigma), dim(S))
    expect_identical(dimnames(Sigma), dimnames(S))
    # Symmetry
    expect_equal(Sigma, t(Sigma), tolerance = 1e-10)
    # Diagonal convention (Cho et al. 2022, Eq. A-11): within composite
    # blocks Sigma reproduces S, so composite indicators have unit implied
    # variances. GSCA-M-treated (effect) indicators have implied variance
    # lambda_j^2 + d_j^2 < 1 instead (eps2 = 0) -- the indicator variance
    # left unexplained by the common and unique parts stays visible as
    # misfit.
    z2 <- z2_indicators(res)
    expect_equal(diag(Sigma)[!z2], diag(S)[!z2], tolerance = 1e-12)
    if(any(z2)) {
      expect_true(all(diag(Sigma)[z2] < diag(S)[z2]), info = i)
      # Implied variance of an effect indicator: lambda_j^2 * V(eta_p) + d_j^2
      # with p the indicator's block (V(eta) is not renormalized to unit
      # variances for GSCA-type estimates).
      V       <- fit(res, .saturated = FALSE, .type_vcv = "construct")
      M       <- res$Information$Model$measurement
      blk     <- rownames(V)[apply(M[rownames(V), colnames(S)] == 1, 2, which)]
      lambda2 <- diag(t(res$Estimates$Loading_estimates) %*%
                        res$Estimates$Loading_estimates)
      d       <- res$Estimates$Unique_loading_estimates[colnames(S)]
      expect_equal(diag(Sigma)[z2], (lambda2 * diag(V)[blk] + d^2)[z2],
                   tolerance = 1e-10)
    }
  }
})

test_that("fit() equals the literal Cho et al. (2022, Eq. A-11) reproduced covariance matrix", {
  for(i in names(res_single)) {
    expect_equal(
      fit(res_single[[i]], .saturated = FALSE, .type_vcv = "indicator"),
      sigma_cho_A11(res_single[[i]]),
      tolerance = 1e-12, info = i
    )
  }
  for(g in names(res_igsca_mg)) {
    expect_equal(
      fit(res_igsca_mg[[g]], .saturated = FALSE, .type_vcv = "indicator"),
      sigma_cho_A11(res_igsca_mg[[g]]),
      tolerance = 1e-12, info = g
    )
  }
})

test_that("Saturated implied covariances of common factor indicators follow lambda_i * lambda_j (GSCA-M)", {
  res    <- res_gscam_ncmp
  Sigma  <- fit(res, .saturated = TRUE, .type_vcv = "indicator")
  Lambda <- res$Estimates$Loading_estimates
  P      <- res$Estimates$Construct_VCV

  # All constructs are common factors; no composite block or correlated-error
  # overwrites apply. Within a block the construct variance is 1, hence the
  # implied covariance of two indicators i, j of block k is lambda_ki * lambda_kj.
  expect_equal(
    Sigma["cei1", "cei2"],
    Lambda["OrgPres", "cei1"] * Lambda["OrgPres", "cei2"],
    tolerance = 1e-10
  )
  expect_equal(
    Sigma["ma1", "ma3"],
    Lambda["OrgIden", "ma1"] * Lambda["OrgIden", "ma3"],
    tolerance = 1e-10
  )
  # Cross-block: lambda_i * P[k, l] * lambda_j
  expect_equal(
    Sigma["cei1", "ma1"],
    Lambda["OrgPres", "cei1"] * P["OrgPres", "OrgIden"] * Lambda["OrgIden", "ma1"],
    tolerance = 1e-10
  )
})

test_that("Estimator invariants that fit() relies on hold for GSCA-type objects", {
  # fit() takes Loading_estimates and Unique_loading_estimates as-is: no CCMP
  # loading substitution and no reordering/NA-cleaning of the unique
  # loadings. That is only sound if the estimators guarantee the following.
  for(i in names(res_single)) {
    res <- res_single[[i]]
    L   <- res$Estimates$Loading_estimates
    W   <- res$Estimates$Weight_estimates
    d   <- res$Estimates$Unique_loading_estimates
    ct  <- res$Information$Model$construct_type
    M   <- res$Information$Model$measurement

    # No loading row of a measured construct is all zero (CCMP rows are
    # filled with the implied loadings W %*% S, masked to the blocks, at
    # estimation time).
    expect_true(all(rowSums(L != 0) > 0 | rowSums(W != 0) == 0), info = i)
    if(i == "GSCA_CCMP") {
      expect_equal(L, (W %*% res$Estimates$Indicator_VCV) * M,
                   tolerance = 1e-12, ignore_attr = TRUE)
    }

    # Unique loadings: NULL for plain GSCA; for GSCAm/IGSCA complete, named
    # in indicator order (= colnames(S)), NA-free, and exactly zero for
    # composite indicators.
    if(any(ct == "Common factor")) {
      expect_identical(names(d), colnames(L), label = i)
      expect_false(anyNA(d))
      comp_ind <- colSums(M[names(ct)[ct == "Composite"], , drop = FALSE]) == 1
      expect_true(all(d[comp_ind] == 0), info = i)
    } else {
      expect_null(d, info = i)
    }
  }
})

test_that("CCMP plain GSCA: implied cross-block covariances are nonzero and composite blocks equal S", {
  res   <- res_gsca_ccmp
  S     <- res$Estimates$Indicator_VCV
  Sigma <- fit(res, .saturated = FALSE, .type_vcv = "indicator")

  # Regression test: CCMP composites must not produce all-zero cross-block
  # covariances (their loading rows are fixed to zero during the ALS
  # iterations; the estimators store the implied loadings S %*% w restricted
  # to the block at output, and fit() uses those as-is).
  block_OrgPres <- c("cei1", "cei2", "cei3")
  block_OrgIden <- c("ma1", "ma2", "ma3")
  expect_true(all(abs(Sigma[block_OrgPres, block_OrgIden]) > 0))

  # Within-composite blocks are overwritten with the empirical S
  expect_equal(Sigma[block_OrgPres, block_OrgPres], S[block_OrgPres, block_OrgPres])
  expect_equal(Sigma[block_OrgIden, block_OrgIden], S[block_OrgIden, block_OrgIden])
})

test_that("fit(.type_vcv = 'construct') works for GSCA-type objects", {
  for(i in names(res_single)) {
    res <- res_single[[i]]
    P   <- res$Estimates$Construct_VCV

    vcv_construct <- fit(res, .saturated = FALSE, .type_vcv = "construct")

    expect_true(is.matrix(vcv_construct), info = i)
    expect_identical(dim(vcv_construct), dim(P))
    expect_true(all(rownames(vcv_construct) %in% rownames(P)))
    expect_equal(vcv_construct, t(vcv_construct), tolerance = 1e-10)
    # For GSCA-type estimates the implied construct VCV is not renormalized
    # to unit variances (Cho et al. 2022, Eq. A-11). In these chain models
    # the implied variances equal 1 whenever the path coefficients are OLS
    # with respect to Construct_VCV; that holds up to the ALS convergence
    # tolerance for plain GSCA/IGSCA but only up to ~1e-2 for GSCAm: the
    # GSCA-M U-update has P criterion-flat directions (the informative
    # subspace (I - P_Gamma) Z D has rank J - P only) that the SVD fills
    # arbitrarily, generically not orthogonal to the constant vector, so the
    # scores (Z - UD) %*% W pick up nonzero column means and the centered
    # cor(scores) reported as Construct_VCV differs from the uncentered
    # metric crossprod(Gamma) in which B is exact OLS. Parameters are
    # invariant to the flat directions; only score-based outputs are
    # affected. Tolerance-invariant; see
    # dev/igsca/updateUD/diagnose_centering.R for the full diagnosis.
    expect_equal(unname(diag(vcv_construct)), rep(1, nrow(P)), tolerance = 1e-2)

    # Saturated case: identical to the estimated construct correlation matrix
    vcv_construct_sat <- fit(res, .saturated = TRUE, .type_vcv = "construct")
    expect_equal(vcv_construct_sat, P)
  }
})

test_that("Estimates$Indicator_VCV and Estimates$Construct_VCV are proper matrices", {
  n_ind <- 12L
  n_con <- 4L
  for(i in names(res_single)) {
    res <- res_single[[i]]
    expect_true(is.matrix(res$Estimates$Indicator_VCV), info = i)
    expect_identical(dim(res$Estimates$Indicator_VCV), c(n_ind, n_ind))
    expect_true(is.matrix(res$Estimates$Construct_VCV), info = i)
    expect_identical(dim(res$Estimates$Construct_VCV), c(n_con, n_con))
  }
  for(g in res_igsca_mg) {
    expect_true(is.matrix(g$Estimates$Indicator_VCV))
    expect_identical(dim(g$Estimates$Indicator_VCV), c(n_ind, n_ind))
    expect_true(is.matrix(g$Estimates$Construct_VCV))
    expect_identical(dim(g$Estimates$Construct_VCV), c(n_con, n_con))
  }
})

test_that("fit() returns a per-group list of matrices for multigroup GSCA-type objects", {
  Sigma_list <- fit(res_igsca_mg, .saturated = FALSE, .type_vcv = "indicator")

  expect_type(Sigma_list, "list")
  expect_length(Sigma_list, length(res_igsca_mg))
  expect_identical(names(Sigma_list), names(res_igsca_mg))

  for(g in names(Sigma_list)) {
    S  <- res_igsca_mg[[g]]$Estimates$Indicator_VCV
    z2 <- z2_indicators(res_igsca_mg[[g]])
    expect_true(is.matrix(Sigma_list[[g]]))
    expect_identical(dim(Sigma_list[[g]]), dim(S))
    # Composite indicators reproduce diag(S); GSCA-M-treated ones fall short
    # of it (Cho et al. 2022, Eq. A-11; see the single-group test above).
    expect_equal(diag(Sigma_list[[g]])[!z2], diag(S)[!z2], tolerance = 1e-12)
    expect_true(all(diag(Sigma_list[[g]])[z2] < diag(S)[z2]))
  }

  vcv_list <- fit(res_igsca_mg, .saturated = FALSE, .type_vcv = "construct")
  expect_type(vcv_list, "list")
  expect_length(vcv_list, length(res_igsca_mg))
})

test_that("fit() results for PLS objects are unchanged (regression guard)", {
  model_pls <- "
  eta1 =~ y11 + y12 + y13
  eta2 =~ y21 + y22 + y23
  eta3 =~ y31 + y32 + y33

  eta2 ~ eta1
  eta3 ~ eta1 + eta2
  "
  res_pls <- csem(.data = threecommonfactors, .model = model_pls)

  F1 <- fit(res_pls, .saturated = FALSE, .type_vcv = "indicator")
  F2 <- fit(res_pls, .saturated = TRUE,  .type_vcv = "indicator")
  F3 <- fit(res_pls, .saturated = FALSE, .type_vcv = "construct")

  # Reference values computed with the implementation before the GSCA
  # branch was added to fit() (commit cb0f0a56). These must stay
  # bit-identical for PLS objects.
  expect_equal(F1["y11", "y21"], 0.22989519320542023, tolerance = 1e-15)
  expect_equal(F1["y12", "y33"], 0.32192169958553024, tolerance = 1e-15)
  expect_equal(F1["y11", "y12"], 0.43051663355925296, tolerance = 1e-15)
  expect_equal(F2["y11", "y21"], 0.22989519320542023, tolerance = 1e-15)
  expect_equal(F2["y23", "y31"], 0.40305034442955345, tolerance = 1e-15)
  expect_equal(F3["eta1", "eta3"], 0.66336491778896245, tolerance = 1e-15)
  expect_equal(sum(F1), 35.523004770798586, tolerance = 1e-15)
  expect_equal(sum(F2), 35.523004770798586, tolerance = 1e-15)
  expect_equal(sum(F3), 6.8953207750622107, tolerance = 1e-15)
})
