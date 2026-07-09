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
    # Unit diagonal convention: implied variances equal the empirical ones
    expect_equal(diag(Sigma), diag(S), tolerance = 1e-8)
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

test_that("CCMP plain GSCA: implied cross-block covariances are nonzero and composite blocks equal S", {
  res   <- res_gsca_ccmp
  S     <- res$Estimates$Indicator_VCV
  Sigma <- fit(res, .saturated = FALSE, .type_vcv = "indicator")

  # Regression test: CCMP composites must not produce all-zero cross-block
  # covariances (loading rows fixed to zero during estimation must be
  # replaced by implied loadings S %*% w restricted to the block).
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
    expect_equal(unname(diag(vcv_construct)), rep(1, nrow(P)), tolerance = 1e-10)

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
    S <- res_igsca_mg[[g]]$Estimates$Indicator_VCV
    expect_true(is.matrix(Sigma_list[[g]]))
    expect_identical(dim(Sigma_list[[g]]), dim(S))
    expect_equal(diag(Sigma_list[[g]]), diag(S), tolerance = 1e-8)
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
