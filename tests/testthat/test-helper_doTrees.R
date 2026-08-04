# Tests for R/helper_doTrees.R. 
model_pls <- "
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33

eta2 ~ eta1
eta3 ~ eta1 + eta2
"

res_pls <- csem(.data = threecommonfactors, .model = model_pls)

test_that("bdiagFit() repeats a single-group implied VCV block-diagonally", {
  for(type in c("indicator", "construct")) {
    Sigma <- fit(res_pls, .saturated = FALSE, .type_vcv = type)
    B3    <- bdiagFit(res_pls, .n_blocks = 3, .type_vcv = type)

    # One comparison covers dimensions, block content, and zero off-blocks
    expect_equal(unname(B3), kronecker(diag(3), unname(Sigma)))
    # Names: block-suffixed and unique
    expect_identical(rownames(B3)[seq_len(nrow(Sigma))],
                     paste0(rownames(Sigma), "_1"))
    expect_identical(anyDuplicated(rownames(B3)), 0L)
  }
})

test_that("bdiagFit() block-diagonalizes multigroup fits in group order", {
  # Model syntax and csem() arguments as in tests/testthat/test-csem_fit.R
  model_IGSCA <- "
  OrgPres =~ cei1 + cei2 + cei3
  OrgIden <~ ma1 + ma2 + ma3
  AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
  AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

  OrgIden ~ OrgPres
  AffLove ~ OrgIden
  AffJoy  ~ OrgIden"

  res_mg <- csem(
    .data = BergamiBagozzi2000,
    .model = model_IGSCA,
    .approach_weights = "GSCA",
    .id = "gender",
    .tolerance = 1e-4,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "NCMP",
    .iter_max = 600
  )

  Sigma_list <- fit(res_mg, .saturated = FALSE, .type_vcv = "indicator")
  B          <- bdiagFit(res_mg, .type_vcv = "indicator")

  expect_identical(dim(B), rep(sum(sapply(Sigma_list, nrow)), 2))
  pos <- 0
  for(g in names(Sigma_list)) {
    idx <- pos + seq_len(nrow(Sigma_list[[g]]))
    expect_equal(B[idx, idx], Sigma_list[[g]], ignore_attr = TRUE)
    pos <- pos + nrow(Sigma_list[[g]])
  }

  # .n_blocks is ignored (with a warning) for multigroup objects
  expect_warning(bdiagFit(res_mg, .n_blocks = 2, .type_vcv = "indicator"))
})

test_that("bdiagFit() validates .n_blocks", {
  expect_error(bdiagFit(res_pls, .n_blocks = 0))
  expect_error(bdiagFit(res_pls, .n_blocks = 1.5))
})
