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

# doTrees node-level helpers ---------------------------------------------
# These call the partition machinery directly, one node at a time. At the tree
# level the root gives a kernel no freedom -- z_true is a 2-level factor, so
# candidate_partitions() returns exactly one candidate -- so this is the only
# layer where cutpoint choice and the scan cache can be observed at all.
load(testthat::test_path("data/igscaTrees.Rdata")) # Creates dat

model_trees <- "# Latent variable model
 eta1 =~ x11 + x12 + x13
 eta2 =~ x21 + x22 + x23

 # Composite model
 eta3 <~ x31 + x32 + x33
 eta4 <~ x41 + x42 + x43 + x44

 # Structural model
 eta4 ~ eta3 + eta1
 eta2 ~ eta1 + eta4 + eta3
 "

# One pooled node fit, packaged the way the trafo hands it to a kernel.
node_fixture <- function(max_cuts = 8L, minbucket = 100L) {
  ind <- parseModel(model_trees)$indicators
  ft <- try_fit(dat[, ind, drop = FALSE], model_trees)
  list(
    model = list(object = ft$fit),
    mf = dat,
    subset = seq_len(nrow(dat)),
    ctrl = list(
      collector = new_collector(),
      model = model_trees,
      indicators = ind,
      max_cuts = max_cuts,
      minbucket = minbucket
    )
  )
}

test_that("partition_stat() records a failed auxiliary fit instead of throwing", {
  local_mocked_bindings(
    try_fit = function(.data, .model, .id = NULL) list(fit = NULL, ok = FALSE)
  )
  coll <- new_collector()
  ctrl <- list(collector = coll, model = model_trees, indicators = "x11")
  goes_left <- rep(c(TRUE, FALSE), length.out = nrow(dat))

  # An unfittable candidate partition must come back as NA, not an exception:
  # anything that escapes here would make SimDesign redraw the whole rep.
  expect_identical(
    partition_stat(
      "DLi", list(object = NULL), dat, seq_len(nrow(dat)), goes_left, ctrl
    ),
    NA_real_
  )
  expect_identical(coll$n_fail_resample, 1L)
})

# NB: this must not be combined with local_mocked_bindings() on the kernels --
# the cache key is closure identity, and a mock would change it.
test_that("argmax_split() reuses a matched scan and rescans a mismatched one", {
  fx <- node_fixture()
  j <- match("noise_1", names(fx$mf))
  coll <- fx$ctrl$collector

  # Arm the cache the way partition_select() does before handing control back
  # to the engine's splitfun.
  sc <- scan_covariate("FITdiff", j, fx$model, fx$mf, fx$subset, fx$ctrl)
  coll$scan <- stats::setNames(list(sc), as.character(j))
  coll$scan_subset <- fx$subset
  coll$scan_splitter <- split_max_fitdiff

  # Matched pair: the selector's argmax comes back untouched, and the kernel is
  # never evaluated -- an unchanged n_split_scan is what "cache hit" means.
  expect_identical(
    argmax_split(
      split_max_fitdiff, coll, fx$model, fx$mf, fx$subset, j, fx$ctrl
    ),
    sc$split
  )
  expect_identical(coll$n_split_scan, 0L)

  # Mismatched pair: a different kernel must rescan with its own statistic ...
  sp_dli <- argmax_split(
    split_max_dli, coll, fx$model, fx$mf, fx$subset, j, fx$ctrl
  )
  expect_s3_class(sp_dli, "partysplit")
  expect_identical(coll$n_split_scan, 1L)
  # ... and DLi genuinely disagrees with FIT about where to cut noise_1.
  expect_false(isTRUE(all.equal(
    partykit::breaks_split(sp_dli),
    partykit::breaks_split(sc$split)
  )))

  # The other half of the key: same kernel, different node.
  coll$n_split_scan <- 0L
  argmax_split(
    split_max_fitdiff, coll, fx$model, fx$mf, fx$subset[-1L], j, fx$ctrl
  )
  expect_identical(coll$n_split_scan, 1L)
})
