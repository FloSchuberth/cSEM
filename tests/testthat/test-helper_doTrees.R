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

# The argument list doTrees() replays at every node, as csem() stores it.
args_trees <- csem(
  .data = dat,
  .model = model_trees,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = "sum_diff_absolute",
  .iter_max = 100,
  .GSCA_modes = "CCMP",
  .tolerance = 0.0001
)$Information$Arguments

# One pooled node fit, packaged the way the trafo hands it to a kernel.
node_fixture <- function(max_cuts = 8L, minbucket = 100L) {
  ind <- parseModel(args_trees$.model)$indicators
  ft <- try_fit(dat[, ind, drop = FALSE], args_trees)
  list(
    model = list(object = ft$fit),
    mf = dat,
    subset = seq_len(nrow(dat)),
    ctrl = list(
      collector = new_collector(),
      args = args_trees,
      indicators = ind,
      max_cuts = max_cuts,
      minbucket = minbucket
    )
  )
}

test_that("fit_csem() replays the arguments of the fit it was given", {
  ind <- parseModel(args_trees$.model)$indicators
  node <- fit_csem(dat[1:500, ind, drop = FALSE], args_trees)
  used <- node$Information$Arguments
  expect_identical(used$.approach_weights, "GSCA")
  expect_identical(used$.GSCA_modes, args_trees$.GSCA_modes)
  expect_identical(used$.tolerance, args_trees$.tolerance)
  # .id is dropped, not passed as NULL, so csem()'s own default applies
  expect_s3_class(node, "cSEMResults_default")

  # The MGA fit the split kernels compare is the same call with .id set
  g <- cbind(
    group = factor(rep(1:2, each = 500L)),
    dat[, ind, drop = FALSE]
  )
  expect_s3_class(fit_csem(g, args_trees, .id = "group"), "cSEMResults_multi")
})

# A kernel that throws is caught by argmax_split() and turned into NA, which
# partykit reads as "no admissible split" -- so at the tree level a broken kernel
# and a genuinely unsplittable node look identical. Called directly, a kernel
# either returns a finite statistic or it does not, and the failure names it.
test_that("every split kernel returns a finite statistic", {
  fx <- node_fixture()
  # A median split on a noise covariate: no group difference is expected here,
  # but every kernel must still put a number on it.
  goes_left <- fx$mf$noise_1 <= stats::median(fx$mf$noise_1)
  for (k in c("FITdiff", "DLi", "DGi")) {
    expect_true(
      is.finite(
        partition_stat(stat_kind = k, 
          model = fx$model, 
          mf = fx$mf,
          subset = fx$subset,
          goes_left = goes_left,
          ctrl = fx$ctrl)
      ),
      info = k
    )
  }
})

test_that("partition_stat() records a failed auxiliary fit instead of throwing", {
  local_mocked_bindings(
    try_fit = function(.data, .args, .id = NULL) list(fit = NULL, ok = FALSE)
  )
  coll <- new_collector()
  ctrl <- list(collector = coll, args = args_trees, indicators = "x11")
  goes_left <- rep(c(TRUE, FALSE), length.out = nrow(dat))

  # An unfittable candidate partition must come back as NA, not an exception:
  # anything that escapes here would make SimDesign redraw the whole rep.
  expect_identical(
    partition_stat(
      stat_kind = "DLi",
      model = list(object = NULL),
      mf = dat,
      subset = seq_len(nrow(dat)),
      goes_left = goes_left,
      ctrl = ctrl
    ),
    NA_real_
  )
  expect_identical(coll$n_fail_candidate, 1L)
})

# This is where cutpoint choice has to be tested: at the tree level the root
# offers no freedom, since z_true is a 2-level factor and candidate_partitions()
# returns exactly one candidate for it. noise_1 is numeric and yields max_cuts
# of them, so a kernel that runs but picks the wrong cutpoint is visible here
# and nowhere else.
test_that("argmax_split() returns an admissible cutpoint the kernels disagree on", {
  fx <- node_fixture()
  j <- match("noise_1", names(fx$mf)) # column number of noise_1 in the data
  z <- fx$mf[[j]]
  # Function for getting the cutpoint in mf$j to split the data on
  cutpoint <- function(kern) {
    sp <- argmax_split(
      splitter = kern,
      collector = fx$ctrl$collector,
      model = fx$model,
      mf = fx$mf,
      subset = fx$subset,
      whichvar = j,
      ctrl = fx$ctrl
    )
    # NULL is what argmax_split() returns when no candidate produced a finite
    # statistic -- the dead kernel this layer exists to localise.
    expect_true(inherits(sp, "partysplit"))
    partykit::breaks_split(sp)
  }

  brs <- vapply(
    list(FIT = split_max_fitdiff, DLi = split_max_dli, DGi = split_max_dgi),
    cutpoint,
    numeric(1)
  )
  for (k in names(brs)) {
    # Strictly inside the observed range: an endpoint sends every row one way,
    # which partykit accepts and grows a degenerate child from.
    expect_true(
      is.finite(brs[[k]]) && brs[[k]] > min(z) && brs[[k]] < max(z),
      info = k
    )
  }

  # The kernels are only worth pairing with a selector if they disagree about
  # where to cut. Asserted for FIT vs DLi only: DLi and DGi are both distances
  # between the same model-implied indicator VCVs and may legitimately land on
  # the same candidate, so distinctness across all three would be an assertion
  # about this fixture rather than about the kernels.
  expect_false(isTRUE(all.equal(brs[["FIT"]], brs[["DLi"]])))
})
