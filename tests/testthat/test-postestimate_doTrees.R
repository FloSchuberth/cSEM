# doTrees ----------------------------------------------------------------
load(testthat::test_path("data/igscaTrees.Rdata")) # Creates dat
covs <- c("z_true", grep("^noise_", names(dat), value = TRUE))

model <- "# Latent variable model
 eta1 =~ x11 + x12 + x13
 eta2 =~ x21 + x22 + x23

 # Composite model
 eta3 <~ x31 + x32 + x33
 eta4 <~ x41 + x42 + x43 + x44

 # Structural model
 eta4 ~ eta3 + eta1
 eta2 ~ eta1 + eta4 + eta3
 "

# doTrees() reads its data, model and estimator settings off this fit and
# replays them at every node, so the covariates ride along in the csem() call --
# non-indicator columns are ignored when fitting.
res <- csem(
  .data = dat,
  .model = model,
  .approach_weights = "GSCA",
  .disattenuate = TRUE,
  .conv_criterion = "sum_diff_absolute",
  .iter_max = 100,
  .GSCA_modes = "CCMP",
  .tolerance = 0.0001
)

# Controls ---------------------------------------------------------------
# The two influence families cost an order of magnitude apart, so they get
# separate budgets. The conditional-test family ("mat"/"vec") inherits COIN
# variable selection from libcoin and is cheap even at the defaults; a mixed
# splitter only re-scans candidate cutpoints, which adds a second or two. The
# partition family ("FIT"/"DLi"/"DGi") runs an R_test permutation test per
# covariate per node and every permutation is a two-group MGA fit, so at the
# defaults a single call takes about an hour and a half.
#
# R_test is inert here: it reaches the conditional-test family only as libcoin's
# nresample, which the default coin_distribution = "asymptotic" never draws. It
# is left in place so the setting still applies if a test below opts back into
# "approximate".
ctl_mixed <- igsca_tree_control(R_test = 50L, maxdepth = 2L, max_cuts = 8L)

# R_test >= 20L is a correctness floor, not just a budget: permutation_pvalue()
# returns (1 + k) / (R + 1), so the smallest attainable p-value is 1 / (R + 1),
# and the split criterion is a strict inequality against alpha = 0.05. At
# R_test = 19L that floor is exactly 0.05 and no node can ever split, whatever
# the data says. 24L keeps headroom at 1/25 = 0.04. maxdepth = 1L costs nothing
# in coverage on this fixture -- measured, the children do not split at
# maxdepth = 2L either, they only double the runtime.
#
# bonferroni is off here against the package default, and that is a budget
# decision rather than a statistical one. The engine applies Šidák across the
# k = 3 covariates scanned at each node, which raises the floor above from
# R_test >= 20L to R_test >= 59L -- 1 - (1 - 1/25)^3 = 0.115 does not clear
# alpha, so every partition test below would assert against a stump. Paying for
# it means 59 two-group MGA refits per covariate per node instead of 24. The
# adjustment is exercised on the conditional-inference path instead, where a
# resample costs a permutation of a precomputed matrix rather than a model fit.
ctl_part <- igsca_tree_control(
  R_test = 24L,
  max_cuts = 3L,
  maxdepth = 1L,
  minbucket = 200L,
  bonferroni = FALSE
)

grow_tree <- function(influence, splitter, control, object = res) {
  doTrees(
    .object = object,
    .covariates = covs,
    .influence = influence,
    .splitter = splitter,
    .control = control
  )
}

# A returned tree is only evidence of a working configuration if it actually
# split. argmax_split() turns a throwing kernel into NA and partykit reads that
# as "no admissible split", so a broken splitter yields a clean-looking stump
# rather than an error -- which is why absence of an exception asserts nothing
# here.
expect_grew <- function(tree) {
  expect_s3_class(tree, "igsca_tree")
  info <- attr(tree, "igsca_info")
  expect_identical(info$n_fail_full, 0L)
  expect_identical(info$n_fail_split, 0L)
  # width() returns a double, so expect_identical(..., 1L) would fail on type.
  expect_gt(partykit::width(partykit::node_party(tree)), 1)
}



# Matrix of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Matrix of Residuals Runs as expected", {
  set.seed(12353)
  # The native path is cheap enough to snapshot at the package defaults.
  trees_mx <- grow_tree("mat", "native", igsca_tree_control())
  expect_grew(trees_mx)
  expect_snapshot(trees_mx)
})

for (sp in c("FIT", "DLi", "DGi")) {
  test_that(paste0("Matrix of Residuals selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("mat", sp, ctl_mixed))
  })
}



# Vector of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Vector of Residuals Runs as expected", {
  set.seed(12353)
  trees_vec <- grow_tree("vec", "native", igsca_tree_control())
  expect_grew(trees_vec)
  expect_snapshot(trees_vec)
})

for (sp in c("FIT", "DLi", "DGi")) {
  test_that(paste0("Vector of Residuals selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("vec", sp, ctl_mixed))
  })
}



# NPT-FIT Split ---------------------------------------------------------------
test_that("IGSCA Trees Variable Selection on FIT Runs as expected", {
  set.seed(12353)
  trees_NPT_FIT <- grow_tree("FIT", "FIT", ctl_part)
  expect_grew(trees_NPT_FIT)
  expect_snapshot(trees_NPT_FIT)
})

for (sp in c("DLi", "DGi")) {
  test_that(paste0("FIT selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("FIT", sp, ctl_part))
  })
}

# dLi Split --------------------------------------------------------------
test_that("IGSCA Trees Conditional Test on Squared-Euclidean Distance Runs as expected", {
  set.seed(12353)
  trees_DLi_FIT <- grow_tree("DLi", "FIT", ctl_part)
  expect_grew(trees_DLi_FIT)
  expect_snapshot(trees_DLi_FIT)
})

for (sp in c("DLi", "DGi")) {
  test_that(paste0("DLi selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("DLi", sp, ctl_part))
  })
}

# dGi Split --------------------------------------------------------------
test_that("IGSCA Trees Conditional Test on Geodesic Distance Runs as expected", {
  set.seed(12353)
  trees_DGi_FIT <- grow_tree("DGi", "FIT", ctl_part)
  expect_grew(trees_DGi_FIT)
  expect_snapshot(trees_DGi_FIT)
})

for (sp in c("DLi", "DGi")) {
  test_that(paste0("DGi selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree("DGi", sp, ctl_part))
  })
}

# Dead split kernel ------------------------------------------------------
# argmax_split() turns a throwing kernel into NA, which partykit reads as "no
# admissible split" -- so a broken kernel yields a clean-looking stump. These
# tests pin the diagnostic that separates that from a genuinely unsplittable
# node.
test_that("a dead split kernel is reported, not silently a stump", {
  local_mocked_bindings(
    split_max_dli = function(model, mf, subset, goes_left, ctrl) {
      stop("kernel is broken")
    }
  )
  set.seed(11)
  expect_warning(
    tr <- grow_tree("mat", "DLi", ctl_mixed),
    "produced no usable statistic"
  )
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_fail_split, 0L)
  expect_identical(info$n_fail_split, info$n_split_scan)
  # The tree really is the stump the diagnostic is warning about
  expect_equal(partykit::width(partykit::node_party(tr)), 1)
})

test_that("a working split kernel records scans but no failures", {
  set.seed(11)
  tr <- expect_no_warning(grow_tree("mat", "DLi", ctl_mixed))
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_split_scan, 0L)
  expect_grew(tr)
})

test_that("the native split path never touches the kernel counters", {
  set.seed(11)
  tr <- grow_tree("mat", "native", ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_split_scan, 0L)
  expect_identical(info$n_fail_split, 0L)
})

# Input contract ---------------------------------------------------------
# doTrees() takes a fitted object rather than data + model, so everything it
# needs is either derivable from that fit or a mistake worth naming.

test_that("doTrees() refits nodes with the estimator the fit used", {
  # Both branches: they build their node fits through separate trafos, and a
  # wrong argument list does not error -- `$<-` coerces, so csem() falls back to
  # its own default estimator and returns a fit that converges and looks fine.
  for (fam in list(c("mat", "native"), c("FIT", "FIT"))) {
    set.seed(11)
    tr <- grow_tree(
      fam[1], fam[2],
      if (fam[1] == "mat") ctl_mixed else ctl_part
    )
    # Leaf refits can legitimately fail on a small node (see n_fail_leaf), so
    # take the first leaf that produced a fit rather than the first leaf.
    fits <- Filter(Negate(is.null), partykit::nodeapply(
      tr,
      ids = partykit::nodeids(tr, terminal = TRUE),
      FUN = function(n) partykit::info_node(n)$object
    ))
    expect_gt(length(fits), 0L)
    used <- fits[[1]]$Information$Arguments
    for (a in c(
      ".approach_weights", ".GSCA_modes", ".conv_criterion",
      ".disattenuate", ".iter_max", ".tolerance"
    )) {
      expect_identical(
        used[[a]], res$Information$Arguments[[a]],
        info = paste(fam[1], fam[2], a)
      )
    }
  }
})

test_that("fit_csem() refuses an argument list that is not one", {
  # The failure this prevents is silent: a model string coerces to a list whose
  # unnamed element matches .model by position, so the fit runs under csem()'s
  # default estimator instead of the object's.
  expect_error(fit_csem(dat, model), "argument list")
})

test_that("doTrees() rejects input it cannot grow a tree from", {
  expect_error(doTrees(dat, covs), "cSEMResults")

  # Grouping is what the tree is for, so a multigroup fit is refused rather
  # than unwrapped.
  res_multi <- csem(
    .data = dat,
    .model = model,
    .id = "z_true",
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .iter_max = 100,
    .GSCA_modes = "CCMP",
    .tolerance = 0.0001
  )
  expect_error(doTrees(res_multi, covs), "single-group")

  expect_error(doTrees(res, c("z_true", "nope")), "nope")
  expect_error(doTrees(res, character(0)), "at least one")
  expect_error(doTrees(res, c("z_true", "x11")), "also indicators")
})

test_that("the GSCA-only statistics refuse a non-GSCA fit", {
  res_pls <- csem(.data = dat, .model = model)
  # calculateGSCAErrors() returns NA rather than erroring off GSCA, and
  # calculateFIT() errors into a silent NA, so both are caught up front.
  expect_error(
    doTrees(res_pls, covs, .influence = "mat", .control = ctl_mixed),
    "needs a GSCA fit"
  )
  expect_error(
    doTrees(res_pls, covs, .influence = "DLi", .splitter = "FIT",
            .control = ctl_part),
    "needs a GSCA fit"
  )
  # The distances read model-implied indicator VCVs and are estimator-agnostic,
  # so this pairing is allowed through and estimates its nodes as PLS-PM. Only
  # that it runs is asserted: whether a PLS fit yields a significant split on
  # this fixture is a question about the data, not about the guard.
  set.seed(11)
  tr <- doTrees(res_pls, covs, .influence = "DLi", .splitter = "DLi",
                .control = ctl_part)
  expect_s3_class(tr, "igsca_tree")
  expect_identical(attr(tr, "igsca_info")$n_fail_full, 0L)
})

# Collector diagnostics --------------------------------------------------
# The counters on attr(tr, "igsca_info") are the port's only window into partial
# failure: every one of them is reached by a path that still returns a
# well-formed tree, so nothing else in the suite would notice them go wrong.

test_that("a root fit failure is a full failure, not a node failure", {
  local_mocked_bindings(
    try_fit = function(.data, .args, .id = NULL) list(fit = NULL, ok = FALSE)
  )
  set.seed(11)
  tr <- grow_tree("mat", "native", ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_fail_full, 1L)
  # The root is counted once and never as a node, which is the whole job of the
  # collector's root_seen flag.
  expect_identical(info$n_fail_node, 0L)
  # No fit means no test ran, so there is nothing to split on and no criteria.
  expect_null(info$root_criteria)
  expect_equal(partykit::width(partykit::node_party(tr)), 1)
})

test_that("a fit failure below the root is a node failure", {
  # Fail on anything smaller than the full sample: the root fit succeeds and
  # every child fit fails. Capturing try_fit first is what avoids recursion.
  real_try_fit <- try_fit
  local_mocked_bindings(
    try_fit = function(.data, .args, .id = NULL) {
      if (nrow(.data) < nrow(dat)) {
        list(fit = NULL, ok = FALSE)
      } else {
        real_try_fit(.data, .args, .id)
      }
    }
  )
  set.seed(11)
  tr <- grow_tree("mat", "native", ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_fail_full, 0L)
  expect_gt(info$n_fail_node, 0L)
  # The root test still ran, so the tree split before the children failed.
  expect_gt(partykit::width(partykit::node_party(tr)), 1)
})

test_that("root_criteria records the test that chose the root split", {
  set.seed(11)
  rc <- attr(grow_tree("mat", "native", ctl_mixed), "igsca_info")$root_criteria
  expect_identical(colnames(rc), covs)
  expect_true(all(c("statistic", "p.value", "criterion") %in% rownames(rc)))
  # z_true is the fixture's true grouping variable, so it must win the root
  # test outright and the noise covariates must not clear alpha.
  expect_identical(colnames(rc)[which.max(rc["criterion", ])], "z_true")
  expect_lt(rc["p.value", "z_true"], 0.05)
  expect_true(all(rc["p.value", c("noise_1", "noise_2")] > 0.05))
})

test_that("a failed leaf refit is counted and surfaces in coef()", {
  set.seed(12353)
  tr <- grow_tree("mat", "native", igsca_tree_control())
  # With 13 indicators the default minbucket = 30L admits leaves the model
  # cannot fit; this tree's n = 68 leaf is one of them.
  expect_identical(attr(tr, "igsca_info")$n_fail_leaf, 1L)

  cf <- coef(tr)
  # One row per leaf, named by node id: rbind() drops NULLs, which is how this
  # method used to collapse a whole tree to NULL.
  expect_identical(
    rownames(cf),
    as.character(partykit::nodeids(tr, terminal = TRUE))
  )
  # The failed leaf is an NA objfun, not a dropped row -- it still reports nobs.
  expect_identical(sum(is.na(cf[, "objfun"])), 1L)
  expect_false(anyNA(cf[, "nobs"]))
})
