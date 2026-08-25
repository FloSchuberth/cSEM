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
# doTrees() has one selector -- COIN variable selection, from libcoin -- and
# four cutpoint rules. The native rule is cheap enough to run at the package
# defaults; a "FIT"/"DLi"/"DGi" rule re-scans candidate cutpoints with a model
# comparison statistic, which adds a second or two.
ctl_mixed <- igsca_tree_control(R_test = 50L, maxdepth = 2L, max_cuts = 8L)

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
  expect_identical(info$n_fail_full, 0L) # TODO: ?
  expect_identical(info$n_fail_split, 0L) # No failures to split
  # The width of the grown tree should be greater than 1 (there should be more than a stub)
  expect_gt(partykit::width(partykit::node_party(tree)), 1)
}

# Matrix of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Matrix of Residuals Runs as expected", {
  set.seed(12353)
  # The native path is cheap enough to snapshot at the package defaults.
  trees_mx <- grow_tree(influence = "mat", splitter = "native", control = igsca_tree_control())
  expect_grew(trees_mx)
  expect_snapshot(trees_mx)
})

## Variable selection on matrix of residuals and split point selection on different criteria
for (sp in c("FIT", "DLi", "DGi")) {
  test_that(paste0("Matrix of Residuals selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree(influence = "mat", splitter = sp, control = ctl_mixed))
  })
}

# Vector of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Vector of Residuals Runs as expected", {
  set.seed(12353)
  trees_vec <- grow_tree(influence = "vec", splitter = "native", control = igsca_tree_control())
  expect_grew(trees_vec)
  expect_snapshot(trees_vec)
})

for (sp in c("FIT", "DLi", "DGi")) {
  test_that(paste0("Vector of Residuals selection splits on ", sp), {
    set.seed(12353)
    expect_grew(grow_tree(influence = "vec", splitter = sp, control = ctl_mixed))
  })
}



# Splitfun contract ------------------------------------------------------
# doTrees() plants its kernel in cc$splitfun and relies on partykit reading it
# back out of the control object. Nothing in ?ctree_control promises that, so it
# is an undocumented contract, and it has already survived one partykit upgrade
# unverified. If a release stops honouring it, every mixed configuration
# silently falls back to partykit's own maxstat scan and still returns a
# perfectly plausible tree -- so counting the kernel's invocations is the only
# way to see it happen.
test_that("doTrees() installs its splitfun into partykit's split search", {
  calls <- new.env(parent = emptyenv())
  # Represents the number of times that split_max_fitdiff is called
  calls$n <- 0L
  # Capturing the real kernel first is what stops the mock recursing into itself.
  real <- split_max_fitdiff
  local_mocked_bindings(
    split_max_fitdiff = function(model, mf, subset, goes_left, ctrl) {
      # `calls`` will be modified even without `<<-` syntax because it is an environment
      calls$n <- calls$n + 1L
      real(model, mf, subset, goes_left, ctrl)
    }
  )
  set.seed(11)
  tr <- grow_tree(influence = "mat", splitter = "FIT", control = ctl_mixed)
  # print(calls$n)
  expect_gt(calls$n, 0L)
  # ... and the tree that came back is one the kernel actually shaped.
  expect_grew(tr)
})

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
    tr <- grow_tree(influence = "mat", splitter = "DLi", control = ctl_mixed),
    "produced no usable statistic"
  )
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_fail_split, 0L)
  expect_identical(info$n_fail_split, info$n_split_scan) # TODO: Why should these be identical?
  # The tree really is the stump the diagnostic is warning about
  expect_equal(partykit::width(partykit::node_party(tr)), 1)
})

test_that("a working split kernel records scans but no failures", {
  set.seed(11)
  tr <- expect_no_warning(grow_tree(influence = "mat", splitter = "DLi", control = ctl_mixed))
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_split_scan, 0L)
  expect_grew(tr)
})

test_that("the native split path never touches the kernel counters", {
  set.seed(11)
  tr <- grow_tree(influence = "mat", splitter = "native", control = ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_split_scan, 0L)
  expect_identical(info$n_fail_split, 0L)
})

# Input contract ---------------------------------------------------------
# doTrees() takes a fitted object rather than data + model, so everything it
# needs is either derivable from that fit or a mistake worth naming.

test_that("doTrees() refits nodes with the estimator the fit used", {
  # A wrong argument list does not error -- `$<-` coerces, so csem() falls
  # back to its own default estimator and returns a fit that converges and
  # looks fine.
  # One trafo now, but two cutpoint routes: the native scan and a kernel that
  # re-fits candidate partitions of its own. Both must reach csem() with the
  # object's argument list.
  for (fam in list(c("mat", "native"), c("mat", "DLi"))) {
    set.seed(11)
    tr <- grow_tree(
      influence = fam[1],
      splitter = fam[2],
      control = ctl_mixed
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

test_that("fit_csem() refuses inappropriate arguments", {
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

test_that("doTrees() no longer offers the partition influence families", {
  # "FIT", "DLi" and "DGi" survive as .splitter values only: they are split
  # kernels, and the selector they are paired with is always the COIN one.
  for (bad in c("FIT", "DLi", "DGi")) {
    expect_error(
      doTrees(res, covs, .influence = bad, .control = ctl_mixed),
      "should be one of",
      fixed = TRUE,
      info = bad
    )
  }
})

test_that("every configuration refuses a non-GSCA fit", {
  res_pls <- csem(.data = dat, .model = model)
  # calculateGSCAErrors() returns NA rather than erroring off a non-GSCA fit,
  # so the node statistic would fail deep inside the trafo. Both influence
  # values read it, so there is no configuration left that a PLS fit can grow.
  for (inf in c("mat", "vec")) {
    expect_error(
      doTrees(res_pls, covs, .influence = inf, .control = ctl_mixed),
      "needs a GSCA fit",
      fixed = TRUE,
      info = inf
    )
  }
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
  tr <- grow_tree(influence = "mat", splitter = "native", control = ctl_mixed)
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
  tr <- grow_tree(influence = "mat", splitter = "native", control = ctl_mixed)
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_fail_full, 0L)
  expect_gt(info$n_fail_node, 0L)
  # The root test still ran, so the tree split before the children failed.
  expect_gt(partykit::width(partykit::node_party(tr)), 1)
})

test_that("root_criteria records the test that chose the root split", {
  set.seed(11)
  rc <- attr(
    grow_tree(influence = "mat", splitter = "native", control = ctl_mixed),
    "igsca_info"
  )$root_criteria
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
  tr <- grow_tree(influence = "mat", splitter = "native", control = igsca_tree_control())
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

test_that("csem_tree_args passes raw data and not standardized data", {
  expect_equal(csem_tree_args(res)$.data, dat)
  local_mocked_bindings(
    csem_tree_args = function(.object) {
      stop("This function uses csem_tree_args")
    }
  )
  expect_error(grow_tree(influence = "mat", splitter = "native", control = igsca_tree_control()), "This function uses csem_tree_args")
})
