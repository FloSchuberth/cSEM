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

# Controls ---------------------------------------------------------------
# The two influence families cost an order of magnitude apart, so they get
# separate budgets. The conditional-test family ("mat"/"vec") inherits COIN
# variable selection from libcoin and is cheap even at the defaults; a mixed
# splitter only re-scans candidate cutpoints, which adds a second or two. The
# partition family ("FIT"/"DLi"/"DGi") runs an R_test permutation test per
# covariate per node and every permutation is a two-group MGA fit, so at the
# defaults a single call takes about an hour and a half.
ctl_mixed <- igsca_tree_control(R_test = 50L, maxdepth = 2L, max_cuts = 8L)

# R_test >= 20L is a correctness floor, not just a budget: permutation_pvalue()
# returns (1 + k) / (R + 1), so the smallest attainable p-value is 1 / (R + 1),
# and the split criterion is a strict inequality against alpha = 0.05. At
# R_test = 19L that floor is exactly 0.05 and no node can ever split, whatever
# the data says. 24L keeps headroom at 1/25 = 0.04. maxdepth = 1L costs nothing
# in coverage on this fixture -- measured, the children do not split at
# maxdepth = 2L either, they only double the runtime.
ctl_part <- igsca_tree_control(
  R_test = 24L,
  max_cuts = 3L,
  maxdepth = 1L,
  minbucket = 200L
)

grow_tree <- function(influence, splitter, control) {
  doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = influence,
    splitter = splitter,
    control = control
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
