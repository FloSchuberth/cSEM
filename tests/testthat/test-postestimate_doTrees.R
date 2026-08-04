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



# Matrix of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Matrix of Residuals Runs as expected", {
  set.seed(12353)
  trees_mx <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "mat",
    splitter = NULL,
    control = igsca_tree_control()
  )
  # Dirty substitute for snapshot
  expect_true(length(trees_mx) == 5)
  expect_snapshot(trees_mx)

  expect_no_error({
    trees_mx_FIT <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = "FIT",
      control = igsca_tree_control()
    )

    trees_mx_DLi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = "DLi",
      control = igsca_tree_control()
    )

    trees_mx_DGi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = "DGi",
      control = igsca_tree_control()
    )
  })
})



# Vector of Residuals ----------------------------------------------------
test_that("IGSCA Trees Conditional Test on Vector of Residuals Runs as expected", {
  set.seed(12353)
  trees_vec <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "vec",
    splitter = NULL,
    control = igsca_tree_control()
  )
  expect_true(length(trees_vec) == 5)
  expect_snapshot(trees_vec)

  expect_no_error({
    trees_vec_FIT <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "vec",
      splitter = "FIT",
      control = igsca_tree_control()
    )

    trees_vec_DLi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "vec",
      splitter = "DLi",
      control = igsca_tree_control()
    )

    trees_vec_DGi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "vec",
      splitter = "DGi",
      control = igsca_tree_control()
    )
  })
})



# NPT-FIT Split ---------------------------------------------------------------
test_that("IGSCA Trees Variable Selection on FIT Runs as expected", {
  set.seed(12353)
  # debug(doTrees)
  trees_NPT_FIT <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "FIT",
    splitter = "FIT",
    control = igsca_tree_control()
  )
  expect_true(length(trees_NPT_FIT) == 5)
  expect_snapshot(trees_NPT_FIT)

  expect_no_error({
    trees_NPT_DLI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "FIT",
      splitter = "DLi",
      control = igsca_tree_control()
    )

    trees_NPT_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "FIT",
      splitter = "DGi",
      control = igsca_tree_control()
    )
  })
  
})

# dLi Split --------------------------------------------------------------
test_that("IGSCA Trees Conditional Test on Squared-Euclidean Distance Runs as expected", {
  set.seed(12353)
  trees_DLi_FIT <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "DLi",
    splitter = "FIT",
    control = igsca_tree_control()
  )
  expect_true(length(trees_DLi_FIT) == 5)
  expect_snapshot(trees_DLi_FIT)

  expect_no_error({
    trees_DLi_DLI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DLi",
      splitter = "DLi",
      control = igsca_tree_control()
    )

    trees_DLi_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DLi",
      splitter = "DGi",
      control = igsca_tree_control()
    )
  })
})

# dGi Split --------------------------------------------------------------
test_that("IGSCA Trees Conditional Test on Geodesic Distance Runs as expected", {
  set.seed(12353)
  trees_DGi_FIT <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "DGi",
    splitter = "FIT",
    control = igsca_tree_control()
  )
  expect_true(length(trees_DGi_FIT) == 5)
  expect_snapshot(trees_DGi_FIT)

  expect_no_error({
    trees_DGi_DLI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DGi",
      splitter = "DLi",
      control = igsca_tree_control()
    )

    trees_DGi_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DGi",
      splitter = 'DGi',
      control = igsca_tree_control()
    )
  })
})

# Dead split kernel ------------------------------------------------------
# argmax_split() turns a throwing kernel into NA, which partykit reads as "no
# admissible split" -- so a broken kernel yields a clean-looking stump. These
# tests pin the diagnostic that separates that from a genuinely unsplittable
# node. Mixed splitters run no permutation test, so a reduced control keeps
# these to a couple of seconds.
ctl_mixed <- igsca_tree_control(R_test = 50L, maxdepth = 2L, max_cuts = 8L)

test_that("a dead split kernel is reported, not silently a stump", {
  local_mocked_bindings(
    split_max_dli = function(model, mf, subset, goes_left, ctrl) {
      stop("kernel is broken")
    }
  )
  set.seed(11)
  expect_warning(
    tr <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = "DLi",
      control = ctl_mixed
    ),
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
  tr <- expect_no_warning(doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "mat",
    splitter = "DLi",
    control = ctl_mixed
  ))
  info <- attr(tr, "igsca_info")
  expect_gt(info$n_split_scan, 0L)
  expect_identical(info$n_fail_split, 0L)
  expect_gt(partykit::width(partykit::node_party(tr)), 1L)
})

test_that("the native split path never touches the kernel counters", {
  set.seed(11)
  tr <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = "mat",
    splitter = "native",
    control = ctl_mixed
  )
  info <- attr(tr, "igsca_info")
  expect_identical(info$n_split_scan, 0L)
  expect_identical(info$n_fail_split, 0L)
})
