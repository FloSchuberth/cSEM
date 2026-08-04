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
      splitter = split_max_fitdiff,
      control = igsca_tree_control()
    )

    trees_mx_DLi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = split_max_dli,
      control = igsca_tree_control()
    )

    trees_mx_DGi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "mat",
      splitter = split_max_dgi,
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
      splitter = split_max_fitdiff,
      control = igsca_tree_control()
    )

    trees_vec_DLi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "vec",
      splitter = split_max_dli,
      control = igsca_tree_control()
    )

    trees_vec_DGi <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "vec",
      splitter = split_max_dgi,
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
    splitter = split_max_fitdiff,
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
      splitter = split_max_dli,
      control = igsca_tree_control()
    )

    trees_NPT_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "FIT",
      splitter = split_max_dgi,
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
      splitter = split_max_dli,
      control = igsca_tree_control()
    )

    trees_DLi_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DLi",
      splitter = split_max_dgi,
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
      splitter = split_max_dli,
      control = igsca_tree_control()
    )

    trees_DGi_DGI <- doTrees(
      data = dat,
      model = model,
      covariates = covs,
      influence = "DGi",
      splitter = split_max_dgi,
      control = igsca_tree_control()
    )
  })
})
