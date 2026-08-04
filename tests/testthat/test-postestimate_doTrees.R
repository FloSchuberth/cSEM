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
  trees_out <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = influence_mat,
    control = igsca_tree_control(),
  )
  expect_snapshot(trees_out)

})



# Vector of Residuals ----------------------------------------------------

test_that("IGSCA Trees Conditional Test on Vector of Residuals Runs as expected", {
  set.seed(12353)
  trees_out <- doTrees(
    data = dat,
    model = model,
    covariates = covs,
    influence = influence_vec,
    control = igsca_tree_control(),
  )
  expect_snapshot(trees_out)
})
