context("csem")

# ==============================================================================
# Tests 
# ==============================================================================
### 2 Common factors -----------------------------------------------------------
## Load data
load(file = "../data/1_DGPs_models_pop_params.RData")

## Draw data
dat <- MASS::mvrnorm(200, rep(0, 6), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)
res_sum <- summarize(res)

### Test

test_that("No data/model provided should give an error", {
  expect_error(
    csem(.model = model)
  )
  expect_error(
    csem(.data = dat)
  )
})

test_that("Linear model: correctly specified models are correctly estimated", {
  expect_equal(res_sum$Estimates$Path_estimates$Estimate, pop_params_Sigma$Path_coefficients)
  expect_equal(res_sum$Estimates$Loading_estimates$Estimate, unname(pop_params_Sigma$Loadings))
})