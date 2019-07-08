context("testMGD")

## Model linear
model_linear <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Model nonlinear
model_nonlinear <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2 + eta1.eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Data
dat <- list(threecommonfactors, 
            threecommonfactors[1:200, ], 
            threecommonfactors[130:250,])

## Estimates (.R is small to save computation time)
res_error <- csem(threecommonfactors, model_linear)
res_nonlinear <- csem(dat, model_nonlinear)
res_linear <-  csem(dat, model_linear, .resample_method = "bootstrap", .R = 80)


# ==============================================================================
# Tests 
# ==============================================================================

test_that("Not providing a cSEMResults_multi object causes an error", {
  expect_error(
    testMGD(res_error)
  )
})

test_that(".approach_mgd = 'Klesel' does not work for linear models", {
  expect_error(
    testMGD(res_nonlinear, .approach_mgd = "Klesel")
  )
  expect_error(
    testMGD(res_nonlinear, .approach_mgd = "all")
  )
})
