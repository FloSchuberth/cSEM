context("csem")

#
model <- "
# Structural model
EXPE ~ IMAG

# Measurement model

IMAG <~ imag1 + imag2
EXPE <~ expe1 + expe2 + expe3
"

test_that("No data provided should give an error", {
  expect_error(
    csem(.model = model)
  )
})

### 2 Common factors ===========================================================

# load(file = "DGP_linear_2commonfactors/1_DGPs.RData") # Loads only Sigma
# load(file = "DGP_linear_2commonfactors/1_DGPs_models_pop_params.RData") # Loads Sigma and models
load(file = "Temp/DGPs/DGP_linear_2commonfactors/1_DGPs_models_pop_params.RData") # Loads Sigma and models

## Draw data

dat <- MASS::mvrnorm(200, rep(0, 6), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate

res <-  csem(dat, model_Sigma)

# Compare
compare(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
compare(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)



test_that("Linear model: correctly specified models are correctly returned", {
  expect_s3_class(parseModel(model4), "cSEMModel")

})