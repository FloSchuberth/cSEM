## Source and set-up test-fixture if necessary
source(testthat::test_path("test-main.R"))
dir.create(testthat::test_path("test_results_exportToExcel"), showWarnings = FALSE)

# ==============================================================================
# General tests 
# ==============================================================================
test_that("No data/model provided causes an error", {
  expect_error(
    csem(.model = model)
  )
  expect_error(
    csem(.data = dat)
  )
})

# ==============================================================================
# DGPs
# ==============================================================================
### DGP_linear_2commonfactors --------------------------------------------------
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_linear_2commonfactors.RData"))

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)

## Test
test_that("DPG_linear_2commonfactors is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
})

### DGP_linear_2composites =====================================================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_linear_2composites.RData"))

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Test
test_that("DPG_linear_2composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})

### DGP_linear_3commonfactors ==================================================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_linear_3commonfactors.RData"))

## Draw data
dat <- MASS::mvrnorm(300, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
for (i in c(
  "PLS-PM",
  "GSCA",
  "SUMCORR",
  "MAXVAR",
  "MINVAR",
  "GENVAR",
  "PCA",
  "unit",
  "bartlett",
  "regression"
)) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularly producing differences
  ##   between estimate and population value larger than 0.01.

  if (i == "PLS-PM") {
    for (j in c("modeA", "modeB")) {
      res <- csem(
        dat,
        model_Sigma,
        .approach_weights = i,
        .PLS_modes = j,
        .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31")
      )

      ## Comparison
      path <- comparecSEM(
        res,
        .what = "Path_estimates",
        pop_params_Sigma$Path_coefficients
      )
      loadings <- comparecSEM(
        res,
        .what = "Loading_estimates",
        pop_params_Sigma$Loadings
      )

      ## Test
      # Note: the tolerance is necessary since Croon correction requires a CFA (i.e. ML estimation)
      # which is only consistent (i.e. will not exactly reproduce the population values for a finite sample size)
      # Similarly for GSCAm.
      test_that(paste("Weighting approach", i, "yields correct results"), {
        expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
        expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.01)
      })
    }
  } else {
    res <- csem(
      dat,
      model_Sigma,
      .approach_weights = i,
      .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31")
    )

    ## Comparison
    path <- comparecSEM(
      res,
      .what = "Path_estimates",
      pop_params_Sigma$Path_coefficients
    )
    loadings <- comparecSEM(
      res,
      .what = "Loading_estimates",
      pop_params_Sigma$Loadings
    )

    ## Test
    # Note: the tolerance is necessary since Croon correction requires a CFA (i.e. ML estimation)
    # which is only consistent (i.e. will not exactly reproduce the population values for a finite sample size)
    # Similarly for GSCAm.
    test_that(paste("Weighting approach", i, "yields correct results"), {
      expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
      expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.01)
    })
  }

  # Export to Excel test
  exportToExcel(
    assess(res),
    .filename = paste0("test_assess_", i, ".xlsx"),
    .path = testthat::test_path("test_results_exportToExcel")
  )
  exportToExcel(
    summarize(res),
    .filename = paste0("test_summarize_", i, ".xlsx"),
    .path = testthat::test_path("test_results_exportToExcel")
  )
  exportToExcel(
    predict(res, .handle_inadmissibles = "ignore", .benchmark = "lm", .disattenuate = FALSE),
    .filename = paste0("test_predict_", i, ".xlsx"),
    .path = testthat::test_path("test_results_exportToExcel")
  )
  exportToExcel(
    testOMF(res, .R = 10),
    .filename = paste0("test_testOMF_", i, ".xlsx"),
    .path = testthat::test_path("test_results_exportToExcel")
  )
}

### DGP_linear_3compostites ====================================================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_linear_3composites.RData"))

## Draw data
dat <- MASS::mvrnorm(300, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
for(i in c("PLS-PM", "GSCA", "SUMCORR", "MAXVAR", "MINVAR", "GENVAR")) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularly producing differences
  ##   between estimate and population value larger than 0.01.
  ## - "PCA" is excluded as weights obtained by PCA are not population weights but 
  ##   simply the first principal component of S_jj.
  ## - "unit" is excluded as unit weights are inconsistent "estimates" for the
  ##   population weights (weights are simply set to 1 and scaled).
  ## - "bartlett" and "regression" are excluded as they are not meaningful 
  ##   for models containing concepts modeled as composites
  res <-  csem(dat, model_Sigma, .approach_weights = i, 
               .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"))
  
  ## Comparison
  path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
  loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
  weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)
  
  ## Test
  test_that(paste("Weighting approach", i, "yields correct results"), {
    expect_equal(path$Estimate, path$Pop_value, tolerance = 0.001)
    expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.001)
    expect_equal(weights$Estimate, weights$Pop_value, tolerance = 0.001)
  })
}

### DGP_2ndorder - Common factor of common factors =============================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_2ndorder_cf_of_cfs.RData"))

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_cf_of_cfs is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
})

# Export to Excel test
exportToExcel(summarize(res), .filename = "test_summarize", .path = testthat::test_path("test_results_exportToExcel"))
exportToExcel(assess(res), .filename = "test_assess", .path = testthat::test_path("test_results_exportToExcel"))
exportToExcel(testOMF(res, .R = 20), .filename = "test_testOMF", .path = testthat::test_path("test_results_exportToExcel"))

### DGP_2ndorder - Common factor of composites =================================
# Loads Sigma, models and population values
load(testthat::test_path("data/DGP_2ndorder_cf_of_composites.RData"))


## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_cf_of_composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})

### DGP_2ndorder - Composite of common factors =================================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_2ndorder_composite_of_cfs.RData"))

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_composites_of_cfs is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})
### DGP_2ndorder - Composite of composites =====================================
# Loads Sigma, models and population values
load(
  file = testthat::test_path("data/DGP_2ndorder_composite_of_composites.RData")
)

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_composites_of_composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})
