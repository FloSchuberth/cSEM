## Source and set-up test-fixture if necessary
source(testthat::test_path("test-main.R"))
dir.create(testthat::test_path("test_results_exportToExcel"), showWarnings = FALSE)

withr::local_seed(442363)

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
for(i in c("PLS-PM", "GSCA", "SUMCORR", "MAXVAR", "MINVAR", "GENVAR", "PCA",
           "unit", "bartlett", "regression")) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularly producing differences
  ##   between estimate and population value larger than 0.01.
  
  if(i == "PLS-PM"){
    for(j in c("modeA","modeB")){
      res <-  csem(dat, model_Sigma, .approach_weights = i,.PLS_modes = j, 
                   .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"))
      
      ## Comparison
      path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
      loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
      
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
      .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"),
      .conv_criterion = if (i == "GSCA") {
              "sum_diff_absolute"
            } else {
              "diff_absolute"
            },
            .iter_max = if (i == "GSCA") {
              1000
            } else {
              100
            }
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

  test_that(paste("tidy method fails with weighting approach", i), {
    expect_no_error(tidy(res))
    expect_s3_class(tidy(res), c("tbl", "tbl_df"))
  })

  test_that(paste("glance method fails with weighting approach", i), {
    expect_no_error(glance(res))
    expect_s3_class(glance(res), c("tbl", "tbl_df"))
  })

  # Export to Excel test
  exportToExcel(assess(res), .filename = paste0("test_assess_", i, ".xlsx"),
                .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
  exportToExcel(summarize(res), .filename = paste0("test_summarize_", i, ".xlsx"),
                .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
  exportToExcel(predict(res, .handle_inadmissibles = 'set_NA', .benchmark = "lm", .disattenuate = FALSE),
                .filename = paste0("test_predict_", i, ".xlsx"),
                .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
  exportToExcel(testOMF(res, .R = 10, .handle_inadmissibles = 'drop'),
                        .filename = paste0("test_testOMF_", i, ".xlsx"),
                        .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
}

### DGP_linear_3compostites ====================================================
# Loads Sigma, models and population values
load(file = testthat::test_path("data/DGP_linear_3composites.RData"))

## Draw data
dat <- MASS::mvrnorm(300, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
for (i in c("PLS-PM", "GSCA", "SUMCORR", "MAXVAR", "MINVAR", "GENVAR")) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularly producing differences
  ##   between estimate and population value larger than 0.01.
  ## - "PCA" is excluded as weights obtained by PCA are not population weights but
  ##   simply the first principal component of S_jj.
  ## - "unit" is excluded as unit weights are inconsistent "estimates" for the
  ##   population weights (weights are simply set to 1 and scaled).
  ## - "bartlett" and "regression" are excluded as they are not meaningful
  ##   for models containing concepts modeled as composites
  res <- csem(
    dat,
    model_Sigma,
    .approach_weights = i,
    .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"),
    .conv_criterion = if (i == "GSCA") {
      "sum_diff_absolute"
    } else {
      "diff_absolute"
    }
  )
  
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

  test_that(paste("tidy method runs without failure with weighting approach", i), {
    expect_no_error(tidy(res))
    expect_s3_class(tidy(res), c("tbl", "tbl_df"))
  })

  test_that(paste("glance method runs without failure with weighting approach", i), {
    expect_no_error(glance(res))
    expect_s3_class(glance(res), c("tbl", "tbl_df"))
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

# Guard against future bugs
# Origin: Errors in printing the results
expect_no_error(capture.output(print(res)))

# Export to Excel test
exportToExcel(summarize(res), .filename = "test_summarize_sole.xlsx", .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
exportToExcel(suppressWarnings(assess(res)), .filename = "test_assess_sole.xlsx", .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)
exportToExcel(testOMF(res, .R = 20), .filename = "test_testOMF_sole.xlsx", .path = testthat::test_path("test_results_exportToExcel"),
              .quiet = TRUE)

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
load(file = testthat::test_path("data/DGP_2ndorder_composite_of_composites.RData"))

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

# GSCA multi-start starting values --------------------------------------------

mstart_model <- "
OrgIden <~ ma1 + ma2 + ma3
AffJoy  <~ orgcmt1 + orgcmt2 + orgcmt3
OrgPres =~ cei1 + cei2 + cei3
OrgIden ~ OrgPres
AffJoy  ~ OrgIden
"

test_that("random-spec starting values run and select a best fit", {
  set.seed(123)
  res <- csem(BergamiBagozzi2000, mstart_model,
              .approach_weights = "GSCA", .GSCA_modes = "NCMP",
              .starting_values = c(n = 4, min = -1, max = 1))
  ms <- res$Information$Weight_info$Multistart
  expect_equal(ms$n_candidates, 4)
  expect_equal(ms$FIT_selected, max(ms$FIT_candidates, na.rm = TRUE))
})

test_that("random-spec multi-start is reproducible under set.seed()", {
  set.seed(99)
  r1 <- csem(BergamiBagozzi2000, mstart_model,
             .approach_weights = "GSCA", .GSCA_modes = "NCMP",
             .starting_values = c(n = 4, min = -1, max = 1))
  set.seed(99)
  r2 <- csem(BergamiBagozzi2000, mstart_model,
             .approach_weights = "GSCA", .GSCA_modes = "NCMP",
             .starting_values = c(n = 4, min = -1, max = 1))
  expect_equal(r1$Estimates$Weight_estimates, r2$Estimates$Weight_estimates)
})

test_that("multi-start is rejected for non-GSCA approaches", {
  set1 <- list(OrgIden = c(ma1 = 1, ma2 = 1, ma3 = 1))
  expect_error(
    csem(BergamiBagozzi2000, mstart_model,
         .approach_weights = "PLS-PM",
         .starting_values = c(n = 4, min = -1, max = 1)))
  expect_error(
    csem(BergamiBagozzi2000, mstart_model,
         .approach_weights = "PLS-PM",
         .starting_values = list(set1, set1)))
})

test_that("multi-start GSCA selects per group for multigroup data", {
  set.seed(7)
  res <- csem(BergamiBagozzi2000, mstart_model,
              .approach_weights = "GSCA", .GSCA_modes = "NCMP",
              .id = "gender",
              .starting_values = c(n = 3, min = -1, max = 1))
  expect_s3_class(res, "cSEMResults_multi")
  expect_true(all(vapply(res, function(g)
    isTRUE(g$Information$Weight_info$Multistart$n_candidates == 3), logical(1))))
})

test_that("multi-start GSCA supports bootstrapping (multi-start per resample)", {
  set1 <- list(OrgIden = c(ma1 = 1, ma2 = 1, ma3 = 1),
               AffJoy  = c(orgcmt1 = 1, orgcmt2 = 1, orgcmt3 = 1),
               OrgPres = c(cei1 = 1, cei2 = 1, cei3 = 1))
  set2 <- list(OrgIden = c(ma1 = 1, ma2 = -1, ma3 = 1),
               AffJoy  = c(orgcmt1 = -1, orgcmt2 = 1, orgcmt3 = 1),
               OrgPres = c(cei1 = 1, cei2 = 1, cei3 = -1))

  res <- csem(BergamiBagozzi2000, mstart_model,
              .approach_weights = "GSCA", .GSCA_modes = "NCMP",
              .starting_values = list(set1, set2),
              .resample_method = "bootstrap", .R = 20, .seed = 1)

  expect_true(inherits(res, "cSEMResults_resampled"))
  expect_true(nrow(res$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled) > 0)
})

