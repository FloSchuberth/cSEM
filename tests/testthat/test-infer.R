source(testthat::test_path("test-main.R"))
#
# # ==============================================================================
# # Tests
# # ==============================================================================

test_that("Run without errors", {
  infer(res_single_linear_boot) |> expect_no_error()

  infer(res_single_linear_boot, .quantity = "mean") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "sd") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "bias") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_standard_z") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_standard_t") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_percentile") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_basic") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_bc") |> expect_no_error()
  infer(res_single_linear_boot, .quantity = "CI_bca") |> expect_no_error()
  
  infer(res_single_2ndorder_boot, .quantity = "mean") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "sd") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "bias") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_standard_z") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_standard_t") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_percentile") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_basic") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_bc") |> expect_no_error()
  infer(res_single_2ndorder_boot, .quantity = "CI_bca") |> expect_no_error()
})

# ==============================================================================
# GSCA-M / IGSCA: bootstrap of Unique_loading_estimates
# ==============================================================================

model_GSCA_boot <- "
# Measurement models
OrgPres <~ cei1 + cei2 + cei3
OrgIden <~ ma1 + ma2 + ma3
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_GSCAM_boot <- "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3
OrgIden =~ ma1 + ma2 + ma3
AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

model_IGSCA_boot <- "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3
OrgIden <~ ma1 + ma2 + ma3
AffJoy <~ orgcmt1 + orgcmt2 + orgcmt3
AffLove  <~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

res_single_GSCA_boot <- csem(
  .data                 = BergamiBagozzi2000,
  .model                = model_GSCA_boot,
  .approach_weights     = "GSCA",
  .GSCA_modes           = "NCMP",
  .tolerance            = 1e-5,
  .conv_criterion       = "sum_diff_absolute",
  .resample_method      = "bootstrap",
  .R                    = 32,
  .seed                 = 1234,
  .handle_inadmissibles = "replace"
)

res_single_GSCAM_boot <- csem(
  .data                 = BergamiBagozzi2000,
  .model                = model_GSCAM_boot,
  .approach_weights     = "GSCA",
  .GSCA_modes           = "NCMP",
  .tolerance            = 1e-5,
  .conv_criterion       = "sum_diff_absolute",
  .resample_method      = "bootstrap",
  .R                    = 32,
  .seed                 = 1234,
  .handle_inadmissibles = "replace"
)

res_single_IGSCA_boot <- csem(
  .data                 = BergamiBagozzi2000,
  .model                = model_IGSCA_boot,
  .approach_weights     = "GSCA",
  .GSCA_modes           = "NCMP",
  .tolerance            = 1e-5,
  .conv_criterion       = "sum_diff_absolute",
  .resample_method      = "bootstrap",
  .R                    = 32,
  .seed                 = 1234,
  .handle_inadmissibles = "replace"
)

res_multi_IGSCA_boot <- csem(
  .data                 = BergamiBagozzi2000,
  .model                = model_IGSCA_boot,
  .approach_weights     = "GSCA",
  .GSCA_modes           = "NCMP",
  .id                   = "gender",
  .tolerance            = 1e-5,
  .conv_criterion       = "sum_diff_absolute",
  .resample_method      = "bootstrap",
  .R                    = 32,
  .seed                 = 1234,
  .handle_inadmissibles = "replace"
)

test_that("infer() bootstraps Unique_loading_estimates for GSCA-M", {
  out <- infer(res_single_GSCAM_boot, .quantity = c("sd", "CI_percentile"))

  expect_true("Unique_loading_estimates" %in% names(out))
  expect_true(all(is.finite(out$Unique_loading_estimates$sd)))
  expect_true(all(out$Unique_loading_estimates$sd > 0))
  expect_true(all(is.finite(out$Unique_loading_estimates$CI_percentile)))
})

test_that("infer() bootstraps Unique_loading_estimates for IGSCA and only includes common-factor indicators", {
  out <- infer(res_single_IGSCA_boot, .quantity = c("sd", "CI_percentile"))

  expect_true("Unique_loading_estimates" %in% names(out))
  expect_true(all(is.finite(out$Unique_loading_estimates$sd)))
  expect_true(all(out$Unique_loading_estimates$sd > 0))

  # Only OrgPres (cei1, cei2, cei3) is modeled as a common factor in model_IGSCA_boot.
  # Names follow the `Name` column format built in postestimate_summarize.R
  # ("ind =~ ind").
  expect_setequal(names(out$Unique_loading_estimates$sd),
                  c("cei1 =~ cei1", "cei2 =~ cei2", "cei3 =~ cei3"))
})

test_that("infer() bootstraps Unique_loading_estimates per group for multigroup IGSCA", {
  out <- infer(res_multi_IGSCA_boot, .quantity = c("sd", "CI_percentile"))

  expect_true(is.list(out))
  for(group_out in out) {
    expect_true("Unique_loading_estimates" %in% names(group_out))
    expect_true(all(is.finite(group_out$Unique_loading_estimates$sd)))
    expect_true(all(group_out$Unique_loading_estimates$sd > 0))
  }
})

test_that("summarize() fills Std_err/CI for Unique_loading_estimates (GSCA-M)", {
  sum_out <- summarize(res_single_GSCAM_boot)
  ul <- sum_out$Estimates$Unique_loading_estimates

  expect_true(nrow(ul) > 0)
  expect_true(all(is.finite(ul$Std_err)))
  expect_true(all(is.finite(ul$t_stat)))
  expect_true(all(is.finite(ul$p_value)))

  # There should be at least one CI column added, mirroring Loading_estimates
  loading_ci_colnames <- setdiff(colnames(sum_out$Estimates$Loading_estimates),
                                 c("Name", "Construct_type",
                                   "Estimate", "Std_err", "t_stat", "p_value"))
  expect_true(length(loading_ci_colnames) > 0)
  for(cn in loading_ci_colnames) {
    expect_true(cn %in% colnames(ul))
    expect_true(all(is.finite(ul[[cn]])))
  }
})

test_that("plain-GSCA and PLS bootstrap objects are unaffected (no Unique_loading_estimates quantity)", {
  ## Plain GSCA (all composite, <~): no common-factor indicators, hence no
  ## Unique_loading_estimates quantity should be resampled/inferred.
  expect_no_error(infer(res_single_GSCA_boot, .quantity = "sd"))
  out_gsca <- infer(res_single_GSCA_boot, .quantity = "sd")
  expect_null(out_gsca$Unique_loading_estimates)

  expect_no_error(summarize(res_single_GSCA_boot))
  sum_gsca <- summarize(res_single_GSCA_boot)
  expect_null(sum_gsca$Estimates$Unique_loading_estimates)

  ## Existing PLS bootstrap fixture (from test-main.R): regression guard
  expect_no_error(infer(res_single_linear_boot, .quantity = "sd"))
  out_pls <- infer(res_single_linear_boot, .quantity = "sd")
  expect_null(out_pls$Unique_loading_estimates)

  expect_no_error(summarize(res_single_linear_boot))
})

