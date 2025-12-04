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

