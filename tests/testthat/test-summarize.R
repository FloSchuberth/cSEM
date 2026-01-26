source(testthat::test_path("test-main.R"))


# ==============================================================================
# Tests
# ==============================================================================
test_that("Run without errors", {
  summarize(res_single_linear) |> expect_no_error()
  summarize(res_single_nonlinear) |> expect_no_error()
  summarize(res_single_2ndorder) |> expect_no_error()
  summarize(res_single_nonlinear_2ndorder) |> expect_no_error()
  summarize(res_single_linear_boot) |> expect_no_error()
  summarize(res_single_linear_boot) |> expect_no_error()
  summarize(res_single_2ndorder_boot) |> expect_no_error()
  summarize(res_single_nonlinear_2ndorder_boot) |> expect_no_error()
})
