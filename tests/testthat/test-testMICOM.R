context("testMGD")

source("test-main.R")
# source("tests/testthat/test-main.R") 
# ==============================================================================
# Tests 
# ==============================================================================
## Errors

test_that("testMICOM() fails for 2ndorder models", {
  expect_error(testMICOM(res_multi_2ndorder, .verbose = FALSE))
  expect_error(testMICOM(res_multi_id_2ndorder, .verbose = FALSE))
  expect_error(testMICOM(res_multi_nonlinear_2ndorder, .verbose = FALSE))
})

test_that("testMICOM() works for linear and nonlinear models (data as list)", {
  expect_output({
    testMICOM(
      .object = res_multi_linear,
      .handle_inadmissibles = "replace",
      .R = 5
    )
  })
  expect_output({
    testMICOM(
      .object = res_multi_nonlinear,
      .handle_inadmissibles = "replace",
      .R = 5
    )
  })
})

test_that("testMICOM() works for linear and nonlinear models (data with id)", {
  expect_output({
    testMICOM(
      .object = res_multi_id_linear,
      .handle_inadmissibles = "replace",
      .R = 5
    )
  })
  expect_output({
    testMICOM(
      .object = res_multi_id_nonlinear,
      .handle_inadmissibles = "replace",
      .R = 5
    )
  })
})