context("testMGD")

source("test-main.R")
# source("tests/testthat/test-main.R")
# ==============================================================================
# Tests 
# ==============================================================================
## Errors

test_that(paste("Supported mgd approaches are:", args_default(.choices = TRUE)$.approach_mgd, collapse = ", "), {
  expect_error(testMGD(res_multi_linear, .approach_mgd = "Hello World"))
})

test_that("Not providing a cSEMResults_multi object causes an error", {
  expect_error(testMGD(res_single_linear))
  expect_error(testMGD(res_single_linear_boot))
})

test_that(".approach_mgd = 'Klesel' does not work for nonlinear models", {
  expect_error(testMGD(res_single_nonlinear, .approach_mgd = "Klesel"))
  expect_error(testMGD(res_single_nonlinear_boot, .approach_mgd = "all"))
})

test_that(".approach_mgd = 'Sarstedt' cant be combined with .handle_inadmissibles = 'drop'", {
  expect_error(
    testMGD(res_single_linear, .approach_mgd = "Sarstedt", .handle_inadmissibles = "drop")
  )
  expect_error(
    testMGD(res_nonlinear, .approach_mgd = "all", .handle_inadmissibles = "drop")
  )
})

## Correct

test_that("All .approach_mgd options work for linear models", {
  expect_output(
    testMGD(
      .object       = res_multi_linear,
      .approach_mgd = "Klesel",
      .R_permutation = 10
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_linear,
      .approach_mgd = "Chin",
      .R_permutation = 10
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_linear,
      .approach_mgd = "Keil",
      .R_permutation = 10
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_linear_boot,
      .approach_mgd = "Sarstedt",
      .R_permutation = 10,
      .handle_inadmissibles = "replace"
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_linear_boot,
      .approach_mgd = "all",
      .R_permutation = 10,
      .handle_inadmissibles = "replace"
    )
  )
})

test_that("Chin, Keil and Sarstedt work for nonlinear models", {
  expect_output(
    testMGD(
      .object       = res_multi_nonlinear,
      .approach_mgd = "Chin",
      .R_permutation = 10
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_nonlinear,
      .approach_mgd = "Keil",
      .R_permutation = 10
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_nonlinear_boot,
      .approach_mgd = "Sarstedt",
      .R_permutation = 10,
      .handle_inadmissibles = "replace"
    )
  )
})

test_that("All .approach_mgd options work for '2ndorder' models", {
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder,
      .approach_mgd = "Klesel",
      .R_permutation = 10,
      .handle_inadmissibles = "replace" # to make sure we have enough admissibles
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder,
      .approach_mgd = "Chin",
      .R_permutation = 10,
      .handle_inadmissibles = "replace" # to make sure we have enough admissibles
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder,
      .approach_mgd = "Keil",
      .R_permutation = 10,
      .handle_inadmissibles = "replace" # to make sure we have enough admissibles
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder_boot,
      .approach_mgd = "Sarstedt",
      .R_permutation = 10,
      .handle_inadmissibles = "replace"
    )
  )
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder_boot,
      .approach_mgd = "all",
      .R_permutation = 10,
      .handle_inadmissibles = "replace"
    )
  )
})
