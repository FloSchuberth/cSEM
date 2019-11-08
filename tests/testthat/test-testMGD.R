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

test_that(".approach_mgd = 'Sarstedt' can not be combined with .handle_inadmissibles = 'drop'", {
  expect_error(
    testMGD(res_single_linear, .approach_mgd = "Sarstedt", .handle_inadmissibles = "drop")
  )
  expect_error(
    testMGD(res_nonlinear, .approach_mgd = "all", .handle_inadmissibles = "drop")
  )
})

## Correct
for(i in args_default(.choices = TRUE)$.approach_mgd) { 
  test_that(paste("testMGD works for .approach_mgd = ", i), {
    expect_output(
      testMGD(
        .object       = res_multi_linear_boot,
        .approach_mgd = i,
        .R_permutation = 5,
        .handle_inadmissibles = "replace"
      )
    )
  })
  if(i %in% c("Chin", "Keil", "Nitzl", "Sarstedt", "Henseler"))
  test_that("Chin, Keil, Nitzl and Sarstedt work for nonlinear models", {
    expect_output(
      testMGD(
        .object       = res_multi_nonlinear,
        .approach_mgd = i,
        .R_permutation = 5,
        .handle_inadmissibles = "replace"
      )
    )
  })
}



test_that("testMGD() works for second order models (.approach_mgd = 'all')", {
  expect_output(
    testMGD(
      .object       = res_multi_2ndorder_boot,
      .approach_mgd = "all",
      .R_permutation = 5,
      .handle_inadmissibles = "replace"
    )
  )
})
