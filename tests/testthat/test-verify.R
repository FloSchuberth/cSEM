context("verify")

source("test-main.R")

test_that("Verify work for all cSEMResults classes and subclasses", {
  expect_identical(class(verify(res_single_linear)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_nonlinear)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_2ndorder)), 
                   c("cSEMVerify", "cSEMVerify_2ndorder"))
  expect_identical(class(verify(res_single_2ndorder$First_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_2ndorder$Second_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_linear_boot)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_nonlinear_boot)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_single_2ndorder_boot)), 
                   c("cSEMVerify", "cSEMVerify_2ndorder"))
  expect_identical(class(verify(res_single_2ndorder_boot$First_stage)), 
                   "cSEMVerify", "cSEMVerify_2ndorder")
  expect_identical(class(verify(res_single_2ndorder_boot$Second_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_linear)), 
                   c("cSEMVerify", "cSEMVerify_multi"))
  expect_identical(class(verify(res_multi_linear$Data_1)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_nonlinear)), 
                   c("cSEMVerify", "cSEMVerify_multi"))
  expect_identical(class(verify(res_multi_nonlinear$Data_1)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_2ndorder$Data_1)), 
                   c("cSEMVerify", "cSEMVerify_2ndorder"))
  expect_identical(class(verify(res_multi_2ndorder$Data_1$First_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_2ndorder$Data_1$Second_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_linear_boot)), 
                   c("cSEMVerify", "cSEMVerify_multi"))
  expect_identical(class(verify(res_multi_linear_boot$Data_1)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_nonlinear_boot$Data_1)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_2ndorder_boot$Data_1)), 
                   c("cSEMVerify", "cSEMVerify_2ndorder"))
  expect_identical(class(verify(res_multi_2ndorder_boot$Data_1$First_stage)), 
                   "cSEMVerify")
  expect_identical(class(verify(res_multi_2ndorder_boot$Data_1$Second_stage)), 
                   "cSEMVerify")
})

test_that("Not providing a cSEMResults object causes an error", {
  expect_error(verify(list(1:3)))
})