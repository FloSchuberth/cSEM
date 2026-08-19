# context("testMGD")

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


# Test whethter testMICOM produces same results for list and ordered datasets with an id varibale. 

set.seed(20240101)
shuffle_idx <- sample(seq_len(nrow(threecommonfactors_id)))
threecommonfactors_id_unordered <- threecommonfactors_id[shuffle_idx, ]
rownames(threecommonfactors_id_unordered) <- NULL

test_that("sanity check: the shuffled id data set is actually unsorted by Group_id", {
  expect_true(is.unsorted(threecommonfactors_id_unordered$Group_id))
})

res_multi_id_linear_unordered <- csem(threecommonfactors_id_unordered, model_linear,
                                      .id = "Group_id")

test_that("csem() per-group estimates are unaffected by row order in the .id data set", {
  expect_equal(
    lapply(res_multi_id_linear, function(x) x$Estimates$Path_estimates),
    lapply(res_multi_id_linear_unordered, function(x) x$Estimates$Path_estimates),
    tolerance = 1e-8
  )
})

test_that("testMICOM() Step 2 test statistic is unaffected by row order in the .id data set", {
  out_ordered <- testMICOM(res_multi_id_linear, .R = 5, .seed = 1,
                           .handle_inadmissibles = "replace", .verbose = FALSE)
  out_unordered <- testMICOM(res_multi_id_linear_unordered, .R = 5, .seed = 1,
                             .handle_inadmissibles = "replace", .verbose = FALSE)
  
  expect_equal(
    unname(out_ordered$Step2$Test_statistic),
    unname(out_unordered$Step2$Test_statistic),
    tolerance = 1e-8
  )
})

test_that("testMICOM() Step 3 test statistics match between a list-built object and an .id-built object on row-shuffled data", {
  out_list <- testMICOM(res_multi_linear, .R = 5, .seed = 1,
                        .handle_inadmissibles = "replace", .verbose = FALSE)
  out_id_unordered <- testMICOM(res_multi_id_linear_unordered, .R = 5, .seed = 1,
                                .handle_inadmissibles = "replace", .verbose = FALSE)
  
  expect_equal(
    unname(out_list$Step3$Mean$Test_statistic),
    unname(out_id_unordered$Step3$Mean$Test_statistic),
    tolerance = 1e-8
  )
  
  expect_equal(
    unname(out_list$Step3$Var$Test_statistic),
    unname(out_id_unordered$Step3$Var$Test_statistic),
    tolerance = 1e-8
  )
})

test_that("testMICOM() Step 3 p-values match between list-built and id-built (shuffled) objects", {
  out_list <- testMICOM(res_multi_linear, .R = 5, .seed = 1,
                        .handle_inadmissibles = "replace", .verbose = FALSE)
  out_id_unordered <- testMICOM(res_multi_id_linear_unordered, .R = 5, .seed = 1,
                                .handle_inadmissibles = "replace", .verbose = FALSE)
  
  expect_equal(out_list$Step3$Mean$P_value, out_id_unordered$Step3$Mean$P_value, tolerance = 1e-8)
  expect_equal(out_list$Step3$Var$P_value, out_id_unordered$Step3$Var$P_value, tolerance = 1e-8)
})
