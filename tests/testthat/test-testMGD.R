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
  expect_error(testMGD(res_multi_nonlinear, .approach_mgd = "Klesel"))
  expect_error(testMGD(res_multi_nonlinear_boot, .approach_mgd = "all"))
  expect_error(testMGD(res_multi_nonlinear_2ndorder, .approach_mgd = "all"))
})

test_that(".approach_mgd = 'Sarstedt' can not be combined with .handle_inadmissibles = 'drop'", {
  expect_error(
    testMGD(res_multi_linear, .approach_mgd = "Sarstedt", .handle_inadmissibles = "drop")
  )
  expect_error(
    testMGD(res_multi_nonlinear, .approach_mgd = "all", .handle_inadmissibles = "drop")
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
        .R_bootstrap = 5,
        .handle_inadmissibles = "replace"
      )
    )
    
    expect_output(
      testMGD(
        .object       = res_multi_nonlinear_2ndorder,
        .approach_mgd = i,
        .R_permutation = 5,
        .R_bootstrap = 5,
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

# test_that("does have one (and only one) of the following types (correction, type_ci, distance", {
#   expect_equal(
#     testMGD(
#       .output_type = "data.frame"
# 
# results$checkNA <- results %>% select(correction, type_ci, distance) %>% 
#   apply(1, function(x) sum(is.na(x))) 
# 
#     )
#   )
# })


# ==============================================================================
# dropNAResamples(): removal of permutation runs containing NA
# ==============================================================================
# Tests for the NA filter used to post-process the permutation
# reference distribution in testMGD(). The bug being guarded against: the filter
# was originally written as Filter(Negate(anyNA), ref_dist), which misses NAs
# nested inside the $Chin and $Klesel elements because anyNA() defaults to
# recursive = FALSE. Such an NA then reached the p-value computation.

# Build a permutation run in the shape produced inside testMGD()'s repeat loop.
make_run <- function(klesel = c("dG" = 0.10, "dL" = 0.20),
                     chin   = c("eta2 ~ eta1" = 0.05, "eta3 ~ eta1" = -0.03),
                     sarstedt = c("F" = 1.5)) {
  list(
    Klesel   = klesel,
    Chin     = list("group1_group2" = chin),   # nested one level deeper
    Sarstedt = sarstedt
  )
}

test_that("dropNAResamples() removes runs with an NA nested inside $Chin", {
  ref_dist <- list(
    run1 = make_run(),
    run2 = make_run(chin = c("eta2 ~ eta1" = NA_real_, "eta3 ~ eta1" = -0.03)),
    run3 = make_run()
  )
  
  out <- cSEM:::dropNAResamples(ref_dist)
  
  expect_length(out, 2)
  expect_named(out, c("run1", "run3"))
  expect_identical(out$run1, ref_dist$run1)
  expect_identical(out$run3, ref_dist$run3)
  expect_false(anyNA(unlist(out)))
})

test_that("the original filter did NOT remove a nested NA (documents the bug)", {
  # Characterisation test. If this ever starts failing, anyNA()'s treatment of
  # lists has changed and the rationale in ?dropNAResamples needs revisiting.
  ref_dist <- list(
    run1 = make_run(),
    run2 = make_run(chin = c("eta2 ~ eta1" = NA_real_, "eta3 ~ eta1" = -0.03)),
    run3 = make_run()
  )
  
  expect_length(Filter(Negate(anyNA), ref_dist), 3)  # old behaviour: kept it
  expect_length(dropNAResamples(ref_dist), 2)        # new behaviour: drops it
})

test_that("dropNAResamples() removes an NA in the multi-element $Klesel vector", {
  # anyNA(recursive = FALSE) only flags *length-one* atomic elements, so a
  # length-two Klesel vector containing NA was missed as well.
  ref_dist <- list(
    run1 = make_run(),
    run2 = make_run(klesel = c("dG" = NA_real_, "dL" = 0.20))
  )
  
  expect_length(Filter(Negate(anyNA), ref_dist), 2)  # old behaviour
  expect_length(dropNAResamples(ref_dist), 1)        # new behaviour
  expect_named(dropNAResamples(ref_dist), "run1")
})

test_that("dropNAResamples() still removes the literal NA sentinel", {
  # Guards against over-fixing: inadmissible runs under
  # .handle_inadmissibles = 'drop' are stored as a bare NA and must still go.
  ref_dist <- list(run1 = make_run(), run2 = NA, run3 = make_run())
  
  expect_length(dropNAResamples(ref_dist), 2)
  expect_named(dropNAResamples(ref_dist), c("run1", "run3"))
})

test_that("dropNAResamples() treats NaN like NA", {
  # Near-singular matrices tend to produce NaN (0/0) rather than NA.
  ref_dist <- list(
    run1 = make_run(),
    run2 = make_run(chin = c("eta2 ~ eta1" = NaN, "eta3 ~ eta1" = -0.03)),
    run3 = make_run(sarstedt = c("F" = NaN))
  )
  
  expect_length(dropNAResamples(ref_dist), 1)
  expect_named(dropNAResamples(ref_dist), "run1")
})

test_that("dropNAResamples() leaves a clean reference distribution untouched", {
  ref_dist <- list(run1 = make_run(), run2 = make_run(), run3 = make_run())
  expect_identical(dropNAResamples(ref_dist), ref_dist)
})

test_that("dropNAResamples() handles the degenerate inputs", {
  expect_identical(dropNAResamples(list()), list())
  expect_length(dropNAResamples(list(a = NA, b = NA)), 0)
  # Runs where only some approaches were requested must survive.
  expect_length(dropNAResamples(list(r1 = list(Chin = list(g = c(a = 0.1))))), 1)
})
