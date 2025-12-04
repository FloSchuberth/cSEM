source(testthat::test_path("test-main.R"))

# ==============================================================================
# Tests 
# ==============================================================================
### Default
test_that("predict() fails for non-linear and 2ndorder models", {
  expect_error(predict(res_single_nonlinear, .cv_folds = 2, .r = 2))
  expect_error(predict(res_single_2ndorder, .cv_folds = 2, .r = 2))
  expect_error(predict(res_single_nonlinear_2ndorder, .cv_folds = 2, .r = 2))
  
  expect_error(predict(res_multi_nonlinear, .cv_folds = 2, .r = 2))
  expect_error(predict(res_multi_2ndorder, .cv_folds = 2, .r = 2))
  
  expect_error(predict(res_multi_id_nonlinear, .cv_folds = 2, .r = 2))
  expect_error(predict(res_multi_id_2ndorder, .cv_folds = 2, .r = 2))
})

test_that("predict() works for linear models", {
  expect_s3_class({
    predict(res_single_linear, .cv_folds = 2, .r = 2, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(res_multi_linear, .cv_folds = 2, .r = 2, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(res_multi_id_linear, .cv_folds = 2, .r = 2, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
})

test_that("predict() works with ordinal categorical indicators", {
  expect_s3_class({
    predict(res_OrdPLS_all_ordinal, .cv_folds = 2, .r = 2, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(res_OrdPLS_ordinal_continuous, .cv_folds = 2, .r = 2, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
})


### Using test data
index <- sample(1:500, 400, replace = FALSE)
train_dat <- cSEM::threecommonfactors[index, ]
test_dat  <- cSEM::threecommonfactors[-index, ]

test_that("predict() works for linear models if test data is given ", {
  expect_s3_class({
    predict(csem(train_dat, model_linear), .test_data = test_dat, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(csem(dat, model_linear), .test_data = test_dat, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(csem(threecommonfactors_id, model_linear, .id = "Group_id"), .test_data = test_dat, .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
  expect_s3_class({
    predict(res_OrdPLS_ordinal_continuous, .test_data = dat_OrdPLS_ordinal_continuous[1:100,], .benchmark = "lm", .disattenuate = FALSE)
  }, "cSEMPredict")
})

### Check the different benchmarks
## Linear
model_linear2 <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

### Test benchmarks
for(i in args_default(.choices = TRUE)$.benchmark) {
  if(i != "PLS-PM") {
    expect_s3_class({
      predict(csem(threecommonfactors, model_linear2), .cv_folds = 2, .r = 1,
              .benchmark = i, .disattenuate = if (i == "lm") FALSE else TRUE) 
    }, "cSEMPredict")
  }
}
