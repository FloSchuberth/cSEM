context("predict")

source("test-main.R")
# source("tests/testthat/test-main.R")

# ==============================================================================
# Tests 
# ==============================================================================
# Default
predict(res_single_linear, .cv_folds = 2, .r = 2)

# Using test data
index <- sample(1:500, 400, replace = FALSE)
train_dat <- threecommonfactors[index, ]
test_dat  <- threecommonfactors[-index, ] 
rownames(test_dat) <- 1:nrow(test_dat)

predict(csem(train_dat, model_linear), .test_data = test_dat)

# Different benchmarks
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

# Test benchmarks
for(i in args_default(.choices = TRUE)$.benchmark) {
  predict(csem(threecommonfactors, model_linear2), .cv_folds = 2, .r = 1,
          .benchmark = i)
}

