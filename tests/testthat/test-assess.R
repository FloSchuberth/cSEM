context("assess")

## Linear
model_linear <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Single data set
res    <- csem(threecommonfactors, model_linear, .PLS_weight_scheme_inner = "factorial")
a <- assess(res, .only_common_factors = FALSE)

test_that("Assess works for all choices of .quality_criterion", {
  expect_identical(class(a), "cSEMAssess")
  expect_equivalent(a$AVE, c(0.4802837, 0.6318118, 0.5557484), tolerance = 1e-07)
  expect_equivalent(a$RhoC, c(0.7339026, 0.8329352, 0.7883384), tolerance = 1e-07)
  expect_equivalent(a$RhoC_mm, c(0.7341254, 0.9446712, 0.7876519), tolerance = 1e-07)
  expect_equivalent(a$RhoC_weighted, c(0.7388285, 1.0000000, 0.7958698 ), tolerance = 1e-07)
  expect_equivalent(a$RhoC_weighted_mm, c(0.7388285, 0.8845458, 0.7958698), tolerance = 1e-07)
  expect_equivalent(a$DG, 0.005499167 , tolerance = 1e-07)
  expect_equivalent(a$DL, 0.01041484, tolerance = 1e-07)
  expect_equivalent(a$DML, 0.02928953, tolerance = 1e-07)
  expect_equivalent(a$Df, 18, tolerance = 1e-07)
  expect_equivalent(c(a$F2), c(0.53061872, 0.34344601, 0,
                               0.06958792, 0.00000000, 0), 
                    tolerance = 1e-07)
  expect_equivalent(a$Chi_square, 14.61547, tolerance = 1e-05)
  expect_equivalent(a$Chi_square_df, 0.8119708, tolerance = 1e-07)
  expect_equivalent(a$CFI, 1, tolerance = 1e-07)
  expect_equivalent(a$GFI, 0.9989305, tolerance = 1e-07)
  expect_equivalent(a$IFI, 1.002361, tolerance = 1e-06)
  expect_equivalent(a$NFI, 0.9899318 , tolerance = 1e-07)
  expect_equivalent(a$NNFI, 1.004782, tolerance = 1e-06)
  expect_equivalent(a$RMSEA, 0, tolerance = 1e-07)
  expect_equivalent(a$RMS_theta, 0.08470806, tolerance = 1e-07)
  expect_equivalent(a$RMS_theta_mi, 0.08470806, tolerance = 1e-07)
  expect_equivalent(a$SRMR, 0.01521318, tolerance = 1e-07)
  expect_equivalent(c(a$`Fornell-Larcker`), c(0.4802837, 0.3466694, 0.4402532, 
                                              0.3466694, 0.6318118, 0.2969352,
                                              0.4402532, 0.2969352, 0.5557484), 
                    tolerance = 1e-07)
  expect_equivalent(a$GoF, 0.4784006 , tolerance = 1e-07)
  expect_equivalent(c(a$HTMT), c(0, 0.6782752, 0.6668841,
                                 0, 0.0000000, 0.6124418,
                                 0, 0.0000000, 0.0000000), tolerance = 1e-07)
  expect_equivalent(a$R2, c(0.3466694, 0.4766706), tolerance = 1e-07)
  expect_equivalent(a$R2_adj, c(0.3453575, 0.4745646), tolerance = 1e-07)
  expect_equivalent(a$RhoT, c(0.7317595, 0.7280738, 0.7859542), tolerance = 1e-07)
  expect_equivalent(a$RhoT_weighted, c(0.7288400, 0.6637134, 0.7821770) , tolerance = 1e-07)
  expect_equivalent(c(a$VIF), c(1.530619, 1.530619), tolerance = 1e-06)
  expect_equivalent(c(a$VIF_modeB$eta2), c(1.859990, 2.621641, 2.623095), tolerance = 1e-06)
})

### Test individual assess's helper functions individually ---------------------
source("test-main.R")

## Calculatef2
test_that("Test that calculatef2() returns the correct output:", {
  expect_true(inherits(calculatef2(res_single_linear), "matrix"))
  expect_true(inherits(calculatef2(res_single_nonlinear), "matrix"))
  expect_true(inherits(calculatef2(res_single_2ndorder), "matrix"))
  
  expect_equal(class(calculatef2(res_multi_linear)), "list")
  expect_equal(class(calculatef2(res_multi_nonlinear)), "list")
  expect_equal(class(calculatef2(res_multi_2ndorder)), "list")
})