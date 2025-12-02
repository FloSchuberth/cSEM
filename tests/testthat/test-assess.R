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
  expect_equal(a$AVE, c(0.4802837, 0.6318118, 0.5557484), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$AIC, c(-209.8370237, -318.7730588), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$AICc, c(292.2113634, 183.3077493), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$AICu, c(-207.8330130, -315.7640227), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$BIC, c(-201.4078075, -306.1292345), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$FPE, c(0.6572610489, 0.5285880026), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$GM, c(511.4292162, 517.6438243), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$HQ, c(-206.5294131, -313.8116428), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$HQc, c(-206.4704807, -313.7009215), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$Mallows_Cp, c(3, 5), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoC, c(0.7339026, 0.8329352, 0.7883384), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoC_mm, c(0.7341254, 0.9446712, 0.7876519), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoC_weighted, c(0.7388285, 1.0000000, 0.7958698 ), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoC_weighted_mm, c(0.7388285, 0.8845458, 0.7958698), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$DG, 0.005499167 , tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$DL, 0.0104148358, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$DML, 0.02928953, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$Df, 22, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(c(a$F2), c(0.53061872, 0.34344601, 0,
                               0.06958792, 0.00000000, 0), 
                    tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$Chi_square, 14.61547, tolerance = 1e-05, ignore_attr = TRUE)
  expect_equal(a$Chi_square_df, 0.6643398, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$CFI, 1, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$GFI, 0.9927503, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$CN, 1159.24462058346, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$IFI, 1.005165, tolerance = 1e-06, ignore_attr = TRUE)
  expect_equal(a$NFI, 0.9899318 , tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$NNFI, 1, tolerance = 1e-06, ignore_attr = TRUE)
  expect_equal(a$RMSEA, 0, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RMS_theta, 0.08470806, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$SRMR, 0.0152131783, tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(c(a$`Fornell-Larcker`), c(0.4802837, 0.3466694, 0.4402532, 
                                              0.3466694, 0.6318118, 0.2969352,
                                              0.4402532, 0.2969352, 0.5557484), 
                    tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$GoF, 0.4784006 , tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(c(a$HTMT$htmts), c(1, 0.6782752, 0.6668841,
                                 0, 1.0000000, 0.6124418,
                                 0, 0.0000000, 1.0000000), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(c(a$HTMT2$htmts), c(1.0000000, 0.6724003, 0.6652760, 
                                         0.0000000, 1.0000000, 0.5958725,
                                         0.0000000, 0.0000000, 1.0000000), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$R2, c(0.3466694, 0.4766706), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$R2_adj, c(0.3453575, 0.4745646), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoT, c(0.7317595, 0.7280738, 0.7859542), tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(a$RhoT_weighted, c(0.7288400, 0.6637134, 0.7821770) , tolerance = 1e-07, ignore_attr = TRUE)
  expect_equal(c(a$VIF), c(1.530619, 1.530619), tolerance = 1e-06, ignore_attr = TRUE)
  expect_equal(c(a$VIF_modeB), c(1.271908, 1.619737, 1.620424), tolerance = 1e-06, ignore_attr = TRUE)
})

## Assess using additional arguments
assess(res,   
            .absolute            = TRUE,
            .alpha               = 0.05,
            .ci                  = "CI_percentile",
            .closed_form_ci      = FALSE,
            .handle_inadmissibles= "drop",
            .inference           = TRUE,
            .null_model          = FALSE,
            .R                   = 199,
            .saturated           = FALSE,
            .seed                = NULL,
            .type_gfi            = "ML",
            .type_vcv            = "indicator"
            )

### Test assess for other classes ----------------------------------------------
source("test-main.R")

test_that("assess() works for list of data", {
  expect_equal(class(assess(res_multi_linear)), "list", ignore_attr = TRUE)  
  expect_error(class(assess(res_multi_nonlinear)))
  expect_equal(class(assess(res_multi_2ndorder)), "list", ignore_attr = TRUE)
})

### Test individual assess's helper functions individually ---------------------

## calculateAVE()
test_that("Test that calculateAVE() returns the correct output:", {
  expect_length(calculateAVE(res_single_linear), 2)
  expect_named(calculateAVE(res_single_linear), c("eta1", "eta3"))
  expect_length(calculateAVE(res_single_nonlinear), 2)
  expect_named(calculateAVE(res_single_nonlinear), c("eta1", "eta3"))
  expect_length(calculateAVE(res_single_2ndorder), 3)
  expect_named(calculateAVE(res_single_2ndorder), c("eta1", "eta2", "c4"))
  
  expect_equal(class(calculateAVE(res_multi_linear)), "list", ignore_attr = TRUE)
  expect_equal(class(calculateAVE(res_multi_nonlinear)), "list", ignore_attr = TRUE)
  expect_equal(class(calculateAVE(res_multi_2ndorder)), "list", ignore_attr = TRUE)
})

## calculateDf()
test_that("Test that calculateDf() returns the correct output:", {
  expect_identical(calculateDf(res_single_linear), 22)
  expect_warning(calculateDf(res_single_nonlinear))
  expect_identical(calculateDf(res_single_2ndorder), 132)
  expect_warning(calculateDf(res_single_nonlinear_2ndorder))
  
  expect_identical(unlist(calculateDf(res_multi_linear)), 
                   c("group1" = 22, "group2" = 22, "group3" = 22))
  expect_warning(calculateDf(res_multi_nonlinear))
  expect_identical(unlist(calculateDf(res_multi_2ndorder)), 
                   c("group1" = 132, "group2" = 132))
})

## Calculatef2
test_that("Test that calculatef2() returns the correct output:", {
  expect_true(inherits(calculatef2(res_single_linear), "matrix"))
  expect_true(inherits(calculatef2(res_single_nonlinear), "matrix"))
  expect_true(inherits(calculatef2(res_single_2ndorder), "matrix"))
  expect_true(inherits(calculatef2(res_single_nonlinear_2ndorder), "matrix"))
  
  expect_equal(class(calculatef2(res_multi_linear)), "list", ignore_attr = TRUE)
  expect_equal(class(calculatef2(res_multi_nonlinear)), "list", ignore_attr = TRUE)
  expect_equal(class(calculatef2(res_multi_2ndorder)), "list", ignore_attr = TRUE)
})

## Fit indices
test_that("Test that calculateGFI returns the correct output:", {
  expect_length(calculateGFI(res_single_linear), 1)
  expect_length(calculateGFI(res_single_linear, .type = "ULS"), 1)
  expect_error(calculateGFI(res_single_nonlinear))
  expect_error(calculateGFI(res_single_nonlinear, .type = "ULS"))
  expect_length(calculateGFI(res_single_2ndorder), 1)
  expect_length(calculateGFI(res_single_2ndorder, .type = "ULS"), 1)
  
})