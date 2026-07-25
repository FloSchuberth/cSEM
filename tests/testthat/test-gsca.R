test_that("GSCA handles correlated composite models", {
  data(threecommonfactors)
  
  model_correlated_composites <- "
  eta1 <~ y11 + y12 + y13
  eta2 <~ y21 + y22 + y23
  eta3 <~ y31 + y32 + y33
  eta1 ~~ eta2
  eta1 ~~ eta3
  eta2 ~~ eta3
  "
  
  res <- csem(
    threecommonfactors,
    model_correlated_composites,
    .approach_weights = "GSCA"
  )
  
  expect_true(res$Information$Weight_info$Convergence_status)
  expect_equal(rownames(res$Estimates$Weight_estimates), c("eta1", "eta2", "eta3"))
})
