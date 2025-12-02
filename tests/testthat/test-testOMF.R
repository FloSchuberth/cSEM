
DGPs <- list.files(testthat::test_path("data/"), pattern = "DGP_")

for(i in DGPs) {
  ## Model and Sigma matrix
  load(testthat::test_path(paste0("data/", i)))
  
  ## Draw data
  dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)
  
  ## Estimate
  res <- csem(dat, model_Sigma) 
  
  ## Test
  test_that(paste("testOMF works for DGP: ", i, "with default values"),  {
      testOMF(
        .object = res,
        .R      = 4,
        .handle_inadmissibles = "replace" # to make sure there are enough admissibles
    )
  })
  
  test_that(paste("All arguments of testOMF work for DGP: ", i),  {
      testOMF(
        .object = res,
        .R      = 4,
        .fit_measures = TRUE,
        .alpha  = c(0.1, 0.05),
        .handle_inadmissibles = "replace", # to make sure there are enough admissibles
        .seed   = 2010
      )
  })
}

## Checks that dont need to be checked for all DGPS:

test_that(paste(".seed in testOMF works corretly"),  {
  # Save .Random.seed before calling testOMF()
  r1 <- .Random.seed
  
  a <- testOMF(
    .object = res,
    .R      = 10,
    .seed   = 1303
  )
  
  # Save after calling testOMF()
  r2 <- .Random.seed
  
  b <- testOMF(
    .object = res,
    .R      = 10,
    .seed   = 1303
  )
  
  # .seed should produce the same results
  expect_equal(a$Information$Bootstrap_values, b$Information$Bootstrap_values)
  # .seed does not affect the global seed
  expect_identical(r1, r2)
})


