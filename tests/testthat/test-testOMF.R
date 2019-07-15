
DGPs <- list.files("../data/")

for(i in DGPs) {
  ## Model and Sigma matrix
  load(paste0("../data/", i))
  
  ## Draw data
  dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = FALSE)
  
  ## Estimate
  res <- csem(dat, model_Sigma) 
  
  ## Test
  test_that(paste("testOMF works for DGP: ", i, "with default values"),  {
    expect_output(
      testOMF(
        .object = res,
        .R      = 50
      )
    )
  })
  
  test_that(paste("All arguments of testOMF work for DGP: ", i),  {
    expect_output(
      testOMF(
        .object = res,
        .R      = 50,
        .alpha  = c(0.1, 0.05),
        .seed   = 2010
      )
    )
  })
}
