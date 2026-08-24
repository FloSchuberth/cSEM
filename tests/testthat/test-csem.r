# context("csem")

## Source 
source("test-main.R")

# ==============================================================================
# General tests 
# ==============================================================================
test_that("No data/model provided causes an error", {
  expect_error(
    csem(.model = model)
  )
  expect_error(
    csem(.data = dat)
  )
})

# ==============================================================================
# DGPs
# ==============================================================================
### DGP_linear_2commonfactors --------------------------------------------------
# Loads Sigma, models and population values
load(file = "../data/DGP_linear_2commonfactors.RData")

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)

## Test
test_that("DPG_linear_2commonfactors is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
})

### DGP_linear_2composites =====================================================
# Loads Sigma, models and population values
load(file = "../data/DGP_linear_2composites.RData")

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Test
test_that("DPG_linear_2composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})

### DGP_linear_3commonfactors ==================================================
# Loads Sigma, models and population values
load(file = "../data/DGP_linear_3commonfactors.RData") 

## Draw data
dat <- MASS::mvrnorm(300, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
for(i in c("PLS-PM", "GSCA", "SUMCORR", "MAXVAR", "MINVAR", "GENVAR", "PCA",
           "unit", "bartlett", "regression")) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularly producing differences
  ##   between estimate and population value larger than 0.01.
  
  if(i == "PLS-PM"){
    for(j in c("modeA","modeB")){
      res <-  csem(dat, model_Sigma, .approach_weights = i,.PLS_modes = j, 
                   .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"))
      
      ## Comparison
      path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
      loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
      
      ## Test
      # Note: the tolerance is necessary since Croon correction requires a CFA (i.e. ML estimation)
      # which is only consistent (i.e. will not exactly reproduce the population values for a finite sample size)
      # Similarly for GSCAm.
      test_that(paste("Weighting approach", i, "yields correct results"), {
        expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
        expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.01)
      })
    }
  }else{
    
    res <-  csem(dat, model_Sigma, .approach_weights = i, 
                 .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"))
    
    ## Comparison
    path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
    loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
    
    ## Test
    # Note: the tolerance is necessary since Croon correction requires a CFA (i.e. ML estimation)
    # which is only consistent (i.e. will not exactly reproduce the population values for a finite sample size)
    # Similarly for GSCAm.
    test_that(paste("Weighting approach", i, "yields correct results"), {
      expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
      expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.01)
    })
  }
  # Export to Excel test
  exportToExcel(assess(res), .filename = paste0("test_assess_", i, ".xlsx"),
                .path = "../test_results_exportToExcel")
  exportToExcel(summarize(res), .filename = paste0("test_summarize_", i, ".xlsx"),
                .path = "../test_results_exportToExcel")
  exportToExcel(predict(res, .handle_inadmissibles = "ignore"), .filename = paste0("test_predict_", i, ".xlsx"),
                .path = "../test_results_exportToExcel")
  exportToExcel(testOMF(res, .R = 10), .filename = paste0("test_testOMF_", i, ".xlsx"),
                .path = "../test_results_exportToExcel")
}

### DGP_linear_3compostites ====================================================
# Loads Sigma, models and population values
load(file = "../data/DGP_linear_3composites.RData") 

## Draw data
dat <- MASS::mvrnorm(300, rep(0, nrow(Sigma$Sigma)), 
                     Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
for(i in c("PLS-PM", "GSCA", "SUMCORR", "MAXVAR", "MINVAR", "GENVAR")) {
  ## - "SSQCOR" is excluded as it is rather unstable, regularily producing differences
  ##   between estimate and population value larger than 0.01.
  ## - "PCA" is excluded as weights obtained by PCA are not population weights but 
  ##   simply the first principal component of S_jj.
  ## - "unit" is excluded as unit weights are inconsitent "estimates" for the
  ##   population weights (weights are simply set to 1 and scaled).
  ## - "bartlett" and "regression" are excluded as they are not meaningful 
  ##   for models containing concepts modeled as composites
  res <-  csem(dat, model_Sigma, .approach_weights = i, 
               .dominant_indicators = c("eta1" = "y11", "eta2" = "y21", "eta3" = "y31"))
  
  ## Comparison
  path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
  loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
  weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)
  
  ## Test
  test_that(paste("Weighting approach", i, "yields correct results"), {
    expect_equal(path$Estimate, path$Pop_value, tolerance = 0.001)
    expect_equal(loadings$Estimate, loadings$Pop_value, tolerance = 0.001)
    expect_equal(weights$Estimate, weights$Pop_value, tolerance = 0.001)
  })
}

### DGP_2ndorder - Common factor of common factors =============================
# Loads Sigma, models and population values
load(file = "../data/DGP_2ndorder_cf_of_cfs.RData")  

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_cf_of_cfs is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
})

# Export to Excel test
exportToExcel(summarize(res), .filename = "test_summarize", .path = "../test_results_exportToExcel")
exportToExcel(assess(res), .filename = "test_assess", .path = "../test_results_exportToExcel")
exportToExcel(testOMF(res, .R = 20), .filename = "test_testOMF", .path = "../test_results_exportToExcel")



### DGP_2ndorder - Common factor of composites =================================
# Loads Sigma, models and population values
load(file = "../data/DGP_2ndorder_cf_of_composites.RData")  

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_cf_of_composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})

### DGP_2ndorder - Composite of common factors =================================
# Loads Sigma, models and population values
load(file = "../data/DGP_2ndorder_composite_of_cfs.RData")  

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_composites_of_cfs is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value, tolerance = 0.01)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})
### DGP_2ndorder - Composite of composites =====================================
# Loads Sigma, models and population values
load(file = "../data/DGP_2ndorder_composite_of_composites.RData")  

## Draw data
dat <- MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)

## Estimate
res <-  csem(dat, model_Sigma)

## Comparison
path     <- comparecSEM(res, .what = "Path_estimates", pop_params_Sigma$Path_coefficients)
loadings <- comparecSEM(res, .what = "Loading_estimates", pop_params_Sigma$Loadings)
weights  <- comparecSEM(res, .what = "Weight_estimates", pop_params_Sigma$Weights)

## Reorder to match estimate and population value
l1 <- loadings[, c("Pop_value", "Pop_value_name")]
l1 <- l1[c(15:17,1:14,18:23), ]
loadings <- cbind(loadings[, 1:2], l1)

## Test
test_that("DPG_2ndorder_composites_of_composites is correctly estimated", {
  expect_equal(path$Estimate, path$Pop_value)
  expect_equal(loadings$Estimate, loadings$Pop_value)
  expect_equal(weights$Estimate, weights$Pop_value)
})

### Non-linear model with cubic term ===========================================
# helper function to generate data
simulate_nonlinear_data <- function(n, seed) {
  set.seed(seed)
  
  rho12 <- -0.30
  gamma <- c(
    gamma1   =  0.50, gamma2   = -0.30, gamma12  = -0.20,
    gamma11  =  0.10, gamma22  = -0.15,
    gamma111 =  0.08, gamma222 = -0.06,
    gamma112 =  0.12, gamma122 = -0.10
  )

  var_signal <- 0.6827232
  sd_zeta <- sqrt(1 - var_signal)  # forces Var(eta3) = 1; R^2 ends up ~0.68
  loading <- 0.80
  err_sd  <- sqrt(1 - loading^2)
  
  L <- chol(matrix(c(1, rho12, rho12, 1), 2, 2))
  eta_exo <- matrix(rnorm(2 * n), n) %*% L
  eta1 <- eta_exo[, 1]
  eta2 <- eta_exo[, 2]
  
  zeta <- rnorm(n, sd = sd_zeta)
  eta3 <- gamma["gamma1"] * eta1 + gamma["gamma2"] * eta2 +
    gamma["gamma12"] * eta1 * eta2 +
    gamma["gamma11"] * eta1^2 + gamma["gamma22"] * eta2^2 +
    gamma["gamma111"] * eta1^3 + gamma["gamma222"] * eta2^3 +
    gamma["gamma112"] * eta1^2 * eta2 + gamma["gamma122"] * eta1 * eta2^2 +
    zeta
  
  make_indicators <- function(eta, prefix) {
    d <- data.frame(
      v1 = loading * eta + rnorm(n, sd = err_sd),
      v2 = loading * eta + rnorm(n, sd = err_sd),
      v3 = loading * eta + rnorm(n, sd = err_sd)
    )
    names(d) <- paste0(prefix, 1:3)
    d
  }
  
  dat <- cbind(
    make_indicators(eta1, "y1_"),
    make_indicators(eta2, "y2_"),
    make_indicators(eta3, "y3_")
  )
  
  list(data = dat, truth = gamma)
}

# helper function to extract results
get_est <- function(name) {
  tab <- s$Estimates$Path_estimates
  row <- tab[tab$Name == name, ]
  if (nrow(row) != 1) {
    stop(
      "Could not find a unique path named '", name, "' in Path_estimates.\n",
      "Available names:\n", paste(tab$Name, collapse = "\n"),
      "\nAdjust get_est()/the Name strings below to match."
    )
  }
  row$Estimate
}

# Allowed relative difference between estimate and true value
tol <- 0.06   

# generate data
data_nonlinear_cubic <- simulate_nonlinear_data(n = 100000, seed = 20260819)

# Estimate model
res_single_nonlinear_cubic <- csem(.data = data_nonlinear_cubic$data, .model = model_nonlinear_cubic)
s   <- summarize(res_single_nonlinear_cubic)


test_that("csem estimates the linear effects correctly", {
  expect_equal(get_est("ETA3 ~ ETA1"), 0.50, tolerance = tol)
  expect_equal(get_est("ETA3 ~ ETA2"), -0.30, tolerance = tol)
})

test_that("csem estimates the quadratic and two-way interaction effects correctly", {
  expect_equal(get_est("ETA3 ~ ETA1.ETA2"), -0.20, tolerance = tol)
  expect_equal(get_est("ETA3 ~ ETA1.ETA1"), 0.10, tolerance = tol)
  expect_equal(get_est("ETA3 ~ ETA2.ETA2"), -0.15, tolerance = tol)
})

test_that("csem recovers the cubic effects (exercises CubicCubic)", {
  expect_equal(get_est("ETA3 ~ ETA1.ETA1.ETA1"), 0.08, tolerance = tol)
  expect_equal(get_est("ETA3 ~ ETA2.ETA2.ETA2"), -0.06, tolerance = tol)
})

test_that("csem recovers the quadratic-times-linear effects", {
  expect_equal(get_est("ETA3 ~ ETA1.ETA1.ETA2"), 0.12, tolerance = tol)
  expect_equal(get_est("ETA3 ~ ETA2.ETA2.ETA1"), -0.10, tolerance = tol)
})
