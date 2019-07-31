################################################################################
#
#   Purpose: central file to be sourced to several test-xxx.R files
#
################################################################################
## Function to compare path and/or loading and/or weight estimates from a cSEM object
## to a vector of population parameters
comparecSEM <- function(.object, .what, .pop_parameters) {
  # .object: cSEM object
  # .what: what to compare
  # .pop_parameters: a vector of population values
  
  x <- cSEM::summarize(.object)
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    x1 <- x$First_stage$Estimates
    x2 <- x$Second_stage$Estimates
  } else {
    x1 <- NULL
    x2 <- x$Estimates
  }
  
  if(.what == "Path_estimates") {
    est <- x2$Path_estimates[, c("Name", "Estimate")]
    
  } else if(.what == "Loading_estimates") {
    est <- rbind(x1$Loading_estimates[, c("Name", "Estimate")],
                 x2$Loading_estimates[, c("Name", "Estimate")])
    
  } else if(.what == "Weight_estimates") {
    ## Compare only weights for composites, since only those weights are
    ## specified when creating the DGP
    x1$Weight_estimates
    est <- rbind(
      x1$Weight_estimates[x1$Weight_estimates$Construct_type == "Composite", 
                          c("Name", "Estimate")],
      x2$Weight_estimates[x2$Weight_estimates$Construct_type == "Composite", 
                          c("Name", "Estimate")])
    
  } else {
    stop("Error") 
  }
  
  data.frame(est, 
             "Pop_value" = unname(.pop_parameters),
             "Pop_value_name" = names(.pop_parameters),
             stringsAsFactors = FALSE)
}


### ----------------------------------------------------------------------------
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

## Nonlinear model
model_nonlinear <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2 + eta1.eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

# Model and Sigma matrix for 2nd order DGP
load("../data/DGP_2ndorder_cf_of_composites.RData")
# load("tests//data/DGP_2ndorder_cf_of_composites.RData")
model_2ndorder <- model_Sigma

## Data
dat <- list(threecommonfactors, 
            threecommonfactors[1:200, ], 
            threecommonfactors[130:250,])

dat2ndorder <- lapply(c(200, 300), function(x) {
  MASS::mvrnorm(x, rep(0, nrow(Sigma$Sigma)), Sigma = Sigma$Sigma, empirical = TRUE)
})

## Estimates (.R is small to save computation time)

## Single data set
res_single_linear    <- csem(threecommonfactors, model_linear)
res_single_nonlinear <- csem(threecommonfactors, model_nonlinear)
res_single_2ndorder  <- csem(dat2ndorder[[1]], model_2ndorder)

## Multiple data sets
res_multi_linear    <- csem(dat, model_linear)
res_multi_nonlinear <- csem(dat, model_nonlinear)
res_multi_2ndorder  <- csem(dat2ndorder, model_2ndorder)

## Single data set including bootstrap 
res_single_linear_boot    <- csem(threecommonfactors, model_linear, 
                                  .resample_method = "bootstrap", .R = 30)
res_single_nonlinear_boot <- csem(threecommonfactors, model_nonlinear, 
                                  .resample_method = "bootstrap", .R = 30)
res_single_2ndorder_boot  <- csem(dat2ndorder[[1]], model_2ndorder, 
                                  .resample_method = "bootstrap", .R = 30)

## Multiple data sets including bootstrap 
res_multi_linear_boot    <- csem(dat, model_linear, 
                                 .resample_method = "bootstrap", .R = 20)
res_multi_nonlinear_boot <- csem(dat, model_nonlinear, 
                                 .resample_method = "bootstrap", .R = 20)
res_multi_2ndorder_boot  <- csem(dat2ndorder, model_2ndorder, 
                                 .resample_method = "bootstrap", .R = 20)
