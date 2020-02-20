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


### Models ---------------------------------------------------------------------
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

## Model and Sigma matrix for 2nd order DGP
load("../data/DGP_2ndorder_cf_of_composites.RData")
# load("tests//data/DGP_2ndorder_cf_of_composites.RData") # uncomment to source
# on local machine
model_2ndorder <- model_Sigma

### Data -----------------------------------------------------------------------
# Add unused columns to threecommonfactors to check if they get removed correctly
# when not part of the model
threecommonfactors <- as.data.frame(threecommonfactors)
threecommonfactors$not_used_numeric <- rnorm(nrow(threecommonfactors))
threecommonfactors$not_used_character <- sample(letters, 
                                                size = (nrow(threecommonfactors)), 
                                                replace = TRUE)

# Suffle columns to make sure columns get sorted correctly
threecommonfactors <- threecommonfactors[, sample(1:ncol(threecommonfactors))]

## List of data without id column
dat <- list(
  group1 = threecommonfactors, 
  # Shuffle again to make sure that this gets correctly detected
  group2 = threecommonfactors[1:200, ][, sample(1:ncol(threecommonfactors))], 
  group3 = threecommonfactors[130:250,])

# Remove one column of group2 to make sure ordering also works when the data
# in the list has different number of columns
dat$group2$not_used_numeric <- NULL

## Data with id column
threecommonfactors_id <- as.data.frame(rbind(threecommonfactors,
                                             threecommonfactors[1:200, ],
                                             threecommonfactors[130:250,]))

threecommonfactors_id$Group_id <- rep(c(1, 2, 3), times = c(nrow(threecommonfactors),
                                                            nrow(threecommonfactors[1:200, ]),
                                                            nrow(threecommonfactors[130:250,])))

## Data for 2ndorder model without id column 
dat2ndorder1_a <- dat2ndorder1_b <- as.data.frame(MASS::mvrnorm(100, rep(0, nrow(Sigma$Sigma)), 
                                                                Sigma = Sigma$Sigma, empirical = TRUE))
dat2ndorder2_a <- dat2ndorder2_b <- as.data.frame(MASS::mvrnorm(200, rep(0, nrow(Sigma$Sigma)), 
                                              Sigma = Sigma$Sigma, empirical = TRUE))

# Add unused columns
dat2ndorder1_a$not_used_character <- sample(letters, 
                                          size = (nrow(dat2ndorder1_a)), 
                                          replace = TRUE)
dat2ndorder2_a$not_used_numeric <- rnorm(nrow(dat2ndorder2_a))

## List of data without id column
dat2ndorder <- list(
  group1 = dat2ndorder1_a, 
  # Shuffle again to make sure that this gets correctly detected
  group2 = dat2ndorder2_a[, sample(1:ncol(dat2ndorder2_a))])


## Dat for 2norder model with id column

dat2ndorder_id <- as.data.frame(rbind(dat2ndorder1_b, dat2ndorder2_b))
dat2ndorder_id$group <- rep(c("A", "B"), times  = c(100, 200))

### Estimates (.R is small to save computation time) ---------------------------

## Single data set
res_single_linear      <- csem(threecommonfactors, model_linear)
res_single_nonlinear   <- csem(threecommonfactors, model_nonlinear)
res_single_2ndorder    <- csem(dat2ndorder1_a, model_2ndorder)

## Multiple data sets using list
res_multi_linear       <- csem(dat, model_linear)
res_multi_nonlinear    <- csem(dat, model_nonlinear)
res_multi_2ndorder     <- csem(dat2ndorder, model_2ndorder)

## Multiple data sets using id
res_multi_id_linear    <- csem(threecommonfactors_id, model_linear, .id = "Group_id")
res_multi_id_nonlinear <- csem(threecommonfactors_id, model_nonlinear, .id = "Group_id")
res_multi_id_2ndorder  <- csem(dat2ndorder_id, model_2ndorder, .id = "group")

## Single data set including bootstrap 
res_single_linear_boot    <- csem(threecommonfactors, model_linear, 
                                  .resample_method = "bootstrap", .R = 6,
                                  .handle_inadmissibles = "replace")
res_single_nonlinear_boot <- csem(threecommonfactors, model_nonlinear, 
                                  .resample_method = "bootstrap", .R = 6,
                                  .handle_inadmissibles = "replace")
res_single_2ndorder_boot  <- csem(dat2ndorder1_a, model_2ndorder, 
                                  .resample_method = "bootstrap", .R = 6,
                                  .handle_inadmissibles = "replace")

## Multiple data sets including bootstrap 
res_multi_linear_boot    <- csem(dat, model_linear, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace")
res_multi_nonlinear_boot <- csem(dat, model_nonlinear, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace")
res_multi_2ndorder_boot  <- csem(dat2ndorder, model_2ndorder, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace")

## Multiple data sets using id including bootstrap 
res_multi_id_linear_boot    <- csem(threecommonfactors_id, model_linear, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace",
                                 .id = "Group_id")
res_multi_id_nonlinear_boot <- csem(threecommonfactors_id, model_nonlinear, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace",
                                 .id = "Group_id")
res_multi_id_2ndorder_boot  <- csem(dat2ndorder_id, model_2ndorder, 
                                 .resample_method = "bootstrap", .R = 6,
                                 .handle_inadmissibles = "replace",
                                 .id = "group")

