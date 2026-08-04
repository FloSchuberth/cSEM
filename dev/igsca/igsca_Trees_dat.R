# Create Data with noise covariates ------------------------------------------------------------

dat <- readRDS("dev/igsca/igscaTrees.RDS")
covs <- c("z_true", grep("^noise_", names(dat), value = TRUE))

model <- "# Latent variable model
 eta1 =~ x11 + x12 + x13
 eta2 =~ x21 + x22 + x23
 
 # Composite model
 eta3 <~ x31 + x32 + x33 
 eta4 <~ x41 + x42 + x43 + x44

 # Structural model
 eta4 ~ eta3 + eta1
 eta2 ~ eta1 + eta4 + eta3
 "

# Fit single group model -------------------------------------------------
pooled <- csem(
    .data = dat,
    .model = model,
    .disattenuate = TRUE,
    .approach_weights = "GSCA",
    .conv_criterion = "sum_diff_absolute", # Default in gsca_m.m and gsca.m
    .iter_max = 100, 
    .GSCA_modes = "CCMP", 
    .tolerance = 0.001 
)



# Apply Do Trees for Variable selection ---------------------------------------------------------

## Matrix of residuals ----------------------------------------------------
debug(doTrees)
doTrees(
    .object = pooled,
    .covariates = covs,
    .model = model,
    # .data = dat,
    .ctree_control = partykit::ctree_control(),
    .igsca_tree_control = igsca_tree_control(),
    influence = influence_mat,
    splitter = NULL
)
undebug(doTrees)


# TODO: Try to reuse the ctree_control() function as much as possible

## Vector of residuals ----------------------------------------------------
debug(doTrees)
doTrees(
    .object = pooled,
    .covariates = covs,
    .model = model,
    # .data = dat,
    .ctree_control = partykit::ctree_control(),
    .igsca_tree_control = igsca_tree_control(),
    influence = influence_vec,
    splitter = NULL
)
undebug(doTrees)


## FIT Diff ---------------------------------------------------------------


## Squared-Euclidean Distance ---------------------------------------------


## Geodesic Distance ------------------------------------------------------


# Apply Do Trees for Split Point Selection ---------------------------------------------------------


## Matrix of residuals ----------------------------------------------------

## Vector of residuals ----------------------------------------------------


## FIT Diff ---------------------------------------------------------------


## Squared-Euclidean Distance ---------------------------------------------


## Geodesic Distance ------------------------------------------------------

