# ===========================================================================
# Using the three common factors dataset
# ===========================================================================
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., modeled as latent variable
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res <- csem(threecommonfactors, model)
a   <- assess(res) # computes all quality criteria (.quality_criterion = "all")
a

## The return value is a named list. Type for example:
a$HTMT

# You may also just compute a subset of the quality criteria
assess(res, .quality_criterion = c("ave", "rho_C", "htmt"))

\donttest{## Resampling ---------------------------------------------------------------
# To resample a given quality criterion use csem()'s .user_funs argument
# Note: The output of the quality criterion needs to be a vector or a matrix.
#       Matrices will be vectorized columnwise.
res <- csem(threecommonfactors, model, 
            .resample_method = "bootstrap", 
            .R               = 40,
            .user_funs       = cSEM:::calculateSRMR
)

## Look at the resamples
res$Estimates$Estimates_resample$Estimates1$User_fun$Resampled[1:4, ]

## Use infer() to compute e.g., the 95% percentile confidence interval
res_infer <- infer(res, .quantity = "CI_percentile")

## The results are saved under the name "User_fun"
res_infer$User_fun 

## Several quality criteria can be resampled simultaneously
res <- csem(threecommonfactors, model, 
            .resample_method = "bootstrap",
            .R               = 40,
            .user_funs       = list(
              "SRMR" = cSEM:::calculateSRMR,
              "RMS_theta" = cSEM:::calculateRMSTheta
            ),
            .tolerance = 1e-04
)
res$Estimates$Estimates_resample$Estimates1$SRMR$Resampled[1:4, ]
res$Estimates$Estimates_resample$Estimates1$RMS_theta$Resampled[1:4]
}