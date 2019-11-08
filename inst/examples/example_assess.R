\donttest{# ===========================================================================
# Using the threecommonfactors dataset
# ===========================================================================
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept os measured by 3 indicators, i.e., modeled as latent variable
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res <- csem(threecommonfactors, model)
a   <- assess(res) # computes all quality criteria (.quality_criterion = "all")
a

## The return value is a named list
str(a)
a$HTMT

# You may also just compute a subset of quality criteria
assess(res, .quality_criterion = c("ave", "rho_C", "htmt"))

## Resampling ---------------------------------------------------------------
# To resample a given quality criterion use csem()'s .user_funs argument
res <- csem(threecommonfactors, model, 
            .resample_method = "bootstrap", 
            .user_funs       = cSEM:::calculateHTMT
)

## Look at the resamples
res$Estimates$Estimates_resample$Estimates1$User_fun$Resampled[1:4, ]

## Use infer() to compute e.g. the 95% percentile confidence interval
res_infer <- infer(res, .quantity = "CI_percentile")
res_infer$User_fun 
}

  
# Note that .user_funs expects a function that returns a vector or a matrix!