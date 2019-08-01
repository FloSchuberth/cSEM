\donttest{# ===========================================================================
# Basic usage
# ===========================================================================
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Estimate using two different estimators
out <- csem(threecommonfactors, model, .approach_weights = "PLS-PM")
out1 <- csem(threecommonfactors, model, .approach_weights = "GSCA")

## Test (.R is small to save time; should be higher in real applications)
testOMF(out, .R = 50)
testOMF(out1, .R = 50)}
