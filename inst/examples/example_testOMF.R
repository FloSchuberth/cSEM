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

## Estimate
out <- csem(threecommonfactors, model, .approach_weights = "PLS-PM")

## Test
testOMF(out, .R = 50, .seed = 320)
}
