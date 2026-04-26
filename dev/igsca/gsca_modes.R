
debugonce(initializeAlsEstimates)
debugonce(calculateWeightsGSCA)

debugonce(updateCBDU)

model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., two are common factors and one is a composite
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

igsca <- csem(threecommonfactors, model, .approach_weights = "GSCA", .disattenuate = TRUE, .GSCA_modes = "CCMP")


model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., two are common factors and one is a composite
eta1 <~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 <~ y31 + y32 + y33
"

GSCA <- csem(
    threecommonfactors,
    model,
    .approach_weights = "GSCA",
    .disattenuate = FALSE
)

GSCA_ncmp <- csem(
    threecommonfactors,
    model,
    .approach_weights = "GSCA",
    .disattenuate = FALSE,
    .GSCA_modes = "NCMP"
)

GSCA_ccmp <- csem(
    threecommonfactors,
    model,
    .approach_weights = "GSCA",
    .disattenuate = FALSE,
    .GSCA_modes = "CCMP"
)

GSCA_mix <- csem(
    threecommonfactors,
    model,
    .approach_weights = "GSCA",
    .disattenuate = FALSE,
    .GSCA_modes = list("eta1"  = "CCMP", "eta2" = "NCMP", "eta3" = "CCMP")
)



# igsca <- csem(
#   threecommonfactors,
#   model,
#   .approach_weights = "GSCA",
#   .starting_values = list("eta1" = c(), "eta2" = c(), "eta3" = c())
# )


# .PLS_modes Reference ---------------------------------------------------
modes <- list("eta1" = "unit", "eta2" = "modeB", "eta3" = "unit")
debugonce(calculateWeightsPLS)
res   <- csem(threecommonfactors, model, .PLS_modes = modes, .approach_weights = "PLS-PM")
