model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., two are common factors and one is a composite
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

igsca <- csem(threecommonfactors, model, .approach_weights = "GSCA")



# igsca <- csem(
#   threecommonfactors,
#   model,
#   .approach_weights = "GSCA",
#   .starting_values = list("eta1" = c(), "eta2" = c(), "eta3" = c())
# )