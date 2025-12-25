if (FALSE) {model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res <- csem(threecommonfactors, model)

glance(res)

threecommonfactors_id <- cbind(
  "id" = sample(1:3, nrow(threecommonfactors), replace = TRUE),
  threecommonfactors
)

res_mg <- csem(
  threecommonfactors_id,
  model
)

glance(res_mg)}