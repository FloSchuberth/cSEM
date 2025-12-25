model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

# Single Group Example
res_boot <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 40)

tidy(res, conf.int = TRUE, conf.level = .95, conf.method = "CI_percentile")


# Multi-Group Example
threecommonfactors_id <- cbind(
  "id" = sample(1:3, nrow(threecommonfactors), replace = TRUE),
  threecommonfactors
)

res_mg_boot <- csem(
  threecommonfactors_id,
  model,
  .resample_method = "bootstrap",
  .R = 40,
  .id = "id"
)

tidy(res_mg_boot)