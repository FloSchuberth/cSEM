## Skript to be run on ubuntu
require(cSEM)

model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Measurement models
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

# res <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 5000,
#             .eval_plan = "multiprocess")

runif(1)

system.time(csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 1000,
                 .eval_plan = "multisession"))
