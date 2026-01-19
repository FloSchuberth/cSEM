# GSCA_M -----------------------------------------------------------------
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., modeled as latent variable
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

# debugonce(csem)
# debugonce(resamplecSEMResults)
debugonce(selectAndVectorize) # TODO: This has to be updated to include D
debugonce(resamplecSEMResultsCore) # TODO: The actual resampling
debugonce(resamplecSEMResults) # TODO: This processes the results after resampling
gscam <- csem(
  threecommonfactors,
  model,
  .approach_weights = "GSCA",
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .resample_method = "bootstrap",
  .R = 5
)

(tidied_gscam <- tidy(gscam, conf.int = TRUE))
# View(tidied_gscam)


# IGSCA ------------------------------------------------------------------


model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Each concept is measured by 3 indicators, i.e., two are common factors and one is a composite
eta1 =~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res <- csem(threecommonfactors, model)


