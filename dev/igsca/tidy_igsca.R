
# Simple Example ---------------------------------------------------------

model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

res_boot <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 40)

tidy(res, conf.int = TRUE, conf.level = .95, conf.method = "CI_percentile")

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

## Glance Examples ------------------------------------------------------

res <- csem(threecommonfactors, model)

glance(res)

glance(res_mg)


# IGSCA ------------------------------------------------------------------



tutorial_igsca_model <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
Numberofjobinterviews <~ Interview1 + Interview2
Numberofjoboffers <~ Offer1 + Offer2

# Reflective Measurement Model
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10

# Structural Model
NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
Numberofjobinterviews ~ NetworkingBehavior
Numberofjoboffers ~ NetworkingBehavior
"

## Compute and tabulate igsca ----------------------------------------------------
igsca_mod <- csem(
  .data = LeDang2022,
  .model = tutorial_igsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)



debugonce(tidy.cSEMResults)
debugonce(summarize)
tidy(igsca_mod, conf.int = TRUE, conf.level = .95)






# Multigroup Things ------------------------------------------------------
igsca_mod_mg <- csem(
  .data = LeDang2022,
  .id = "Gender",
  .model = tutorial_igsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

debugonce(tidy.cSEMResults)
debugonce(summarize)
tidy(igsca_mod_mg)

# debugonce(assess)
debugonce(glance.cSEMResults)
glance(igsca_mod_mg)
