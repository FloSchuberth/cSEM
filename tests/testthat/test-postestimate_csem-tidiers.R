# Partly based on broom's lavaan test suite
library(modeltests)
test_that("cSEM tidier arguments", {
  expect_true(modeltests::check_arguments(tidy.cSEMResults))
  expect_true(modeltests::check_arguments(glance.cSEMResults))
})


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

igsca_mod_mg <- csem(
  .data = LeDang2022,
  .id = "Gender",
  .model = tutorial_igsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

gsca_model <- "
# Measurement Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
Numberofjobinterviews <~ Interview1 + Interview2
Numberofjoboffers <~ Offer1 + Offer2
HonestyHumility <~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
Emotionality <~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
Extraversion <~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
Agreeableness <~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
Conscientiousness <~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
OpennesstoExperience <~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10

# Structural Model
NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
Numberofjobinterviews ~ NetworkingBehavior
Numberofjoboffers ~ NetworkingBehavior
"

gsca_mod <- csem(
  .data = LeDang2022,
  gsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

# gsca_mod$Estimates$Loading_estimates |> View()

gsca_m_model <- "
# Measurement Model
NetworkingBehavior =~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
Numberofjobinterviews =~ Interview1 + Interview2
Numberofjoboffers =~ Offer1 + Offer2
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

gsca_m_mod <- csem(
  .data = LeDang2022,
  gsca_m_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

test_that("tidy.cSEMResults for GSCA models", {
  tidy_igsca <- tidy(igsca_mod)
  tidy_igsca_mg <- tidy(igsca_mod_mg)
  tidy_gsca <- tidy(gsca_mod)
  tidy_gsca_m <- tidy(gsca_m_mod)

  expect_true(modeltests::check_tidy_output(tidy_igsca))
  expect_true(modeltests::check_tidy_output(tidy_igsca_mg))
  expect_true(modeltests::check_tidy_output(tidy_gsca))
  expect_true(modeltests::check_tidy_output(tidy_gsca_m))

  modeltests::check_dims(tidy_igsca, 238, 7)
  modeltests::check_dims(tidy_igsca_mg, 476, 8)
  modeltests::check_dims(tidy_gsca, 187, 7)
  modeltests::check_dims(tidy_gsca_m, 249, 7)
})


test_that("glance.cSEMResults for GSCA", {
  glance_igsca <- glance(igsca_mod)
  glance_igsca_mg <- glance(igsca_mod_mg)
  glance_gsca <- glance(gsca_mod)
  glance_gsca_m <- glance(gsca_m_mod)

  check_glance_outputs(glance_igsca)

  expect_error(check_glance_outputs(glance_igsca_mg),
  regexp = "Output column names not in the column glossary: group"
  )
  
  check_glance_outputs(glance_gsca)
  check_glance_outputs(glance_gsca_m)
})


test_that("Confidence interval functionality works", {
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
  res_boot <- csem(
    threecommonfactors,
    model,
    .resample_method = "bootstrap",
    .R = 40
  )

  expect_identical(
    tidy(
      res_boot,
      conf.int = TRUE,
      conf.level = .95,
      conf.type = "CI_percentile"
    ),
    tidy(res_boot, conf.int = TRUE, conf.level = .95, conf.type = NULL)
  )

  tidy_boot <- tidy(res_boot, conf.int = TRUE, conf.level = .95, conf.type = NULL)

  expect_true(modeltests::check_tidy_output(tidy_boot))
  modeltests::check_dims(tidy_boot, 28, 10)

  
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

  expect_identical(
    tidy(res_mg_boot, conf.int = TRUE, conf.type = "CI_percentile"),
    tidy(res_mg_boot, conf.int = TRUE, conf.type = NULL)
  )
  tidy_boot_mg <- tidy(res_mg_boot, conf.int = TRUE, conf.type = "CI_percentile")
  expect_true(modeltests::check_tidy_output(tidy_boot_mg))
  modeltests::check_dims(tidy_boot_mg, 84, 11)
})

