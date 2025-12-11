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
  # FIXME: Write the tests that I'd like for the method here

  # FIXME: Test multigroup models glancing + tidy
  # FIXME: conf.int level TRUE/FALSE
  # FIXME: conf.level supported
  tidy_igsca <- tidy(igsca_mod)
  tidy_gsca <- tidy(gsca_mod)
  tidy_gsca_m <- tidy(gsca_m_mod)

  expect_true(modeltests::check_tidy_output(tidy_igsca))
  expect_true(modeltests::check_tidy_output(tidy_gsca))
  expect_true(modeltests::check_tidy_output(tidy_gsca_m))

  modeltests::check_dims(tidy_igsca, 249, 6)
  modeltests::check_dims(tidy_gsca, 187, 6)
  modeltests::check_dims(tidy_gsca_m, 249, 6)
})


test_that("glance.cSEMResults for GSCA", {

  glance_igsca <- glance(igsca_mod)
  glance_gsca <- glance(gsca_mod)
  glance_gsca_m <- glance(gsca_m_mod)

  # Claude wrote the regexp for me
  expect_error(
    check_glance_outputs(glance_igsca),
    regexp = "Expected `length\\(unacceptable\\) == 0` to be TRUE"
  )
  expect_error(
    check_glance_outputs(glance_gsca),
    regexp = "Expected `length\\(unacceptable\\) == 0` to be TRUE"
  )
  expect_error(
    check_glance_outputs(glance_gsca_m),
    regexp = "Expected `length\\(unacceptable\\) == 0` to be TRUE"
  )
})