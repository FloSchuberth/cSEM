# General Pre-Test -------------------------------------------------------------

## Model Specification and Load Data ---------------------------------------
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

data("LeDang2022")


igsca_in <- extract_parseModel(model = tutorial_igsca_model,
                                   data = LeDang2022,
                                   ind_domi_as_first = TRUE)

## Testing Input Matrices From extract_parseModel -------------------------------

testthat::expect_identical(igsca_in$W0, igsca_in$C0,
                           label = "All indicators for composite and factorial LVs should have loadings in I-GSCA")

## Compute and tabulate igsca ----------------------------------------------------

igsca_out <- with(igsca_in, igsca(z0 = z0,
                           W0 = W0,
                           C0 = C0,
                           B0 = B0,
                           lv_type = lv_type,
                           ov_type = ov_type,
                           ind_domi = ind_domi,
                           nbt = 0))

igsca_r_table <- with(
  igsca_out,
  get_lavaan_table_igsca_matrix(
    model = tutorial_igsca_model,
    weights = Weights,
    loadings = Loadings,
    uniqueD = `Uniqueness Terms`,
    paths = `Path Coefficients`
  )
)

# Comparison Against GSCAPro and igsca_sim.m ------------------------------



## Pre-Test ----------------------------------------------------------------

### igsca_sim.m -------------------------------------------------------------

# TODO: Write matlab results to a csv and write separate reproducible code for how to generate these matlab results.

## Use custom parser for loading GSCAPro Results
gscapro <- parse_GSCAPro_FullResults()

gscapro_tabulated <-
  get_lavaan_table_igsca_gscapro(gscapro_in = gscapro, model = tutorial_igsca_model)


## GSCAPro and Matlab ------------------------------------------------------
try(testthat::expect_equal(object = end_comparisons_table[["matlab"]],
                       expected = gscapro_tabulated))

# See https://r-pkgs.org/testing-basics.html
waldo::compare(end_comparisons_table[["matlab"]], gscapro_tabulated,
               max_diffs = Inf)

all.equal(end_comparisons_table[["matlab"]], gscapro_tabulated)

## GSCAPro and R ---------------------------------------------------
compared_R_gscapro <-
  try(testthat::expect_equal(igsca_r_table, gscapro_tabulated)
  )

waldo::compare(igsca_out, gscapro_tabulated,
               max_diffs = Inf)

all.equal(end_comparisons_table[["noswap"]], gscapro_tabulated)
