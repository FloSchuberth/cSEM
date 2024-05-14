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
                               data = LeDang2022)

## Testing Input Matrices From extract_parseModel -------------------------------

testthat::expect_identical(igsca_in$W0, igsca_in$C0,
                           label = "All indicators for composite and factorial LVs should have loadings in I-GSCA")

## Compute and tabulate igsca ----------------------------------------------------

igsca_out <- with(igsca_in, igsca(Z0 = Z0,
                           W0 = W0,
                           C0 = C0,
                           B0 = B0,
                           con_type = con_type,
                           indicator_type = indicator_type,
                           .dominant_indicators =  NULL))

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

# Comparisons between cSEM::igsca(), GSCAPro and igsca_sim.-----------------

## Load Matlab Results -----------------------------------------------------
# Loads into igsca_sim_m_table
load(testthat::test_path("data", "igsca_matlab.RData"))



## Load GSCAPro Results ----------------------------------------------------
# Loads into igsca_gscapro
load(testthat::test_path("data", "igsca_gscapro.RData"))



## Compare Matlab and cSEM::igsca()------------------------------------------
# See https://r-pkgs.org/testing-basics.html
testthat::expect_equal(object = igsca_sim_m_table,
                           expected = igsca_r_table)

# waldo::compare(igsca_sim_m_table, igsca_r_table, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_r_table)

## GSCAPro and R ---------------------------------------------------
testthat::expect_equal(igsca_r_table,
                      igsca_gscapro ,
                      tolerance = 10 ^ -2
                      )


# waldo::compare(igsca_r_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_r_table, igsca_gscapro)


## Compare GSCAPro and Matlab ----------------------------------------------
testthat::expect_equal(object = igsca_sim_m_table,
                       expected = igsca_gscapro ,
                       tolerance = 10 ^ -2
                       )

# waldo::compare(igsca_sim_m_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_gscapro)







