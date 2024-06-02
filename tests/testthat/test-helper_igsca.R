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

## Compute and tabulate igsca ----------------------------------------------------
mod <- csem(.data = LeDang2022,
                  tutorial_igsca_model,
                  .approach_weights = "IGSCA",
                  .dominant_indicators = NULL,
                  .tolerance = 0.0001,
                  .conv_criterion = "sum_diff_absolute")


### Custom Function for Organizing IGSCA Results ----------------------------

#' Converts Output of igsca functions into a table to facilitate comparisons
#' 
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#' 
#' I chose not to use the tidy method for the tests because it would be a lot of
#' work to make the GSCAPro results fit in the tidy format.
#' 
#' @param weights Weights matrix
#' @param loadings Loadings matrix
#' @param uniqueD Vector of Uniqueness for each indicator of a common factor
#' @param paths Path coefficients matrix
#' @importFrom lavaan lavaanify
#' @author Michael S. Truong
#' @return Table of Weights, Loadings, Path-Coefficients and Uniqueness terms from i-gsca algorithms in Matlab or R.
#'
get_lavaan_table_igsca_matrix <- function(model, weights, loadings, uniqueD, paths) {
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"),]
  # Pre-allocate Columns
  table <-
    cbind(table, list(
      "weights" = 0,
      "loadings" = 0,
      "uniqueD" = 0,
      "paths" = 0
    ))
  
  # Slide in weights
  for (indicator in rownames(weights)) {
    for (lv in colnames(weights)) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "weights"] <- weights[indicator, lv]
    }
  }
  
  # Slide in loadings
  for (indicator in rownames(loadings)) {
    for (lv in colnames(loadings)) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "loadings"] <- loadings[indicator, lv]
    }
  }
  
  # Slide in uniqueD
  for (indicator in names(uniqueD)) {
    # This assumes that every indicator only loads onto one factor
    # Cross-factor loadings will not work with this
    table[((table$rhs == indicator) &
             (table$op == "=~")), "uniqueD"] <- uniqueD[indicator]
  }
  
  # Slide in Paths
  for (lv_from in rownames(paths)) {
    for (lv_to in colnames(paths)) {
      table[((table$rhs == lv_from &
                table$lhs == lv_to) &
               table$op == "~"), "paths"] <- paths[lv_from, lv_to]
    }
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  return(table)
  
}


### Fetching igsca_r_table Results ------------------------------------------
igsca_r_table <- with(
  mod$Estimates,
  get_lavaan_table_igsca_matrix(
    model = tutorial_igsca_model,
    weights = t(Weight_estimates),
    loadings = t(Loading_estimates),
    uniqueD = D2,
    paths = t(Path_estimates)
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
testthat::expect_equal(object = igsca_r_table,
                           expected = igsca_sim_m_table)

testthat::expect_failure(
  testthat::expect_identical(object = igsca_r_table, expected = igsca_sim_m_table),
  info = "Matlab and R versions should be very similar, but not identical"
)

# waldo::compare(igsca_sim_m_table, igsca_r_table, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_r_table)

## GSCAPro and R ---------------------------------------------------
testthat::expect_failure(testthat::expect_equal(igsca_r_table,
                                                igsca_gscapro))

testthat::expect_success(testthat::expect_equal(igsca_r_table,
                                                igsca_gscapro,
                                                tolerance = .1))

# waldo::compare(igsca_r_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_r_table, igsca_gscapro)


## Compare GSCAPro and Matlab ----------------------------------------------
testthat::expect_failure(testthat::expect_equal(igsca_sim_m_table,
                                                igsca_gscapro))

testthat::expect_success(testthat::expect_equal(igsca_sim_m_table,
                                                igsca_gscapro,
                                                tolerance = .1))

# waldo::compare(igsca_sim_m_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_gscapro)


# Comparing Different Ways of Fitting Group Models ------------------------
# TODO: Is fitting multiple models separately for each sub-group the same as
# fitting a 'global model'? Wasn't sure if the number of iterations + exit
# condition might make it so that a global model might exit earlier than a model
# for each sub-group?

