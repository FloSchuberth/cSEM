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
                  .model = tutorial_igsca_model,
                  .approach_weights = "GSCA",
                  .dominant_indicators = NULL,
                  .tolerance = 0.0001,
                  .conv_criterion = "sum_diff_absolute")



### Compute GSCA_M and GSCA for Reference ----------------------------------
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

gsca_mod <- csem(.data = LeDang2022,
                  gsca_model,
                  .approach_weights = "GSCA",
                  .dominant_indicators = NULL,
                  .tolerance = 0.0001,
                  .conv_criterion = "sum_diff_absolute")

test_that("GSCA estimates are nominal", {
  expect_true(all(verify(gsca_mod) == FALSE))
  expect_true(all(gsca_mod$Estimates$Reliabilities == 1))
})

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

gsca_m_mod <- csem(.data = LeDang2022,
                  gsca_m_model,
                  .approach_weights = "GSCA",
                  .dominant_indicators = NULL,
                  .tolerance = 0.0001,
                  .conv_criterion = "sum_diff_absolute")

test_that("GSCA-M estimates are nominal", {
  expect_true(all(verify(gsca_m_mod) == FALSE))
  expect_true(all(gsca_m_mod$Estimates$Reliabilities <= 1))
})

# gsca_m_mod$Estimates$Loading_estimates |> View()


# Parameter Estimates Only on Allowed Parameters -------------------------
test_that("Only specified path-coefficients are non-zero", {
  # Zero path coefficients where expected
  expect_true(all(
    c(mod$Estimate$Path_estimates)[c(mod$Information$Model$structural == 0)] ==
      0
  ))
  # Non-zero path coefficients where expected
  expect_true(all(
    c(mod$Estimate$Path_estimates)[c(mod$Information$Model$structural == 1)] !=
      0
  ))
})

test_that("Only specified loadings are non-zero", {
  # Zero loadings where expected
  expect_true(all(
    c(mod$Estimate$Loading_estimates)[c(mod$Information$Model$measurement == 0)] ==
      0
  ))
  # Non-zero loadings where expected
  expect_true(all(
    c(mod$Estimate$Loading_estimates)[c(mod$Information$Model$measurement == 1)] !=
      0
  ))
})

test_that("Only specified weights are non-zero", {
  # Note that we can use `mod$information$Model$measurement` as a substitute because for GSCA,
  # all constructs will always have both weights and loadings.

  # Zero weights where expected
  expect_true(all(
    c(mod$Estimate$Weight_estimates)[c(
      mod$Information$Model$measurement == 0
    )] ==
      0
  ))
  # Non-zero weights where expected
  expect_true(all(
    c(mod$Estimate$Weight_estimates)[c(
      mod$Information$Model$measurement == 1
    )] !=
      0
  ))
})

test_that("Only indicators of common factors have uniqueness scores and unique loadings", {
  names_cf <- names(mod$Information$Model$construct_type[
    mod$Information$Model$construct_type == "Common factor"
  ])

  names_c <- names(mod$Information$Model$construct_type[
    mod$Information$Model$construct_type == "Composite"
  ])

  indicator_cf <- apply(
    mod$Information$Model$measurement[names_cf, ],
    2,
    function(col) as.logical(col) |> any()
  )

  indicator_c <- apply(
    mod$Information$Model$measurement[names_c, ],
    2,
    function(col) as.logical(col) |> any()
  )

  absolute_sum_U <- colSums(abs(mod$Estimate$Unique_scores))

  # Zero uniqueness scores and unique loadings where expected
  expect_true(all(absolute_sum_U[indicator_c] == 0))
  expect_true(all(mod$Estimate$Unique_loading_estimates[indicator_c] == 0))

  # Non-zero uniqueness scores and unique loadings where expected
  expect_true(all(absolute_sum_U[indicator_cf] != 0))
  expect_true(all(mod$Estimate$Unique_loading_estimates[indicator_cf] != 0))
})

# Valid Reliabilities ----------------------------------------------------
test_that("Reliabilities are correctly estimated", {
  expect_true(all(mod$Estimates$Reliabilities <= 1))
})

test_that("Model estimation passes standards", {
  # Counter-intuitively, FALSE means that convergence was OK
  expect_true(all(!verify(mod))) 
})

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
    uniqueD = Unique_loading_estimates,
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
testthat::expect_failure(testthat::expect_equal(object = igsca_r_table,
                                                expected = igsca_sim_m_table,
                                                tolerance= .0341))

testthat::expect_equal(object = igsca_r_table,
                       expected = igsca_sim_m_table, tolerance = .0342)

# If the kronecker bypass for Loadings is done:
# testthat::expect_success(testthat::expect_equal(object = igsca_r_table,
#                                                 expected = igsca_sim_m_table, tolerance = .037))

testthat::expect_failure(
  testthat::expect_identical(object = igsca_r_table, expected = igsca_sim_m_table),
  info = "Matlab and R versions should be very similar, but not identical"
)

# waldo::compare(igsca_sim_m_table, igsca_r_table, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_r_table)

## GSCAPro and R ---------------------------------------------------
testthat::expect_success(testthat::expect_equal(
  igsca_r_table,
  igsca_gscapro,
  tolerance = 0.0335
))

# waldo::compare(igsca_r_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_r_table, igsca_gscapro)

## Compare GSCAPro and Matlab ----------------------------------------------
testthat::expect_success(testthat::expect_equal(
  igsca_sim_m_table,
  igsca_gscapro,
  tolerance = 0.0335
))

# waldo::compare(igsca_sim_m_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_gscapro)

# Comparing Different Ways of Fitting Group Models ------------------------
# TODO: Is fitting multiple models separately for each sub-group the same as
# fitting a 'global model'? Wasn't sure if the number of iterations + exit
# condition might make it so that a global model might exit earlier than a model
# for each sub-group?

# Accessory Functions ----------------------------------------------------
# TODO: Write these tests for getIgscaInputs

test_that("getIgscaInputs output is of the correct dimensions", {
  out <- getIgscaInputs(
    .data = processData(
      .data = LeDang2022,
      .model = tutorial_igsca_model,
      .instruments = NULL
    ),
    .model = tutorial_igsca_model
  )
  Z0 <- out$Z0
  W0 <- out$W0
  B0 <- out$B0
  C0 <- out$C0
  con_type <- out$con_type
  indicator_type <- out$indicator_type

  # Data matrix (Z0) does not correctly correspond with Weights matrix (W0)
  expect_equal(colnames(Z0), rownames(W0))
  # Weights matrix (W0) does not correctly correspond with Structural matrix (B0)
  expect_equal(colnames(W0), colnames(B0))
  # Construct indicator does not correctly correspond with Weights Matrix (W0)
  expect_equal(names(con_type), colnames(W0))
  # All indicators for composite and common factor should have loadings in I-GSCA
  expect_equal(W0, t(C0))
  # Every indicator should only have loadings from one construct
  expect_equal(unique(colSums(C0)), 1)
  # Every indicator should only have weights to one construct
  expect_equal(unique(rowSums(W0)), 1)
})
