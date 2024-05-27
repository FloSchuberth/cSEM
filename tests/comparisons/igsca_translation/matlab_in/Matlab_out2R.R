#' Converts Output of igsca functions into a table to facilitate comparisons
#' 
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#' @param weights Weights matrix
#' @param loadings Loadings matrix
#' @param uniqueD Vector of Uniqueness for each indicator of a common factor
#' @param paths Path coefficients matrix
#' @importFrom lavaan lavaanify
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

mat_end_extract_order <- c("W", "C", "B", "uniqueD")

mat_end <-
  R.matlab::readMat(
    "~/Documents/RStudio/cSEM/tests/comparisons/igsca_translation/matlab_out/end_results.MAT" # Hard-Coded value
  ) |>
  {
    \(mat_env) mat_env[mat_end_extract_order]
  }() 

mat_end$uniqueD <- as.vector(mat_end$uniqueD) # Matlab vectors are matrices when loaded into R, so we convert to a vector here

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

mod <- csem(.data = LeDang2022,
            tutorial_igsca_model,
            .approach_weights = "IGSCA",
            .dominant_indicators = NULL,
            .tolerance = 0.0001,
            .conv_criterion = "sum_diff_absolute")

igsca_out <- mod[["Estimates"]][c("Weight_estimates", "Loading_estimates", "Path_estimates", "D2")]

if (identical(names(igsca_out), c("Weight_estimates", "Loading_estimates", "Path_estimates", "D2"))) {
  names(igsca_out) <- c("Weights", "Loadings", "Path Coefficients", "Uniqueness Terms")
} else {
  stop("Renaming won't work -- names of matrices have changed")
}

if(identical(names(mat_end), c("W", "C", "B", "uniqueD"))) {
  if(identical(names(igsca_out), c("Weights", "Loadings", "Path Coefficients", "Uniqueness Terms"))) {
    names(mat_end) <- names(igsca_out)
    # The Weights and Path Coefficients Matrices  need to be transposed from Matlab version
    mat_end$Weights<-t(mat_end$Weights)
    mat_end$`Path Coefficients`<-t(mat_end$`Path Coefficients`)
  }
} else {
  stop("Can't rename the Matlab end results matrices based on cSEM::igsca() because they don't matchup")
}


for (i in names(mat_end)) {
  if (isTRUE(sapply(mat_end, is.matrix)[i])) {
    rownames(mat_end[[i]]) <-
      rownames(igsca_out[[i]])
    colnames(mat_end[[i]]) <-
      colnames(igsca_out[[i]])
  } else if (isTRUE(sapply(mat_end, is.vector)[i])) {
    names(mat_end[[i]]) <- names(igsca_out[[i]])
  } else {
    stop("One of the compared objects is neither a vector nor a matrix")
  }
}

(igsca_sim_m_table <-
  with(
    mat_end,
    get_lavaan_table_igsca_matrix(
      model = tutorial_igsca_model,
      weights = t(Weights),
      loadings = t(Loadings),
      uniqueD = `Uniqueness Terms`,
      paths = t(`Path Coefficients`)
    )
  ))


save(igsca_sim_m_table, file = "tests/testthat/data/igsca_matlab.RData")


testthat::expect_equal(mat_end$Weights, igsca_out$Weights)
testthat::expect_equal(mat_end$Loadings, igsca_out$Loadings)
testthat::expect_equal(mat_end$`Path Coefficients`, igsca_out$`Path Coefficients`)
testthat::expect_equal(mat_end$`Uniqueness Terms`, igsca_out$`Uniqueness Terms`)
cat("On a matrix level and not looking at the individual elements selected for the lavaan table, R and matlab implementations are equivalent")