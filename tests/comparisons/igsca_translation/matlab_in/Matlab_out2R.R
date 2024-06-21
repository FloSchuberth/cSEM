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

igsca_in <- extract_parseModel(model = tutorial_igsca_model,
                               data = LeDang2022)

igsca_out <- with(igsca_in, igsca(Z0 = Z0,
                                  W0 = W0,
                                  C0 = C0,
                                  B0 = B0,
                                  con_type = con_type,
                                  indicator_type = indicator_type,
                                  .dominant_indicators =  NULL))

if(identical(names(mat_end), c("W", "C", "B", "uniqueD"))) {
  if(identical(names(igsca_out), c("Weights", "Loadings", "Path Coefficients", "Uniqueness Terms"))) {
    names(mat_end) <- names(igsca_out)
    # The Loadings need to be transposed from Matlab version of C
    mat_end$Loadings<-t(mat_end$Loadings)
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
      weights = Weights,
      loadings = Loadings,
      uniqueD = `Uniqueness Terms`,
      paths = `Path Coefficients`
    )
  ))


save(igsca_sim_m_table, file = "tests/testthat/data/igsca_matlab.RData")


testthat::expect_equal(mat_end$Weights, igsca_out$Weights)
testthat::expect_equal(mat_end$Loadings, igsca_out$Loadings)
testthat::expect_equal(mat_end$`Path Coefficients`, igsca_out$`Path Coefficients`)
testthat::expect_equal(mat_end$`Uniqueness Terms`, igsca_out$`Uniqueness Terms`)
cat("On a matrix level and not looking at the individual elements selected for the lavaan table, R and matlab implementations are equivalent")