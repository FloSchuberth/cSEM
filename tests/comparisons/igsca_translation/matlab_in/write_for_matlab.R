#' Writes the extracted matrices to run with igsca_sim_test.m
#'
#' @param extracted_matrices Object returned by extract_parseModel
#'
#' @return Writes the matrices from extract_parseModel() to the appropriate test directory for the Matlab implementation of igsca_sim to run.
#' @importFrom here here
#' @examples
#' 
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' data("LeDang2022")
#' 
#' write_for_matlab(extract_parseModel(model = tutorial_igsca_model, data = LeDang2022, ind_domi_as_first = TRUE))
write_for_matlab <- function(extracted_matrices) {
  indir <- list("tests", "comparisons", "igsca_translation", "matlab_in")
  extracted_matrices$lv_type <-
    as.numeric(extracted_matrices$lv_type)
  extracted_matrices$ov_type <-
    as.numeric(extracted_matrices$ov_type)
  
  mapply(
    write.csv,
    x = extracted_matrices,
    file = paste0(here::here(indir), "/", names(extracted_matrices), ".csv"),
    row.names = FALSE
  )
  
  invisible()
}