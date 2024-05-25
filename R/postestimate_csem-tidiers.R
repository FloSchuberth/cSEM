#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @importFrom generics tidy
#' @export
#'
#' @examples
#' # Specify the model according to GSCA Pro's example
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
#' data(LeDang2022)
#' 
#' mod <- csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "IGSCA",
#' .dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion = "sum_diff_absolute")
#' tidy(mod)
tidy.cSEMResults <- function(x,
                             format = c("list", "tibble"),
                             what = c(
                               "all",
                               "Path_estimates",
                               "Loadings_estimates",
                               "Weight_estimates",
                               "D2",
                               "Effect_estimates",
                               "Indicator_correlation",
                               "Exo_construct_correlation"
                             ),
                             ...) {
  
  summarized_cSEMResults <- summarize(x)
  
  if (what == "all") { 
    general_estimates <- summarized_cSEMResults[["Estimates"]][c("Path_estimates", "Loading_estimates", "Weight_estimates", "D2", "Indicator_correlation", "Exo_construct_correlation")]
    effect_estimates <- summarized_cSEMResults[["Estimates"]][["Effect_estimates"]][c("Direct_effect", "Indirect_effect", "Total_effect")]
    out <- c(general_estimates, effect_estimates) # TODO: Test that this works -- I want it flat
    
  } else if (what != "Effect_estimates") {
    general_estimates <- summarized_cSEMResults[["Estimates"]][what] # TODO: Test that this works
    out <- general_estimates
    
  } else if (what == "Effect_estimates") {
    effect_estimates <- summarized_cSEMResults[["Estimates"]][c("Direct_effect", "Indirect_effect", "Total_effect")]
    out <- effect_estimates
    
  } else {
    stop("An unsupported result of summarize was selected.")
    
  }
  
  
  # TODO: Implement whether to return as single list the way summarize does or a tibble like tidy.lavaan
  if (format == "list") {
    return(out)
    
  } else if (format == "tibble") {
    # TODO: Gotta make some remark if Path_estimates is not equal to Direct effect estimates
    # TODO: Flatten and return as a tibble in the fashion of lavaan-tidiers
  } else {
    stop("An unsupported return format of tidy.cSEMResults was selected.")
    
  }
    
  
  
  
}

#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @importFrom generics glance
#' @export
#'
#' @examples
glance.cSEMResults <- function(x, ...) {
  # TODO: May want to look at assess
  # cSEM:::assess()
}

#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @importFrom generics augment
#' @export
#'
#' @examples
augment.cSEMResults <- function(x, ...) {
  # TODO: Look at predict and maybe also get the construct scores
  # cSEM:::predict()
}