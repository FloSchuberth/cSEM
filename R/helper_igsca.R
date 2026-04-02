#' Integrated Generalized Structured Component Analysis
#'
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#'
#' In the example section, the specified model is based on the tutorial I-GSCA model associated with GSCA Pro \insertCite{hwangetal2023StructuralEquationModelingAMultidisciplinaryJournal}{cSEM}.
#'
#' **Note**: Here, we assume that there is only one unique loading per indicator.
#'
#' @param .Z0 Data matrix of N cases (measurements) x J indicators with named
#'   columns, unstandardized.
#' @param .W0 Indicator matrix of weights: J indicators (rows) and their
#'   corresponding Gamma construct variables (columns).
#' @param .C0 Indicator matrix of loadings: J indicators (rows) and their
#'   corresponding Gamma construct variable (columns).
#' @param .B0 Square indicator matrix of path coefficients:
#'   from-construct-variable (rows) and to-construct-variable (columns). The
#'   order of Gamma construct variables should match the order in C0 and W0.
#' @param .con_type A vector that denotes whether each construct variable
#'   (columns in W0 and C0) is a common factor or composite. Its length should
#'   be equal to the number of columns of W0 and C0.
#' @param .indicator_type An indicator vector that indices whether a j indicator
#'   (rows of W0 and C0) corresponds to a common factor variable (1) or a
#'   composite variable (0). This vector is important for computing the
#'   uniqueness terms (D) because it zeros the entries for composite indicators.
#' @inheritParams csem
#' @inheritParams csem_arguments
#'
#' @author Michael S. Truong
#' @return List of 4 matrices that make up a fitted I-GSCA Model:
#' * (1) Weights
#' * (2) Loadings
#' * (3) Squared Unique Loadings D^2
#' * (4) Path Coefficients.
#' * (5) Unique Component of Indicators (for common factors) DU
#'
#' @importFrom MASS ginv
#'
#' @references
#'   \insertAllCited{}
#' @examples
#' \dontrun{
#' # Specify the model according to GSCA Pro's example
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 +
#'                       Behavior7 + Behavior8 +  Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2
#'
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 +
#'                     Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 +
#'                 Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 +
#'                 Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 +
#'                  Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 +
#'                      Conscientious6 + Conscientious7 + Conscientious8 +
#'                      Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 +
#'                         Openness7 + Openness8 + Openness9 + Openness10
#'
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion +
#'                      Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#'
#' data(LeDang2022)
#'
#' csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "GSCA",
#' .dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion =
#' "sum_diff_absolute")
#' }
igsca <-
  function(
    .X = args_default()$.X,
    .S = args_default()$.S,
    .csem_model = args_default()$.csem_model,
    .conv_criterion = args_default()$.conv_criterion,
    .GSCA_modes = args_default()$.GSCA_modes,
    .iter_max = args_default()$.iter_max,
    .starting_values = args_default()$.starting_values,
    .tolerance = args_default()$.tolerance,
    .Z0,
    .W0,
    .C0,
    .B0,
    .con_type,
    .indicator_type
  ) {
    
  }

