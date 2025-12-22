# Example ----------------------------------------------------------------
# Specify the model according to GSCA Pro's example
tutorial_igsca_model <- "
# Composite (Formative Measurement) Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + 
                      Behavior7 + Behavior8 +  Behavior9
Numberofjobinterviews <~ Interview1 + Interview2
Numberofjoboffers <~ Offer1 + Offer2 

# Latent Variable (Reflective Measurement) Model
HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + 
                    Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + 
                Emotion5 + Emotion6 + Emotion8 + Emotion10
Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + 
                Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 +
                 Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + 
                     Conscientious6 + Conscientious7 + Conscientious8 + 
                     Conscientious9 + Conscientious10
OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + 
                        Openness7 + Openness8 + Openness9 + Openness10

# Structural Model
NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + 
                     Agreeableness + Conscientiousness + OpennesstoExperience
Numberofjobinterviews ~ NetworkingBehavior
Numberofjoboffers ~ NetworkingBehavior
"

data(LeDang2022)


# What I iterate through
#' devtools::load_all(".")
#' devtools::test_active_file()

# debugonce(csem)
# debugonce(foreman)
# debugonce(calculateWeightsIGSCA)
# debugonce(igsca)
# debugonce(initializeAlsEstimates)
# debugonce(initializeIgscaEstimates)
# debugonce(updateCBDU)
# debugonce(calculateReliabilities)
# debugonce(scaleWeights)
debugonce(flipGammaCBSigns)
(debug_mod <- csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "GSCA",
.dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion =
"sum_diff_absolute")
)
