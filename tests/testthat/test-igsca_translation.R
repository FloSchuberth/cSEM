## Libraries
require(R.matlab) # TODO: Should check whether one should load libraries in testthat
require(testthat) # TODO: This certainly doesn't show up in the cSEM test files...
require(here)
require(readxl)
require(future)
require(future.apply)

# Generate R Outputs by Swap ----------------------------------------------

## Model Specification and Load Data ---------------------------------------

tutorial_igsca_model <- "
# Composite Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
NumberofJobInterviews <~ Interview1 + Interview2
NumberofJobOffers <~ Offer1 + Offer2

# Reflective Measurement Model
Honesty_Humility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
Openness_to_Experience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10

# Structural Model
NetworkingBehavior ~ Honesty_Humility + Emotionality + Extraversion + Agreeableness + Conscientiousness + Openness_to_Experience
NumberofJobInterviews ~ NetworkingBehavior
NumberofJobOffers ~ NetworkingBehavior
"

dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx"))


igsca_sim_in <- extract_parseModel(model = tutorial_igsca_model,
                                   data = dat,
                                   ind_domi_as_first = TRUE)


## Generate Each Swap ------------------------------------------------------

# Pre-make Folder Directories
here::here(eval(formals(igsca_sim)$devdir),
           "R_out",
           eval(formals(igsca_sim)$swap_step)) |>
  lapply(FUN = \(input_path) dir.create(
    path = input_path,
    recursive = TRUE,
    showWarnings = FALSE
  )) |>
  invisible()

# Generate Swaps
# TODO: In principle, we could parallelize this. But this function and its tests are unlikely to reach the end
generate_swaps <- formals(igsca_sim)$swap_step |>
  eval() |>
  as.list() |>
  lapply(FUN = \(sw_step) try(igsca_sim(
    z0 = igsca_sim_in$z0,
    W0 = igsca_sim_in$W0,
    C0 = igsca_sim_in$C0,
    B0 = igsca_sim_in$B0,
    lv_type = igsca_sim_in$lv_type,
    ov_type = igsca_sim_in$ov_type,
    ind_domi = igsca_sim_in$ind_domi,
    nbt = 0,
    devmode = TRUE,
    swap_step = sw_step
  ))
  ) 

names(generate_swaps) <- formals(igsca_sim)$swap_step |>
  eval() |>
  as.list()

lapply(generate_swaps, is, 'try-error')

# TODO: Set some testthat expectation that there are no errors here


# Compare Matlab and R Object Names ---------------------------------------
igsca_sim(z0 = igsca_sim_in$z0,
          W0 = igsca_sim_in$W0,
          C0 = igsca_sim_in$C0,
          B0 = igsca_sim_in$B0,
          lv_type = igsca_sim_in$lv_type,
          ov_type = igsca_sim_in$ov_type,
          ind_domi = igsca_sim_in$ind_domi,
          nbt = 0,
          devmode = TRUE,
          swap_step = "noswap",
          devmode_checkobj = TRUE)

# TODO: Set testthat expectation of what the current names of the objects are...

 

# Compare GSCAPro and Matlab ----------------------------------------------

# TODO: Compare GSCAPro and matlab
 


# Compare Matlab and R ----------------------------------------------------

# TODO: Compare matlab and R


# Compare GSCAPro and R ---------------------------------------------------

# TODO: Compare GSCAPro and R
 







