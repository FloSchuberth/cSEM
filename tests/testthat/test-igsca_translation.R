## Libraries
require(testthat) # TODO: This certainly doesn't show up in the cSEM test files...
require(R.matlab) # TODO: Should check whether one should load libraries in testthat
require(here)
require(readxl)

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
generate_swaps <- formals(igsca_sim)$swap_step |>
  eval() |>
  as.list() |>
  lapply(FUN = \(sw_step) testthat::expect_no_error(igsca_sim(
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

# Compare Matlab and R Object Names ---------------------------------------
matlab_r_compared_names <- igsca_sim(z0 = igsca_sim_in$z0,
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

# Check to see if the list of acceptable non-overlapping names between R and Matlab is the same as before
testthat::test_that("Matlab_R_Compared_Names", {
  testthat::expect_snapshot(igsca_sim(z0 = igsca_sim_in$z0,
                            W0 = igsca_sim_in$W0,
                            C0 = igsca_sim_in$C0,
                            B0 = igsca_sim_in$B0,
                            lv_type = igsca_sim_in$lv_type,
                            ov_type = igsca_sim_in$ov_type,
                            ind_domi = igsca_sim_in$ind_domi,
                            nbt = 0,
                            devmode = TRUE,
                            swap_step = "noswap",
                            devmode_checkobj = TRUE))
})
 


# Compare Matlab and R ----------------------------------------------------

# TODO: Compare matlab and R


# Compare GSCAPro and R ---------------------------------------------------

# TODO: Compare GSCAPro and R
 

# Compare GSCAPro and Matlab ----------------------------------------------


# TODO: Need to make a custom parser for the full results
GSCAPro_in <- list("dev", "Notes", "data", "GSCAPro_1_2_1Output") |>
  {\(x) here::here(x, list.files(here::here(x)))}()
GSCAPro <- vector(mode = "list", length = length(GSCAPro_in))
names(GSCAPro) <- list("dev", "Notes", "data", "GSCAPro_1_2_1Output") |>
  {\(x) list.files(here::here(x))}()

for (i in seq_len(length(GSCAPro))) {
  if (is(try(read.csv(file = GSCAPro_in[[i]]), silent = TRUE)
         , 'try-error')) {
    try({
      x <- read.csv(file = GSCAPro_in[[i]], skip = 1)
    })
  } else {
    x <- read.csv(file = GSCAPro_in[[i]], skip = 0)
  }
  GSCAPro[[i]] <- x
}

# Rename parts for facility in interaction
names(GSCAPro) <-
  gsub(
    pattern = "Tutorial_IGSCA_model1_",
    replacement = "",
    x = names(GSCAPro),
    fixed = TRUE
  )

names(GSCAPro) <-
  gsub(
    pattern = "_ver0.csv",
    replacement = "",
    x = names(GSCAPro),
    fixed = TRUE
  )



# TODO: Compare GSCAPro and matlab
# TODO: Figure out how to compare the values of GSCAPro and Matlab