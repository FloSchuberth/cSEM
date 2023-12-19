# Libraries ---------------------------------------------------------------
require(testthat) # TODO: This certainly doesn't show up in the cSEM test files...
require(R.matlab) # TODO: Should check whether one should load libraries in testthat
require(here)
require(readxl)
require(future)
require(future.apply)
require(withr)

# future Specifications ---------------------------------------------------
## Future Plans copied from csem_resample.R
oplan <- future::plan()
on.exit(future::plan(oplan), add = TRUE)
future::plan("multisession", workers = 4)

# Model Specification and Load Data ---------------------------------------
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

dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx"))


igsca_sim_in <- extract_parseModel(model = tutorial_igsca_model,
                                   data = dat,
                                   ind_domi_as_first = TRUE)



# Testing Input Matrices From extract_parseModel -------------------------------

expect_identical(igsca_sim_in$W0, igsca_sim_in$C0,
                 label = "All indicators for composite and factorial LVs should have loadings in I-GSCA")

# Generate R Outputs by Swap ----------------------------------------------


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
# TODO: For some reason future.apply::future_lapply() doesn't work here because it can't find igsca_sim?
generated_swaps_end_results <- formals(igsca_sim)$swap_step |>
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

names(generated_swaps_end_results) <- eval(formals(igsca_sim)$swap_step)

# R-Swaps Should be different From Each Other -----------------------------

# Load R end_results from each swap step

end_results <- here::here(eval(formals(igsca_sim)$devdir),
                          "R_out",
                          eval(formals(igsca_sim)$swap_step),
                          "end_results.rds") |>
  lapply(FUN = readRDS)

names(end_results) <- eval(formals(igsca_sim)$swap_step)
end_comparisons <-
  lapply(end_results,
         FUN = \(computations) computations[c("Weights", "Loadings",
                                              "PathCoefficients", "uniqueD")]) 

end_comparisons_table <- lapply(
  end_comparisons,
  FUN = function(swap_results)
    with(
      swap_results,
      get_lavaan_table_igsca_matrix(
        model = tutorial_igsca_model,
        weights = Weights,
        loadings = Loadings,
        uniqueD = uniqueD,
        paths = PathCoefficients
      )
    )
)
## Compare the R end_results across all different unique combination --------

unique_combn_swaps <-
  combn(names(end_comparisons_table[!names(end_comparisons_table) %in% "matlab"]),
        m = 2,
        simplify = TRUE)

compared_R_swaps <- mapply(
  FUN = function(obj, expect)
    identical(end_comparisons_table[[obj]],
              end_comparisons_table[[expect]]),
  obj = unique_combn_swaps[1, ],
  expect = unique_combn_swaps[2,],
  SIMPLIFY = TRUE
)
## Not expecting them to be identical due to differences between R and matlab, so expect_false

testthat::expect_false(all(compared_R_swaps))

## R Swaps Comparison ------------------------------------------

if (identical(names(compared_R_swaps), unique_combn_swaps[1, ])) {
  compared_R_swaps_named <-
    rbind(compared_R_swaps, unique_combn_swaps[2, ]) |>
    t()
} else {
  stop(
    "It is unsafe to rbind the compared_R_swaps and unique_combn_swaps[2,] because the correct correspondence is not there"
  )
}

# Interpret Table of identicalness
compared_R_swaps_named


# The only two times any two swaps are identical to each other are:
# noswap and first_iteration_update_C_B_D_uniqueD_est
# TODO: Difficult to interpret...
# prepare_for_ALS and first_factor_update
# - This makes sense to me because the first iteration in this model is on a factor (Honesty_Humility),
#   so I don't expect much to change in-between loading for prepare_for_ALS and the first_factor_update


# Compare Matlab and R Object Names ---------------------------------------
# Check to see if the list of acceptable non-overlapping names between R and Matlab is the same as before
testthat::test_that("Matlab_R_Compared_Names", {
  testthat::expect_snapshot(with(
    igsca_sim_in,
    igsca_sim(
      z0 = z0,
      W0 = W0,
      C0 = C0,
      B0 = B0,
      lv_type = lv_type,
      ov_type = ov_type,
      ind_domi = ind_domi,
      nbt = 0,
      devmode = TRUE,
      swap_step = "noswap",
      devmode_checkobj = TRUE
    )
  ))
})
 


# Compare Matlab and R ----------------------------------------------------

## Comparing End Results ---------------------------------------------------

# Load Matlab end_results
matlab_end_results <-
  here::here(list("dev", "Notes", "data", "matlab_out", "end_results.MAT")) |>
  R.matlab::readMat(fixNames = FALSE) |>
  {
    \(mat_env) mat_env[c("W", "C", "B", "uniqueD")]
  }() |>
  convert_matlab2R(vectorness = unlist(lapply(end_comparisons[[1]], is.vector)))

if (with(matlab_end_results, exists("C"))) {
  # C needs to be transposed
  matlab_end_results$C <- t(matlab_end_results$C)
}


### Add the names back to the matlab objects and expect matching dimensions --------

## Generate Dimension Sizes

end_comparisons_dims <-
  lapply(end_comparisons$noswap, \(result) try(dim(result))
  )

end_comparisons_dims[unlist(lapply(end_comparisons_dims, is.null))] <- 
  lapply(end_comparisons$noswap[sapply(end_comparisons_dims, is.null)], length)

matlab_end_results_dims <-
  lapply(matlab_end_results, \(result) try(dim(result))
  )

matlab_end_results_dims[unlist(lapply(matlab_end_results_dims, is.null))] <- 
  lapply(matlab_end_results[sapply(matlab_end_results_dims, is.null)], length)


## Expect the dimensions of the extracted objects to match
mapply(FUN = testthat::expect_identical, object = end_comparisons_dims, expected = matlab_end_results_dims)

## Add back the names
### Name of each list
names(matlab_end_results) <- names(end_comparisons$noswap)


### Name of the columns, rows and vectors
for (i in names(matlab_end_results)) {
  if (isTRUE(sapply(matlab_end_results, is.matrix)[i])) {
    rownames(matlab_end_results[[i]]) <-
      rownames(end_comparisons$noswap[[i]])
    colnames(matlab_end_results[[i]]) <-
      colnames(end_comparisons$noswap[[i]])
  } else if (isTRUE(sapply(matlab_end_results, is.vector)[i])) {
    names(matlab_end_results[[i]]) <- names(end_comparisons$noswap[[i]])
  } else {
    stop("One of the compared objects is neither a vector nor a matrix")
  }
}

# Tabulate both R and Matlab Output
end_comparisons_table <- lapply(
  c(end_comparisons, list("matlab" = matlab_end_results)),
  FUN = function(swap_results)
    with(
      swap_results,
      get_lavaan_table_igsca_matrix(
        model = tutorial_igsca_model,
        weights = Weights,
        loadings = Loadings,
        uniqueD = uniqueD,
        paths = PathCoefficients
      )
    )
)

## Compares each R with the Matlab end result ------------------------------
compared_R_matlab <-
  lapply(end_comparisons_table[!(names(end_comparisons_table) == "matlab")],
         FUN = \(selected_R_table) try(
           testthat::expect_equal(object = selected_R_table,
                                  expected = end_comparisons_table[["matlab"]])
           )
  )

# Check which ones failed to be equivalent
# Only flip_signs_ind_domi should match matlab because no further computations are performed


testthat::expect_equal(end_comparisons_table$noswap, end_comparisons_table$matlab)

# compared_R_matlab[!names(compared_R_matlab) %in% "flip_signs_ind_domi"] |>
#   sapply(is, 'try-error') |>
#   all() |>
#   testthat::expect_true()
# 
# compared_R_matlab[names(compared_R_matlab) %in% "flip_signs_ind_domi"] |>
#   sapply(is, 'try-error') |>
#   all() |>
#   testthat::expect_false()


 

# Compare GSCAPro and Matlab ----------------------------------------------

## Use custom parser for loading GSCAPro Results
gscapro <- parse_GSCAPro_FullResults()

gscapro_tabulated <-
  get_lavaan_table_igsca_gscapro(gscapro_in = gscapro, model = tutorial_igsca_model)

testthat::expect_equal(object = end_comparisons_table[["matlab"]],
                       expected = gscapro_tabulated)

# See https://r-pkgs.org/testing-basics.html
withr::with_options(list(width = 20),
                    waldo::compare(end_comparisons_table[["matlab"]], gscapro_tabulated))

# Compare GSCAPro and R ---------------------------------------------------

# TODO: Compare GSCAPro and R

compared_R_gscapro <-
  lapply(end_comparisons_table[!(names(end_comparisons_table) == "matlab")],
         FUN = \(selected_R_table) try(
           testthat::expect_equal(object = selected_R_table,
                                  expected = gscapro_tabulated)
         )
  )

testthat::expect_equal(
  end_comparisons_table[["noswap"]],
  end_comparisons_table[["matlab"]]
  )

try(testthat::expect_equal(end_comparisons_table[["noswap"]], gscapro_tabulated))
all.equal(end_comparisons_table[["noswap"]], gscapro_tabulated)


withr::with_options(list(width = 20),
                    waldo::compare(end_comparisons_table[["noswap"]], gscapro_tabulated))


## TODO: Interpret
## 
## TODO: See what happens when I increase the maximum number of iterations and the tolerance