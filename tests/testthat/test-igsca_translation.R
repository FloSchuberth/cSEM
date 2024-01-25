if (exists(run)) { # To disable for now
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

dat <- readxl::read_excel(here::here("dev", "Notes", "data", "mmc1.xlsx"))


igsca_sim_in <- extract_parseModel(model = tutorial_igsca_model,
                                   data = dat,
                                   ind_domi_as_first = TRUE)



## Testing Input Matrices From extract_parseModel -------------------------------

testthat::expect_identical(igsca_sim_in$W0, igsca_sim_in$C0,
                           label = "All indicators for composite and factorial LVs should have loadings in I-GSCA")

## Generate R Outputs by Swap ----------------------------------------------


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
# FIXME: For some reason future.apply::future_lapply() doesn't work here because it can't find igsca_sim?

generate_swaps <- FALSE
if (isTRUE(generate_swaps)) {
  generated_swaps_end_results <- formals(igsca_sim)$swap_step |>
    eval() |>
    as.list() |>
    lapply(FUN = \(sw_step) testthat::expect_no_error(
      igsca_sim(
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
      )
    ))
  
  names(generated_swaps_end_results) <-
    eval(formals(igsca_sim)$swap_step)
} else if (isFALSE(generate_swaps)) {
  warning("Swaps were not generated, comparisons with swaps may be out-of-date.")
}

# Comparing R-Swaps-----------------------------
## Comparing R-Swaps against each other

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
  FUN = function(end_comparisons_iter)
    with(
      end_comparisons_iter,
      get_lavaan_table_igsca_matrix(
        model = tutorial_igsca_model,
        weights = Weights,
        loadings = Loadings,
        uniqueD = uniqueD,
        paths = PathCoefficients
      )
    )
)

## Compare All Unique Combinations of R end_results Based on Swap --------

# Generate Combinations
unique_combn_swaps <-
  combn(names(end_comparisons_table[!names(end_comparisons_table) %in% "matlab"]),
        m = 2,
        simplify = TRUE)


### Rough Equivalence -------------------------------------------------------
# Compare for rough equivalence
compared_R_swaps <- mapply(
  FUN = function(obj, expect)
    try(testthat::expect_equal(end_comparisons_table[[obj]],
              end_comparisons_table[[expect]])),
  obj = unique_combn_swaps[1, ],
  expect = unique_combn_swaps[2,],
  SIMPLIFY = FALSE
)

## Expect that all of the R-swaps to be roughly equivalent 
compared_R_swaps |>
  sapply(is, 'try-error') |>
  any() |>
  testthat::expect_false()

### Exact Identicalness -----------------------------------------------------
# Compare for exact equivalence
compared_R_swaps_identical <- mapply(
  FUN = function(obj, expect)
    try(testthat::expect_identical(end_comparisons_table[[obj]],
                 end_comparisons_table[[expect]])),
  obj = unique_combn_swaps[1, ],
  expect = unique_combn_swaps[2,],
  SIMPLIFY = TRUE
)
## Not expecting them to be identical due to differences between R and matlab, so expect_false
compared_R_swaps_identical |>
  sapply(is, 'try-error') |>
  all() |>
  testthat::expect_false()

#### Interpretation of Identicalness ----------------------------
if (identical(names(compared_R_swaps_identical), unique_combn_swaps[1, ])) {
  compared_R_swaps_identical_named <-
    rbind(compared_R_swaps_identical, unique_combn_swaps[2, ]) |>
    t()
} else {
  stop(
    "It is unsafe to rbind the compared_R_swaps and unique_combn_swaps[2,] because the correct correspondence is not there"
  )
}

compared_R_swaps_identical_named

testthat::expect_identical(end_comparisons_table[["noswap"]], end_comparisons_table[["first_iteration_update_C_B_D_uniqueD_est"]])

# None of the swaps are precisely identical except for when Matlab is substituted in after the first iteration... This makes the result identical to noswap
## TODO: This seems to require a cup of coffee to interpret


# Compare Matlab and R ----------------------------------------------------

## Pre-Test ----------------------------------------------------------------

### Compare Matlab and R Object Names ---------------------------------------
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


### Load Matlab end_results -------------------------------------------------

# Load Matlab end_results
matlab_end_results <-
  here::here(list(
    "tests",
    "comparisons",
    "igsca_translation",
    "matlab_out",
    "end_results.MAT"
  )) |> 
  R.matlab::readMat(fixNames = FALSE) |>
  {
    \(mat_env) mat_env[c("W", "C", "B", "uniqueD")]
  }() |>
  convert_matlab2R(vectorness = unlist(lapply(end_comparisons[[1]], is.vector)))

if (with(matlab_end_results, exists("C"))) {
  # C-matrix needs to be transposed
  matlab_end_results$C <- t(matlab_end_results$C)
}


#### Add Names to Matlab Objects ---------------------------------------------
# Doing this requires that the dimensions between the R and matlab objects 
# are matching. So first, we compare the dimension sizes.

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
mapply(
  FUN = testthat::expect_identical,
  object = end_comparisons_dims,
  expected = matlab_end_results_dims
)

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


### Tabulate both R and Matlab Output ---------------------------------------
# This facilitates the comparison, so that we are comparing a two tables
#  against each other instead of 4 matrices on each side.

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

## Compare each R_swap with the Matlab end result ------------------------------
### Rough Equivalence -------------------------------------------------------
compared_R_matlab <-
  lapply(end_comparisons_table[!(names(end_comparisons_table) == "matlab")],
         FUN = \(selected_R_table)
         try(testthat::expect_equal(object = selected_R_table,
                                    expected = end_comparisons_table[["matlab"]]))
  )

# Tests that all the different R swaps are equivalent to the matlab version within error
compared_R_matlab |>
  sapply(is, 'try-error') |>
  any() |>
  testthat::expect_false()

### Exact Identicalness -----------------------------------------------------
compared_R_matlab_exact <-
  lapply(
    end_comparisons_table[!(names(end_comparisons_table) == "matlab")],
    FUN = \(selected_R_table)
    try(testthat::expect_identical(object = selected_R_table,
                           expected = end_comparisons_table[["matlab"]]))    
  )

# Expect TRUE because these should all swaps except for `flip_signs_ind_domi`
# should not be identical to matlab
compared_R_matlab_exact[
  !names(compared_R_matlab_exact) %in% "flip_signs_ind_domi"
  ] |>
  sapply(is, 'try-error') |>
  all() |>
  testthat::expect_true()

# Expect FALSE because flip_signs_ind_domi should be identical to matlab, so no issue of try-error
compared_R_matlab_exact[
  names(compared_R_matlab_exact) %in% "flip_signs_ind_domi"
  ] |>
  sapply(is, 'try-error') |>
  testthat::expect_false()

# Comparison Against GSCAPro ----------------------------------------------

## Pre-Test ----------------------------------------------------------------
## Use custom parser for lo ading GSCAPro Results
gscapro <- parse_GSCAPro_FullResults()

gscapro_tabulated <-
  get_lavaan_table_igsca_gscapro(gscapro_in = gscapro, model = tutorial_igsca_model)


## GSCAPro and Matlab ------------------------------------------------------
try(testthat::expect_equal(object = end_comparisons_table[["matlab"]],
                       expected = gscapro_tabulated))

# See https://r-pkgs.org/testing-basics.html
waldo::compare(end_comparisons_table[["matlab"]], gscapro_tabulated,
               max_diffs = Inf)

all.equal(end_comparisons_table[["matlab"]], gscapro_tabulated)

## GSCAPro and R ---------------------------------------------------
compared_R_gscapro <-
  try(testthat::expect_equal(end_comparisons_table$noswap, gscapro_tabulated)
  )

waldo::compare(end_comparisons_table[["noswap"]], gscapro_tabulated,
               max_diffs = Inf)

all.equal(end_comparisons_table[["noswap"]], gscapro_tabulated)
}