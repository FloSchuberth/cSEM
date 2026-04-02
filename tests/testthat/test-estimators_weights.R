
# GSCA, GSCAM and IGSCA ------------------------------------------------------------------

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

### Compute and tabulate igsca ----------------------------------------------------
mod <- csem(
  .data = LeDang2022,
  .model = tutorial_igsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)


#### Compute GSCA_M and GSCA for Reference ----------------------------------
gsca_model <- "
# Measurement Model
NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
Numberofjobinterviews <~ Interview1 + Interview2
Numberofjoboffers <~ Offer1 + Offer2
HonestyHumility <~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
Emotionality <~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
Extraversion <~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
Agreeableness <~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
Conscientiousness <~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
OpennesstoExperience <~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10

# Structural Model
NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
Numberofjobinterviews ~ NetworkingBehavior
Numberofjoboffers ~ NetworkingBehavior
"

gsca_mod <- csem(
  .data = LeDang2022,
  gsca_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute",
  .GSCA_modes = "NCMP"
)

test_that("GSCA estimates are nominal", {
  expect_true(all(verify(gsca_mod) == FALSE))
  expect_true(all(gsca_mod$Estimates$Reliabilities == 1))
})

# gsca_mod$Estimates$Loading_estimates |> View()

gsca_m_model <- "
# Measurement Model
NetworkingBehavior =~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
Numberofjobinterviews =~ Interview1 + Interview2
Numberofjoboffers =~ Offer1 + Offer2
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

gsca_m_mod <- csem(
  .data = LeDang2022,
  gsca_m_model,
  .approach_weights = "GSCA",
  .dominant_indicators = NULL,
  .tolerance = 0.0001,
  .conv_criterion = "sum_diff_absolute"
)

test_that("GSCA-M estimates are nominal", {
  expect_true(all(verify(gsca_m_mod) == FALSE))
  expect_true(all(gsca_m_mod$Estimates$Reliabilities <= 1))
})

##### .conv_criterion for GSCAM ----------------------------------------------
test_that(".conv_criterion affects GSCAM", {
  HolzingerSwineford1939 <- lavaan::HolzingerSwineford1939

  HS.model_csem <- ' visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9 
                visual ~~ textual
                speed ~~ textual
                visual ~~ speed'
  
  sum_diff_absolute_mod <- cSEM::csem(
    .data = HolzingerSwineford1939,
    .model = HS.model_csem,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = 'sum_diff_absolute'
  )
  
  mean_diff_absolute_mod <- cSEM::csem(
    .data = HolzingerSwineford1939,
    .model = HS.model_csem,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = 'mean_diff_absolute'
  )
  expect_failure(expect_equal(sum_diff_absolute_mod, mean_diff_absolute_mod))
})

# gsca_m_mod$Estimates$Loading_estimates |> View()

## Parameter Estimates Only on Allowed Parameters -------------------------
test_that("Only specified path-coefficients are non-zero", {
  # Zero path coefficients where expected
  expect_true(all(
    c(mod$Estimate$Path_estimates)[c(mod$Information$Model$structural == 0)] ==
      0
  ))
  # Non-zero path coefficients where expected
  expect_true(all(
    c(mod$Estimate$Path_estimates)[c(mod$Information$Model$structural == 1)] !=
      0
  ))
})

test_that("Only specified loadings are non-zero", {
  # Zero loadings where expected
  expect_true(all(
    c(mod$Estimate$Loading_estimates)[c(
      mod$Information$Model$measurement == 0
    )] ==
      0
  ))
  # Non-zero loadings where expected
  expect_true(all(
    c(mod$Estimate$Loading_estimates)[c(
      mod$Information$Model$measurement == 1
    )] !=
      0
  ))
})

test_that("Only specified weights are non-zero", {
  # Note that we can use `mod$information$Model$measurement` as a substitute because for GSCA,
  # all constructs will always have both weights and loadings.

  # Zero weights where expected
  expect_true(all(
    c(mod$Estimate$Weight_estimates)[c(
      mod$Information$Model$measurement == 0
    )] ==
      0
  ))
  # Non-zero weights where expected
  expect_true(all(
    c(mod$Estimate$Weight_estimates)[c(
      mod$Information$Model$measurement == 1
    )] !=
      0
  ))
})

test_that("Only indicators of common factors have uniqueness scores and unique loadings", {
  names_cf <- names(mod$Information$Model$construct_type[
    mod$Information$Model$construct_type == "Common factor"
  ])

  names_c <- names(mod$Information$Model$construct_type[
    mod$Information$Model$construct_type == "Composite"
  ])

  indicator_cf <- apply(
    mod$Information$Model$measurement[names_cf, , drop = FALSE],
    2,
    function(col) as.logical(col) |> any()
  )

  indicator_c <- apply(
    mod$Information$Model$measurement[names_c, , drop = FALSE],
    2,
    function(col) as.logical(col) |> any()
  )

  absolute_sum_U <- colSums(abs(mod$Estimate$Unique_scores))

  # Zero uniqueness scores and unique loadings where expected
  expect_true(all(absolute_sum_U[indicator_c] == 0))
  expect_true(all(mod$Estimate$Unique_loading_estimates[indicator_c] == 0))

  # Non-zero uniqueness scores and unique loadings where expected
  expect_false(is.null(mod$Estimate$Unique_scores))
  expect_false(is.null(mod$Estimate$Unique_loading_estimates))
  expect_true(all(absolute_sum_U[indicator_cf] != 0))
  expect_true(all(mod$Estimate$Unique_loading_estimates[indicator_cf] != 0))
})


## Standardized Construct, Unique and Data --------------------------------

test_that("Returned data, construct scores and unique scores are standardized, as opposed to normalized", {
  # This test assumes that the normalization factor is sqrt(n_case - 1)
  # https://github.com/FloSchuberth/cSEM/issues/581

  normalization_factor <- sqrt(nrow(LeDang2022) - 1)

  normed_sd <- 1 / normalization_factor

  cf_indicator_names <- apply(mod$Estimates$Unique_scores, 2, sum) |>
    vapply(function(x) x != 0, TRUE)

  cf_indicator_names <- names(cf_indicator_names)[cf_indicator_names == TRUE]

  standardized_data_scores <- Reduce(
    cbind,
    list(
      mod$Information$Data,
      mod$Estimates$Construct_scores,
      mod$Estimates$Unique_scores[, cf_indicator_names],
      gsca_mod$Estimates$Construct_scores,
      gsca_mod$Information$Data,
      gsca_m_mod$Estimates$Construct_scores,
      gsca_m_mod$Estimates$Unique_scores,
      gsca_m_mod$Information$Data
    )
  ) |>
    apply(2, sd)

  # If everything is standardized then the standard deviation should be closer to 1 than the normed_sd
  expect_true(
    all(
      abs(standardized_data_scores - 1) <
        abs(standardized_data_scores - normed_sd)
    )
  )
})


## Valid Reliabilities ----------------------------------------------------
test_that("Reliabilities are correctly estimated", {
  expect_true(all(mod$Estimates$Reliabilities <= 1))
})

test_that("Model estimation passes standards", {
  # Counter-intuitively, FALSE means that convergence was OK
  expect_true(all(!verify(mod)))
})


## DGP Data ---------------------------------------------------------------
test_that("Tests pass use of cSEM.DGP", {
  skip_if_not_installed(c("cSEM.DGP", "purrr", "dplyr"))

  #' Build a data frame of expected population values for a single model
  #'
  #' Constructs a reference data frame mapping model terms to their known
  #' population values. Terms are built from the names of the supplied vectors,
  #' so all population value vectors must be named.
  #'
  #' @param mod_name Character. Model identifier for the `mod` column.
  #' @param weights Named list of named numeric vectors. Each element represents
  #'   a composite construct (list name = construct name). Vector names are
  #'   indicator names, values are population weights.
  #'   Example: `list(xi1 = c(x11 = 0.64, x12 = 0.48, x13 = 0.32))`
  #' @param loadings Named list of named numeric vectors. Same structure as
  #'   `weights` but for common factor loadings (uses the `=~` operator).
  #' @param paths Named numeric vector. Names are path terms in `"lhs ~ rhs"`
  #'   format, values are population path coefficients.
  #'   Example: `c("xi2 ~ xi1" = 0.5)`
  #'
  #' @return A data.frame with columns: `mod`, `term`, `pop_value`.
  build_expected_popvalues <- function(
    mod_name,
    weights = NULL,
    loadings = NULL,
    paths = NULL
  ) {
    rows <- list()

    if (!is.null(weights)) {
      for (construct in names(weights)) {
        w <- weights[[construct]]
        stopifnot("weights must be a named vector" = !is.null(names(w)))
        rows <- c(
          rows,
          list(data.frame(
            mod = mod_name,
            term = paste(construct, "<~", names(w)),
            pop_value = unname(w)
          ))
        )
      }
    }

    if (!is.null(loadings)) {
      for (construct in names(loadings)) {
        l <- loadings[[construct]]
        stopifnot("loadings must be a named vector" = !is.null(names(l)))
        rows <- c(
          rows,
          list(data.frame(
            mod = mod_name,
            term = paste(construct, "=~", names(l)),
            pop_value = unname(l)
          ))
        )
      }
    }

    if (!is.null(paths)) {
      stopifnot("paths must be a named vector" = !is.null(names(paths)))
      rows <- c(
        rows,
        list(data.frame(
          mod = mod_name,
          term = names(paths),
          pop_value = unname(paths)
        ))
      )
    }

    return(Reduce(rbind, rows))
  }

  #' Generate expected population values from model name conventions
  #'
  #' Parses model names to automatically determine which constructs get weights
  #' vs loadings, and whether they use single or triple indicators.
  #'
  #' Name format: `"{count}[Type]_{count}[Type]"` where the first part maps to
  #' xi1 and the second to xi2.
  #' - count: `"uni"` (1 indicator) or `"tri"` (3 indicators)
  #'
  #' @param mod_names Character vector of model names (e.g., `"uniC_triF"`).
  #' @param paths Named numeric vector of path coefficients.
  #' @param tri_weights Named list (`xi1`, `xi2`) of population weight vectors
  #'   for triple-indicator composites.
  #' @param tri_loadings Named list (`xi1`, `xi2`) of population loading vectors
  #'   for triple-indicator factors.
  #'
  #' @return A data.frame with columns: `mod`, `term`, `pop_value`.
  make_expected_from_names <- function(
    mod_names,
    paths,
    tri_weights = list(xi1 = xi1_tri_cmp_weights, xi2 = xi2_tri_cmp_weights),
    tri_loadings = list(xi1 = xi1_tri_fct_loadings, xi2 = xi2_tri_fct_loadings)
  ) {
    constructs <- c("xi1", "xi2")

    parse_part <- function(part) {
      count <- sub("[CF]$", "", part)
      type <- sub(".*([CF])$", "\\1", part)
      return(list(count = count, type = type))
    }

    configs <- lapply(mod_names, function(mn) {
      parts <- strsplit(mn, "_")[[1]]
      stopifnot(length(parts) == 2)

      cfg <- list(mod_name = mn, paths = paths)

      for (i in seq_along(parts)) {
        parsed <- parse_part(parts[i])
        xi <- constructs[i]
        uni_name <- paste0("x", i, "1")

        if (parsed$type == "C") {
          val <- if (parsed$count == "uni") {
            setNames(1, uni_name)
          } else {
            tri_weights[[xi]]
          }
          cfg$weights <- c(cfg$weights, setNames(list(val), xi))
        } else {
          val <- if (parsed$count == "uni") {
            setNames(1, uni_name)
          } else {
            tri_loadings[[xi]]
          }
          cfg$loadings <- c(cfg$loadings, setNames(list(val), xi))
        }
      }

      return(cfg)
    })

    configs |>
      lapply(\(cfg) do.call(build_expected_popvalues, cfg)) |>
      (\(dfs) Reduce(rbind, dfs))()
  }

  # General Population Values ----------------------------------------------

  # Note: cSEM.DGP assumes canonical weights

  # These are the actual weights used to construct the population correlation matrix,
  # from which the data is sampled
  Sxi1 <- matrix(c(1, 0.4, -0.3,
                 0.4, 1, 0.4,
                 -0.3, 0.4, 1),
                  3, 3)
  w1 <- c(0.4, 0.3, 0.2)

  w1 %*% Sxi1 %*% w1
  # 0.386

  (cmp_canonical_weights <- w1 / c(sqrt(w1 %*% Sxi1 %*% w1)))
  # [1] 0.6438228 0.4828671 0.3219114

  (single_indicator_weight <- .4 / sqrt(.4 %*% 1 %*% .4))
  # 1

  # Named population vectors for tests
  xi1_tri_cmp_weights <- setNames(cmp_canonical_weights, c("x11", "x12", "x13"))
  xi2_tri_cmp_weights <- setNames(cmp_canonical_weights, c("x21", "x22", "x23"))
  xi1_tri_fct_loadings <- c(x11 = 0.6, x12 = 0.8, x13 = 0.7)
  xi2_tri_fct_loadings <- c(x21 = 0.6, x22 = 0.8, x23 = 0.7)
  path_xi2_xi1 <- c("xi2 ~ xi1" = 0.5)

  # IGSCA ------------------------------------------------------------------

  igsca_pop <- list(
    uniC_uniF = 'xi1 <~ .4*x11

               xi2 =~ 1*x21

               xi2 ~ .5*xi1',
    uniC_triF = 'xi1 <~ .4*x11

               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23

               xi2 ~ .5*xi1',
    triC_uniF = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    triC_triF = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23
               
               xi2 ~ 0.5*xi1',
    uniF_uniC = 'xi1 =~ 1*x11

               xi2 <~ .4*x21

               xi2 ~ .5*xi1',
    uniF_triC = 'xi1 =~ 1*x11

               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23

               xi2 ~ .5*xi1',
    triF_uniC = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    triF_triC = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23
               
               xi2 ~ 0.5*xi1'
  )

  igsca_datapop <- lapply(igsca_pop, cSEM.DGP::generateData, .empirical = TRUE)

  igsca_model_spec <- list(
    uniC_uniF = 'xi1 <~ x11

               xi2 =~ x21

               xi2 ~ xi1',
    uniC_triF = 'xi1 <~ x11

               xi2=~x21 + x22 + x23

               xi2 ~ xi1',
    triC_uniF = 'xi1<~x11 + x12 + x13
               
               xi2 =~ x21
               
               xi2 ~ xi1',
    triC_triF = 'xi1<~x11 + x12 + x13
               
               xi2=~x21 + x22 + x23
               
               xi2 ~ xi1',
    uniF_uniC = 'xi1 =~ x11

               xi2 <~ x21

               xi2 ~ xi1',
    uniF_triC = 'xi1 =~ x11

               xi2<~x21 + x22 + x23

               xi2 ~ xi1',
    triF_uniC = 'xi1=~x11 + x12 + x13
               
               xi2 <~ x21
               
               xi2 ~ xi1',
    triF_triC = 'xi1=~x11 + x12 + x13
               
               xi2<~x21 + x22 + x23
               
               xi2 ~ xi1'
  )

  igsca_mods <- mapply(
    cSEM::csem,
    .data = igsca_datapop,
    .model = igsca_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.001,
    .iter_max = 1000
  )
  # Note that the uniC_uniF can sometimes create false positive convergence failures despite appearing to produce correct parameter estimates
  tidied_igsca_mods <- lapply(igsca_mods, function(x) {
    tidy(x) |>
      dplyr::filter(op %in% c('=~', '~', '<~')) |>
      dplyr::select(term, estimate)
  }) |>
    purrr::list_rbind(names_to = 'mod') |>
    dplyr::filter(
      !((grepl('xi2 <~', term, fixed = TRUE) &
        grepl(pattern = '_...F', x = mod, fixed = FALSE)) &
        (grepl('xi1 <~', term, fixed = TRUE) &
          grepl(pattern = '...F_', x = mod, fixed = FALSE)))
    )

  ## Test IGSCA parameter recovery ------------------------------------------
  # Note: There seems to sometimes be convergence problems when using factors
  # with only one indicator, but not always.

  igsca_expected <- make_expected_from_names(
    c(
      "uniC_uniF",
      "uniC_triF",
      "triC_uniF",
      "triC_triF",
      "uniF_uniC",
      "uniF_triC",
      "triF_uniC",
      "triF_triC"
    ),
    paths = path_xi2_xi1
  )

  igsca_joined <- merge(
    tidied_igsca_mods,
    igsca_expected,
    by = c("mod", "term")
  )

  testthat::test_that("IGSCA recovers population weights, loadings, and path coefficients", {
    testthat::expect_equal(
      nrow(igsca_joined),
      nrow(igsca_expected),
      info = "All expected terms should be present in tidy output"
    )
    is_loading <- grepl("=~", igsca_joined$term, fixed = TRUE)
    is_weight <- grepl("<~", igsca_joined$term, fixed = TRUE)
    is_path <- !is_loading & !is_weight
    testthat::expect_equal(
      igsca_joined$estimate[is_weight],
      igsca_joined$pop_value[is_weight] #,
      # tolerance = 0.0001
    )
    testthat::expect_equal(
      igsca_joined$estimate[is_path],
      igsca_joined$pop_value[is_path],
      tolerance = 0.01
    )
    testthat::expect_equal(
      igsca_joined$estimate[is_loading],
      igsca_joined$pop_value[is_loading],
      tolerance = 0.02
    )
  })

  # View(igsca_joined)

  # GSCA -------------------------------------------------------------------

  gsca_pops <- list(
    uniC_uniC = 'xi1 <~ .4*x11

               xi2 <~ .3*x21

               xi2 ~ .5*xi1',
    uniC_triC = 'xi1 <~ .4*x11

               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23

               xi2 ~ .5*xi1',
    triC_uniC = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2 <~ .3*x21
               
               xi2 ~ .5*xi1',
    triC_triC = 'xi1<~0.4*x11 + 0.3*x12 + 0.2*x13
               x11~~0.4*x12 + -0.3*x13
               x12~~0.4*x13
               
               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23
               
               xi2 ~ 0.5*xi1'
  )

  gsca_datapop <- lapply(gsca_pops, cSEM.DGP::generateData, .empirical = TRUE)

  gsca_model_spec <- list(
    uniC_uniC = 'xi1 <~ x11

               xi2 <~ x21

               xi2 ~ xi1',
    uniC_triC = 'xi1 <~ x11

               xi2<~x21 + x22 + x23

               xi2 ~ xi1',
    triC_uniC = 'xi1<~x11 + x12 + x13
               
               xi2 <~ x21
               
               xi2 ~ xi1',
    triC_triC = 'xi1<~x11 + x12 + x13
               
               xi2<~x21 + x22 + x23
               
               xi2 ~ xi1'
  )

  gsca_mods <- mapply(
    cSEM::csem,
    .data = gsca_datapop,
    .model = gsca_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = FALSE,
    .GSCA_modes = "CCMP"
  )

  tidied_gsca_mods <- lapply(gsca_mods, function(x) {
    tidy(x) |>
      dplyr::filter(op %in% c('<~', '~')) |>
      dplyr::select(term, estimate)
  }) |>
    purrr::list_rbind(names_to = 'mod')

  ## Test GSCA parameter recovery -------------------------------------------

  gsca_expected <- make_expected_from_names(
    c(
      "uniC_uniC",
      "uniC_triC",
      "triC_uniC",
      "triC_triC"
    ),
    paths = path_xi2_xi1
  )

  gsca_joined <- merge(tidied_gsca_mods, gsca_expected, by = c("mod", "term"))

  testthat::test_that("GSCA recovers population weights and path coefficients", {
    testthat::expect_equal(
      nrow(gsca_joined),
      nrow(gsca_expected),
      info = "All expected terms should be present in tidy output"
    )
    testthat::expect_equal(
      gsca_joined$estimate,
      gsca_joined$pop_value
    )
  })

  # GSCA M -----------------------------------------------------------------
  set.seed(3883)
  gscam_pop <- list(
    uniF_uniF = 'xi1 =~ 1*x11

               xi2 =~ 1*x21

               xi2 ~ .5*xi1',
    uniF_triF = 'xi1 =~ 1*x11

               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23

               xi2 ~ .5*xi1',
    triF_uniF = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2 =~ 1*x21
               
               xi2 ~ .5*xi1',
    triF_triF = 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2=~0.6*x21 + 0.8*x22 + 0.7*x23
               
               xi2 ~ 0.5*xi1'
  )

  gscam_datapop <- lapply(gscam_pop, cSEM.DGP::generateData, .empirical = TRUE)

  gscam_model_spec <- list(
    uniF_uniF = 'xi1 =~ x11

               xi2 =~ x21

               xi2 ~ xi1',
    uniF_triF = 'xi1 =~ x11

               xi2=~x21 + x22 + x23

               xi2 ~ xi1',
    triF_uniF = 'xi1=~x11 + x12 + x13
               
               xi2 =~ x21
               
               xi2 ~ xi1',
    triF_triF = 'xi1=~x11 + x12 + x13
               
               xi2=~x21 + x22 + x23
               
               xi2 ~ xi1'
  )

  gscam_mods <- mapply(
    cSEM::csem,
    .data = gscam_datapop,
    .model = gscam_model_spec,
    SIMPLIFY = FALSE,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute"
  )

  tidied_gscam_mods <- lapply(gscam_mods, function(x) {
    tidy(x) |>
      dplyr::filter(op %in% c('=~', '~')) |>
      dplyr::select(term, estimate)
  }) |>
    purrr::list_rbind(names_to = 'mod')

  ## Test GSCAM parameter recovery ------------------------------------------
  # Note: The current GSCA_M may sometimes have trouble in computing D2 and U,
  # especially when there's only one indicator for a common factor with a small loading

  gscam_expected <- make_expected_from_names(
    c(
      "uniF_uniF",
      "uniF_triF",
      "triF_uniF",
      "triF_triF"
    ),
    paths = path_xi2_xi1,
  )

  gscam_joined <- merge(
    tidied_gscam_mods,
    gscam_expected,
    by = c("mod", "term")
  )

  testthat::test_that("GSCAM recovers population loadings and path coefficients", {
    testthat::expect_equal(
      nrow(gscam_joined),
      nrow(gscam_expected),
      info = "All expected terms should be present in tidy output"
    )
    is_loading <- grepl("=~", gscam_joined$term, fixed = TRUE)
    is_path <- !is_loading
    testthat::expect_equal(
      gscam_joined$estimate[is_path],
      gscam_joined$pop_value[is_path],
      tolerance = 0.02
    )
    testthat::expect_equal(
      gscam_joined$estimate[is_loading],
      gscam_joined$pop_value[is_loading],
      tolerance = 0.01
    )
  })
})


#### Custom Function for Organizing IGSCA Results ----------------------------

#' Converts Output of igsca functions into a table to facilitate comparisons
#'
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#'
#' I chose not to use the tidy method for the tests because it would be a lot of
#' work to make the GSCAPro results fit in the tidy format.
#'
#' @param weights Weights matrix
#' @param loadings Loadings matrix
#' @param uniqueD Vector of Uniqueness for each indicator of a common factor
#' @param paths Path coefficients matrix
#' @importFrom lavaan lavaanify
#' @author Michael S. Truong
#' @return Table of Weights, Loadings, Path-Coefficients and Uniqueness terms from i-gsca algorithms in Matlab or R.
#'
get_lavaan_table_igsca_matrix <- function(
  model,
  weights,
  loadings,
  uniqueD,
  paths
) {
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"), ]
  # Pre-allocate Columns
  table <-
    cbind(
      table,
      list(
        "weights" = 0,
        "loadings" = 0,
        "uniqueD" = 0,
        "paths" = 0
      )
    )

  # Slide in weights
  for (indicator in rownames(weights)) {
    for (lv in colnames(weights)) {
      table[
        ((table$lhs == lv &
          table$rhs == indicator) &
          table$op %in% c("<~", "=~")),
        "weights"
      ] <- weights[indicator, lv]
    }
  }

  # Slide in loadings
  for (indicator in rownames(loadings)) {
    for (lv in colnames(loadings)) {
      table[
        ((table$lhs == lv &
          table$rhs == indicator) &
          table$op %in% c("<~", "=~")),
        "loadings"
      ] <- loadings[indicator, lv]
    }
  }

  # Slide in uniqueD
  for (indicator in names(uniqueD)) {
    # This assumes that every indicator only loads onto one factor
    # Cross-factor loadings will not work with this
    table[
      ((table$rhs == indicator) &
        (table$op == "=~")),
      "uniqueD"
    ] <- uniqueD[indicator]
  }

  # Slide in Paths
  for (lv_from in rownames(paths)) {
    for (lv_to in colnames(paths)) {
      table[
        ((table$rhs == lv_from &
          table$lhs == lv_to) &
          table$op == "~"),
        "paths"
      ] <- paths[lv_from, lv_to]
    }
  }

  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA

  return(table)
}


#### Fetching igsca_r_table Results ------------------------------------------
igsca_r_table <- with(
  mod$Estimates,
  get_lavaan_table_igsca_matrix(
    model = tutorial_igsca_model,
    weights = t(Weight_estimates),
    loadings = t(Loading_estimates),
    uniqueD = Unique_loading_estimates,
    paths = t(Path_estimates)
  )
)

## Comparisons between cSEM::igsca(), GSCAPro and igsca_sim.-----------------

### Load Matlab Results -----------------------------------------------------
# Loads into igsca_sim_m_table
load(testthat::test_path("data", "igsca_matlab.RData"))
# The original matlab code squares the unique loadings
igsca_sim_m_table$uniqueD <- sqrt(igsca_sim_m_table$uniqueD)

### Load GSCAPro Results ----------------------------------------------------
# Loads into igsca_gscapro
load(testthat::test_path("data", "igsca_gscapro.RData"))
# Presumably, the original GSCA Pro reports the squared unique loadings
igsca_gscapro$uniqueD <- sqrt(igsca_gscapro$uniqueD)


### Compare Matlab and cSEM::igsca()------------------------------------------
testthat::expect_failure(testthat::expect_equal(
  object = igsca_r_table,
  expected = igsca_sim_m_table,
  tolerance = .00034
))

testthat::expect_equal(
  object = igsca_r_table,
  expected = igsca_sim_m_table,
  tolerance = .006
)

# If the kronecker bypass for Loadings is done:
# testthat::expect_success(testthat::expect_equal(object = igsca_r_table,
#                                                 expected = igsca_sim_m_table, tolerance = .037))

testthat::expect_failure(
  testthat::expect_identical(
    object = igsca_r_table,
    expected = igsca_sim_m_table
  ),
  info = "Matlab and R versions should be very similar, but not identical"
)

# waldo::compare(igsca_sim_m_table, igsca_r_table, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_r_table)

### GSCAPro and R ---------------------------------------------------
testthat::expect_success(testthat::expect_equal(
  igsca_r_table,
  igsca_gscapro,
  tolerance = 0.0335
))

# waldo::compare(igsca_r_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_r_table, igsca_gscapro)

### Compare GSCAPro and Matlab ----------------------------------------------
testthat::expect_success(testthat::expect_equal(
  igsca_sim_m_table,
  igsca_gscapro,
  tolerance = 0.0335
))

# waldo::compare(igsca_sim_m_table, igsca_gscapro, max_diffs = Inf)

# all.equal(igsca_sim_m_table, igsca_gscapro)

## Comparing Different Ways of Fitting Group Models ------------------------

## Starting Values --------------------------------------------------------
test_that("Starting Values run with IGSCA and IGSCA Primer results are replicated", {
  # Using IGSCA  Primer results
  data(corp_rep_data, package = "cSEM")

  # Code from https://osf.io/9tm2y/files/wk5vz?view_only=59d9792bab994892a735e4efc763b511
  # Schamberger et al. (2025)
  x_clean <- corp_rep_data |>
    dplyr::mutate(across(everything(), ~ ifelse(.x == -99, NA, .x))) |>
    na.omit()

  model <- '
Quality <~ qual_1 + qual_2 + qual_3 + qual_4 + qual_5 + qual_6 + qual_7 + qual_8
Performance <~ perf_1 + perf_2 + perf_3 + perf_4 + perf_5
CorpSocResp <~ csor_1 + csor_2 + csor_3 + csor_4 + csor_5
Attractiveness <~ attr_1 + attr_2 + attr_3

Competence =~ comp_1 + comp_2 + comp_3
CustomerLoyality =~ cusl_1 + cusl_2 + cusl_3
Likeability =~ like_1 + like_2 + like_3
CustomerSatisfaction =~ cusa

Competence ~ Quality + Performance + CorpSocResp + Attractiveness
Likeability ~ Quality + Performance + CorpSocResp + Attractiveness
CustomerSatisfaction ~ Competence + Likeability
CustomerLoyality ~ CustomerSatisfaction + Competence + Likeability'

  igsca <- csem(
    x_clean,
    model,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .dominant_indicators = NULL,
    .tolerance = 0.001,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "NCMP",
    .starting_values = list(
      "Quality" = c(
        "qual_1" = 0.173,
        "qual_2" = 0.143,
        "qual_3" = 0.184,
        "qual_4" = 0.161,
        "qual_5" = 0.179,
        "qual_6" = 0.193,
        "qual_7" = 0.177,
        "qual_8" = 0.142
      ),
      "Performance" = c(
        "perf_1" = 0.301,
        "perf_2" = 0.313,
        "perf_3" = 0.243,
        "perf_4" = 0.285,
        "perf_5" = 0.268
      ),
      "CorpSocResp" = c(
        "csor_1" = 0.252,
        "csor_2" = 0.252,
        "csor_3" = 0.282,
        "csor_4" = 0.258,
        "csor_5" = 0.27
      ),
      "Attractiveness" = c("attr_1" = 0.469, "attr_2" = 0.404, "attr_3" = 0.463)
    )
  )
  tidied_igsca <- tidy(igsca)

  igscaPrimer <- read.csv(testthat::test_path(
    "data",
    "corp_rep_igscaPrimer.csv"
  )) |>
    subset(select = c(term, estimate))

  tidied_igsca <- subset(
    tidied_igsca,
    term %in% igscaPrimer$term & (op %in% c("~", "=~", "<~")),
    select = c(term, estimate)
  )

  tidied_igsca <- tidied_igsca[order(tidied_igsca$term), ]

  igscaPrimer <- igscaPrimer[order(igscaPrimer$term), ]

  testthat::expect_equal(
    object = tidied_igsca,
    expected = igscaPrimer,
    tolerance = .012,
    ignore_attr = TRUE
  )
})