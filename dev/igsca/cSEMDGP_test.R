# Test IGSCA
library(cSEM.DGP)
library(lavaan)
library(dplyr)
library(purrr)
# library(cSEM)
devtools::load_all()
source("dev/igsca/helpers.R")


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

tidied_igsca_mods <- lapply(igsca_mods, function(x) {
    tidy(x) |>
        dplyr::filter(op %in% c('=~', '~', '<~')) |>
        dplyr::select(term, estimate)
}) |>
    list_rbind(names_to = 'mod') |>
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

igsca_joined <- merge(tidied_igsca_mods, igsca_expected, by = c("mod", "term"))

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
    tolerance = 0.001
  )
  testthat::expect_equal(
    igsca_joined$estimate[is_loading],
    igsca_joined$pop_value[is_loading],
    tolerance = 0.001
  )
})

View(igsca_joined)

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
  list_rbind(names_to = 'mod')

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
set.seed(1234)
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
  list_rbind(names_to = 'mod')

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

gscam_joined <- merge(tidied_gscam_mods, gscam_expected, by = c("mod", "term"))

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
    tolerance = 0.001
  )
  testthat::expect_equal(
    gscam_joined$estimate[is_loading],
    gscam_joined$pop_value[is_loading],
    tolerance = 0.01
  )
})
