#' Calculate the difference in FIT
#'
#'
#' @param data Data passed by boot
#' @param idx Row-indices passed by boot
#' @inheritParams csem
#' @return Named vector of the delta between the FIT of the multi-group model
#'   versus the FIT of the single-group model--positive values favor
#'   multi-group, negative, single-group; the FIT of the single-group model; and
#'   the FIT of the multi-group model.
#' @keywords internal
#'
calculateFITForSplit <- function(data,
                                 idx,
                                 .model,
                                 .id,
                                 .approach_weights,
                                 .conv_criterion,
                                 .disattenuate,
                                 .dominant_indicators,
                                 .iter_max,
                                 .tolerance) {
  
  # Robust Boot by-pass by Thomas on June 11,2013 https://stackoverflow.com/a/17040580
  ret <- tryCatch({
    sg_mod <- csem(
      .data = data[idx, ],
      .model = .model,
      .approach_weights = .approach_weights,
      .tolerance = .tolerance,
      .conv_criterion = .conv_criterion
    )
    sg_fit <- calculateFIT(sg_mod)
    
    if (length(.id) == 1) {
      mg_mod <- csem(
        .data = data[idx, ],
        .id = .id,
        .model = .model,
        .approach_weights = .approach_weights,
        .tolerance = .tolerance,
        .conv_criterion = .conv_criterion
      )
      
      mg_fit <- bdiagGSCA(mg_mod) |>
        calculateFIT()
      
      names(mg_fit) <- .id
      
      return(c("sg_fit" = sg_fit, mg_fit))
      
    } else if (length(.id) > 1) {
      mg_fits <- lapply(.id, function(.id_iter) {
        mg_mod <- csem(
          .data = data[idx, ],
          .id = .id_iter,
          .model = .model,
          .approach_weights = .approach_weights,
          .tolerance = .tolerance,
          .conv_criterion = .conv_criterion
        )
        
        mg_fit <- bdiagGSCA(mg_mod) |>
          calculateFIT()
        names(mg_fit) <- .id_iter
        
        return(mg_fit)
      })
      
      return(c("sg_fit" = sg_fit, unlist(mg_fits)))
      
    } else {
      stop("Inappropriate length or type of .splitvars passed")
      
    }
  }, error = function(error) {
    return(c("sg_fit" = NA, "mg_fit" = NA))
  })

  return(ret)

}


#' Block-diagonalized model-implied correlation matrices
#'
#' Make it easy to compare the matrix distance between the model-implied
#' correlation matrix of a single-group model and the model-implied matrices
#' of a multigroup model: the single-group model's implied matrix is
#' duplicated block-diagonally a user-chosen number of times, while for a
#' multigroup model each group's implied matrix is placed along the block
#' diagonal. The block-diagonal matrices of a single-group model (repeated
#' as many times as the multigroup model has groups) and of the multigroup
#' model then have identical dimensions and can be compared with any matrix
#' distance.
#'
#' @usage bdiagFit(
#'   .object    = NULL,
#'   .n_blocks  = args_default()$.n_blocks,
#'   .type_vcv  = args_default()$.type_vcv,
#'   .saturated = args_default()$.saturated
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A block-diagonal `matrix`.
#' @seealso [fit()], [csem()]
#' @export
#' @importFrom Matrix bdiag
bdiagFit <- function(.object    = NULL,
                     .n_blocks  = args_default()$.n_blocks,
                     .type_vcv  = args_default()$.type_vcv,
                     .saturated = args_default()$.saturated) {

## Warning Checks ----------------------------------------------------------
  stopifnot(
    '`.n_blocks` must be a single positive whole number' =
      length(.n_blocks) == 1 && is.numeric(.n_blocks) &&
      !is.na(.n_blocks) && .n_blocks > 0 && .n_blocks == round(.n_blocks)
  )

## Multigroup objects -------------------------------------------------------
  if(inherits(.object, "cSEMResults_multi")) {

    if(.n_blocks != 1) {
      warning(paste0(".n_blocks is ignored for multigroup objects: the ",
                     "number of blocks is set to the number of groups."))
    }

    Sigma_list <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)

    Sigma_bdiag <- Matrix::bdiag(Sigma_list) |>
      as.matrix()

    new_names <- unlist(lapply(names(Sigma_list), function(g) {
      paste0(rownames(Sigma_list[[g]]), "_", g)
    }))

    dimnames(Sigma_bdiag) <- list(new_names, new_names)

    return(Sigma_bdiag)
  }

## Single-group objects (including 2nd-order objects) -----------------------
  Sigma <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)

  Sigma_bdiag <- Matrix::bdiag(replicate(.n_blocks, Sigma, simplify = FALSE)) |>
    as.matrix()

  new_names <- unlist(lapply(seq_len(.n_blocks), function(i) {
    paste0(rownames(Sigma), "_", i)
  }))

  dimnames(Sigma_bdiag) <- list(new_names, new_names)

  return(Sigma_bdiag)
}


igsca_tree_control <- function(alpha = 0.05,
                               bonferroni = TRUE,
                               minbucket = 30L,
                               minsplit = 2L * minbucket,
                               maxdepth = 3L,
                               max_cuts = 20L,
                               R_test = 500L,
                               coin_distribution = c("approximate", "asymptotic")) {
  list(
    alpha = alpha,
    bonferroni = bonferroni,
    minbucket = as.integer(minbucket),
    minsplit = as.integer(minsplit),
    maxdepth = as.integer(maxdepth),
    max_cuts = as.integer(max_cuts),
    R_test = as.integer(R_test),
    coin_distribution = match.arg(coin_distribution)
  )
}

#' Casewise sum of squared GSCA residuals (n x 1) -- COIN_ssr influence.
influence_ssr <- function(E) {
  matrix(rowSums(E^2), ncol = 1L)
}


#' Casewise squared-residual matrix (n x q, multivariate) -- COIN_mat.
influence_mat <- function(E) {
  E^2
}

#' `.id = NULL` yields the pooled single-group fit; `.id = "group"` the MGA fit.
fit_csem <- function(.data, .model, .id = NULL) {
  csem(
    .data = .data,
    .model = .model,
    .id = .id,
    .approach_weights = "GSCA",
    .disattenuate = TRUE, # to get igsca
    .conv_criterion = "sum_diff_absolute", # Default in gsca_m.m and gsca.m
    .iter_max = 100, # Default in igsca_sim.m
    .GSCA_modes = "CCMP", # Unknown (to me) if it affects convergence
    .tolerance = 0.0001 # Default in gsca_m.m and gsca.m
  )
}

try_fit <- function(.data, .model, .id = NULL) {
  fit <- suppressWarnings(tryCatch(
    fit_csem(.data, .model, .id),
    error = function(e) NULL
  ))
  ok <- !is.null(fit) &&
    isTRUE(tryCatch(sum(verify(fit)), error = function(e) FALSE))
  list(fit = fit, ok = ok)
}

