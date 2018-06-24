#' Test measurement invariance of composites
#'
#' This functions performs the test for measurement invariance of composites
#' proposed by Henseler et al. (2016).
#'
#' The test is only meaningful for composite models.
#'
#' If more than two groups are to be compared issus related to multiple testing
#' should be taken into account.
#'
#' The number of permutations is defaults to 100 for performance reasons.
#' According to Henseler et al. (2016) the number of permutations should be at least 5000 for
#' assessement to be reliable.
#'
#' @usage testMICOM(
#'   .data             = NULL,
#'   .model            = NULL,
#'   .group_var        = NULL,
#'   .approach_weights = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR",
#'                         "GSCA", "fixed", "unit")
#'   .permutations     = 100,
#'   .alpha            = c(0.1, 0.05, 0.01),
#'   .verbose          = TRUE,
#'   ...)
#'
#' @inheritParams csem_arguments
#'
#' @return An object of class "testMICOM"
#'
#' An object of class "testMICOM" is a named list containing the following list elemement:
#' \describe{
#' \item{**Step2**}{A list containing the results of the test for composistional invariance (Step 2).}
#' \item{**Step3**}{A list containing the results of the test for mean and variance equality (Step 3).}
#' \item{**Meta_information**}{A list of additional information on the test.}
#' }
#'
#' @export
#'
#' @examples
#'
#' still to implement
#'

testMICOM <- function(.data             = NULL,
                      .model            = NULL,
                      .group_var        = NULL,
                      .approach_weights = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR",
                                            "GSCA", "fixed", "unit"),
                      .PLS_weight_scheme_inner = c("centroid", "factorial", "path"),
                      .permutations      = 100,
                      .alpha             = c(0.1, 0.05, 0.01),
                      .verbose           = TRUE,
                      ...) {

  # Implementation and notation is based on:
  # Henseler et al. (2016) - Testing measurement invariance of composites using
  #                          partial least squares

  ### Information ==============================================================
  if(.verbose) {
    cat(cli::rule(center = "Test for measurement invariance based on Henseler et al. (2016)",
                  line = "bar3"), "\n\n")
    cat(
      "Note:\t", "Due to the permutation runs the test procedure is computationally intensive.\n\t",
      "Depending on the complexity of the model, the number of permutations\n\t",
      "choosen, and the available hardware the procedure may take a while.\n\n",
      sep = "")
  }

  ### Preparation ==============================================================
  ## Parse model and capture arguments
  csem_model <- parseModel(.model)

  ## Warnings, messages, and errors
  if(!(.group_var %in% names(.data))) {
    stop("No grouping variable provided.")
  }

  if(all(csem_model$construct_type$Type == "Common factor")) {
    warning("\tAll constructs are modelled as common factors.\n",
            "\tTest results are only meaningful for composite models!",
            call. = FALSE)
  }

  if(length(unique(.data[, .group_var])) > 2) {
    warning("\tComparing multiple groups inflates the familywise error rate.\n",
            "\tInterpret statistical significance with caution.\n",
            "\n\tFuture versions of the package will likely include\n",
            "\tappropriate correction options.",
            call = FALSE)
  }

  approach_weights        <- match.arg(.approach_weights)
  if(approach_weights == "PLS") {
    PLS_weight_scheme_inner <- match.arg(.PLS_weight_scheme_inner)
  } else {
    PLS_weight_scheme_inner <- NULL
  }

  ## Definitions
  .data <- .data[, c(colnames(csem_model$measurement), .group_var)]
  data_no_group_var <- .data[, -which(names(.data) == .group_var)]
  data_split        <- split(x = .data, f = .data[, .group_var])
  data_ordered      <- .data[order(.data[, .group_var]), ]

  no_obs            <- data.frame("x" = c("Total", unique(data_ordered[, .group_var])),
                                  "n" = c(nrow(.data), table(data_ordered[, .group_var])),
                                  row.names = NULL,
                                  stringsAsFactors = FALSE)
  ### Step 1 - Configural invariance ===========================================

  # Has to be assessed by the user prior to using the testMICOM function. See
  # the original paper for details.

  ### Step 2 - Compositional invariance ========================================
  ## Compute weights for each group and use these to compute proxies/scores using
  # the pooled data

  scores <- lapply(data_split, function(x, ...) {
    W <- workhorse(
      .data                    = x[, -which(names(x) == .group_var)],
      .model                   = csem_model,
      .approach_cf             = "dist_euclid",
      .approach_nl             = NULL,
      .approach_paths          = NULL,
      .approach_weights        = approach_weights,
      .disattenuate            = FALSE,
      .estimate_structural     = FALSE,
      .ignore_structural_model = FALSE,
      .iter_max                = 100,
      .normal                  = TRUE,
      .PLS_mode                = NULL,
      .PLS_weight_scheme_inner = "centroid",
      .tolerance               = 1e-06
    )$Estimates$Weight_estimates

    # W <- workhorse(.data = x[, -which(names(x) == .group_var)],
    #                .model = csem_model,
    #                .approach_weights = approach_weights,
    #                .PLS_weight_scheme_inner = PLS_weight_scheme_inner,
    #                .estimate_structural = FALSE,
    #                ...)$Estimates$Weight_estimates

    H <- scale(as.matrix(data_no_group_var)) %*% t(W)

  })

  ## Compute the correlation of the scores for all group combinations
  # Get the scores for all group combinations
  temp <- combn(scores, 2, simplify = FALSE)

  # Compute the correlation c for each group combination
  c <- lapply(temp, function(x) diag(cor(x[[1]], x[[2]])))

  # Set the names for each group combination
  names(c) <- combn(names(scores), 2, FUN = paste0, collapse = "_", simplify = FALSE)

  ## Permutation ---------------------------------------------------------------
  # Procedure: see page 414 and 415 of the paper
  if(.verbose) {
    cat("Progress:\n")
  }

  # Run all permuations
  perm1 <- lapply(1:.permutations, function(x) {

    if(.verbose) {
      cat(x, " / ", .permutations, "\n", sep = "")
    }

    d <- data_no_group_var[sample(x = nrow(.data)), ]
    temp <- split(x = d, rep(unique(data_ordered[, .group_var]),
                             times = table(data_ordered[, .group_var])))

    scores <- lapply(temp, function(y) {

      W <- workhorse(
        .data                    = y,
        .model                   = csem_model,
        .approach_cf             = "dist_euclid",
        .approach_nl             = NULL,
        .approach_paths          = NULL,
        .approach_weights        = approach_weights,
        .disattenuate            = FALSE,
        .estimate_structural     = FALSE,
        .ignore_structural_model = FALSE,
        .iter_max                = 100,
        .normal                  = TRUE,
        .PLS_mode                = NULL,
        .PLS_weight_scheme_inner = "centroid",
        .tolerance               = 1e-06
      )$Estimates$Weight_estimates

      H <- scale(as.matrix(d)) %*% t(W)

    })
    # Make the combinations of list elements
    temp <- combn(scores, 2, simplify = FALSE)

    # Compute correlation
    c_u <- lapply(temp, function(x) diag(cor(x[[1]], x[[2]])))

    # Get the names
    names(c_u) <- combn(names(scores), 2, FUN = paste0, collapse = "_", simplify = FALSE)

    return(c_u)
  })

  ## Bring data to form and compute quantiles
  temp <- do.call(rbind, lapply(perm1, function(x) do.call(rbind, x)))
  temp <- split(as.data.frame(temp), rownames(temp))

  step2_out <- lapply(temp, function(x) as.data.frame(t(apply(x, 2, quantile, probs = .alpha))))
  step2_out <- mapply(function(x, y) cbind("c" = y, x),
                      x = step2_out,
                      y = c,
                      SIMPLIFY = FALSE)

  if(.verbose) {
    cat("\nCompositional invariance test finished.\n",
        "Proceding to test equality of composite mean values and variances.\n", sep = "")
  }

  ### Step 3 - Equal mean values and variances==================================

  H <- workhorse(
    .data                    = data_ordered[, -which(names(data_ordered) == .group_var)],
    .model                   = csem_model,
    .approach_cf             = "dist_euclid",
    .approach_nl             = NULL,
    .approach_paths          = NULL,
    .approach_weights        = approach_weights,
    .disattenuate            = FALSE,
    .estimate_structural     = FALSE,
    .ignore_structural_model = FALSE,
    .iter_max                = 100,
    .normal                  = TRUE,
    .PLS_mode                = NULL,
    .PLS_weight_scheme_inner = "centroid",
    .tolerance               = 1e-06
  )$Estimates$Proxies

  # H <- workhorse(.data  = data_ordered[, -which(names(data_ordered) == .group_var)],
  #                .model = csem_model,
  #                .approach_weights    = approach_weights,
  #                .PLS_weight_scheme_inner = PLS_weight_scheme_inner,
  #                .estimate_structural = FALSE,
  #                ...)$Estimates$Proxies

  temp <- split(x = as.data.frame(H), rep(unique(data_ordered[, .group_var]),
                                          times = table(data_ordered[, .group_var])))
  temp_m <- lapply(temp, colMeans)
  temp_v <- lapply(temp, function(x) {
    a <- matrixStats::colVars(as.matrix(x))
    names(a) <- colnames(x)
    a
  })

  # Make the combinations of list elements
  temp_m <- combn(temp_m, 2, simplify = FALSE)
  temp_v <- combn(temp_v, 2, simplify = FALSE)

  # Compute difference for each group combination
  m <- lapply(temp_m, function(x) x[[1]] - x[[2]])
  v <- lapply(temp_v, function(x) x[[1]] - x[[2]])

  # Get the names
  names(m) <- names(v) <- combn(names(temp), 2, FUN = paste0, collapse = "_", simplify = FALSE)

  ### Permutation --------------------------------------------------------------
  if(.verbose) {
    cat("\nProgress:\n", sep = "")
  }

  perm2 <- lapply(1:.permutations, function(x) {

    if(.verbose & x %in% seq(0, .permutations, by = 10)) {
      cat(x, " / ", .permutations, "\n", sep = "")
    }

    H <- workhorse(
      .data                    = data_no_group_var[sample(x = nrow(data_no_group_var)), ],
      .model                   = csem_model,
      .approach_cf             = "dist_euclid",
      .approach_nl             = NULL,
      .approach_paths          = NULL,
      .approach_weights        = approach_weights,
      .disattenuate            = FALSE,
      .estimate_structural     = FALSE,
      .ignore_structural_model = FALSE,
      .iter_max                = 100,
      .normal                  = TRUE,
      .PLS_mode                = NULL,
      .PLS_weight_scheme_inner = "centroid",
      .tolerance               = 1e-06
    )$Estimates$Proxies

    temp <- split(x = as.data.frame(H), rep(unique(data_ordered[, .group_var]),
                                            times = table(data_ordered[, .group_var])))
    temp_m <- lapply(temp, colMeans)
    temp_v <- lapply(temp, function(x) {
      a <- matrixStats::colVars(as.matrix(x))
      names(a) <- colnames(x)
      a
    })

    # Make the combinations of list elements
    temp_m <- combn(temp_m, 2, simplify = FALSE)
    temp_v <- combn(temp_v, 2, simplify = FALSE)

    # Compute difference for each group combination
    m <- lapply(temp_m, function(x) x[[1]] - x[[2]])
    v <- lapply(temp_v, function(x) log(x[[1]]) - log(x[[2]]))

    # Get the names
    names(m) <- names(v) <- combn(names(temp), 2, FUN = paste0, collapse = "_", simplify = FALSE)
    list("m" = m, "v" = v)
  })

  temp_mv <- lapply(perm2, function(x) lapply(x, function(y) do.call(rbind, y)))
  temp_mv <- lapply(purrr::transpose(temp_mv), function(x) do.call(rbind, x))

  temp_m <- split(as.data.frame(temp_mv[["m"]]), names(m))
  temp_v <- split(as.data.frame(temp_mv[["v"]]), names(v))

  step3_out_m <- lapply(temp_m, function(x) as.data.frame(t(apply(x, 2, quantile, probs = c(.alpha/2, 1 - .alpha/2)))))
  step3_out_m <- mapply(function(x, y) cbind("Diff_mean" = y, x),
                      x = step3_out_m,
                      y = m,
                      SIMPLIFY = FALSE)

  step3_out_v <- lapply(temp_v, function(x) as.data.frame(t(apply(x, 2, quantile, probs = c(.alpha/2, 1 - .alpha/2)))))
  step3_out_v <- mapply(function(x, y) cbind("Diff_log_var" = y, x),
                        x = step3_out_v,
                        y = v,
                        SIMPLIFY = FALSE)
  ### Results ==================================================================
  out <- list("Step2"            = step2_out,
              "Step3"            = list("Mean_diff" = step3_out_m,
                                        "Var_diff"  = step3_out_v),
              "Meta_information" = list("Number_of_observations" = no_obs,
                                        "Number_of_Groups"       = nrow(no_obs) - 1,
                                        "Grouping_variable"      = .group_var)
              )
  ## Set class
  class(out) <- "testMICOM"
  return(out)
}

#' @export
print.testMICOM <- function(x, ...) {

  cat(cli::rule(), "\n")
  cat(cli::rule(center = "Overview", line = "bar3"), "\n\n",
      crayon::col_align("\tNumber of Observations", 25), "= ", x$Meta_information$Number_of_observations[1, 2], "\n", sep = "")
  for(i in 2:nrow(x$Meta_information$Number_of_observations)) {
    cat("\t\t", x$Meta_information$Number_of_observations[i, "x"], " : ",
    x$Meta_information$Number_of_observations[i, "n"], "\n")
  }
  cat(
      crayon::col_align("\tNumber of Groups", 25), "= ", x$Meta_information$Number_of_Groups, "\n",
      crayon::col_align("\tGrouping Variable", 25), "= ", x$Meta_information$Grouping_variable, "\n\n",
      sep = "")
  cat(cli::rule(center = "Details", line = "bar3"), "\n")
  cat(cli::rule(center = "Step 1 - Configural invariance", line = 2), "\n\n",
      "\tConfigural invariance is a precondition for step 2 and 3.\n",
      "\tDo not proceed to interpret results unless\n",
      "\tconfigural invariance has been established.\n\n",
      sep = "")
  cat(cli::rule(center = "Step 2 - Compositional invariance", line = 2), "\n\n",
      cli::boxx("H0: Compositional measurement invariance holds", float = "center"), "\n\n",
      sep = "")

  l <- max(nchar(c("Construct", rownames(x$Step2[[1]]))))

    for(i in seq_along(x$Step2)) {

    cat("Groups: ", names(x$Step2)[i], "\n\t",
        crayon::col_align("", width = l + 2, align = "center"),
        crayon::col_align("", 10, align = "center"), "\t",
        crayon::col_align("Critical Value(s)", 8*(ncol(x$Step2[[1]]) - 1), align = "center"), "\n\t",
        crayon::col_align("Construct", width = l + 2, align = "center"),
        crayon::col_align("c", 10, align = "center"), "\t",
        sep = "")
    for(j in colnames(x$Step2[[1]])[-1]) {
      cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
    }
    cat("\n\t")

    for(j in 1:nrow(x$Step2[[i]])) {
      cat(crayon::col_align(row.names(x$Step2[[i]])[j], l + 2), ": ",
          sprintf("%7.4f", x$Step2[[i]][j, "c"]) , "\t", sep = "")
      for(k in 2:ncol(x$Step2[[i]])) {
        cat(sprintf("%7.4f", x$Step2[[i]][j, k]), "\t",
            sep = "")
      }
      cat("\n\t")
    }
    cat("\n")
  }

  cat(cli::rule(center = "Step 3 - Equality of the mean values and variances", line = 2), "\n\n",
      cli::boxx(c("1. H0: Difference between group means is zero",
                  "2. H0: Log of the ratio of the group variances is zero"),
                float = "center"),
      sep = "")

  cat("\n\nEquality of the means:\n", "______________________", sep = "")
  for(i in seq_along(x$Step3$Mean_diff)) {

    cat("\n\nGroups: ", names(x$Step3$Mean_diff)[i], "\n\t",
        crayon::col_align("", width = l + 2, align = "center"),
        crayon::col_align("", 10, align = "center"), "\t",
        crayon::col_align("Critical Value(s)", 8*(ncol(x$Step3$Mean_diff[[1]]) - 1), align = "center"), "\n\t",
        crayon::col_align("Construct", width = l + 2, align = "center"),
        crayon::col_align("Mean diff.", 11, align = "center"), "\t",
        sep = "")

    for(j in colnames(x$Step3$Mean_diff[[1]][-1])) {
      cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
    }
    cat("\n\t")

    for(j in 1:nrow(x$Step3$Mean_diff[[i]])) {
      cat(crayon::col_align(row.names(x$Step3$Mean_diff[[i]])[j], l + 2), ": ",
          crayon::col_align(sprintf("%7.4f", x$Step3$Mean_diff[[i]][j, "Diff_mean"]), 11), sep = "")
      for(k in 2:ncol(x$Step3$Mean_diff[[i]])) {
        cat(sprintf("%7.4f", x$Step3$Mean_diff[[i]][j, k]), "\t",
            sep = "")
      }
      cat("\n\t")
    }
    cat("\n")
  }

  cat("\n\nEquality of the variances:\n", "__________________________", sep = "")
  for(i in seq_along(x$Step3$Var_diff)) {

    cat("\n\nGroups: ", names(x$Step3$Var_diff)[i], "\n\t",
        crayon::col_align("", width = l + 2, align = "center"),
        crayon::col_align("", 10, align = "center"), "\t",
        crayon::col_align("Critical Value(s)", 8*(ncol(x$Step3$Var_diff[[1]]) - 1), align = "center"), "\n\t",
        crayon::col_align("Construct", width = l + 2, align = "center"),
        crayon::col_align("Var diff.", 11, align = "center"), "\t",
        sep = "")

    for(j in colnames(x$Step3$Var_diff[[1]][-1])) {
      cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
    }
    cat("\n\t")

    for(j in 1:nrow(x$Step3$Var_diff[[i]])) {
      cat(crayon::col_align(row.names(x$Step3$Var_diff[[i]])[j], l + 2), ": ",
          crayon::col_align(sprintf("%7.4f", x$Step3$Var_diff[[i]][j, "Diff_log_var"]), 11), sep = "")
      for(k in 2:ncol(x$Step3$Var_diff[[i]])) {
        cat(sprintf("%7.4f", x$Step3$Var_diff[[i]][j, k]), "\t",
            sep = "")
      }
      cat("\n\t")
    }
    cat("\n")
  }
}

