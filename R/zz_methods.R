## Methods

#' @export
print.cSEMResults <- function(x, ...) {

  cat(cli::rule(line = "bar2"), "\n")
  cat(cli::rule(center = "Overview"), "\n\n")
  cat("Estimation was successful. The result is a named list of class " %+% bold("cSEMResults") %+%" \n",
      "containing the following list elements (:\n\n\t",
      "- ", crayon::green("Estimates\n\t"),
      "- ", crayon::green("Tests\n\t"),
      "- ", crayon::green("Meta_information\n\n"), sep = "")
  cat("To get an overview or help type: \n\n\t",
      "- ", crayon::magenta("?"), crayon::cyan("cSEMResults"),"\n\t",
      "- ", crayon::magenta("str"), "(", crayon::cyan("<listname>"), ")\n\t",
      "- ", crayon::magenta("summary"), "(", crayon::cyan("<listname>"), ") or \n\t",
      "- ", crayon::magenta("listviewer"), crayon::yellow("::"), crayon::magenta("jsondedit"),
      "(", crayon::cyan("<listname>"), ", ", crayon::red("mode"), " = ", crayon::cyan("'view'"), ")\n\n", sep = "")
  cat("If you wish to access the list elements directly type e.g. \n\t",
      "- ", crayon::cyan("<listname>"), crayon::yellow("$"), crayon::green("Estimates"), "\n", sep = "")
  cat(cli::rule(line = "bar2"), "\n")
}

#' @export
print.cSEMResultssummary <- function(x, ...) {
  cat(cli::rule(line = "bar2"), "\n",
      cli::rule(center = "General Information"), "\n\n\t", sep = "")

  cat(crayon::col_align("Number of Observations", 25), "= ", x$Number_of_observations, "\n\t",
      crayon::col_align("Weight estimator", 25), "= ", x$Weight_estimator, "\n\t", sep = "")
  if(x$Weight_estimator == "PLS") {
  cat(crayon::col_align("Inner Weighting Scheme ", 25), "= ", x$PLS_weight_scheme_inner, "\n\t", sep = "")
  }
  cat(
      crayon::col_align("Path estimator", 25), "= ", x$Path_estimator, "\n\t",
      crayon::col_align("Convergence Status", 25), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("Overall model fit", 25), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("Degrees of Freedom", 25), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("Computation Time", 25), "= ", c("not yet implemented"), "\n\n\t",
      sep = "")

  cat("Construct Types:\n\t","----------------","\n\t", sep = "")

  for(i in seq_along(x$Construct_types$Name)) {
    cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$Construct_types$Type[i],"\n\t", sep = "")
  }
  cat("\n")

  if(x$Weight_estimator == "PLS") {
    cat("\tPLS Modes:\n\t","----------------","\n\t", sep = "")

    for(i in seq_along(x$Construct_type$Name)) {
      cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$PLS_modes[i],"\n\t", sep = "")
    }
    cat("\n")
  }

  cat(cli::rule(center = "Estimates"), "\n\n", sep = "")

  cat("Estimated Path Coefficients:\n============================\n", sep = "")
  print(x$Path_estimates, row.names = FALSE)


  cat("\nEstimated Loadings:\n===================\n", sep = "")
  print(x$Loading_estimates, row.names = FALSE)

  cat("\nEstimated Weights:\n==================\n", sep = "")
  print(x$Weight_estimates, row.names = FALSE)

  if(x$Weight_estimator == "PLS") {
    cat("\nEstimated Correction Factors:\n=============================\n", sep = "")
    print(x$Correction_factors)
  }

  cat("\n\n", cli::rule(center = "Other output"), "\n\n\t", sep = "")

  cat("<not yet implemented>")

  cat("\n\n", cli::rule(center = "Fit Indices"), "\n\n\t", sep = "")

  cat(crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\n",
      sep = "")

  cat(cli::rule(line = "bar2"))
}

#' @export
summary.cSEMResults <- function(x, ...) {

  ## Structure loadings output
  temp <- x$Estimates$Loading_estimates
  names_loadings <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                           " =~ ", colnames(temp))
  loading_estimates <- data.frame("Loading" = names_loadings,
                                  "Estimate" = unlist(t(temp)[t(temp) != 0 ]))

  ## Structure weights output
  temp <- x$Estimates$Weight_estimates
  names_weights <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                          " -- ", colnames(temp))
  weight_estimates <- data.frame("Weights" = names_weights,
                                 "Estimate" = unlist(t(temp)[t(temp) != 0 ]))

  ## Create summary list
  summary_out <- list(
    "Construct_types"        = x$Meta_information$Construct_types,
    "Correction_factors"     = x$Estimates$Correction_factors,
    "Loading_estimates"      = loading_estimates,
    "Names_endogenous_var"   = x$Meta_information$modellist$vars_endo,
    "Number_of_observations" = x$Meta_information$Number_of_observations,
    "Path_estimates"         = x$Estimates$Path_estimates,
    "Path_estimator"         = x$Meta_information$Path_approach,
    "PLS_modes"              = x$Meta_information$PLS_Modes,
    "PLS_weight_scheme_inner"= x$Meta_information$PLS_Inner_Weightning_scheme,
    "Weight_estimates"       = weight_estimates,
    "Weight_estimator"       = x$Meta_information$Weight_approach
  )

  ## Set class for printing and return
  class(summary_out) <- "cSEMResultssummary"
  return(summary_out)
}

