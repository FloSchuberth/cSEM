#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' tidy Method for cSEMResults Objects
#'
#' This method and documentation is based on the `tidy` method for `lavaan`, as found in `broom:::tidy.lavaan()`.
#'
#' **Cautionary Note**: As with all `tidy()` methods, mis-specified arguments may be silently ignored.
#'
#' The attribute `verification` shows whether the model is admissible.
#'
#' @param x Fitted `csem` model with class `cSEMResults` as returned from [cSEM::csem()]
#' @param parameters Character string. Selects the desired parameter estimates as outputted from [cSEM::summarize()]. Defaults to `all`.
#' @param conf.int Boolean. Controls whether confidence intervals are returned or not. Defaults to `FALSE`.
#' @param conf.level Controls the size of the confidence interval returned. Passed to the `.alpha` arguments of [cSEM::summarize()] and [cSEM::infer()] as `1 - conf.level`.
#' @param conf.type Type of confidence interval returned. See  the `.ci` argument of [cSEM::summarize()]. Unlike [cSEM::summarize()], only one type of confidence interval at-a-time is supported.
#' @inheritDotParams summarize
#' @return A `data.frame` with one row for each estimated parameter. The
#'   description of the column names are as follows:
#'
#'   \item{term}{The result of `paste(lhs, op, rhs)`}
#'   \item{op}{The operator in the model syntax (e.g., `~~` for correlation, `<~` for weights, `=~` for loadings, and `~` for path estimates). Indicators of common factors that load onto themselves denote the unique loading of the common factor (e.g., cf_ind_1 =~ cf_ind_1).}
#'   \item{group}{The group as specified by `.id` in the [cSEM::csem()].}
#'   \item{estimate}{The parameter estimate}
#'   \item{construct.type}{Whether the modelled construct is a common factor or composite variable}
#'   \item{std.error}{}
#'   \item{statistic}{The t-statistic as returned from [cSEM::summarize()] and `cSEM:::addInfer()`}
#'   \item{p_value}{}
#'
#' @author Michael S. Truong
#'
#' @seealso [cSEM::summarize()]
#' @importFrom generics tidy
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#'model <- "
#'# Structural model
#'eta2 ~ eta1
#'eta3 ~ eta1 + eta2
#'
#'# (Reflective) measurement model
#'eta1 =~ y11 + y12 + y13
#'eta2 =~ y21 + y22 + y23
#'eta3 =~ y31 + y32 + y33
#'"
#'
#'# Single Group Example
#'res_boot <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 40)
#'
#'tidy(res_boot, conf.int = TRUE, conf.level = .95, conf.type = "CI_percentile")
#'tidy(res_boot, conf.int = TRUE, conf.level = .95, conf.type = NULL)
#'# Multi-Group Example
#'threecommonfactors_id <- cbind(
#'  "id" = sample(1:3, nrow(threecommonfactors), replace = TRUE),
#'  threecommonfactors
#')
#'
#'res_mg_boot <- csem(
#'  threecommonfactors_id,
#'  model,
#'  .resample_method = "bootstrap",
#'  .R = 40,
#'  .id = "id"
#')
#'
#'tidy(res_mg_boot, conf.int = TRUE)
#' }
#' @export
#'
tidy.cSEMResults <- function(
  x,
  parameters = c(
    "all",
    "Path_estimates",
    "Loadings_estimates",
    "Weight_estimates",
    "Unique_loading_estimates",
    "Effect_estimates",
    "Exo_construct_correlation"
  ),
  conf.int = FALSE,
  conf.level = 0.95,
  conf.type = NULL,
  ...
) {
  parameters <- match.arg(
    parameters,
    c(
      "all",
      "Path_estimates",
      "Loadings_estimates",
      "Weight_estimates",
      "Unique_loading_estimates",
      "Effect_estimates",
      "Exo_construct_correlation"
    )
  )

  match.arg(conf.type, args_default(.choices = TRUE)$.ci, several.ok = FALSE)

  # Create Table -----------------------------------------------------------

  if (inherits(x, "cSEMResults_default")) {
    # Single Group Summary
    summarized_cSEMResults_list <- list(summarize(
      x,
      .alpha = 1 - conf.level,
      .ci = conf.type,
      ...
    ))
  } else if (inherits(x, "cSEMResults_multi")) {
    # Multi Group Summary
    summarized_cSEMResults_list <- lapply(
      x,
      summarize,
      .alpha = 1 - conf.level,
      .ci = conf.type,
      ...
    )
  } else {
    stop("Unsupported object passed to tidy.cSEMResults()")
  }

  out <- lapply(
    summarized_cSEMResults_list,
    function(summarized_cSEMResults, parameters) {
      ## Selection ---------------------------------------------------------------
      if (parameters == "all") {
        if (
          is.null(summarized_cSEMResults$Estimates$Unique_loading_estimates)
        ) {
          general_estimates <- summarized_cSEMResults[["Estimates"]][c(
            "Path_estimates",
            "Loading_estimates",
            "Weight_estimates",
            "Exo_construct_correlation"
          )]
        } else {
          general_estimates <- summarized_cSEMResults[["Estimates"]][c(
            "Path_estimates",
            "Loading_estimates",
            "Weight_estimates",
            "Unique_loading_estimates",
            "Exo_construct_correlation"
          )]
        }

        effect_estimates <- summarized_cSEMResults[["Estimates"]][[
          "Effect_estimates"
        ]][c("Indirect_effect", "Total_effect")]

        # Add type of effect as column
        # See user1317221_G on April 26/2015 https://stackoverflow.com/a/29878732

        effect_estimates <- mapply(
          function(df, df_name) {
            # If there are no indirect effects, then the matrix will be NULL
            try(df$op <- rep(df_name, nrow(df)), silent = TRUE)
            return(df)
          },
          effect_estimates,
          names(effect_estimates),
          SIMPLIFY = FALSE
        )

        out <- c(general_estimates, effect_estimates)

        out <- out[!(unlist(lapply(out, is.null)))]
      } else if (parameters != "Effect_estimates") {
        out <- summarized_cSEMResults[["Estimates"]][parameters]
      } else if (parameters == "Effect_estimates") {
        effect_estimates <- summarized_cSEMResults[["Estimates"]][[
          "Effect_estimates"
        ]][c("Indirect_effect", "Total_effect")]

        # Add type of effect as column
        # See user1317221_G on April 26/2015 https://stackoverflow.com/a/29878732
        effect_estimates <- mapply(
          function(df, df_name) {
            try(df$op <- rep(df_name, nrow(df)), silent = TRUE)
            return(df)
          },
          effect_estimates,
          names(effect_estimates),
          SIMPLIFY = FALSE
        )

        out <- effect_estimates[!(unlist(lapply(effect_estimates, is.null)))]
      }

      # Take out the empty data.frames to facilitate the later concatenation
      out <- out[sapply(out, function(x) nrow(x) > 0)]
      # Add operand as op column
      out <- mapply(
        function(df, df_name) {
          if (
            (df_name %in% c("Indirect_effect", "Total_effect"))
          ) {
            return(df)
          } else {
            # Test for what operand is consistent for the iterated dataframe
            ## This needs spaces so that we don't conflate ~~ with ~
            operands <- c(" =~ ", " <~ ", " ~~ ", " ~ ")
            names(operands) <- c("=~", "<~", "~~", "~")
            df_type <- sapply(operands, function(x) grepl(x, df$Name[1]))

            if (length(names(df_type)[df_type]) == 1) {
              df$op <- rep(names(df_type)[df_type], nrow(df))
            } else {
              stop(
                "There is more than one type of allowed operand in the iterated data.frame. Cannot determine what to assign to op column."
              )
            }

            return(df)
          }
        },
        out,
        names(out),
        SIMPLIFY = FALSE
      )

      ## Conversion to broom Style -------------------------------------------------------------
      # Merge list of data.frames into one
      # See response by Charles on November 11/2011 https://stackoverflow.com/a/8097519
      out <- Reduce(function(...) merge(..., all = TRUE, sort = FALSE), out)

      # Rename columns to be consistent with tidy glossary standard
      # See Zon Shi Wu Jie from May 5/2014 https://stackoverflow.com/a/23475492

      if (isTRUE(conf.int)) {
        # We can't pass conf.type directly because summarize sets
        #  the conf.type to a non-NULL string when NULL is passed as an argument
        # Claude Sonnet 4.5 helped with string extraction
        out$conf.type <- sub(
          "\\.95%U.*$",
          "",
          grep(".95%U", colnames(out), value = TRUE)
        )
      }
      colnames(out)[colnames(out) == "Name"] <- 'term'
      colnames(out)[colnames(out) == "Construct_type"] <- 'type'
      colnames(out)[colnames(out) == "Estimate"] <- 'estimate'
      colnames(out)[colnames(out) == "Std_err"] <- 'std.error'
      colnames(out)[colnames(out) == "t_stat"] <- 'statistic'
      colnames(out)[colnames(out) == "p_value"] <- 'p.value'
      colnames(out)[grepl(".95%L", colnames(out))] <- 'conf.low'
      colnames(out)[grepl(".95%U", colnames(out))] <- 'conf.high'

      # Column reordering
      cols <- c(
        "term",
        "op",
        "type",
        "estimate",
        "std.error",
        "statistic",
        "p.value"
      )

      if (isTRUE(conf.int)) {
        cols <- c(cols, 'conf.low', 'conf.high', 'conf.type')
      }

      out <- out[, cols[cols %in% colnames(out)]]

      return(out)
    },
    parameters = parameters
  )

  if (length(out) == 1) {
    out <- out[[1]]
  } else if (length(out) > 1) {
    out <- mapply(
      FUN = function(summary, name) {
        summary$group <- name
        cols <- c(
          "term",
          "op",
          "type",
          "group",
          "estimate",
          "std.error",
          "statistic",
          "p.value"
        )
        if (isTRUE(conf.int)) {
          cols <- c(cols, 'conf.low', 'conf.high', 'conf.type')
        }
        summary <- summary[, cols[cols %in% colnames(summary)]]
        return(summary)
      },
      summary = out,
      name = names(out),
      SIMPLIFY = FALSE
    )

    out <- Reduce(function(...) merge(..., all = TRUE, sort = FALSE), out)
  }

  # Format output ----------------------------------------------------------
  out$type <- ifelse(
    out$op %in% c("~", "Indirect_effect", "Total_effect"),
    yes = NA,
    no = out$type
  )
  if (is.null(conf.type)) {
    conf.type <- attributes(summarized_cSEMResults_list[[1]])$.ci
  }
  out$conf.type <- conf.type
  out <- as_tibble(out)
  attr(out, "verification") <- verify(x)

  if (any(unlist(attr(out, "verification")))) {
    warning2(attr(x, "verification"))
  }

  return(out)
}

#' glance Method for cSEMResults Objects
#'
#'
#' This method and documentation is based on the `glance` method for `lavaan`, as found in `broom:::glance.lavaan()`.
#'
#' Note: Only fit statistics that can always be returned as part of a one-row tibble are included.
#'
#' @param x Fitted `csem` model with class `cSEMResults` as returned from [cSEM::csem()]
#' @inheritDotParams assess
#'
#' @return A one-row `data.frame` with columns:
#' @importFrom generics glance
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' 
#'
#' @seealso [cSEM::assess()]
#' @author Michael S. Truong
#' 
#' @examples 
#' \dontrun{
#'   model <- "
#' # Structural model
#' eta2 ~ eta1
#' eta3 ~ eta1 + eta2
#' # (Reflective) measurement model
#' eta1 =~ y11 + y12 + y13
#' eta2 =~ y21 + y22 + y23
#' eta3 =~ y31 + y32 + y33
#' "
#'   res <- csem(threecommonfactors, model)
#'   glance(res)
#'   threecommonfactors_id <- cbind(
#'     "id" = sample(1:3, nrow(threecommonfactors), replace = TRUE),
#'     threecommonfactors
#'   )
#'   res_mg <- csem(
#'     threecommonfactors_id,
#'     model
#'   )
#'   glance(res_mg)
#' }
#' 
#' @export
glance.cSEMResults <- function(
  x,
  ...
) {

  fits <- assess(
    x,
    .quality_criterion = c(
      'dg',
      'dl',
      'dml',
      'df',
      'chi_square',
      'chi_square_df',
      'cfi',
      'gfi',
      'cn',
      'ifi',
      'nfi',
      'nnfi',
      'rmsea',
      'rms_theta',
      'srmr',
      'FIT',
      'FIT_m',
      'FIT_s',
      'gof'
    )
  )

  if(inherits(x, "cSEMResults_default")) {
    fits <- list(fits)
  }

  tabulatable_fits <- c(
    'DG',
    'DL',
    'DML',
    'Df',
    'Chi_square',
    'Chi_square_df',
    'CFI',
    'GFI',
    'CN',
    'IFI',
    'NFI',
    'NNFI',
    'RMSEA',
    'RMS_theta',
    'SRMR',
    'FIT',
    'FIT_m',
    'FIT_s',
    'GoF'
  )

  # FIXME: Single line fit statistics for multigroup: https://github.com/FloSchuberth/cSEM/issues/587
  out <- lapply(fits, function(.fits, .tabulatable_fits = tabulatable_fits) {
    as_tibble(t(unlist(.fits[names(.fits)[
      names(.fits) %in% .tabulatable_fits
    ]])))
  })

  if (is.list(out) & (length(out) > 1)) {
    group_names <- names(out)
    out <- Reduce(rbind, out)
    out$group <- group_names
  } else {
    out <- out[[1]]
  }
    
  # Either rename here or rename throughout package to be compliant with tidymodels glossary https://www.tidymodels.org/learn/develop/broom/#glossaries
  names(out)[names(out) == "Chi_square_df"] <- "chi.square.df"
  names(out)[names(out) == "Chi_square"] <- "chi.squared"
  names(out)[names(out) == "RMS_theta"] <- "rms.theta"
  names(out)[names(out) == "FIT_m"] <- "fit.m"
  names(out)[names(out) == "FIT_s"] <- "fit.s"
  names(out) <- tolower(names(out))

  return(out)
}