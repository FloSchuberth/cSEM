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

  Reduce(rbind, rows)
}
