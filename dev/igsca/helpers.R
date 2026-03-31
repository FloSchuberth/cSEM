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
        val <- if (parsed$count == "uni") setNames(1, uni_name) else tri_weights[[xi]]
        cfg$weights <- c(cfg$weights, setNames(list(val), xi))
      } else {
        val <- if (parsed$count == "uni") setNames(1, uni_name) else tri_loadings[[xi]]
        cfg$loadings <- c(cfg$loadings, setNames(list(val), xi))
      }
    }

    return(cfg)
  })

  configs |>
    lapply(\(cfg) do.call(build_expected_popvalues, cfg)) |>
    (\(dfs) Reduce(rbind, dfs))()
}
