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


#' Title
#'
#' @param alpha
#' @param bonferroni
#' @param minbucket
#' @param minsplit
#' @param maxdepth
#' @param max_cuts
#' @param R_test
#' @param coin_distribution
#' @param influence
#' @param splitter
#'
#' @returns
#'
#' @export
#' @examples
igsca_tree_control <- function(
  alpha = 0.05,
  bonferroni = TRUE,
  minbucket = 30L,
  minsplit = 2L * minbucket,
  maxdepth = 3L,
  max_cuts = 20L,
  R_test = 500L, # TODO: Consider reusing ctree_control()$nresample
  coin_distribution = c("approximate", "asymptotic"),
  influence = influence_vec,
  splitter = NULL
) {
  # TODO: consider doing something like utils::modifyList(ctree_control, list(...)) instead to reduce the number of arguments
  list(
    alpha = alpha,
    bonferroni = bonferroni,
    minbucket = as.integer(minbucket),
    minsplit = as.integer(minsplit),
    maxdepth = as.integer(maxdepth),
    max_cuts = as.integer(max_cuts),
    R_test = as.integer(R_test),
    coin_distribution = match.arg(coin_distribution),
    influence = influence,
    splitter = splitter
  )
}

#' Casewise sum of squared GSCA residuals (n x 1) -- COIN_ssr influence.
influence_vec <- function(E) {
  matrix(rowSums(E^2), ncol = 1L)
}


#' Casewise squared-residual matrix (n x q, multivariate) -- COIN_mat.
influence_mat <- function(E) {
  E^2
}

#' `.id = NULL` yields the pooled single-group fit; `.id = "group"` the MGA fit.
fit_csem <- function(.data, model, .id = NULL) {
  csem(
    .data = .data,
    model = model,
    .id = .id,
    .approach_weights = "GSCA",
    .disattenuate = TRUE, # to get igsca
    .conv_criterion = "sum_diff_absolute", # Default in gsca_m.m and gsca.m
    .iter_max = 100, # Default in igsca_sim.m
    .GSCA_modes = "CCMP", # Unknown (to me) if it affects convergence
    .tolerance = 0.0001 # Default in gsca_m.m and gsca.m
  )
}

try_fit <- function(.data, model, .id = NULL) {
  fit <- suppressWarnings(tryCatch(
    fit_csem(.data, model, .id),
    error = function(e) NULL
  ))
  ok <- !is.null(fit) &&
    isTRUE(tryCatch(sum(verify(fit)), error = function(e) FALSE))
  list(fit = fit, ok = ok)
}



argmax_split <- function(
  splitter,
  collector,
  model,
  mf,
  subset,
  whichvar,
  ctrl
) {
  for (j in whichvar) {
    hit <- identical(collector$scan_subset, subset) &&
      identical(collector$scan_splitter, splitter) &&
      !is.null(collector$scan[[as.character(j)]])
    if (hit) {
      return(collector$scan[[as.character(j)]]$split)
    }
    cands <- candidate_partitions(
      j,
      mf[[j]],
      mf[[j]][subset],
      ctrl$max_cuts,
      ctrl$minbucket
    )
    if (!length(cands)) {
      next
    }
    stats <- vapply(
      cands,
      function(cc) {
        tryCatch(
          splitter(model, mf, subset, cc$goes_left, ctrl),
          error = function(e) NA_real_
        )
      },
      numeric(1)
    )
    if (!any(is.finite(stats))) {
      next
    }
    return(cands[[which.max(stats)]]$split)
  }
  NULL
}


candidate_partitions <- function(j, z, zs, max_cuts, minbucket) {
  keep_min <- function(cands) {
    Filter(
      function(cc) {
        nl <- sum(cc$goes_left)
        nl >= minbucket && (length(zs) - nl) >= minbucket
      },
      cands
    )
  }
  if (is.numeric(z) && !is.factor(z)) {
    uz <- sort(unique(zs))
    if (length(uz) < 2L) {
      return(list())
    }
    mids <- (uz[-1L] + uz[-length(uz)]) / 2
    if (length(mids) > max_cuts) {
      mids <- unique(stats::quantile(
        zs,
        probs = seq_len(max_cuts) / (max_cuts + 1),
        type = 7,
        names = FALSE
      ))
    }
    keep_min(lapply(mids, function(ct) {
      list(
        goes_left = zs < ct,
        split = partykit::partysplit(
          as.integer(j),
          breaks = as.double(ct),
          index = 1L:2L
        )
      )
    }))
  } else if (is.ordered(z)) {
    levs <- levels(droplevels(zs))
    K <- length(levs)
    if (K < 2L) {
      return(list())
    }
    keep_min(lapply(1L:(K - 1L), function(i) {
      list(
        goes_left = zs %in% levs[1L:i],
        split = partykit::partysplit(
          as.integer(j),
          breaks = as.integer(match(levs[i], levels(z)))
        )
      )
    }))
  } else if (is.factor(z)) {
    levs <- levels(droplevels(zs))
    K <- length(levs)
    if (K < 2L) {
      return(list())
    }
    olev <- levels(z)
    keep_min(lapply(1L:(2L^(K - 1L) - 1L), function(m) {
      g <- levs[as.logical(intToBits(m))[1L:K]]
      idx <- rep(NA_integer_, length(olev))
      idx[match(levs, olev)] <- ifelse(levs %in% g, 1L, 2L)
      list(
        goes_left = zs %in% g,
        split = partykit::partysplit(as.integer(j), index = idx)
      )
    }))
  } else list()
}

new_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$root_seen <- FALSE      # first trafo call = root fit (full vs node fail)
  e$n_fail_full <- 0L
  e$n_fail_node <- 0L
  e$n_fail_resample <- 0L
  e$scan <- list()
  e$scan_subset <- NULL
  e$scan_splitter <- NULL
  e
}


#' Root-node criteria matrix as stored by extree/ctree (saveinfo = TRUE):
#' columns = tested covariates (all-NA columns dropped by extree), rows
#' "statistic" / "p.value" (both raw scale) / "criterion" (log(1 - p) with
#' extree's own statistic-rank tie-break already added). NULL when the root
#' trafo failed (no test ran).
root_criteria <- function(tree) {
  info <- partykit::info_node(partykit::node_party(tree))
  if (is.null(info)) {
    return(NULL)
  }
  info$criterion
}

## Methods: print falls through to constparty; coef/plot follow
coef.igsca_tree <- function(object, node = NULL, drop = TRUE, ...) {
  if (is.null(node)) {
    node <- partykit::nodeids(object, terminal = TRUE)
  }
  cf <- do.call(
    "rbind",
    partykit::nodeapply(object, ids = node, FUN = function(n) {
      i <- partykit::info_node(n)
      c(nobs = i$nobs, objfun = i$objfun)
    })
  )
  if (drop) {
    drop(cf)
  } else {
    cf
  }
}

plot.igsca_tree <- function(x, terminal_panel = NULL, FUN = NULL, tp_args = NULL, ...) {
  if (is.null(terminal_panel)) {
    if (is.null(FUN)) {
      FUN <- function(info) {
        c(
          sprintf("n = %s", info$nobs),
          sprintf("SSR = %s", format(info$objfun, digits = 4L))
        )
      }
    }
    terminal_panel <- do.call(
      partykit::node_terminal,
      c(list(obj = x, FUN = FUN), tp_args)
    )
    tp_args <- NULL
  }
  partykit::plot.party(x, terminal_panel = terminal_panel, tp_args = tp_args, ...)
}
