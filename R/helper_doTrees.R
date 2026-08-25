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


#' Pull the csem() argument list out of the object doTrees() was given
#'
#' The tree replays these arguments at every node, so the object has to be one
#' whose arguments describe a single fit. A multigroup `cSEMResults_multi` keeps
#' its `Information` per group and has no top-level `Arguments` at all -- and
#' grouping is what the tree is for, so it is refused rather than unwrapped.
#' @noRd
csem_tree_args <- function(.object) {
  if (!inherits(.object, "cSEMResults")) {
    stop2(
      "`.object` must be a cSEMResults object as returned by csem(), not ",
      class(.object)[1], "."
    )
  }
  if (inherits(.object, "cSEMResults_multi")) {
    stop2(
      "`.object` must be a single-group fit: doTrees() finds the groups. ",
      "Refit without `.id` and pass the pooled result."
    )
  }
  args <- .object$Information$Arguments
  if (is.null(args)) {
    stop2("`.object` carries no $Information$Arguments to refit nodes with.")
  }
  args
}


#' Check the tree can be grown from this fit
#'
#' One thing partykit and csem() would otherwise report from far away:
#' `calculateGSCAErrors()` returns `NA` (not an error) off a non-GSCA fit, so the
#' node statistic would fail deep inside the trafo. Every `.influence` value
#' reads it, so every fit `doTrees()` is given must be a GSCA one.
#' @noRd
validate_tree_input <- function(data, indicators, covariates, influence,
                                splitter, args) {
  if (!length(covariates)) {
    stop2("`.covariates` must name at least one column to partition on.")
  }
  missing_cov <- setdiff(covariates, names(data))
  if (length(missing_cov)) {
    stop2(
      "The following covariates are not columns of the data `.object` was ",
      "fitted on: ", paste(missing_cov, collapse = ", "),
      ". Include them in the csem() call -- non-indicator columns are ignored ",
      "when fitting."
    )
  }
  overlap <- intersect(covariates, indicators)
  if (length(overlap)) {
    stop2(
      "The following covariates are also indicators of the model: ",
      paste(overlap, collapse = ", "), "."
    )
  }
  ## Both surviving `.influence` values read calculateGSCAErrors(), which
  ## returns NA rather than erroring off a non-GSCA fit -- so without this the
  ## failure would surface deep inside the trafo as a node that would not fit.
  if (!identical(args$.approach_weights, "GSCA")) {
    stop2(
      "`.influence = \"", influence, "\"` needs a GSCA fit, but `.object` was ",
      "fitted with .approach_weights = \"", args$.approach_weights,
      "\". The casewise residual statistics are GSCA-specific."
    )
  }
  invisible(TRUE)
}


#' Tuning parameters for [doTrees()]
#'
#' `doTrees()` runs one algorithm --
#' conditional-inference variable selection from libcoin, followed by a
#' cutpoint rule chosen with `.splitter` -- so every parameter below applies to
#' every configuration, except `max_cuts`, which only the non-native cutpoint
#' rules consult.
#'
#' @param alpha Significance level a covariate must clear to split a node.
#' @param bonferroni Adjust the node's p-values for the number of covariates
#'   tested. The adjustment partykit applies is Šidák, \eqn{1 - (1 - p)^k}, which
#'   agrees with Bonferroni to first order. `k` is the number of covariates that
#'   returned a p-value at that node, not the number requested. The engine
#'   applies it to whatever the selector returns. Defaults to `TRUE`, matching
#'   [partykit::ctree_control()]: every node scans all covariates, so the
#'   unadjusted p-value is a multiple-comparison statistic reported as if it
#'   were a single test.
#'
#'   Note that "Šidák" is our name for it, not partykit's -- the word appears
#'   nowhere in that package's documentation, which calls the option
#'   "Bonferroni" throughout. `?partykit::ctree_control` documents only the
#'   choice ("`Bonferroni` and `Univariate` relate to p-values from the
#'   asymptotic distribution (adjusted or unadjusted)"), and
#'   `?partykit::ctree` only that the criterion is multiplicity adjusted. The
#'   exponential form is visible solely in the source of partykit's unexported
#'   `.extree_node()`, which multiplies the stored \eqn{\log(1 - p)} by `k` --
#'   i.e. \eqn{\log((1 - p)^k)}. The primary reference for the framework is
#'   Hothorn, Hornik and Zeileis (2006), listed under `?partykit::ctree`.
#' @param minbucket Minimum number of observations in a terminal node. With many
#'   indicators the default admits leaves the model may not fit; see the
#'   `n_fail_leaf` counter on the returned tree.
#' @param minsplit Minimum number of observations required to attempt a split.
#' @param minprob Minimum size of a terminal node as a proportion of its
#'   parent's. partykit raises `minbucket` per node rather than replacing it --
#'   `.extree_node()` computes \eqn{\max(\code{minbucket}, \lceil n_{node}
#'   \times \code{minprob} \rceil)} and hands the raised value to the variable
#'   selector and the cutpoint rule alike, so it constrains the non-native
#'   splitters through `candidate_partitions()` as well as partykit's own scan.
#'   Note that the proportion is of the *node*, not of the sample: at
#'   `minprob = 0.25` a depth-2 leaf can hold a sixteenth of the data. Defaults
#'   to `0.01`, matching [partykit::ctree_control()]; with the default
#'   `minbucket = 30L` it does not bind below n = 3000.
#' @param maxdepth Maximum depth of the tree, the root counting as 0.
#' @param max_cuts Maximum number of candidate cutpoints scanned per numeric
#'   covariate. Only the non-native cutpoint rules use this; partykit's own
#'   scan considers every distinct value.
#' @param R_test Number of permutations libcoin draws for the
#'   conditional-inference variable-selection test. Drawn only under
#'   `coin_distribution = "approximate"`; ignored entirely under
#'   `"asymptotic"`. Defaults to `9999L`, matching
#'   [partykit::ctree_control()]'s `nresample`.
#'
#'   This is the tree's only resampler, whatever `.splitter` is set to. The
#'   cutpoint search never permutes: partykit's own scan sets `nresample = 0L`
#'   unconditionally in the `SPLITONLY` branch of its `.ctree_test_internal()`,
#'   and the `"FIT"`, `"DLi"` and `"DGi"` kernels take a deterministic argmax
#'   over at most `max_cuts` candidate partitions. Making a whole run
#'   permutation-free is therefore `coin_distribution = "asymptotic"`, on any
#'   splitter, and not a property of the splitter itself.
#'
#'   The Monte Carlo estimate carries a standard error of
#'   \eqn{\sqrt{p (1 - p) / R}}. At the default and \eqn{p \approx \alpha =
#'   0.05} that is about 0.0022; at `R_test = 500L` it would be about 0.0097,
#'   a fifth of `alpha` and wide enough to flip a borderline split between
#'   runs.
#' @param coin_distribution How the conditional-inference family evaluates the
#'   null distribution -- not *which* null, which is the permutation null
#'   either way. Both settings test the same hypothesis by the same Strasser
#'   and Weber (1999) framework; they differ only in how its tail is obtained.
#'   "approximate" (the default) estimates the tail by drawing `R_test`
#'   permutations of the null distribution. "asymptotic" draws no permutations at all: libcoin
#'   computes
#'   the linear statistic's exact conditional expectation and covariance under
#'   the permutation null in closed form, standardises by them, and reads the
#'   p-value off the limiting chi-squared distribution -- the large-sample
#'   shape of that same permutation distribution.
#'
#'   So "asymptotic" is still a permutation test; it is the resampling that is
#'   approximated away, not the conditioning. `?libcoin::LinStatExpCov` states
#'   the split directly -- the code "computes the linear statistic, its
#'   expectation and covariance and, *optionally*, `nresample` samples from its
#'   permutation distribution" -- and Strasser and Weber's title, "The
#'   Asymptotic Theory of Permutation Statistics", is the whole point.
#'
#'   Prefer `"asymptotic"` when exact reproducibility matters more than the
#'   large-sample approximation: nothing is drawn, so its p-values do not
#'   depend on the seed. The default trades that for not relying on the
#'   chi-squared limit; see `R_test` for the Monte Carlo standard error it
#'   costs.
#'
#' @returns A named `list` of tuning parameters.
#' @seealso [doTrees()]
#' @export
igsca_tree_control <- function(alpha = 0.05,
                               bonferroni = TRUE,
                               minbucket = 30L,
                               minsplit = 2L * minbucket,
                               minprob = 0.01,
                               maxdepth = 3L,
                               max_cuts = 20L,
                               R_test = 9999L,
                               coin_distribution = c("approximate", "asymptotic")) {
  stopifnot(
    "`minprob` must be a single number in [0, 1]" =
      length(minprob) == 1 && is.numeric(minprob) && !is.na(minprob) &&
      minprob >= 0 && minprob <= 1,
    "`R_test` must be a single positive whole number" =
      length(R_test) == 1 && is.numeric(R_test) && !is.na(R_test) &&
      R_test >= 1 && R_test == round(R_test)
  )
  list(
    alpha = alpha,
    bonferroni = bonferroni,
    minbucket = as.integer(minbucket),
    minsplit = as.integer(minsplit),
    minprob = minprob,
    maxdepth = as.integer(maxdepth),
    max_cuts = as.integer(max_cuts),
    R_test = as.integer(R_test),
    coin_distribution = match.arg(coin_distribution)
  )
}


#' Casewise sum of squared GSCA residuals (n x 1) -- the "vec" influence.
#'
#' 
#'
#' @noRd
influence_vec <- function(E) {
  matrix(rowSums(E^2), ncol = 1L)
}



#' Casewise GSCA squared-residual matrix (N x (J + P), multivariate) -- COIN_mat.
#'
#' Elementwise square of matrix of N cases by (J indicators + P constructs) residuals.
#'
#' @noRd
influence_mat <- function(E) {
  E^2
}


#' Refit the tree's model on a subset of the data
#'
#' Replays the argument list of the [csem()] call [doTrees()] was given, with
#' `.data` swapped for the node's rows, so that every node in the tree is
#' estimated exactly the way its root was -- estimator, modes, convergence
#' criterion and all. `.args` is `.object$Information$Arguments`; its `.model`
#' entry is the parsed `cSEMModel`, which [csem()] accepts unchanged.
#'
#' `.id = NULL` yields the pooled single-group fit; `.id = "group"` the MGA fit;
#' assigning `NULL` drops the entry, so the pooled case falls back to [csem()]'s
#' own default.
#'
#' Do NOT reach for `modifyList()` here. It recurses whenever the old and new
#' values are both lists, and a `data.frame` is a list -- so it merges the node's
#' columns into the stored data instead of replacing it, and `[[<-.data.frame`
#' *recycles* rather than errors whenever the node size divides the original.
#' A node of n/2 then fits happily on a duplicated copy of itself.
#'
#' @param .data Data for this node.
#' @param .args Argument list of the original [csem()] call.
#' @param .id Grouping column, or `NULL` for a pooled fit.
#'
#' @returns A `cSEMResults` object.
#' @noRd
fit_csem <- function(.data, .args, .id = NULL) {
  ## Not paranoia: `$<-` coerces an atomic to a list, so passing a model string
  ## here (as every call site used to) yields a list whose unnamed first element
  ## matches csem()'s `.model` by position -- a fit that runs, converges, and
  ## silently uses csem()'s default estimator instead of the object's.
  if (!is.list(.args)) {
    stop2(
      "`.args` must be the argument list of a csem() call, not ",
      class(.args)[1], "."
    )
  }
  .args$.data <- .data
  .args$.id <- .id
  do.call(csem, .args)
}

#' Fit that reports failure instead of throwing
#'
#' @param .data Data for this node.
#' @param .args Argument list of the original [csem()] call.
#' @param .id Grouping column, or `NULL` for a pooled fit.
#'
#' @returns List of the fit (`NULL` on failure) and whether it converged.
#' @noRd
try_fit <- function(.data, .args, .id = NULL) {
  fit <- suppressWarnings(tryCatch(
    fit_csem(.data, .args, .id),
    error = function(e) NULL
  ))
  ok <- !is.null(fit) &&
    isTRUE(tryCatch(csem_converged(fit), error = function(e) FALSE)) 
  list(fit = fit, ok = ok)
}

#' Did every group of this fit pass verify()?
#' @noRd
csem_converged <- function(fit) {
  if (inherits(fit, "cSEMResults_multi")) {
    all(vapply(fit, function(x) sum(verify(x)) == 0, logical(1)))
  } else {
    sum(verify(fit)) == 0
  }
}



#' Best cutpoint for a covariate under one split kernel
#'
#' Scans `whichvar` in order and returns the first covariate's argmax candidate.
#' Counts scans that ran but produced nothing finite
#' into `n_fail_split`: that is the signature of a broken kernel, which the
#' `tryCatch` below would otherwise make indistinguishable from a node with no
#' admissible partition.
#' @noRd
argmax_split <- function(
  splitter,
  collector,
  model,
  mf,
  subset,
  whichvar,
  ctrl
) {
  scanned <- FALSE
  for (j in whichvar) {
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
    ## The kernel was actually evaluated on this covariate, so the invocation
    ## counts as a scan regardless of what came back.
    scanned <- TRUE
    if (!any(is.finite(stats))) {
      next
    }
    collector$n_split_scan <- collector$n_split_scan + 1L
    return(cands[[which.max(stats)]]$split)
  }
  ## Falling through after a scan means the kernel ran and never returned a
  ## finite statistic -- the signature of a broken kernel, which tryCatch above
  ## turns into NA and partykit would otherwise read as "no admissible split".
  ## A node with no admissible partition never sets `scanned`, so the two cases
  ## stay distinguishable in the collector.
  if (scanned) {
    collector$n_split_scan <- collector$n_split_scan + 1L
    collector$n_fail_split <- collector$n_fail_split + 1L
  }
  NULL
}

#' Warn when a requested split kernel never produced a usable statistic
#'
#' A kernel that throws is caught inside `argmax_split()` and turned into `NA`,
#' which partykit reads as "no admissible split" -- indistinguishable from a
#' genuinely unsplittable node unless the failed scans are counted. Warn when
#' every scan the kernel was asked for came back empty.
#' @noRd
warn_dead_splitter <- function(collector, splitter) {
  n <- collector$n_split_scan
  if (splitter == "native" || n == 0L || collector$n_fail_split < n) {
    return(invisible(NULL))
  }
  warning2(paste0(
    "The '", splitter, "' split kernel produced no usable statistic at any of ",
    "the ", n, " node(s) where it was scanned, so no split point could be ",
    "chosen from it and the tree was grown as if no split were admissible. ",
    "Check that the kernel runs standalone via partition_stat()."
  ))
  invisible(NULL)
}


#' Admissible binary partitions of one covariate
#'
#' Returns at most `max_cuts` candidates, each a `partysplit` plus the logical
#' vector saying which of the node's rows go left. Candidates that would leave
#' a child smaller than `minbucket` are dropped, so an empty list means the
#' covariate offers no admissible split -- not that a kernel failed.
#' @noRd
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
        goes_left = zs <= ct, # TODO: Revisit whether this should be <= or <. Was originally <.
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

#' @returns List
#'
new_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$root_seen <- FALSE      # first trafo call = root fit (full vs node fail)
  e$n_fail_full <- 0L
  e$n_fail_node <- 0L
  e$n_fail_candidate <- 0L
  # Split-kernel scans: how many argmax_split() calls actually evaluated the
  # kernel, and how many of those produced nothing usable. See argmax_split().
  e$n_split_scan <- 0L
  e$n_fail_split <- 0L
  # Leaf refits: failed IGSCA fits at terminal nodes. See attach_leaf_fits().
  e$n_fail_leaf <- 0L
  e
}




#' Fit the IGSCA model at every terminal node and attach it to the tree
#'
#' partykit only calls the trafo at nodes it *attempts to split*, so leaves come
#' back with `info = NULL` and every leaf-reading method sees nothing:
#' [coef.igsca_tree()] returned `NULL` for any tree with at least one split, and
#' [plot.igsca_tree()] drew unlabelled terminal boxes. More importantly the tree
#' carried no per-leaf IGSCA fit at all, which is the main thing a SEM tree is
#' for. Refitting once after the tree is grown fills both gaps.
#'
#' The fits are written into the terminal nodes' own `info`, so the methods can
#' read them through the ordinary [partykit::info_node()] route rather than a
#' side table. Printed output is unaffected: partykit's `formatinfo_node()`
#' falls back to `"*"` for list-valued info, verified to hold even when the
#' info carries a full `cSEMResults` object.
#'
#' Costs one IGSCA fit per leaf. Failures are counted into
#' `collector$n_fail_leaf` and leave `converged = FALSE` on that node.
#' @noRd
attach_leaf_fits <- function(tree, mf, args, indicators, collector) {
  ids <- partykit::nodeids(tree, terminal = TRUE)
  ## Use the party object's own fitted node ids: they are guaranteed to be in
  ## the same row order as `mf`, which re-deriving them would not be.
  leaf <- tree$fitted[["(fitted)"]]

  fits <- lapply(ids, function(id) {
    rows <- which(leaf == id)
    ft <- try_fit(mf[rows, indicators, drop = FALSE], args)
    if (!ft$ok) {
      collector$n_fail_leaf <- collector$n_fail_leaf + 1L
      return(list(
        nobs = length(rows),
        converged = FALSE,
        objfun = NA_real_,
        object = NULL
      ))
    }
    ## Sum of squared casewise GSCA residuals, matching the objfun the ctree
    ## trafo stores on inner nodes. calculateGSCAErrors() can throw on an
    ## otherwise converged fit, so the fit is kept either way.
    E <- tryCatch(calculateGSCAErrors(ft$fit), error = function(e) NULL)
    list(
      nobs = length(rows),
      converged = TRUE,
      objfun = if (is.null(E)) NA_real_ else sum(E^2),
      object = ft$fit
    )
  })
  names(fits) <- as.character(ids)

  ## as.list()/as.partynode() is partykit's own round-trip and is identity-
  ## preserving when the info is left alone, so only the leaves change.
  nd <- as.list(partykit::node_party(tree))
  for (i in seq_along(nd)) {
    key <- as.character(nd[[i]]$id)
    if (is.null(nd[[i]]$kids) && !is.null(fits[[key]])) {
      nd[[i]]$info <- fits[[key]]
    }
  }
  tree$node <- partykit::as.partynode(nd)
  tree
}

#' Strip the fitted cSEMResults objects from a grown tree's inner nodes
#'
#' Every inner node's info carries a full `cSEMResults` object, because
#' `saveinfo = TRUE` persists whatever the trafo returned and the trafo must
#' return it (see the note in [doTrees()]). A `maxdepth = 3` tree can therefore
#' retain up to seven of them. This drops those payloads after growing, leaving
#' the criteria, `nobs` and `objfun` intact. Terminal nodes are untouched, so
#' the leaf fits attached by `attach_leaf_fits()` survive.
#'
#' Not called by default -- see the commented-out line in [doTrees()].
#' @noRd
drop_inner_node_objects <- function(tree) {
  nd <- as.list(partykit::node_party(tree))
  for (i in seq_along(nd)) {
    if (!is.null(nd[[i]]$kids) && !is.null(nd[[i]]$info$object)) {
      nd[[i]]$info$object <- NULL
    }
  }
  tree$node <- partykit::as.partynode(nd)
  tree
}

#' Get root-node criteria of igsca tree
#'
#' Root-node criteria matrix as stored by extree/ctree (saveinfo = TRUE):
#' columns = tested covariates (all-NA columns dropped by extree), rows
#' "statistic" / "p.value" (both raw scale) / "criterion" (log(1 - p) with
#' extree's own statistic-rank tie-break already added). Under
#' `bonferroni = TRUE` the reported p-values are the Šidák-adjusted ones -- the
#' engine adjusts before it stores, so this is the quantity that was actually
#' compared against `alpha`, not the per-covariate p-value. NULL when the root
#' trafo failed (no test ran).
#'
#' @param tree A tree returned by [doTrees()].
#'
#' @returns The root node's criteria `matrix`, or `NULL` if no test ran.
#' @seealso [doTrees()]
#' @importFrom partykit info_node node_party
#' @export
root_criteria <- function(tree) {
  info <- partykit::info_node(partykit::node_party(tree))
  if (is.null(info)) {
    return(NULL)
  }
  info$criterion
}

#' Per-node fit summaries of an igsca tree
#'
#' Reports the number of observations and the objective function value of each
#' node's fit. Terminal nodes carry the fits [doTrees()] attaches after growing;
#' a node whose fit failed returns an `NA` row rather than being dropped, so the
#' rows always correspond to `node`.
#'
#' @param object A tree returned by [doTrees()].
#' @param node Node ids to report. Defaults to the terminal nodes.
#' @param drop If `TRUE` (the default), a single-node result is simplified from
#'   a one-row matrix to a named vector.
#' @param ... Ignored.
#'
#' @returns A `matrix` with columns `nobs` and `objfun`, one row per node, row
#'   names being the node ids; a named vector when one node is requested and
#'   `drop = TRUE`.
#' @seealso [doTrees()]
#' @export
#' @importFrom partykit nodeids info_node nodeapply
coef.igsca_tree <- function(object, node = NULL, drop = TRUE, ...) {
  if (is.null(node)) {
    node <- partykit::nodeids(object, terminal = TRUE)
  }
  cf <- do.call(
    "rbind",
    partykit::nodeapply(object, ids = node, FUN = function(n) {
      i <- partykit::info_node(n)
      ## A node whose fit failed, or one never visited by the trafo, has no
      ## info. Return an NA row rather than NULL: rbind() silently drops NULLs,
      ## which is how this method used to return NULL for a whole tree.
      c(
        nobs = if (is.null(i$nobs)) NA_real_ else as.numeric(i$nobs),
        objfun = if (is.null(i$objfun)) NA_real_ else as.numeric(i$objfun)
      )
    })
  )
  rownames(cf) <- as.character(node)
  if (drop) {
    drop(cf)
  } else {
    cf
  }
}

#' Plot an igsca tree
#'
#' Draws the tree with a terminal panel reporting each leaf's `nobs` and
#' `objfun`. A leaf whose fit failed prints `?` rather than an empty box.
#'
#' @param x A tree returned by [doTrees()].
#' @param terminal_panel Panel function for the terminal nodes. Defaults to the
#'   `nobs`/`objfun` panel described above.
#' @param FUN Formatting function applied to each leaf's info before printing.
#' @param tp_args Arguments passed on to `terminal_panel`.
#' @param ... Passed to [partykit::plot.party()].
#'
#' @returns Called for the plot it draws. The return value is
#'   [partykit::plot.party()]'s and is not meaningful.
#' @seealso [doTrees()]
#' @export
plot.igsca_tree <- function(x, terminal_panel = NULL, FUN = NULL, tp_args = NULL, ...) {
  if (is.null(terminal_panel)) {
    if (is.null(FUN)) {
      FUN <- function(info) {
        ## info is NULL for a node the trafo never reached and objfun is NA
        ## when the leaf refit failed; sprintf() on NULL yields character(0),
        ## which silently draws an empty box, so both are shown explicitly.
        c(
          sprintf("n = %s", if (is.null(info$nobs)) "?" else info$nobs),
          sprintf(
            "SSR = %s",
            if (is.null(info$objfun)) "?" else format(info$objfun, digits = 4L)
          )
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



# Split kernels for doTrees() ----------------------------------------------
#
# doTrees() has one selector -- partykit's own COIN variable selection, driven
# by the influence_mat/influence_vec statistic -- and four cutpoint rules. This
# section holds the three that are not partykit's own:
#
#   split_max_fitdiff / split_max_dli / split_max_dgi
#
# Splitter contract (a statistic kernel argmaxed by argmax_split(), which
# doTrees() installs as the engine's splitfun):
#   split_*(model, mf, subset, goes_left, ctrl) -> scalar observed statistic
#
# Each is evaluated at every admissible candidate partition of the selected
# covariate (one two-group MGA fit each; the node's pooled fit comes free from
# the trafo as model$object) and the argmax wins. No permutation is involved:
# see the note under `R_test` in igsca_tree_control().



#'Node data for one candidate partition: children as `group` levels 1/2.
#' 
#'
#' @noRd
node_group_data <- function(mf, subset, indicators, goes_left) {
  d <- mf[subset, indicators, drop = FALSE]
  cbind(group = factor(ifelse(goes_left, 1L, 2L), levels = c(1L, 2L)), d)
}



#' Pooled model-implied VCVs replicated to 2 blocks (constant per node;
#' cached in the collector under $ndt_pools keyed by the subset).
#' 
#'
#' @noRd
ndt_pools <- function(single_fit) {
  list(
    Sc = bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "construct"),
    Si = bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "indicator")
  )
}

#' 
#' Observed statistic for one candidate partition. `stat_kind` is "FITdiff"
#' (NPT) or an ndt_dists() distance name -- "DGi"/"DLi" only: Study 1 works
#' with the model-implied indicator VCV, never the construct VCV (the shared
#' helper also computes "DGc"/"DLc", which stay unused here). Counts every
#' candidate partition whose statistic could not be computed into
#' collector$n_fail_candidate.
#' 
#'
#' @noRd
partition_stat <- function(stat_kind, model, mf, subset, goes_left, ctrl) {
  coll <- ctrl$collector
  d <- node_group_data(mf, subset, ctrl$indicators, goes_left)
  mga <- try_fit(d, ctrl$args, .id = "group")
  if (!mga$ok) { coll$n_fail_candidate <- coll$n_fail_candidate + 1L; return(NA_real_) }
  if (stat_kind == "FITdiff") {
    if (is.null(model$object)) {
      return(NA_real_)
    }
    val <- tryCatch(
      calculateFIT(mga$fit) - calculateFIT(model$object),
      error = function(e) NULL
    )
    if (is.null(val)) {
      coll$n_fail_candidate <- coll$n_fail_candidate + 1L
      return(NA_real_)
    }
    return(val)
  }
  ## NDT distances need the node's pooled 2-block VCVs (computed once/node).
  ## cSEM internals (bdiagFit/calculateDG/calculateDL) can throw; guard so no
  ## exception escapes the selector (would make SimDesign redraw the rep).
  ## Invariant: within one tree-growing call, `subset` uniquely identifies a
  ## node and model$object (the pooled fit) is invariant per node, so this
  ## single-slot cache keyed on subset alone cannot collide across nodes.
  ds <- tryCatch(
    {
      if (
        is.null(coll$ndt_pools_subset) ||
          !identical(coll$ndt_pools_subset, subset)
      ) {
        if (is.null(model$object)) {
          return(NA_real_)
        }
        coll$ndt_pools <- ndt_pools(model$object)
        coll$ndt_pools_subset <- subset
      }
      ndt_dists(coll$ndt_pools$Sc, coll$ndt_pools$Si, mga$fit)
    },
    error = function(e) NULL
  )
  if (is.null(ds)) { coll$n_fail_candidate <- coll$n_fail_candidate + 1L; return(NA_real_) }
  unname(ds[stat_kind])
}



## The three split kernels: deliberately repetitive plain functions (no
## factories). Each is the bare observed-statistic kernel for one candidate
## partition, argmaxed by argmax_split() when doTrees() is given a .splitter
## other than "native". Distances are model-implied INDICATOR-VCV only
## (dGi/dLi) -- Study 1 never works with the construct VCV.


#' Split kernel: FIT difference at one candidate partition.
#' @noRd
split_max_fitdiff <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("FITdiff", model, mf, subset, goes_left, ctrl)
}

#' Split kernel: geodesic indicator-VCV distance at one candidate partition.
#' @noRd
split_max_dgi <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DGi", model, mf, subset, goes_left, ctrl)
}

#' Split kernel: squared-Euclidean indicator-VCV distance at one candidate partition.
#' @noRd
split_max_dli <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DLi", model, mf, subset, goes_left, ctrl)
}




#' The four pooled-vs-MGA distances for one MGA fit, given the precomputed
#' replicated-block pooled VCVs.
#' 
#'
#' @noRd
ndt_dists <- function(Sc_pool, Si_pool, mga_fit) {
  Sc_mga <- bdiagFit(mga_fit, .type_vcv = "construct")
  Si_mga <- bdiagFit(mga_fit, .type_vcv = "indicator")
  c(
    DGc = calculateDG(.matrix1 = Sc_pool, .matrix2 = Sc_mga),  # geodesic,  construct
    DGi = calculateDG(.matrix1 = Si_pool, .matrix2 = Si_mga),  # geodesic,  indicator
    DLc = calculateDL(.matrix1 = Sc_pool, .matrix2 = Sc_mga),  # sq-Euclid, construct
    DLi = calculateDL(.matrix1 = Si_pool, .matrix2 = Si_mga)   # sq-Euclid, indicator
  )
}