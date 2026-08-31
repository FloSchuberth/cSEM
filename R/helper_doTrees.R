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
  ## partykit scans numeric, ordered and unordered covariates and nothing else.
  ## A logical dies as "object 'X' not found" inside .ctree_test_1d() and a
  ## character as "cannot handle objects of class 'character'" inside inum(),
  ## both far from this call. A non-native kernel is quieter and worse:
  ## candidate_partitions() offers no candidate and the covariate goes unsplit
  ## with no diagnostic at all.
  bad <- covariates[!vapply(
    data[covariates],
    function(z) is.numeric(z) || is.factor(z),
    logical(1)
  )]
  if (length(bad)) {
    stop2(
      "`.covariates` must name numeric, ordered or unordered factor columns. ",
      "The following are neither: ", paste(bad, collapse = ", "),
      ". Convert them first -- a logical or character covariate is a factor."
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
#' every configuration.
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
#' @param R_test Number of permutations libcoin draws for the
#'   conditional-inference variable-selection test. Drawn only under
#'   `coin_distribution = "approximate"`; ignored entirely under
#'   `"asymptotic"`.
#'
#' @param coin_distribution How the conditional-inference family evaluates the
#'   null distribution. See `?partykit::ctree_control()` and `?libcoin::LinStatExpCov` for more information.
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
#'
#' This stands in for `partykit:::.split()`, which `doTrees()` bypasses by
#' replacing `ctrl$splitfun` outright rather than plugging a kernel into
#' `.split()`'s own `FUN` slot -- the clean route, but `.split()` is internal.
#' `.split()` does two things around the cutpoint search, and both are
#' reproduced here: it takes unordered factors multiway without consulting the
#' cutpoint rule at all, and under `lookahead` it abandons a covariate whose
#' chosen split leaves a kid the model cannot be refitted in. Neither setting is
#' exposed by [igsca_tree_control()], so both arrive only through `ctrl`.
#' @noRd
argmax_split <- function(
  splitter,
  collector,
  model,
  trafo,
  mf,
  subset,
  whichvar,
  ctrl,
  weights = integer(0)
) {
  scanned <- FALSE
  got_stat <- FALSE
  ret <- NULL
  for (j in whichvar) {
    z <- mf[[j]]
    if (
      isTRUE(ctrl$multiway) && is.factor(z) && !is.ordered(z) &&
        identical(as.numeric(ctrl$maxsurrogate), 0) &&
        nlevels(droplevels(z[subset])) > 1L
    ) {
      ## Structural, not a statistic: partykit picks one kid per level and
      ## never calls the cutpoint rule, so no kernel can override it.
      ret <- multiway_split(j, z, z[subset], ctrl$minbucket)
    } else {
      ## partykit tests a covariate on its non-missing rows alone (`subsetNArm`
      ## in partykit:::.ctree_test()) and routes the missing ones afterwards
      ## through the chosen split's `prob`. Left in, `zs <= ct` is NA, so
      ## `sum(goes_left)` is NA, Filter() drops every candidate, and the
      ## covariate goes unsplit without even counting as a scan.
      sub_j <- subset[!is.na(z[subset])]
      cands <- candidate_partitions(
        j,
        z,
        z[sub_j],
        ctrl$minbucket,
        isTRUE(ctrl$intersplit)
      )
      if (!length(cands)) {
        next
      }
      stats <- vapply(
        cands, # tryCatch is within vapply, so we take the best split of those associated with convergence success.
        function(cc) {
          tryCatch(
            splitter(model, mf, sub_j, cc$goes_left, ctrl),
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
      got_stat <- TRUE
      ret <- cands[[which.max(stats)]]$split
    }
    if (is.null(ret)) {
      next
    }
    ## partykit moves on to the next covariate rather than to the next
    ## cutpoint, so a rejected split takes its whole covariate with it.
    if (isTRUE(ctrl$lookahead) && !kids_converge(ret, trafo, mf, subset, weights)) {
      ret <- NULL
      next
    }
    break
  }
  ## Falling through after a scan that never produced a finite statistic is the
  ## signature of a broken kernel, which tryCatch above turns into NA and
  ## partykit would otherwise read as "no admissible split". A node with no
  ## admissible partition never sets `scanned`, so the two cases stay
  ## distinguishable in the collector -- and a statistic that lookahead then
  ## discarded sets `got_stat`, so the kernel is not blamed for it either.
  if (scanned) {
    collector$n_split_scan <- collector$n_split_scan + 1L
    if (!got_stat) {
      collector$n_fail_split <- collector$n_fail_split + 1L
    }
  }
  ret
}

#' One kid per level of an unordered factor
#'
#' partykit's multiway rule, transcribed from `partykit:::.split()`: levels
#' absent from the node get no kid (`NA`), levels present but smaller than
#' `minbucket` are lumped together into one extra kid rather than dropped, and
#' the resulting ids are re-coded to be consecutive. `NULL` when that leaves a
#' single kid, which is no split at all.
#' @noRd
multiway_split <- function(j, z, zs, minbucket) {
  K <- nlevels(z)
  index <- seq_len(K)
  xt <- tabulate(zs, nbins = K)
  index[xt == 0L] <- NA_integer_
  index[xt > 0L & xt < minbucket] <- K + 1L
  if (length(unique(index)) == 1L) {
    return(NULL)
  }
  partykit::partysplit(
    as.integer(j),
    index = as.integer(unclass(factor(index)))
  )
}

#' Does the node model refit in every kid this split would create?
#'
#' partykit's `lookahead` check, transcribed from `partykit:::.split()`. The
#' trafo is the same one that fits the node, so this costs one IGSCA fit per
#' kid on top of the scan that chose the split -- and a kid that fails to fit
#' counts into `collector$n_fail_node` like any other node fit, which is what
#' partykit's own lookahead would do with this trafo.
#' @noRd
kids_converge <- function(split, trafo, mf, subset, weights) {
  sp <- partykit::kidids_split(split, mf, obs = subset)
  all(vapply(
    unique(sp[!is.na(sp)]),
    function(i) {
      isTRUE(trafo(subset[!is.na(sp) & sp == i], weights = weights)$converged)
    },
    logical(1)
  ))
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
#' Returns every candidate, each a `partysplit` plus the logical vector saying
#' which of the node's rows go left. Candidates that would leave a child
#' smaller than `minbucket` are dropped, so an empty list means the covariate
#' offers no admissible split -- not that a kernel failed. `zs` must already
#' have the covariate's missing rows removed; see `argmax_split()`.
#'
#' The scan is exhaustive: every cut between consecutive distinct values of a
#' numeric covariate, every one of the `K - 1` cuts of an ordered factor, and
#' all `2^(K - 1) - 1` bipartitions of an unordered one. That is the same set
#' of partitions partykit's own scan considers, and the emitted `partysplit`
#' follows the same conventions (`z <= breaks` goes left; an ordered break is
#' an index into `levels(z)`; an unordered `index` is 1/2 per level of `z`,
#' `NA` for levels absent from the node).
#'
#' The costs are not the same, though: the non-native splitters pay a two-group
#' IGSCA fit per candidate where partykit pays C arithmetic, which is why the
#' unordered scan is capped far below the ">= 31 levels" at which libcoin
#' itself refuses. Whether the numeric grid should be binned first, and how, is
#' left open -- partykit would do it through `nmax`, which `doTrees()` pins at
#' `Inf`.
#'
#' @param j Column index of the covariate in the model frame.
#' @param z The covariate over all rows, which fixes the factor levels and
#'   level order the emitted `partysplit` is expressed in.
#' @param zs The covariate over this node's non-missing rows.
#' @param minbucket Smallest child either kid may be left with.
#' @param intersplit Report a numeric break as the midpoint of the gap rather
#'   than the observed value below it. `FALSE` matches ctree's default.
#' @noRd
candidate_partitions <- function(j, z, zs, minbucket, intersplit = FALSE) {
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
    ## `<=` is partykit's convention, not a choice: partysplit() defaults to
    ## right = TRUE, so kidids_split() bins on (-Inf, breaks] and sends
    ## `z <= breaks` left.
    cuts <- uz[-length(uz)]
    ## Which rows go left is fixed by the gap; only the number reported as the
    ## break differs, and it matters for rows that land inside the gap later.
    ## ctree reports the observed value below the gap and switches to the
    ## midpoint under intersplit = TRUE -- see the tail of
    ## partykit:::.ctree_test_internal(). Taking the midpoint unconditionally
    ## also risks rounding up onto uz[i + 1] when the two are adjacent doubles,
    ## which duplicates the next candidate and drops this one.
    brks <- if (intersplit) (uz[-1L] + cuts) / 2 else cuts
    keep_min(Map(
      function(ct, br) {
        list(
          goes_left = zs <= ct,
          split = partykit::partysplit(
            as.integer(j),
            breaks = as.double(br),
            index = 1L:2L
          )
        )
      },
      cuts,
      brks
    ))
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
          breaks = as.integer(match(levs[i], levels(z))),
          index = 1L:2L
        )
      )
    }))
  } else if (is.factor(z)) {
    levs <- levels(droplevels(zs))
    K <- length(levs)
    if (K < 2L) {
      return(list())
    }
    ## Every bipartition costs a two-group IGSCA fit here where partykit pays
    ## C arithmetic, so this scan has to stop well short of the ">= 31 levels"
    ## at which libcoin itself gives up. 11 levels is already 1023 fits for one
    ## covariate at one node.
    if (K > 11L) {
      stop2(
        "Scanning the unordered factor in column ", j, " would take 2^", K - 1L,
        " - 1 candidate partitions, each costing a two-group IGSCA fit. ",
        "Collapse it to 11 levels or fewer, make it an ordered factor, or use ",
        "`.splitter = \"native\"`."
      )
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
#' Costs one IGSCA fit per leaf that does not already have one. Failures are
#' counted into `collector$n_fail_leaf` and leave `converged = FALSE` on that
#' node.
#'
#' Not every leaf has `info = NULL` to begin with: a node partykit attempted to
#' split and then left terminal keeps the criteria of that test. The fits are
#' therefore merged into whatever is already there rather than replacing it, so
#' those criteria survive for [root_criteria()] and for node inspection.
#' @noRd
attach_leaf_fits <- function(tree, mf, args, indicators, collector) {
  ids <- partykit::nodeids(tree, terminal = TRUE)
  ## Use the party object's own fitted node ids: they are guaranteed to be in
  ## the same row order as `mf`, which re-deriving them would not be.
  leaf <- tree$fitted[["(fitted)"]]

  ## as.list()/as.partynode() is partykit's own round-trip and is identity-
  ## preserving when the info is left alone, so only the leaves change.
  nd <- as.list(partykit::node_party(tree))
  names(nd) <- vapply(nd, function(n) as.character(n$id), character(1))

  fits <- lapply(ids, function(id) {
    ## A leaf partykit attempted to split already carries the trafo's fit for
    ## exactly this node's rows, in exactly the fields written below and from
    ## the same `sum(E^2)`, so refitting recomputes numbers we already have.
    ## Reuse costs nothing and saves the most expensive fit in the run on a
    ## stump, where the leaf is the root and its rows are the whole sample. It
    ## also keeps a redundant refit of an already-converged node from being
    ## able to count into `n_fail_leaf`. A node that converged but whose
    ## error calculation threw is stored by the trafo as `converged = FALSE`
    ## with no object, so requiring both here still refits it.
    have <- nd[[as.character(id)]]$info
    if (isTRUE(have$converged) && !is.null(have$object)) {
      return(NULL)
    }
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
  ## A NULL entry is a leaf whose fit was reused; the write loop's own
  ## !is.null() guard then leaves that node's info exactly as it stands.
  names(fits) <- as.character(ids)

  for (i in seq_along(nd)) {
    key <- as.character(nd[[i]]$id)
    if (is.null(nd[[i]]$kids) && !is.null(fits[[key]])) {
      info <- as.list(nd[[i]]$info)
      info[names(fits[[key]])] <- fits[[key]]
      nd[[i]]$info <- info
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
#' compared against `alpha`, not the per-covariate p-value.
#'
#' Present whether or not the root went on to split: a root that tested every
#' covariate and rejected none still reports the criteria that made it stop.
#' `NULL` only when no test ran at all, which means the root trafo failed and
#' `n_fail_full` is 1.
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
#' with the model-implied indicator VCV, never the construct VCV. It is handed
#' straight to ndt_dists() as the one distance to compute, so a "DLi" run never
#' reaches calculateDG() and never builds the construct-VCV block matrix.
#' Counts every candidate partition whose statistic could not be computed into
#' collector$n_fail_candidate.
#' 
#'
#' @noRd
partition_stat <- function(stat_kind, model, mf, subset, goes_left, ctrl) {
  ## Deliberately ahead of the tryCatch below, and ahead of the refit: an
  ## unknown stat_kind used to fall out of ds[stat_kind] as a silent NA, which
  ## partykit reads as "no admissible split", so a typo and a genuinely
  ## unsplittable node were indistinguishable. Raised here it is visible;
  ## raised inside the tryCatch it would be swallowed into an NA and a bogus
  ## n_fail_candidate increment. Exact matching, no regex.
  stat_kinds <- c("FITdiff", "DGc", "DGi", "DLc", "DLi")
  if (
    !(length(stat_kind) == 1L && is.character(stat_kind) &&
        stat_kind %in% stat_kinds)
  ) {
    stop2(
      "The following error occured in the partition_stat() function:\n",
      "`stat_kind` must be one of ", paste(stat_kinds, collapse = ", "), "."
    )
  }
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
      ndt_dists(coll$ndt_pools$Sc, coll$ndt_pools$Si, mga$fit, dists = stat_kind)
    },
    error = function(e) NULL
  )
  if (is.null(ds)) { coll$n_fail_candidate <- coll$n_fail_candidate + 1L; return(NA_real_) }
  val <- unname(ds[stat_kind])
  ## real_scalar() has already turned a complex or NaN distance into NA_real_,
  ## so a non-finite value here means this candidate's statistic could not be
  ## computed -- the same thing a throwing bdiagFit() means, and counted the
  ## same way. Without this the counter reads 0 failures on a DGi tree that
  ## stumped precisely because every geodesic distance came back complex.
  if (!is.finite(val)) {
    coll$n_fail_candidate <- coll$n_fail_candidate + 1L
    return(NA_real_)
  }
  val
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




#' A distance that is not one finite real number is not a statistic.
#'
#' calculateDG() eigen-decomposes solve(S) %*% Sigma_hat, a NON-symmetric
#' product, so its eigenvalues are not guaranteed real: it can return a COMPLEX
#' scalar with no error and no warning, and log() of a negative eigenvalue
#' returns NaN. A complex value is the dangerous one. Written into the double
#' vector ndt_dists() returns it coerces the WHOLE vector to complex, so a
#' distance calculateDL() computed perfectly cleanly comes back complex too;
#' argmax_split()'s vapply(..., numeric(1)) then rejects the kernel's return
#' AFTER FUN has returned -- i.e. outside its own per-candidate tryCatch -- and
#' the type error escapes into partykit and kills the tree instead of degrading
#' to "no admissible split".
#'
#' Anything that is not a length-1 finite double becomes NA_real_, which every
#' consumer (is.finite() in argmax_split(), which.max()) already handles.
#'
#' @noRd
real_scalar <- function(x) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    return(NA_real_)
  }
  as.double(x)
}



#' The pooled-vs-MGA distances for one MGA fit, given the precomputed
#' replicated-block pooled VCVs.
#'
#' `dists` names the distances actually wanted. The others are not computed at
#' all -- not merely discarded -- and the block-diagonal MGA matrix a skipped
#' distance would have read is not built either: bdiagFit() runs once per
#' *candidate partition* and is the dominant cost after the refit itself.
#'
#' Isolating the calls is also what keeps a "DLi" run clean. Computing all four
#' unconditionally into one c() meant a complex value out of calculateDG()
#' coerced the whole vector to complex, so `ds["DLi"]` in partition_stat() came
#' back complex even though calculateDL() had never misbehaved.
#'
#' The return is always the same named length-4 double, unrequested slots
#' NA_real_: partition_stat() subscripts it by name, and NA_real_ rather than a
#' bare (logical) NA is what keeps argmax_split()'s vapply(..., numeric(1))
#' satisfied.
#'
#' @noRd
ndt_dists <- function(Sc_pool, Si_pool, mga_fit,
                      dists = c("DGc", "DGi", "DLc", "DLi")) {
  ## Exact matching, no regex. On the production route partition_stat() has
  ## already validated stat_kind ahead of its own tryCatch, so this can only
  ## fire for a direct caller -- where a typo would otherwise return an all-NA
  ## vector that looks like a node nothing could be computed at.
  all_dists <- c("DGc", "DGi", "DLc", "DLi")
  stopifnot(
    "`dists` must be a non-empty subset of DGc/DGi/DLc/DLi" =
      length(dists) > 0L && is.character(dists) && all(dists %in% all_dists)
  )

  out <- c(DGc = NA_real_, DGi = NA_real_, DLc = NA_real_, DLi = NA_real_)

  ## One bdiagFit() per VCV type, and only if some requested distance reads it:
  ## DGc/DLc are construct-VCV distances, DGi/DLi indicator-VCV ones.
  Sc_mga <- if (any(c("DGc", "DLc") %in% dists)) {
    bdiagFit(mga_fit, .type_vcv = "construct")
  }
  Si_mga <- if (any(c("DGi", "DLi") %in% dists)) {
    bdiagFit(mga_fit, .type_vcv = "indicator")
  }

  ## real_scalar() is what keeps `out` a plain double: assigning a complex or a
  ## NaN geodesic distance straight in would re-poison every other slot.
  if ("DGc" %in% dists) {                                        # geodesic,  construct
    out[["DGc"]] <- real_scalar(calculateDG(.matrix1 = Sc_pool, .matrix2 = Sc_mga))
  }
  if ("DGi" %in% dists) {                                        # geodesic,  indicator
    out[["DGi"]] <- real_scalar(calculateDG(.matrix1 = Si_pool, .matrix2 = Si_mga))
  }
  if ("DLc" %in% dists) {                                        # sq-Euclid, construct
    out[["DLc"]] <- real_scalar(calculateDL(.matrix1 = Sc_pool, .matrix2 = Sc_mga))
  }
  if ("DLi" %in% dists) {                                        # sq-Euclid, indicator
    out[["DLi"]] <- real_scalar(calculateDL(.matrix1 = Si_pool, .matrix2 = Si_mga))
  }

  out
}