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



#' Best cutpoint for a (list of) covariate(s) based on the type of splitter statistic
#'
#' 'Best cutpoint' is defined as the partition of data along the covariate
#' that maximizes the difference in the statistic between the multigroup and pooled
#' csem model.
#'
#' This function is based off of `partykit:::.split()`.
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
  ret <- NULL # looping over whichvar is the same idiom as in partykit:::.split
  if (length(whichvar) == 0) return(ret) 
  for (j in whichvar) { # Presumably whichvar is given in the order of the statistically significant covariates
    z <- mf[[j]]
    if ( # same idiom as in partykit:::.split
      isTRUE(ctrl$multiway) && is.factor(z) && !is.ordered(z) &&
        identical(as.numeric(ctrl$maxsurrogate), 0) &&
        nlevels(droplevels(z[subset])) > 1L
    ) {
      # Placeholder for multiway_split funcitonality (not currently supported)
      ret <- multiway_split(j, z, z[subset], ctrl$minbucket)
    } else {
      sub_j <- subset[!is.na(z[subset])] # Further subset the covariate data on only the non-missing covariate values
      cands <- candidate_partitions(
        j = j,
        z = z,
        zs = z[sub_j],
        minbucket = ctrl$minbucket
      )
      if (!length(cands)) {
        next
      }
      stats <- vapply(
        cands, # tryCatch is within vapply, so we take the best split of those associated with convergence success.
        function(cc) {
          tryCatch(
            splitter(
              model = model,
              mf = mf,
              subset = sub_j,
              goes_left = cc$goes_left,
              ctrl = ctrl
            ),
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
      ret <- cands[[which.max(stats)]]$split # Here is where the maximum is taken. cands is a nested list
    }
    if (!is.null(ret)) {
      (break)()
    }
  }
  if (scanned) {
    collector$n_split_scan <- collector$n_split_scan + 1L # Successful scan of covariate
    if (!got_stat) {
      collector$n_fail_split <- collector$n_fail_split + 1L # Failure to split on this covariate
    }
  }
  return(ret)
}

#' One kid per level of an unordered factor (incomplete)
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
#'
#' **The `partysplit` encoding is ctree's**, from
#' `partykit:::.ctree_test_internal()` (`R/ctree.R`): numeric and ordered
#' covariates get `.partysplit(j, breaks = sp, index = 1L:2L)`, unordered ones
#' `.partysplit(j, index = as.integer(sp) + 1L)`. Driving that function on a
#' node directly confirms all three conventions this code relies on -- an
#' ordered break is a position in `levels(z)`, not in the node's own levels;
#' an unordered `index` is 1/2 over `levels(z)` with `NA` for levels the node
#' never sees; a numeric break is an observed value with `zs <= break` going
#' to kid 1.
#'
#'
#' * **Numeric.** mob takes `uz <- sort(unique(zselect))` and forms
#'   `zs <- zselect <= uz[i]`; `.objfun_test()` does the same over the
#'   integer-coded covariate, `sleft <- subset[ix[subset] <= u]` for `u` in
#'   `which(ixtab > 0)`. The cut dropped here up front -- `zs <= max(uz)`,
#'   which sends every row left -- is one both generate and then reject on
#'   `minbucket`. Only the observed-value convention is implemented, that is
#'   ctree's `intersplit = FALSE` and mob's `numsplit = "left"`; the
#'   interpolated-midpoint variants are not supported.
#' * **Ordered.** The cumulative level sets come from the ordered branch of
#'   `partykit:::.mob_grow_getlevels()` (`R/utils.R`) or, inside
#'   `.objfun_test()`, from the same `ix[subset] <= u` loop the numeric case
#'   uses (`ORDERED <- is.ordered(x) || is.numeric(x)`). `match(levs[i],
#'   levels(z))` is mob's `match(levels(zselect)[which.min(dev)], olevels)`.
#' * **Unordered.** `.mob_grow_getlevels()` builds the same `2^(K - 1) - 1`
#'   bit-pattern matrix `intToBits(m)` builds here, in the same order, over the
#'   levels surviving `factor(zselect)` (that is, `droplevels(zs)`); membership
#'   is `zselect %in% levels(zselect)[w]`. Re-expanding to `z`'s original
#'   levels with `NA` for absent ones is mob's
#'   `ix <- structure(rep.int(NA_integer_, length(olevels)), names = olevels)`
#'   then `ix[colnames(al)] <- !al[which.min(dev), ]; as.integer(ix) + 1L` --
#'   that negation is why the selected group is kid `1L` here.
#'
#' @param j Column index of the covariate in the model frame.
#' @param z The covariate over all rows, which fixes the factor levels and
#'   level order the emitted `partysplit` is expressed in.
#' @param zs The covariate over this node's non-missing rows.
#' @param minbucket Smallest child either kid may be left with.
#' @return A list of candidate partitions, one per admissible split, each a
#'   list with `$split` (a `partysplit`) and `$goes_left` (the logical vector
#'   over `zs` saying which of the node's rows go to the first kid). Empty when
#'   the covariate is constant in the node or every cut violates `minbucket`.
#'
#' @examples
#' ## Numeric: one cut per distinct value except the largest, `zs <= ct`.
#' zn <- c(1, 2, 2, 3, 5, 8)
#' candidate_partitions(1L, zn, zn, minbucket = 1L)
#'
#' ## Ordered: K - 1 cumulative level sets; the break is an index into
#' ## levels(z), so it stays comparable across nodes.
#' zo <- factor(c("lo", "lo", "mid", "hi", "hi"),
#'              levels = c("lo", "mid", "hi"), ordered = TRUE)
#' candidate_partitions(2L, zo, zo, minbucket = 1L)
#'
#' ## Unordered: all 2^(K - 1) - 1 bipartitions, carried in `index` rather
#' ## than `breaks`; 1 = first kid, 2 = second.
#' zf <- factor(c("a", "a", "b", "b", "c", "c"))
#' candidate_partitions(3L, zf, zf, minbucket = 1L)
#'
#' ## Nothing to split on.
#' candidate_partitions(1L, zn, rep(2, 4), minbucket = 1L)
#' @noRd
candidate_partitions <- function(j, z, zs, minbucket) {
  # browser()
  # Returns only candidate partitions whose child nodes would have a sample size that is greater than or equal to the minbucket.
  keep_min <- function(cands) {
    Filter( 
      function(cc) {
        nl <- sum(cc$goes_left)
        nl >= minbucket && (length(zs) - nl) >= minbucket
      },
      cands
    )
  }
  # Create list of break points based on the metric of z
  if (is.numeric(z) && !is.factor(z)) {
    uz <- sort(unique(zs))
    if (length(uz) < 2L) {
      return(list())
    }
    # Vector of cut points in uz
    cuts <- uz[-length(uz)] # Equivalent to uz[1:(length(uz)-1)]

    keep_min(lapply(cuts, function(ct) {
      list(
        goes_left = zs <= ct, # mob_grow_findsplit(): `zs <- zselect <= uz[i]`
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
        goes_left = zs %in% levs[1L:i], # .mob_grow_getlevels(), ordered branch
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
    # Maximum number of levels in a nominal covariate. 
    # 11 levels is 1023 fits for one covariate at one node.
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
        goes_left = zs %in% g, # mob_grow_findsplit(): `zselect %in% levels(zselect)[w]`
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
  e$n_split_scan <- 0L
  e$n_fail_split <- 0L
  e$n_fail_leaf <- 0L
  e
}




#' Fit the IGSCA model at every terminal node and attach it to the tree
#'
#' partykit only calls the trafo at nodes it *attempts to split*, so leaves come
#' back with `info = NULL` and every leaf-reading method sees nothing.
#'
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
  
  leaf <- tree$fitted[["(fitted)"]]

  ## as.list()/as.partynode() is partykit's own round-trip and is identity-
  ## preserving when the info is left alone, so only the leaves change.
  nd <- as.list(partykit::node_party(tree)) # This is the idiom to go
  names(nd) <- vapply(nd, function(n) as.character(n$id), character(1))

  fits <- lapply(ids, function(id) {
    
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
    E <- tryCatch(calculateGSCAErrors(ft$fit), error = function(e) NULL)
    list(
      nobs = length(rows),
      converged = TRUE,
      objfun = if (is.null(E)) NA_real_ else sum(E^2),
      object = ft$fit
    )
  })
  
  names(fits) <- as.character(ids)

  for (i in seq_along(nd)) {
    key <- as.character(nd[[i]]$id)
    if (is.null(nd[[i]]$kids) && !is.null(fits[[key]])) { # Makes sure that we're not over-writing already fitted models or inner nodes
      info <- as.list(nd[[i]]$info)
      info[names(fits[[key]])] <- fits[[key]] # Rewrites the entire info, nobs, converged, objfun and object
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

#' Per-node parameter estimates of an igsca tree
#'
#' Reports the number of observations and the parameter estimates of each
#' node's fit, one column per estimate. The estimates are those of [tidy()]
#' applied to the node's `cSEMResults`, named by that node's `term`. 
#'
#'
#' @param object A tree returned by [doTrees()].
#' @param node Node ids to report. Defaults to the terminal nodes.
#' @param parameters Which parameter estimates to report. Passed on to
#'   [tidy()], which validates it. Defaults to `"all"`. The effect estimates
#'   are not available here; see below.
#' @param drop If `TRUE` (the default), a single-node result is simplified from
#'   a one-row matrix to a named vector.
#' @param ... Ignored.
#'
#' @returns A `matrix` whose first column is `nobs` and whose remaining columns
#'   are one parameter estimate each, one row per node, row names being the
#'   node ids; a named vector when one node is requested and `drop = TRUE`.
#' @seealso [doTrees()], [tidy()]
#' @export
#' @importFrom partykit nodeids info_node nodeapply
coef.igsca_tree <- function(object, node = NULL, parameters = "all",
                            drop = TRUE, ...) {
  if (identical(parameters, "Effect_estimates")) {
    stop2(
      "`parameters = \"Effect_estimates\"` is not supported by coef(): an ",
      "indirect or total effect carries the same term as the path of the ",
      "same name, so its estimates cannot be told apart by column name. ",
      "Use tidy() on a node's fit for those."
    )
  }
  if (is.null(node)) {
    node <- partykit::nodeids(object, terminal = TRUE)
  }
  rows <- partykit::nodeapply(object, ids = node, FUN = function(n) {
    i <- partykit::info_node(n)
    nobs <- if (is.null(i$nobs)) NA_real_ else as.numeric(i$nobs)
    
    if (is.null(i$object)) {
      return(c(nobs = nobs))
    }
    est <- tidy(i$object, parameters = parameters)
    
    est <- est[!est$op %in% c("Indirect_effect", "Total_effect"), ]
    c(
      nobs = nobs,
      stats::setNames(as.numeric(est$estimate), est$term)
    )
  })
  
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  cf <- do.call(
    "rbind",
    lapply(rows, function(r) {
      out <- stats::setNames(rep(NA_real_, length(cols)), cols)
      out[names(r)] <- r
      out
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

#' Reserved name of the two-level grouping column partition_stat() adds to a node's data
#' and hands csem() as `.id`.
#'
#' @noRd
TREE_GROUP_COL <- "TREETEMPGROUP"

#'
#' Observed statistic for one candidate partition. `stat_kind` is "FITdiff"
#' (NPT) or an ndt_dists() distance name -- "DGi"/"DLi" only.
#'
#' @noRd
partition_stat <- function(stat_kind, model, mf, subset, goes_left, ctrl) {
  
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
  ## Two groups, always: the statistic is a pooled-vs-2-group comparison
  d <- mf[subset, ctrl$indicators, drop = FALSE]
  d[[TREE_GROUP_COL]] <- factor(ifelse(goes_left, 1L, 2L), levels = c(1L, 2L))
  mga <- try_fit(d, ctrl$args, .id = TREE_GROUP_COL)
  if (!mga$ok) {
    coll$n_fail_candidate <- coll$n_fail_candidate + 1L
    return(NA_real_)
  }
  if (stat_kind == "FITdiff") {
    if (is.null(model$object)) {
      return(NA_real_)
    }
    val <- tryCatch(
      calculateFIT(mga$fit) - calculateFIT(model$object), # model$object is the parent
      error = function(e) NULL
    )
    if (is.null(val)) {
      coll$n_fail_candidate <- coll$n_fail_candidate + 1L
      return(NA_real_)
    }
    return(val)
  }
  
  ds <- tryCatch(
    {
      if (is.null(model$object)) {
        return(NA_real_)
      } else {
        ndt_dists(
          Sc_pool = bdiagFit(model$object, .n_blocks = 2L, .type_vcv = "construct"),
          Si_pool = bdiagFit(model$object, .n_blocks = 2L, .type_vcv = "indicator"),
          mga_fit = mga$fit,
          dists = stat_kind
        )
      }
    },
    error = function(e) NULL
  )
  if (is.null(ds)) {
    coll$n_fail_candidate <- coll$n_fail_candidate + 1L
    return(NA_real_)
  }
  val <- unname(ds[stat_kind])
  ## A non-finite value here means this candidate's statistic could not be computed.
  if (!is.finite(val)) {
    coll$n_fail_candidate <- coll$n_fail_candidate + 1L
    return(NA_real_)
  }
  val
}

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


#' The pooled-vs-MGA distances for one MGA fit, given the precomputed
#' replicated-block pooled VCVs.
#'
#' `dists` names the distances actually wanted.
#'
#' Isolating the calls is also what keeps a "DLi" run clean. Computing all four
#' unconditionally into one c() meant a complex value out of calculateDG()
#' coerced the whole vector to complex, so `ds["DLi"]` in partition_stat() came
#' back complex even though calculateDL() had never misbehaved.
#'
#'
#' @noRd
ndt_dists <- function(Sc_pool, Si_pool, mga_fit,
                      dists = c("DGc", "DGi", "DLc", "DLi")) {
  
  all_dists <- c("DGc", "DGi", "DLc", "DLi")
  stopifnot(
    "`dists` must be a non-empty subset of DGc/DGi/DLc/DLi" =
      length(dists) > 0L && is.character(dists) && all(dists %in% all_dists)
  )

  out <- c(DGc = NA_real_, DGi = NA_real_, DLc = NA_real_, DLi = NA_real_)

  Sc_mga <- if (any(c("DGc", "DLc") %in% dists)) {
    bdiagFit(mga_fit, .type_vcv = "construct")
  }
  Si_mga <- if (any(c("DGi", "DLi") %in% dists)) {
    bdiagFit(mga_fit, .type_vcv = "indicator")
  }

  ## real_scalar() is forces a plain double. This is necessary because calculateDG() can
  ## return a complex scalar, which would silently convert the type of data of every other 
  ## distance. 
  ##
  ## This happens because calculateDG() eigen-decomposes solve(S) %*% Sigma_hat, a NON-symmetric
  ## product, so its eigenvalues are not guaranteed real: it can return a COMPLEX
  ## scalar with no error and no warning, and log() of a negative eigenvalue returns NaN.
  real_scalar <- function(x) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x)) {
    return(NA_real_)
  }
    as.double(x)
  }
  if ("DGc" %in% dists) { # geodesic,  construct
    out[["DGc"]] <- calculateDG(.matrix1 = Sc_pool, .matrix2 = Sc_mga) |> 
      real_scalar()
  }
  if ("DGi" %in% dists) { # geodesic,  indicator
    out[["DGi"]] <- calculateDG(.matrix1 = Si_pool, .matrix2 = Si_mga) |> 
      real_scalar()
  }
  if ("DLc" %in% dists) { # sq-Euclid, construct
    out[["DLc"]] <- calculateDL(.matrix1 = Sc_pool, .matrix2 = Sc_mga) |> 
      real_scalar()
  }
  if ("DLi" %in% dists) { # sq-Euclid, indicator
    out[["DLi"]] <- calculateDL(.matrix1 = Si_pool, .matrix2 = Si_mga) |> 
      real_scalar()
  }

  return(out)
}