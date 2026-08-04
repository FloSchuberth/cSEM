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
#'
#' @returns
#'
#' @export
#' @examples
igsca_tree_control <- function(alpha = 0.05,
                               bonferroni = FALSE,
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
#'
#' @param E
#'
#' @returns
#'
#' @export
#' @examples
influence_vec <- function(E) {
  matrix(rowSums(E^2), ncol = 1L)
}



#' Casewise squared-residual matrix (n x q, multivariate) -- COIN_mat.
#'
#' @param E
#'
#' @returns
#'
#' @export
#' @examples
influence_mat <- function(E) {
  E^2
}


#' Title
#' 
#' `.id = NULL` yields the pooled single-group fit; `.id = "group"` the MGA fit.
#' 
#' @param .data
#' @param .model
#' @param .id
#'
#' @returns
#'
#' @export
#' @examples
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

#' Title
#'
#' @param .data
#' @param .model
#' @param .id
#'
#' @returns
#'
#' @export
#' @examples
try_fit <- function(.data, .model, .id = NULL) {
  fit <- suppressWarnings(tryCatch(
    fit_csem(.data, .model, .id),
    error = function(e) NULL
  ))
  ok <- !is.null(fit) &&
    isTRUE(tryCatch(csem_converged(fit), error = function(e) FALSE)) 
  list(fit = fit, ok = ok)
}

#' Title
#'
#' @param fit
#'
#' @returns
#'
#' @export
#' @examples
csem_converged <- function(fit) {
  if (inherits(fit, "cSEMResults_multi")) {
    all(vapply(fit, function(x) sum(verify(x)) == 0, logical(1)))
  } else {
    sum(verify(fit)) == 0
  }
}


#' Title
#'
#' Generic extree_fit() driver for the partition family. `d` is an
#' extree_data object; `selector` already satisfies extree_fit's selectfun
#' contract and is passed through VERBATIM; `splitter` is a statistic kernel
#' argmaxed by argmax_split() via the thin splitfun closure below. ctrl must
#' carry $collector (and, for IGSCA use, $model + $indicators for the
#' testers). Exposed separately from igsca_tree() so
#' ctree_equivalence_proof.R can drive it with Gaussian fixtures (no IGSCA
#' fitting). Note: on this path a trafo's estfun (if any) is node-local
#' (rows = subset); only the native ctree path needs the full-n scatter.
#' 
#' @param d
#' @param trafo
#' @param selector
#' @param splitter
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @importFrom partykit ctree_control
grow_extree <- function(d, trafo, selector, splitter, ctrl) {
  collector <- ctrl$collector
  splitfun <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
    argmax_split(
      splitter,
      collector,
      model,
      model.frame(data),
      subset,
      whichvar,
      ctrl
    )
  }

  ## ctree_control() is only a scaffold: it supplies the extree knobs we do
  ## not own (criterion, splittry, saveinfo, ...). Everything igsca-specific
  ## is overridden from ctrl below; the trailing loop also carries the tester
  ## config (model, indicators, collector, R_test, max_cuts, ...) into the
  ## ctrl the callbacks receive.
  ectrl <- partykit::ctree_control()
  ectrl$update <- TRUE
  ectrl$logmincriterion <- log(1 - ctrl$alpha)
  ectrl$maxsurrogate <- 0L
  ectrl$mtry <- Inf
  for (nm in names(ctrl)) {
    ectrl[[nm]] <- ctrl[[nm]]
  }

  n <- nrow(model.frame(d))
  fit <- partykit::extree_fit(
    data = d,
    trafo = trafo,
    converged = TRUE,
    selectfun = selector,
    splitfun = splitfun,
    svselectfun = selector,
    svsplitfun = splitfun,
    partyvars = d$variables$z,
    subset = seq_len(n),
    weights = rep.int(1L, n),
    ctrl = ectrl
  )
  fit$nodes
}

#' Title
#'
#' @param splitter
#' @param collector
#' @param model
#' @param mf
#' @param subset
#' @param whichvar
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
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
#' A kernel that throws is caught inside [argmax_split()] and turned into `NA`,
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


#' Title
#'
#' @param j
#' @param z
#' @param zs
#' @param max_cuts
#' @param minbucket
#'
#' @returns
#'
#' @export
#' @examples
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

#' Title
#'
#' @returns List
#'
new_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$root_seen <- FALSE      # first trafo call = root fit (full vs node fail)
  e$n_fail_full <- 0L
  e$n_fail_node <- 0L
  e$n_fail_resample <- 0L
  # Split-kernel scans: how many argmax_split() calls actually evaluated the
  # kernel, and how many of those produced nothing usable. See argmax_split().
  e$n_split_scan <- 0L
  e$n_fail_split <- 0L
  e$scan <- list()
  e$scan_subset <- NULL
  e$scan_splitter <- NULL
  e
}




#' Get root-node criteria of igsca tree
#'
#' Root-node criteria matrix as stored by extree/ctree (saveinfo = TRUE):
#' columns = tested covariates (all-NA columns dropped by extree), rows
#' "statistic" / "p.value" (both raw scale) / "criterion" (log(1 - p) with
#' extree's own statistic-rank tie-break already added). NULL when the root
#' trafo failed (no test ran).
#' 
#' @param tree
#'
#' @returns Vector
#'
#' @importFrom partykit info_node node_party
#' @export
root_criteria <- function(tree) {
  info <- partykit::info_node(partykit::node_party(tree))
  if (is.null(info)) {
    return(NULL)
  }
  info$criterion
}

#' Get coefficients of tree object
#'
#' @param object
#' @param node
#' @param drop
#' @param ...
#'
#' @returns List of Vectors
#'
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
      c(nobs = i$nobs, objfun = i$objfun)
    })
  )
  if (drop) {
    drop(cf)
  } else {
    cf
  }
}

#' Title
#'
#' @param x
#' @param terminal_panel
#' @param FUN
#' @param tp_args
#' @param ...
#'
#' @returns Graphic
#'
#' @export
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



# Partition-family selectors & splitters for igsca_tree() -------------------
#
# The five study methods (registry cut 2026-07-14):
#   COIN_ssr / COIN_mat     -> igsca_ctree() + influence_ssr/influence_mat
#                              (native partykit::ctree; nothing lives here)
#   NPT / NDT_DGi / NDT_DLi -> the plain selector/splitter functions below.
#
# Selector contract (= extree_fit's selectfun, passed through verbatim; the
# formal names below are checked by extree_fit and must not change):
#   select_*(model, trafo, data, subset, weights, whichvar, ctrl) ->
#     list(criteria = 2 x p matrix, rows "statistic" [log scale] and
#          "p.value" [log1p(-p) scale], columns named like model.frame(data))
# Splitter contract (statistic kernel argmaxed by the engine's splitfun):
#   split_*(model, mf, subset, goes_left, ctrl) -> scalar observed statistic
#
# Statistic-first scan (spec section 3): observed statistic at every
# candidate (one 2-group MGA fit each; the node's pooled fit comes free from
# the trafo as model$object), argmax per covariate, permutation p-value ONLY
# at the argmax. The within-covariate selection effect is deliberately
# uncorrected -- that bias is Study 1's measurand. Scan results are cached in
# ctrl$collector so a matched splitter costs nothing; a mixed pair (any
# select_* with any split_*) simply re-scans with its own kernel.
#



#' Title
#' 
#' Node data for one candidate partition: children as `group` levels 1/2.
#' 
#' @param mf
#' @param subset
#' @param indicators
#' @param goes_left
#'
#' @returns
#'
#' @export
#' @examples
node_group_data <- function(mf, subset, indicators, goes_left) {
  d <- mf[subset, indicators, drop = FALSE]
  cbind(group = factor(ifelse(goes_left, 1L, 2L), levels = c(1L, 2L)), d)
}



#' Title
#'
#' Pooled model-implied VCVs replicated to 2 blocks (constant per node;
#' cached in the collector under $ndt_pools keyed by the subset).
#' 
#' @param single_fit
#'
#' @returns
#'
#' @export
#' @examples
ndt_pools <- function(single_fit) {
  list(
    Sc = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "construct"),
    Si = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "indicator")
  )
}

#' Title
#'
#' 
#' Observed statistic for one candidate partition. `stat_kind` is "FITdiff"
#' (NPT) or an ndt_dists() distance name -- "DGi"/"DLi" only: Study 1 works
#' with the model-implied indicator VCV, never the construct VCV (the shared
#' helper also computes "DGc"/"DLc", which stay unused here). Counts every
#' failed auxiliary fit into collector$n_fail_resample (recorded, never
#' redrawn).
#' 
#' @param stat_kind
#' @param model
#' @param mf
#' @param subset
#' @param goes_left
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
partition_stat <- function(stat_kind, model, mf, subset, goes_left, ctrl) {
  coll <- ctrl$collector
  d <- node_group_data(mf, subset, ctrl$indicators, goes_left)
  mga <- try_fit(d, ctrl$model, .id = "group")
  if (!mga$ok) { coll$n_fail_resample <- coll$n_fail_resample + 1L; return(NA_real_) }
  if (stat_kind == "FITdiff") {
    if (is.null(model$object)) {
      return(NA_real_)
    }
    val <- tryCatch(
      cSEM::calculateFIT(mga$fit) - cSEM::calculateFIT(model$object),
      error = function(e) NULL
    )
    if (is.null(val)) {
      coll$n_fail_resample <- coll$n_fail_resample + 1L
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
  if (is.null(ds)) { coll$n_fail_resample <- coll$n_fail_resample + 1L; return(NA_real_) }
  unname(ds[stat_kind])
}



#' Title
#'
#' 
#' Scan all candidates of covariate j; returns NULL when nothing admissible.
#' 
#' @param stat_kind
#' @param j
#' @param model
#' @param mf
#' @param subset
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
scan_covariate <- function(stat_kind, j, model, mf, subset, ctrl) {
  cands <- candidate_partitions(j, mf[[j]], mf[[j]][subset],
                                ctrl$max_cuts, ctrl$minbucket)
  if (!length(cands)) return(NULL)
  stats <- vapply(
    cands,
    function(cc) {
      partition_stat(stat_kind, model, mf, subset, cc$goes_left, ctrl)
    },
    numeric(1)
  )
  if (!any(is.finite(stats))) return(NULL)
  k <- which.max(stats)
  list(stat = stats[k], split = cands[[k]]$split, goes_left = cands[[k]]$goes_left)
}

#' Title
#'
#' One-sided permutation p-value at the argmax partition (post-BootBalVerify
#' Study 0 conventions, 2026-07-13): permute the child labels, null >= obs
#' without abs(). Failed resamples are dropped from the null (finite counts)
#' and counted by partition_stat.
#' 
#' @param stat_kind
#' @param obs
#' @param model
#' @param mf
#' @param subset
#' @param goes_left
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
permutation_pvalue <- function(stat_kind, obs, model, mf, subset, goes_left, ctrl) {
  n <- length(subset)
  R <- ctrl$R_test
  perm <- idx_permutation(n, R)
  nul <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    nul[r] <- partition_stat(
      stat_kind,
      model,
      mf,
      subset,
      goes_left[perm[r, ]],
      ctrl
    )
  }
  (1 + sum(nul >= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1)
}


#' Title
#'
#' Shared worker behind the three partition selectors: statistic-first scan
#' per covariate, cache the argmax for the splitfun, permutation p at the
#' argmax only. `matched_split_fn` is the splitter function whose kernel
#' equals this selector's scan statistic -- the cache key grow_extree's
#' splitfun compares against with identical().
#' 
#' @param stat_kind
#' @param matched_split_fn
#' @param model
#' @param data
#' @param subset
#' @param whichvar
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
partition_select <- function(stat_kind, matched_split_fn,
                             model, data, subset, whichvar, ctrl) {
  coll <- ctrl$collector
  mf <- model.frame(data)
  crit <- matrix(
    NA_real_,
    2L,
    ncol(mf),
    dimnames = list(c("statistic", "p.value"), names(mf))
  )
  ## (re)arm the node-local scan cache for the engine's splitfun
  coll$scan <- list()
  coll$scan_subset <- subset
  coll$scan_splitter <- matched_split_fn
  if (is.null(model$object)) return(list(criteria = crit))  # node fit failed
  for (j in whichvar) {
    sc <- scan_covariate(stat_kind, j, model, mf, subset, ctrl)
    if (is.null(sc)) {
      next
    }
    coll$scan[[as.character(j)]] <- sc
    p <- permutation_pvalue(
      stat_kind,
      sc$stat,
      model,
      mf,
      subset,
      sc$goes_left,
      ctrl
    )
    if (!is.finite(p)) {
      next
    }
    ## One-sided world: negative observed statistics must not out-rank
    ## positive ones on the log scale -- floor them at eps instead of
    ## folding sign via abs().
    crit["statistic", j] <- log(max(sc$stat, .Machine$double.eps))
    crit["p.value", j] <- log1p(-min(p, 1 - 1e-12))
  }
  list(criteria = crit)
}

## The three partition selectors and the splitter kernels: deliberately
## repetitive plain functions (no factories). Each selector pins its scan
## statistic and names its matched splitter; each splitter is the bare
## observed-statistic kernel. Any selector can be paired with any splitter
## in igsca_tree() -- an unmatched pair just re-scans (cache miss) -- and
## every kernel can also serve as igsca_ctree()'s mixed-pair `splitter`
## (COIN variable selection + model-comparison split point, 2026-07-16).
## Distances are model-implied INDICATOR-VCV only (dGi/dLi) -- Study 1 never
## works with the construct VCV (user correction 2026-07-16).

#' Title
#'
#' @param model
#' @param trafo
#' @param data
#' @param subset
#' @param weights
#' @param whichvar
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
select_npt <- function(model, trafo, data, subset, weights, whichvar, ctrl) {
  partition_select(
    "FITdiff",
    split_max_fitdiff,
    model,
    data,
    subset,
    whichvar,
    ctrl
  )
}

#' Title
#'
#' @param model
#' @param trafo
#' @param data
#' @param subset
#' @param weights
#' @param whichvar
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
select_ndt_dgi <- function(
  model,
  trafo,
  data,
  subset,
  weights,
  whichvar,
  ctrl
) {
  partition_select("DGi", split_max_dgi, model, data, subset, whichvar, ctrl)
}

#' Title
#'
#' @param model
#' @param trafo
#' @param data
#' @param subset
#' @param weights
#' @param whichvar
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
select_ndt_dli <- function(
  model,
  trafo,
  data,
  subset,
  weights,
  whichvar,
  ctrl
) {
  partition_select("DLi", split_max_dli, model, data, subset, whichvar, ctrl)
}

#' Title
#'
#' @param model
#' @param mf
#' @param subset
#' @param goes_left
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
split_max_fitdiff <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("FITdiff", model, mf, subset, goes_left, ctrl)
}

#' Title
#'
#' @param model
#' @param mf
#' @param subset
#' @param goes_left
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
split_max_dgi <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DGi", model, mf, subset, goes_left, ctrl)
}

#' Title
#'
#' @param model
#' @param mf
#' @param subset
#' @param goes_left
#' @param ctrl
#'
#' @returns
#'
#' @export
#' @examples
split_max_dli <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DLi", model, mf, subset, goes_left, ctrl)
}


#' Title
#' 
#' Permutation resampling: sample WITHOUT replacement, within each stratum, so
#' every replicate is a permutation of its stratum (each original row used
#' exactly once per replicate). Mirrors/Near-Copy of boot:::permutation.array.
#' 
#' @param n
#' @param R
#'
#' @returns
#'
#' @export
#' @examples
idx_permutation <- function(n, R) {
  # Argument removed
  # , strata = NULL 

  # row r starts as the identity 1:n
  out <- matrix(rep(seq_len(n), each = R), nrow = R, ncol = n)
  # if (is.null(strata)) {
    for (r in seq_len(R)) out[r, ] <- sample.int(n)
  # } else {
  #   strata <- as.integer(as.factor(strata))
  #   for (s in unique(strata)) {
  #     gp <- which(strata == s)
  #     if (length(gp) > 1L) {
  #       for (r in seq_len(R)) out[r, gp] <- sample(gp)
  #     }
  #   }
  # }
  out
}


#' Title
#'
#' The four pooled-vs-MGA distances for one MGA fit, given the precomputed
#' replicated-block pooled VCVs.
#' 
#' @param Sc_pool
#' @param Si_pool
#' @param mga_fit
#'
#' @returns
#'
#' @export
#' @examples
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