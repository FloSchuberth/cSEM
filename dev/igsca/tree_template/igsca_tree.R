# IGSCA trees on partykit ---------------------------------------------------
#
# Two constructors, one per test framework (spec: docs/superpowers/specs/
# 2026-07-12-study1-igsca-tree-unbiased-selection-design.md; refactor:
# docs/superpowers/plans/2026-07-14-igsca-tree-refactor.md):
#
#   igsca_ctree()  COIN family (COIN_ssr / COIN_mat). Variable selection is
#     partykit::ctree()'s own libcoin machinery; the only IGSCA-specific code
#     is the ytrafo (node IGSCA fit -> casewise residual influence via
#     influence_ssr()/influence_mat()). teststat = "quadratic" with testtype
#     = "MonteCarlo" is the same Strasser-Weber permutation test the previous
#     implementation hand-rolled through coin::independence_test(
#     distribution = approximate(R_test)). Split-point selection defaults to
#     ctree's native maxstat scan (splitstat = "quadratic"); alternatively a
#     partition kernel (split_max_fitdiff / split_max_dgi / split_max_dli --
#     model-implied indicator-VCV distances or FIT, never the construct VCV)
#     can be passed as `splitter`, mixing COIN variable selection with a
#     model-comparison split point (2026-07-16).
#
#   igsca_tree()   partition family (NPT / NDT_*). These statistics need a
#     fresh 2-group MGA refit per candidate partition -- not expressible as a
#     libcoin linear statistic -- so this path drives partykit::extree_fit()
#     with a custom selectfun/splitfun. `selector` and `splitter` are PLAIN
#     functions (defined in igsca_testers.R):
#       selector = extree_fit's own selectfun, passed through verbatim:
#         function(model, trafo, data, subset, weights, whichvar, ctrl) ->
#         list(criteria = 2 x p matrix, rows "statistic" [log scale] and
#         "p.value" [log1p(-p) scale], columns named like model.frame(data)).
#         extree_fit itself applies ctrl$bonferroni, ranks the variables, and
#         stores the raw-scale criteria matrix in every node's info.
#       splitter = observed-statistic kernel for one candidate partition:
#         function(model, mf, subset, goes_left, ctrl) -> scalar. The
#         engine's splitfun argmaxes it over candidate_partitions().
#     On this path the trafo returns no estfun; the `model` reaching the
#     selector/splitter carries the pooled node fit in model$object.
#
# What is a "trafo"? partykit's extree framework separates WHAT is tested
# from HOW the tree machinery tests it, and the trafo ("transformation") is
# the per-node bridge between the two. Each time the tree considers a node,
# the trafo receives the node's rows (`subset`), fits the node's model to
# just those rows, and returns whatever the framework may need -- most
# importantly `estfun`: an n x q matrix of casewise influence scores h_i,
# intuitively "how much, and in what way, does observation i deviate from
# this node's fitted model?". The conditional-inference test then asks
# whether those deviations are systematically associated with any candidate
# covariate; if yes, the node is split and the trafo is called again on each
# child (update = TRUE). For IGSCA the influence comes from the casewise
# GSCA residuals (cSEM:::calculateGSCAErrors). The partition family has no
# casewise influence at all, so its trafo instead hands over the pooled node
# fit (model$object) for the testers to refit candidate 2-group MGAs
# against. Returning `converged = FALSE` from a trafo turns the node into a
# leaf -- the no-throw failure channel.
#
# Definitions only; requires partykit. try_fit() comes from
# R/MGA/csem_test_helpers.R (source it before this file).

igsca_tree_control <- function(alpha = 0.05,
                               # User decision 2026-07-13: no Bonferroni in
                               # the simulations. ctree_equivalence_proof.R
                               # sets TRUE explicitly to match ctree's
                               # testtype = "Bonferroni".
                               bonferroni = FALSE,
                               minbucket = 30L,
                               minsplit = 2L * minbucket,
                               maxdepth = 3L,
                               # Scan-grid cap for partition testers and for
                               # igsca_ctree() mixed-pair splitters; the
                               # native COIN split (splitter = NULL) scans
                               # all cutpoints inside libcoin instead.
                               max_cuts = 20L,
                               R_test = 500L,
                               coin_distribution = c("approximate", "asymptotic")) {
  list(alpha = alpha, bonferroni = bonferroni, minbucket = as.integer(minbucket),
       minsplit = as.integer(minsplit), maxdepth = as.integer(maxdepth),
       max_cuts = as.integer(max_cuts), R_test = as.integer(R_test),
       coin_distribution = match.arg(coin_distribution))
}

#' Failure/diagnostic collector for one tree call. An environment so the
#' trafo and the extree callbacks can write to it despite partykit's fixed
#' signatures. The scan_* fields implement the node-local scan cache: a
#' partition selector has already computed every candidate's observed
#' statistic, so a matched splitter reuses the cached argmax instead of
#' refitting (see partition_select / grow_extree's splitfun).
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

#' Candidate binary partitions of covariate j within the node (rows `subset`).
#'
#' Used wherever a split statistic needs a 2-group MGA refit per candidate
#' (the partition family, and igsca_ctree()'s mixed-pair splitters) -- such
#' scans cannot ride partykit's own split search, which evaluates libcoin
#' linear statistics of a fixed influence and never refits. The native COIN
#' split path (igsca_ctree with splitter = NULL) uses partykit's scan.
#'
#' Numeric: midpoints between consecutive unique values, thinned to the
#' `max_cuts` interior quantiles when there are more than max_cuts of them
#' (keeps the many-cutpoints-vs-binary asymmetry while bounding compute).
#'
#' Unordered factor: all 2^(K-1)-1 binary groupings.
#'
#' Ordered factor: order-respecting cuts. Candidates violating minbucket
#' are dropped here so every downstream consumer sees only admissible
#' partitions. Adapted from R/tree/HowTree/lavaan.R candidate_partitions().
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

#' Engine split-point search shared by grow_extree() and igsca_ctree()'s
#' mixed-pair path: over the criterion-ranked variables in `whichvar`,
#' return the argmax partysplit of the `splitter` kernel across
#' candidate_partitions(), or NULL when nothing is admissible. Reuses a
#' selector's cached scan when the SAME kernel already scanned this node
#' (matched pairs cost nothing; the cache key is the splitter function
#' itself, compared with identical()). Kernel errors become NA candidates --
#' never an exception into partykit.
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

#' Generic extree_fit() driver for the partition family. `d` is an
#' extree_data object; `selector` already satisfies extree_fit's selectfun
#' contract and is passed through VERBATIM; `splitter` is a statistic kernel
#' argmaxed by argmax_split() via the thin splitfun closure below. ctrl must
#' carry $collector (and, for IGSCA use, $model + $indicators for the
#' testers). Exposed separately from igsca_tree() so
#' ctree_equivalence_proof.R can drive it with Gaussian fixtures (no IGSCA
#' fitting). Note: on this path a trafo's estfun (if any) is node-local
#' (rows = subset); only the native ctree path needs the full-n scatter.
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

#' Grow one partition-family IGSCA tree (NPT / NDT_*). `data` must contain
#' exactly the model's indicator columns plus the partitioning covariates
#' named in `covariates`. `selector`/`splitter` are the plain functions from
#' igsca_testers.R (contracts in this file's header).
igsca_tree <- function(data, model, covariates, selector, splitter,
                       control = igsca_tree_control()) {
  indicators <- setdiff(names(data), covariates)
  fml <- stats::as.formula(paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(covariates, collapse = " + ")
  ))
  d <- partykit::extree_data(
    fml,
    data = data,
    yx = "none",
    nmax = c(yx = Inf, z = Inf)
  ) # no binning
  mf <- model.frame(d)
  collector <- new_collector()

  ## Node trafo: the pooled single-group IGSCA fit. The partition testers
  ## work from model$object (each candidate partition refits a 2-group MGA);
  ## no casewise estfun exists on this path, so objfun (casewise SSR) is NA.
  ## A failed node fit makes the node a leaf and is counted -- recorded,
  ## never redrawn.
  trafo <- function(subset, weights, info = NULL, estfun = TRUE, object = TRUE, ...) {
    was_root <- !collector$root_seen
    collector$root_seen <- TRUE
    ft <- try_fit(mf[subset, indicators, drop = FALSE], model)
    if (!ft$ok) {
      if (was_root) {
        collector$n_fail_full <- 1L
      } else {
        collector$n_fail_node <- collector$n_fail_node + 1L
      }
      return(list(
        coefficients = NA_real_,
        objfun = Inf,
        object = NULL,
        estfun = NULL,
        converged = FALSE,
        nobs = length(subset)
      ))
    }
    list(
      coefficients = NA_real_,
      objfun = NA_real_,
      object = ft$fit,
      estfun = NULL,
      converged = TRUE,
      nobs = length(subset)
    )
  }

  ctrl <- control
  ctrl$model <- model
  ctrl$indicators <- indicators
  ctrl$collector <- collector

  nodes <- grow_extree(d, trafo, selector, splitter, ctrl)
  fitted <- data.frame(
    "(fitted)" = partykit::fitted_node(nodes, mf),
    "(weights)" = rep.int(1L, nrow(mf)),
    check.names = FALSE
  )
  ret <- partykit::party(
    nodes,
    data = mf,
    fitted = fitted,
    info = list(model = model, covariates = covariates)
  )
  ret$terms <- d$terms$all
  class(ret) <- c("igsca_tree", "constparty", class(ret))
  attr(ret, "igsca_info") <- list(
    n_fail_full = collector$n_fail_full,
    n_fail_node = collector$n_fail_node,
    n_fail_resample = collector$n_fail_resample,
    root_criteria = root_criteria(ret)
  )
  ret
}

# ---- COIN family: native partykit::ctree() --------------------------------

#' Casewise sum of squared GSCA residuals (n x 1) -- COIN_ssr influence.
influence_ssr <- function(E) matrix(rowSums(E^2), ncol = 1L)

#' Casewise squared-residual matrix (n x q, multivariate) -- COIN_mat.
influence_mat <- function(E) E^2

#' Grow one COIN-family IGSCA tree: partykit::ctree() with an IGSCA ytrafo.
#' `influence` is a plain function E -> numeric matrix (one row per node
#' observation): influence_ssr or influence_mat. Variable selection is
#' always partykit's (libcoin quadratic linear statistics; permutation
#' p-values under coin_distribution = "approximate" via testtype MonteCarlo).
#'
#' Split-point selection: with `splitter = NULL` (default) it is ctree's
#' native maxstat scan over the selected variable. Passing a partition
#' kernel from igsca_testers.R (split_max_fitdiff, split_max_dgi,
#' split_max_dli -- FIT or a model-implied indicator-VCV distance; the
#' construct VCV is never used in Study 1) instead places the cutpoint at
#' the maximum of that model-comparison statistic -- COIN variable
#' selection mixed with a partition split point. The mixed path refits one
#' 2-group MGA per candidate (numeric grids capped at control$max_cuts,
#' minbucket enforced by candidate_partitions) and requires igsca_testers.R
#' to be sourced; failed candidate fits are counted into n_fail_resample.
igsca_ctree <- function(
  data,
  model,
  covariates,
  influence = influence_ssr,
  splitter = NULL,
  control = igsca_tree_control()
) {
  indicators <- setdiff(names(data), covariates)
  fml <- stats::as.formula(paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(covariates, collapse = " + ")
  ))
  collector <- new_collector()

  ## ctree's function-form ytrafo contract: the outer formals must be named
  ## (data, weights, control) -- they shadow igsca_ctree's own data/control,
  ## which are not needed inside; model/indicators/influence/collector are
  ## captured. estfun must span ALL rows (libcoin indexes it by the full-data
  ## subset), so the node influence is scattered into an all-zero matrix.
  ## calculateGSCAErrors can throw even on a converged fit; either failure
  ## makes the node a leaf and is counted -- recorded, never redrawn.
  ytrafo <- function(data, weights, control) {
    mf <- model.frame(data)
    function(subset, weights, info = NULL, estfun = TRUE, object = TRUE, ...) {
      was_root <- !collector$root_seen
      collector$root_seen <- TRUE
      ft <- try_fit(mf[subset, indicators, drop = FALSE], model)
      E <- if (ft$ok) {
        tryCatch(cSEM:::calculateGSCAErrors(ft$fit), error = function(e) NULL)
      } else {
        NULL
      }
      if (is.null(E)) {
        if (was_root) {
          collector$n_fail_full <- 1L
        } else {
          collector$n_fail_node <- collector$n_fail_node + 1L
        }
        return(list(
          estfun = NULL,
          converged = FALSE,
          objfun = Inf,
          object = NULL,
          nobs = length(subset)
        ))
      }
      h <- influence(E)
      ef <- matrix(0, nrow = nrow(mf), ncol = ncol(h))
      ef[subset, ] <- h
      ## object is returned unconditionally: a mixed-pair splitter's kernel
      ## needs the pooled node fit (partition_stat reads model$object), and
      ## extree does not promise object = TRUE on the split path.
      list(
        estfun = ef,
        converged = TRUE,
        objfun = sum(E^2),
        object = ft$fit,
        nobs = length(subset)
      )
    }
  }

  ## igsca_tree_control() -> ctree_control() mapping. ctree derives its
  ## bonferroni flag from testtype ("Bonferroni" %in% testtype), but
  ## .extree_node reads ctrl$bonferroni independently of how the p-values
  ## were computed. ctree_control(testtype = "MonteCarlo")$bonferroni already
  ## defaults to FALSE, so under the registry's default config this override
  ## is a no-op; it only matters for the (currently unexercised)
  ## bonferroni = TRUE + coin_distribution = "approximate" combination.
  testtype <- if (control$coin_distribution == "approximate") {
    "MonteCarlo"
  } else if (isTRUE(control$bonferroni)) {
    "Bonferroni"
  } else {
    "Univariate"
  }
  cc <- partykit::ctree_control(
    teststat = "quadratic",
    splitstat = "quadratic",
    testtype = testtype,
    nresample = control$R_test,
    alpha = control$alpha,
    minsplit = control$minsplit,
    minbucket = control$minbucket,
    maxdepth = control$maxdepth,
    maxsurrogate = 0L,
    nmax = c(yx = Inf, z = Inf),
    saveinfo = TRUE
  )
  cc$bonferroni <- isTRUE(control$bonferroni)

  ## Mixed pair: keep native selection, replace only the split-point search.
  ## ctree() hands extree_fit its splitfun via the control object
  ## (extree_fit(splitfun = ctrl$splitfun, ...)), so planting ours here is
  ## all that is needed (verified against partykit 1.2.28). The kernel's
  ## config rides along in the same object, exactly as in grow_extree(); the
  ## scan cache never arms on this path (no partition selector), so
  ## argmax_split always runs its own scan.
  if (!is.null(splitter)) {
    cc$model <- model
    cc$indicators <- indicators
    cc$collector <- collector
    cc$max_cuts <- control$max_cuts
    cc$splitfun <- function(
      model,
      trafo,
      data,
      subset,
      weights,
      whichvar,
      ctrl
    ) {
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
    cc$svsplitfun <- cc$splitfun # never called (maxsurrogate = 0)
  }

  ret <- partykit::ctree(fml, data = data, ytrafo = ytrafo, control = cc)
  class(ret) <- c("igsca_tree", class(ret))
  attr(ret, "igsca_info") <- list(
    n_fail_full = collector$n_fail_full,
    n_fail_node = collector$n_fail_node,
    ## 0 on the native path (COIN resampling lives inside libcoin); counts
    ## failed candidate MGA fits when a mixed-pair splitter is in play.
    n_fail_resample = collector$n_fail_resample,
    root_criteria = root_criteria(ret)
  )
  ret
}

# ---- Shared readers / metrics ---------------------------------------------

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

#' Split variables used anywhere in the tree (column names of tree$data).
used_split_vars <- function(tree) {
  ids <- partykit::nodeids(tree)
  sv <- unlist(partykit::nodeapply(tree, ids = ids, FUN = function(nd) {
    s <- partykit::split_node(nd)
    if (is.null(s)) NA_integer_ else partykit::varid_split(s)
  }))
  names(tree$data)[unique(stats::na.omit(sv))]
}

#' Per-replication tree metrics (spec section 5). `covariates` and `true_var`
#' name columns of the data the tree was grown on; noise covariates are all
#' covariates except true_var. root_true/root_noise encode the forced-stump
#' reading from the root criteria matrix extree stores in node info: the
#' stored "criterion" row already contains extree's rank tie-break, so
#' which.max reproduces the engine's own pick exactly. ARI follows
#' Study1_rpart's adj_rand (NA when truth degenerate).
tree_metrics <- function(tree, covariates, true_var = "z_true") {
  info <- attr(tree, "igsca_info")
  if (info$n_fail_full == 1L) {
    return(c(
      any_split = NA_real_,
      true_sel = NA_real_,
      noise_sel = NA_real_,
      root_true = NA_real_,
      root_noise = NA_real_,
      n_leaves = NA_real_,
      depth = NA_real_,
      ari = NA_real_,
      n_fail_full = 1,
      n_fail_node = info$n_fail_node,
      n_fail_resample = info$n_fail_resample
    ))
  }
  noise_vars <- setdiff(covariates, true_var)
  used <- used_split_vars(tree)
  leaves <- stats::predict(tree, type = "node")
  truth <- tree$data[[true_var]]

  crit <- info$root_criteria
  root_pick <- if (is.null(crit) || ncol(crit) == 0L ||
                   !any(is.finite(crit["criterion", ]))) {
    NA_character_
  } else {
    colnames(crit)[which.max(crit["criterion", ])]
  }

  ari <- if (length(unique(truth)) < 2L) {
    NA_real_
  } else {
    mclust::adjustedRandIndex(truth, leaves)
  }

  c(
    any_split = as.numeric(length(unique(leaves)) > 1L),
    true_sel = as.numeric(true_var %in% used),
    noise_sel = as.numeric(any(noise_vars %in% used)),
    root_true = if (is.na(root_pick)) {
      NA_real_
    } else {
      as.numeric(root_pick == true_var)
    },
    root_noise = if (is.na(root_pick)) {
      NA_real_
    } else {
      as.numeric(root_pick %in% noise_vars)
    },
    n_leaves = as.numeric(partykit::width(partykit::node_party(tree))),
    # partykit does not export the `depth` generic (it lives in grid and
    # partykit only registers depth.partynode); call it on the root node.
    depth = as.numeric(grid::depth(partykit::node_party(tree))),
    ari = ari,
    n_fail_full = 0,
    n_fail_node = info$n_fail_node,
    n_fail_resample = info$n_fail_resample
  )
}

## Methods: print falls through to constparty; coef/plot follow
## lavaan.R's my_semtree, reading the trafo payload from info_node().
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
