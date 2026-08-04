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
# Requires igsca_tree.R sourced first (candidate_partitions) and
# R/MGA/csem_test_helpers.R (try_fit, idx_permutation, ndt_dists).

#' Node data for one candidate partition: children as `group` levels 1/2.
node_group_data <- function(mf, subset, indicators, goes_left) {
  d <- mf[subset, indicators, drop = FALSE]
  cbind(group = factor(ifelse(goes_left, 1L, 2L), levels = c(1L, 2L)), d)
}

#' Pooled model-implied VCVs replicated to 2 blocks (constant per node;
#' cached in the collector under $ndt_pools keyed by the subset).
ndt_pools <- function(single_fit) {
  list(
    Sc = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "construct"),
    Si = cSEM:::bdiagFit(single_fit, .n_blocks = 2L, .type_vcv = "indicator")
  )
}

#' Observed statistic for one candidate partition. `stat_kind` is "FITdiff"
#' (NPT) or an ndt_dists() distance name -- "DGi"/"DLi" only: Study 1 works
#' with the model-implied indicator VCV, never the construct VCV (the shared
#' helper also computes "DGc"/"DLc", which stay unused here). Counts every
#' failed auxiliary fit into collector$n_fail_resample (recorded, never
#' redrawn).
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

#' Scan all candidates of covariate j; returns NULL when nothing admissible.
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

#' One-sided permutation p-value at the argmax partition (post-BootBalVerify
#' Study 0 conventions, 2026-07-13): permute the child labels, null >= obs
#' without abs(). Failed resamples are dropped from the null (finite counts)
#' and counted by partition_stat.
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

#' Shared worker behind the three partition selectors: statistic-first scan
#' per covariate, cache the argmax for the splitfun, permutation p at the
#' argmax only. `matched_split_fn` is the splitter function whose kernel
#' equals this selector's scan statistic -- the cache key grow_extree's
#' splitfun compares against with identical().
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

split_max_fitdiff <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("FITdiff", model, mf, subset, goes_left, ctrl)
}

split_max_dgi <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DGi", model, mf, subset, goes_left, ctrl)
}

split_max_dli <- function(model, mf, subset, goes_left, ctrl) {
  partition_stat("DLi", model, mf, subset, goes_left, ctrl)
}
