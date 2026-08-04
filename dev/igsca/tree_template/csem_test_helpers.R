# Shared cSEM test machinery ----------------------------------------------
#
# Single source of truth for the model-fitting chokepoint, the resampling
# index generators, and the node-level test statistics shared by:
#   * R/MGA/BootBal_Defaults/  (Study 0 -- sources this file from its Sim.R)
#   * R/tree/igscaTree/        (IGSCA tree engine + testers, Study 1/2)
# Definitions only -- no library() calls, no side effects. All external
# functions are namespace-qualified (cSEM::, SimDesign::); stats/sd/pt come
# from base R. Index-generator equivalence to boot is proved in
# R/MGA/BootBal_Defaults/resample_index_proof.R.
# Extracted verbatim from BootBal_Defaults_Sim.R (see
# docs/superpowers/specs/2026-07-12-study1-igsca-tree-unbiased-selection-design.md §1).

# cSEM fitting helpers ---------------------------------------------------

#' Single chokepoint replacing every direct cSEM::csem() call in this study.
#' `.id = NULL` yields the pooled single-group fit; `.id = "group"` the MGA fit.
fit_csem <- function(.data, .model, .id = NULL) {
  cSEM::csem(
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

#' TRUE when every admissibility check passes. For a multi-group (MGA) fit cSEM
#' returns a list of per-group results; verify() must pass for all of them.
#' Reuses the sum(cSEM::verify(x)) == 0 pattern from
#' R/MGA/kleselSource/klesel2022/Simulation.R:237.
csem_converged <- function(fit) {
  if (inherits(fit, "cSEMResults_multi")) {
    all(vapply(fit, function(x) sum(cSEM::verify(x)) == 0, logical(1)))
  } else {
    sum(cSEM::verify(fit)) == 0
  }
}

#' Never throws, never warns to the caller; returns the fit and an `ok` flag.
#' cSEM's inadmissibility warnings are suppressed here so they can never be
#' promoted to errors by an ambient `options(warn = 2)` / rlang::global_handle().
try_fit <- function(.data, .model, .id = NULL) {
  fit <- suppressWarnings(tryCatch(
    fit_csem(.data, .model, .id),
    error = function(e) NULL
  ))
  ok <- !is.null(fit) &&
    isTRUE(tryCatch(csem_converged(fit), error = function(e) FALSE)) # That parameter estimates are proper solutions is the same bar that Hwang et al. (2021) IGSCA paper put
  list(fit = fit, ok = ok)
}

# Resampling index generators --------------------------------------------
#
# Each generator returns an R x n integer matrix. Row r holds the row indices
# drawn for resample r, so the orientation matches
# boot::boot.array(fit, indices = TRUE). Generators are pure functions of
# (n, R, strata) and the RNG state and carry no boot dependency; their
# equivalence to boot's ordinary.array / balanced.array / permutation.array is
# proved standalone in resample_index_proof.R.
#
# `strata` (when supplied) is any vector of length n whose distinct values mark
# the groups within which resampling is confined; NULL means a single stratum.
# 
#' Ordinary resampling: sample WITH replacement, within each stratum.
#' Near verbatim copy of boot:::ordinary.array: the columns belonging to a stratum are filled
#' with draws (with replacement) from that stratum's own rows, so each replicate
#' contains exactly the stratum's original count of rows from every stratum.
#' 
#' debug(idx_ordinary)
#' idx_ordinary(10, 5, strata = NULL) # Ordinary bootstrap
#' idx_ordinary(10, 5, strata = rep(c("A", "B"), times = c(7,3))) # Ordinary + Stratified bootstrap
#' undebug(idx_ordinary)
idx_ordinary <- function(n, R, strata = NULL) {
  out <- matrix(NA_integer_, nrow = R, ncol = n)
  if (is.null(strata)) {
    out[] <- sample.int(n, n * R, replace = TRUE)
  } else {
    strata <- as.integer(as.factor(strata))
    for (s in unique(strata)) {
      gp <- which(strata == s)
      out[, gp] <- if (length(gp) == 1L) {
        rep.int(gp, R)
      } else {
        matrix(gp[sample.int(length(gp), R * length(gp), replace = TRUE)],
               nrow = R, ncol = length(gp))
      }
    }
  }
  out
}

#' First-order balanced resampling: each original row is used EXACTLY R times
#' in total across the R replicates. Mirrors/Near copy of boot:::balanced.array: within each
#' stratum, R copies of the stratum's rows are permuted once and dealt out
#' across replicates. With a single stratum this is the classic
#' "permute R copies of 1:n, chop into R blocks of length n".
#' 
#' debug(idx_balanced)
#' idx_balanced(10, 5, strata = NULL) # Balanced bootstrap
#' idx_balanced(10, 5, strata = rep(c("A", "B"), times = c(7,3))) # Balanced + Stratified bootstrap
#' undebug(idx_balanced)
idx_balanced <- function(n, R, strata = NULL) {
  if (is.null(strata)) strata <- rep.int(1L, n)
  strata <- as.integer(as.factor(strata))
  # column c (n x R layout) holds replicate c's indices; identity to start
  # This is what is done in Section 9.2 of Davidson & Hinkley
  out <- matrix(rep.int(seq_len(n), R), nrow = n, ncol = R)
  for (s in unique(strata)) {
    gp <- which(strata == s)
    if (length(gp) > 1L) {
      out[gp, ] <- matrix(sample(rep.int(gp, R)), nrow = length(gp), ncol = R)
    }
  }
  t(out)
}

#' Permutation resampling: sample WITHOUT replacement, within each stratum, so
#' every replicate is a permutation of its stratum (each original row used
#' exactly once per replicate). Mirrors/Near-Copy of boot:::permutation.array.
idx_permutation <- function(n, R, strata = NULL) {
  # row r starts as the identity 1:n
  out <- matrix(rep(seq_len(n), each = R), nrow = R, ncol = n)
  if (is.null(strata)) {
    for (r in seq_len(R)) out[r, ] <- sample.int(n)
  } else {
    strata <- as.integer(as.factor(strata))
    for (s in unique(strata)) {
      gp <- which(strata == s)
      if (length(gp) > 1L) {
        for (r in seq_len(R)) out[r, gp] <- sample(gp)
      }
    }
  }
  out
}

#' Full-data / per-resample HT statistic vector. Returns NULL when either the
#' pooled or the MGA fit fails to converge (so the caller can record a failure
#' and drop this fit from the distribution), else the four (A)FIT values.
ht_stat <- function(.data, .model) {
  single <- try_fit(.data, .model)
  mga    <- try_fit(.data, .model, .id = 'group')
  if (!single$ok || !mga$ok) return(NULL)
  c(
    singleFIT  = cSEM::calculateFIT(single$fit),
    singleAFIT = cSEM::calculateAFIT(single$fit),
    groupFIT   = cSEM::calculateFIT(mga$fit),
    groupAFIT  = cSEM::calculateAFIT(mga$fit)
  )
}

#' Internal resampling loop for the HT test. `idx` is an R x n index matrix from
#' one of the generators. Returns the R x k replicate-statistic matrix (rows for
#' non-converged resamples left as NA) and the resample-failure count.
#'
#' `replace_failures` / `max_replace` are the designed-but-not-enabled
#' replacement accommodation (design doc section 5). Only the FALSE path is
#' implemented: a failed resample is recorded as NA and counted. To enable
#' replacement later, on a failure draw a fresh index row from `regen` and refit
#' up to `max_replace` times, tracking an `n_replaced` counter -- localized here,
#' no other code changes.
ht_resample <- function(.data, .model, idx, template,
                        replace_failures = FALSE, max_replace = 0L,
                        regen = NULL) {
  R <- nrow(idx)
  t_ast <- matrix(NA_real_, nrow = R, ncol = length(template),
              dimnames = list(NULL, names(template)))
  n_fail_resample <- 0L
  for (r in seq_len(R)) {
    s <- ht_stat(.data[idx[r, ], , drop = FALSE], .model)
    if (is.null(s)) {
      if (isTRUE(replace_failures)) {
        stop("replace_failures = TRUE is a design stub and is not implemented")
      }
      n_fail_resample <- n_fail_resample + 1L
    } else {
      t_ast[r, ] <- s
    }
  }
  list(t = t_ast, n_fail_resample = n_fail_resample)
}

#' The four pooled-vs-MGA distances for one MGA fit, given the precomputed
#' replicated-block pooled VCVs.
ndt_dists <- function(Sc_pool, Si_pool, mga_fit) {
  Sc_mga <- cSEM:::bdiagFit(mga_fit, .type_vcv = "construct")
  Si_mga <- cSEM:::bdiagFit(mga_fit, .type_vcv = "indicator")
  c(
    DGc = cSEM:::calculateDG(.matrix1 = Sc_pool, .matrix2 = Sc_mga),  # geodesic,  construct
    DGi = cSEM:::calculateDG(.matrix1 = Si_pool, .matrix2 = Si_mga),  # geodesic,  indicator
    DLc = cSEM:::calculateDL(.matrix1 = Sc_pool, .matrix2 = Sc_mga),  # sq-Euclid, construct
    DLi = cSEM:::calculateDL(.matrix1 = Si_pool, .matrix2 = Si_mga)   # sq-Euclid, indicator
  )
}
