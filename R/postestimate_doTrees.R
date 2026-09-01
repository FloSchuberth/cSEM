#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#' Grows a tree that partitions the rows `.object` was fitted on into subgroups
#' whose model estimates differ. Every node is refit by replaying `.object`'s own
#' [csem()] arguments, so the estimator, modes and convergence settings of the
#' tree are those of the fit it was given.
#'
#' `.influence` chooses the statistic that will be permuted (the
#' conditional-inference or COIN procedure) to find which covariate is
#' significantly associated with the transformed 'Y' variable. Currently only
#' `"mat"` is supported: the transformed 'Y' is the casewise GSCA
#' squared-residual matrix \insertCite{Hwang2021a}{cSEM}.
#'
#' `.splitter` then chooses the cutpoint within the selected covariate.
#' `"native"` is the [partykit::ctree()] default, which maximizes the influence
#' statistic. `"FIT"`, `"DLi"` and `"DGi"` instead refit a two-group GSCA model
#' on every candidate partition and maximize either the FIT *difference*, or
#' the `"DLi"`/`"DGi"` matrix *distance*, between the pooled and multigroup
#' model-implied covariance matrices.
#' 
#' @param .object A single-group `cSEMResults` object, as returned by [csem()].
#'   Its data must contain the `.covariates` columns; [csem()] ignores
#'   non-indicator columns, so they can simply ride along in the original call.
#' @param .covariates Character vector of columns of `.object`'s data to
#'   partition on.
#' @param .influence Node statistic driving variable selection.
#' @param .splitter Cutpoint rule. One of "native", "FIT", "DLi" or "DGi".
#' @param .control Tuning parameters, see [igsca_tree_control()].
#'
#' @return A tree of class `c("igsca_tree", "constparty", "party")`, carrying the
#'   per-leaf fits (see [coef.igsca_tree()]) and an `igsca_info` attribute of
#'   failure counters and root criteria.
#' @export
#' @importFrom partykit ctree_control ctree
#' @references
#'   \insertAllCited{}
doTrees <- function(
  .object,
  .covariates,
  .influence = "mat",
  .splitter = c("native", "FIT", "DLi", "DGi"),
  .control = igsca_tree_control()
) {
  influence <- match.arg(.influence)
  splitter <- match.arg(.splitter)
  control <- .control

  args <- csem_tree_args(.object) # Safely fetches .object$Information$Arguments
  data <- as.data.frame(args$.data)
  # partition_stat() adds a column named TREE_GROUP_COL ("TREETEMPGROUP") to each
  # candidate node's data and hands it to csem() as `.id`. A column of that name
  # already in the data would be silently overwritten there, so refuse it here
  # where the user can still see which column is in the way.
  stopifnot(
    "`.object` was fitted on data containing a column named TREETEMPGROUP. doTrees() reserves that name for the grouping column it builds at every candidate split; please rename it." =
      !(TREE_GROUP_COL %in% names(data))
  )
  indicators <- parseModel(args$.model)$indicators
  validate_tree_input(data, indicators, .covariates, influence, splitter, args)
  fml <- paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(.covariates, collapse = " + ")
  ) |>
    stats::as.formula() # See partykit::ctree_control()$splitfun

  # Counters for the convergence failures hit while growing the tree
  collector <- new_collector()

  influence_fn <- switch(
    influence,
    mat = influence_mat #,
    # vec = influence_vec  # Extension point: a second influence function plugs in here
  )
  
  ytrafo <- function(data, weights, control) {
    mf <- model.frame(data)
    function(
      subset,
      weights,
      info = NULL,
      estfun = TRUE,
      object = TRUE,
      ...
    ) {
      # To see what ctree hands a trafo:
      #   unclass(partykit::ctree(Ozone ~ ., subset(airquality, !is.na(Ozone))))$trafo
      # partykit:::.y2infl() is the influence function it uses by default.
      stopifnot("case weights are not supported by doTrees(). partykit has changed since initial doTrees development, please report to developers." = length(weights) == 0L)
      was_root <- !collector$root_seen
      collector$root_seen <- TRUE
      ft <- try_fit(mf[subset, indicators, drop = FALSE], args) # Returns a fitted csem model (ft$fit) and whether the model converged (ft$ok)
      E <- if (ft$ok) {
        tryCatch(calculateGSCAErrors(ft$fit), error = function(e) NULL)
      } else {
        NULL
      }
      if (is.null(E)) {
        if (was_root) {
          collector$n_fail_full <- 1L
        } else {
          collector$n_fail_node <- collector$n_fail_node + 1L
        }
        ## partykit drops this node's info entirely, so attach_leaf_fits() would
        ## otherwise repeat the identical failing fit and count it a second time.
        record_failed_subset(collector, subset)
        return(list(
          estfun = NULL,
          converged = FALSE,
          objfun = Inf,
          object = NULL,
          nobs = length(subset)
        ))
      }
      h <- influence_fn(E)
      ef <- matrix(0, nrow = nrow(mf), ncol = ncol(h))
      ef[subset, ] <- h # ?ctree's own example trafo also writes only the subset rows into a full-length estfun
      list(
        estfun = ef, # The one element libcoin::LinStatExpCov() consumes downstream
        converged = TRUE,
        objfun = sum(E^2), # Eq. 9 of Hwang et al. (2021, IGSCA). Reported only; nothing downstream reads it.
        object = ft$fit,
        nobs = length(subset)
      )
    }
  }
  
  # See ?partykit::ctree_control() for the `testtype` argument and how it's internally handled for why Bonferonni is paired with Monte Carlo
  testtype <- c( 
    if (isTRUE(control$bonferroni)) "Bonferroni", # Bonferonni corrected p-values
    if (control$coin_distribution == "approximate") "MonteCarlo" # Use permutation resamples to generate null distribution
  )
  if (is.null(testtype) && control$coin_distribution == "asymptotic") { # testtype is null if the above don't fire, meaning testtype <- c()
    testtype <- "Univariate" # This means that you use the asymptotic null distribution of the permutation LinStatExpCov statistic
  }
  cc <- partykit::ctree_control(
    teststat = "quadratic", # partykit::ctree_control() default
    splitstat = "quadratic", # partykit::ctree_control() default
    testtype = testtype, 
    nresample = control$R_test, # Number of resamples for the the permutation test
    alpha = control$alpha,
    minsplit = control$minsplit,
    minbucket = control$minbucket,
    minprob = control$minprob,
    maxdepth = control$maxdepth,
    maxsurrogate = 0L, # Surrogate splits for missing covariate values; unsupported here (cf. cc$svsplitfun below)
    nmax = c(yx = Inf, z = Inf), # Default: no binning of covariates or influence values
    saveinfo = TRUE,
    update = TRUE, # TRUE by default because ytrafo is a function, but does not necessarily refit to terminal nodes
    lookahead = FALSE, # Not supported by IGSCA trees
    intersplit = FALSE # Not supported by IGSCA trees; candidate_partitions() always breaks at the observed value
  )

  split_fn <- switch(
    splitter, # From .splitter
    native = NULL,
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )
  # Native path: the search is written in C. Section 4.2 "Splitting criteria" in
  # vignette("ctree", package = "partykit") is the readable account -- the
  # standardized quadratic linear statistic (c_quad) from LinStatExpCov is
  # maximized over cutpoints and libcoin returns the index to split at.
  # To watch it: debug(partykit:::.split); partykit::ctree(Ozone ~ ., subset(airquality, !is.na(Ozone)))
  # Inside, the chain is FUN -> .ctree_test -> .ctree_test_1d -> .ctree_test_internal -> doTest.
  
  # Over-write our ctree_control object even further
  if (!is.null(split_fn)) {
    # The FIT/DLi/DGi kernels score a candidate by refitting it as a TWO-group
    # MGA (partition_stat() builds a grouping column with levels 1/2), so they
    # can only ever put a number on a binary split. partykit's multiway rule
    # asks for one kid per factor level instead, and argmax_split() would answer
    # with multiway_split(), which no kernel ever scores.
    stopifnot(
      "multiway splitting is not supported by the 'FIT', 'DLi' and 'DGi' splitters, which score two-group refits only. partykit has changed since initial doTrees development, please report to developers." =
        !isTRUE(cc$multiway)
    )
    cc$args <- args # How to refit the cSEM models
    cc$indicators <- indicators
    cc$collector <- collector
    cc$splitfun <- function( # See partykit::ctree_control()$splitfun -> partykit:::.split() + partykit:::.ctree_test() for what this function needs to accept and return.
      model,
      trafo,
      data,
      subset,
      weights,
      whichvar,
      ctrl
    ) {
      stopifnot("case weights are not supported by doTrees(). partykit has changed since initial doTrees development, please report to developers." = length(weights) == 0L)
      argmax_split(
        splitter = split_fn,
        collector = collector,
        model = model,        
        trafo = trafo, # Accepted for partykit's splitfun contract but unread: argmax_split() dropped it with lookahead support.
        mf = model.frame(data),
        subset = subset,
        whichvar = whichvar,
        ctrl = ctrl,
        weights = weights
      )
    }
    cc$svsplitfun <- cc$splitfun # never called because (maxsurrogate = 0)
  }
  # The main workhorse
  ret <- partykit::ctree(formula = fml, data = data, ytrafo = ytrafo, control = cc)


  # Return output ----------------------------------------------------------

  class(ret) <- c("igsca_tree", class(ret))

  warn_dead_splitter(collector, splitter) # For FIT, DLi or DGi paths. Just throws a warning if all the scans at the candidate split points don't work
  
  # partykit only calls ytrafo at nodes it attempts to split, so a node stopped
  # by minbucket/maxdepth carries no csem object, while one stopped by a failed
  # test does. Refit only the leaves that need it.
  ret <- attach_leaf_fits(ret, ret$data, args, indicators, collector)

  ## Not wired up: drop_inner_node_objects(ret) frees the inner-node
  ## cSEMResults (RAM, forests) at the cost of node inspection.
  ##   ret <- drop_inner_node_objects(ret)
  
  attr(ret, "igsca_info") <- list(
    n_fail_full = collector$n_fail_full, # Convergence failure on root
    n_fail_node = collector$n_fail_node, # Convergence failure on nodes below the root. Similar to n_fail_leaf, but includes inner nodes.
    n_fail_candidate = collector$n_fail_candidate, # Relevant to non-native splitters (FIT, DLi, DGi)
    n_split_scan = collector$n_split_scan, # Number of scans in a selected covariate
    n_fail_split = collector$n_fail_split,  # Failed to split on non-native splitters (FIT, DLi, DGi)
    n_fail_leaf = collector$n_fail_leaf, # Convergence failure on a terminal node, fitted via attach_leaf_fits
    root_criteria = root_criteria(ret) # Solely on the root node: The different permutation statistics, their p-value and the criterion that was used to decide whether a split should occur.
  )
  return(ret)
}