#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#' Grows a tree that partitions the rows `.object` was fitted on into subgroups
#' whose model estimates differ. Every node is refit by replaying `.object`'s own
#' [csem()] arguments, so the estimator, modes and convergence settings of the
#' tree are those of the fit it was given.
#'
#' `.influence` selects the node statistic that conditional-inference variable
#' selection runs on: `"mat"` the casewise GSCA squared-residual matrix.
#' `.splitter` then chooses the cutpoint within the selected covariate -- `"native"`
#'  is partykit's own maxstat scan, while `"FIT"`, `"DLi"` and `"DGi"` replace it with
#'  a deterministic argmax over candidate partitions of a model-comparison statistic.
#'
#' @param .object A single-group `cSEMResults` object, as returned by [csem()].
#'   Its data must contain the `.covariates` columns; [csem()] ignores
#'   non-indicator columns, so they can simply ride along in the original call.
#' @param .covariates Character vector of columns of `.object`'s data to
#'   partition on.
#' @param .influence Node statistic driving variable selection. Currently only "mat" is supported.
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
  # Preparation
  influence <- match.arg(.influence)
  splitter <- match.arg(.splitter)
  control <- .control

  args <- csem_tree_args(.object)  # Safely fetches .object$Information$Arguments
  data <- as.data.frame(args$.data)
  indicators <- parseModel(args$.model)$indicators
  validate_tree_input(data, indicators, .covariates, influence, splitter, args)
  fml <- paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(.covariates, collapse = " + ")
  ) |>
    stats::as.formula()

  # Environment for collecting  metrics on different kinds of convergence failures found while fitting the tree
  collector <- new_collector()

  # debug(libcoin::LinStatExpCov)
  influence_fn <- switch(influence, # from .influence character string
     mat = influence_mat
     # Commented out because it's not as powerful as the matrix approach
     # , vec = influence_vec
    ) 
  # TODO: Walk through how ctree, extree_fit and .extree_node modify ytrafo
  #  How `subset` reaches the inner closure, i.e. what partykit does to what we hand it:
  #   ctree()        - our formals are (data, weights, control), so it calls us as a factory
  #   extree_fit()   - the returned closure has all of (subset, weights, info, estfun,
  #                    object), so it is used verbatim rather than wrapped
  #   .extree_node() - calls it once per node as trafo(subset = <that node's row indices>,
  #                    weights = <global, unchanged>, ...), recursing on subset[kidids == k]
  # `subset` is therefore the only argument that identifies the node; `weights` is the same
  # vector at every depth. The `weights > 0` idiom in ?ctree belongs to the older
  # (y, x, offset, weights, start) interface, which extree_fit() *does* wrap, re-encoding
  # node membership as a weight vector via libcoin::ctabs().
  ytrafo <- function(data, weights, control) {
    # Compare against `partykit::ctree_control()$selectfun`
    # TODO: Uncomment and see in browser
    # browser()
    mf <- model.frame(data)
    function(
      subset,
      weights,
      info = NULL,
      estfun = TRUE,
      object = TRUE,
      ...
    ) {
      # TODO: Uncomment and see in browser
      # browser()
      # We never pass weights to ctree(), so this is a guard against a future change in partykit:
      # try_fit() below ignores them, which would silently drop case weights.
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
      ef[subset, ] <- h # In case I have forgotten why subset is here, note that the ?ctree does the same in its example function
      list(
        estfun = ef, # Yes, this is what downstream `libcoin::LinStatExpCov()` will work with
        converged = TRUE,
        objfun = sum(E^2), # This is equivalent to the objective function of Equation 9 of Hwang et al. (2021, IGSCA): sum(calculateGSCAErrors(ft$fit)^2). AFAIK, this is just used later for prrinting, it doesn't affect the functionality of anything.
        object = ft$fit,
        nobs = length(subset)
      )
    }
  }
  
  # See ?partykit::ctree_control() for the `testtype` argument and how it's internally handled for why Bonferonni is paired with Monte Carlo
  testtype <- c( # TODO: Needs to be compared against ctree and ctree_control for how they handle this. 
    if (isTRUE(control$bonferroni)) "Bonferroni", # Bonferonni corrected p-values
    if (control$coin_distribution == "approximate") "MonteCarlo" # Use permutation resamples to generate null distribution
  )
  if (is.null(testtype) && control$coin_distribution == "asymptotic") { # testtype is null if the above don't fire, meaning testtype <- c()
    testtype <- "Univariate" # This means that you use the asymptotic null distribution of the permutation LinStatExpCov statistic
  }
  cc <- partykit::ctree_control(
    teststat = "quadratic", # TODO: Consider maximum? 
    splitstat = "quadratic", # TODO: Consider maximum?
    testtype = testtype, # TODO: Document
    nresample = control$R_test, # Number of resamples for the the permutation test
    alpha = control$alpha,
    minsplit = control$minsplit,
    minbucket = control$minbucket,
    minprob = control$minprob,
    maxdepth = control$maxdepth,
    maxsurrogate = 0L,
    nmax = c(yx = Inf, z = Inf),
    saveinfo = TRUE
  )

  split_fn <- switch(
    splitter, # From .splitter
    native = NULL, # See partykit::ctree_control()$splitfun
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )

  if (!is.null(split_fn)) {
    # Sets the splitter to one of the three non-native functions that we're looking for.
    cc$args <- args
    cc$indicators <- indicators
    cc$collector <- collector
    cc$splitfun <- function(
      model,
      trafo,
      data,
      subset,
      weights,
      whichvar,
      ctrl
    ) {
      # .extree_node() hands `subset` and `weights` to splitfun just as it does to ytrafo;
      # argmax_split() partitions on `subset` alone, so weighted rows would be ignored.
      stopifnot("case weights are not supported by doTrees(). partykit has changed since initial doTrees development, please report to developers." = length(weights) == 0L)
      argmax_split(
        split_fn,
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
  # The main workhors
  ret <- partykit::ctree(formula = fml, data = data, ytrafo = ytrafo, control = cc)


# Return output ----------------------------------------------------------


  class(ret) <- c("igsca_tree", class(ret))
  warn_dead_splitter(collector, splitter)

  
  ## TODO: Maybe clarify this comment depending on whether it's about every node or just the terminal node?
  #  partykit never calls the trafo at a node it does not try to split, so
  ## leaves arrive with info = NULL and carry no IGSCA fit. Refit them.
  # TODO: Double check that this is efficient and only fits the nodes that it needs to
  ret <- attach_leaf_fits(ret, ret$data, args, indicators, collector)

  ## TODO: Consider deleting this comment
  ##  MEMORY OPTIMIZATION (opt-in): each inner node's info holds a full
  ## cSEMResults object, so a maxdepth = 3 tree keeps up to 7 of them.
  ## Uncomment the next line to drop them; leaf fits are unaffected.
  ##   ret <- drop_inner_node_objects(ret)
  ## Do NOT instead delete `object = ft$fit` from the ytrafo above: a
  ## non-native splitter reads model$object while the tree is still growing,
  ## so removing it there silently disables every non-native splitter.

  attr(ret, "igsca_info") <- list(
    n_fail_full = collector$n_fail_full,
    n_fail_node = collector$n_fail_node,
    ## 0 under .splitter = "native"; counts candidate partitions whose
    ## statistic could not be computed under any other splitter.
    n_fail_candidate = collector$n_fail_candidate,
    ## Split-kernel scans: n_fail_split counts the scans that produced no
    ## finite statistic. Both stay 0 when splitter = "native".
    n_split_scan = collector$n_split_scan,
    n_fail_split = collector$n_fail_split,
    ## Failed IGSCA refits at terminal nodes (attach_leaf_fits).
    n_fail_leaf = collector$n_fail_leaf,
    root_criteria = root_criteria(ret) # TODO: What is this?
  )
  return(ret)
}