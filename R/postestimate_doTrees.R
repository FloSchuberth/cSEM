#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#' Grows a tree that partitions the rows `.object` was fitted on into subgroups
#' whose model estimates differ. Every node is refit by replaying `.object`'s own
#' [csem()] arguments, so the estimator, modes and convergence settings of the
#' tree are those of the fit it was given.
#'
#' `.influence` chooses the statistic that will be permuted (the conditional-inference or COIN procedure) 
#' to find which covariate is significantly associated with the transformed 'Y' variable. Currently, only
#' `"mat"` is supported, which means that the the transformed 'Y' variable is the casewise GSCA squared-residual matrix. 
#' 
#' `.splitter` then chooses the cutpoint within the selected covariate. `"native"` is the `partykit::ctree()` default that
#'  chooses the data partitions that maximize the influence statistic. `"FIT"`, `"DLi"` and `"DGi"` behave similarly,
#'  but fit multigroup GSCA models on every possible data partition and maximize the difference in either `"FIT"`, `"DLi"` or `"DGi"`
#'  between the pooled vs multigroup GSCA models.
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
    stats::as.formula()# See partykit::ctree_control()$splitfun

  # Environment for collecting  metrics on different kinds of convergence failures found while fitting the tree
  collector <- new_collector()

  influence_fn <- switch(influence, # from .influence character string
     mat = influence_mat #,
     # vec = influence_vec # Uncomment to see the vector influence function. But it is not that great. Code left here to show how future influence functions may be added
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
      # You can get a sense of ctree's trafo by running the following (from ?ctree)
      # airq <- subset(airquality, !is.na(Ozone))
      # airct <- ctree(Ozone ~ ., data = airq) 
      # airct_list<-unclass(airct)
      # airct_list$trafo
      # Relatedly, you can get a sense of the influence function for by looking-at/debugging-through partykit:::.y2infl()
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
        estfun = ef, # This is what downstream `libcoin::LinStatExpCov()` will work with. You can verify that this specific return item (estfun) is what's important by looking at airct_list$trafo (see the above commented code)
        converged = TRUE,
        objfun = sum(E^2), # This is equivalent to the objective function of Equation 9 of Hwang et al. (2021, IGSCA): sum(calculateGSCAErrors(ft$fit)^2). AFAIK, this is just used later for prrinting, it doesn't affect the functionality of anything.
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
    maxsurrogate = 0L, # In the case of missing data in a covariate this is relevant. This is not handled in our function
    nmax = c(yx = Inf, z = Inf), # Set to default, but this would affect whether or not the covariates or influence function are binned to lower computation costs.
    saveinfo = TRUE,
    update = TRUE # TRUE by default because ytrafo is a function, but does not necessarily refit to terminal nodes
  )

  split_fn <- switch(
    splitter, # From .splitter
    native = NULL, 
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )
  # The relevant code for the native split_fn path is quite deep and written in C. To get a glimpse of it you'd have to run the following code.
  # It's better to consult Section 4.2 "Splitting criteria" in `vignette("ctree", package = "partykit")`. Basically, by-default, the standardized quadratic linear statistic (c_quad) (returned by LinStatExpCov) is computed for all possible subsets of the data,  and libcoin will give you the index for where to split along the covariate in-order to maximize c_quad
  # debug(partykit:::.split)
  # airq <- subset(airquality, !is.na(Ozone))
  # airct <- ctree(Ozone ~ ., data = airq)
  # When you're inside .split run
  # debug(FUN)
    # debug(.ctree_test)
      # debug(.ctree_test_1d)
        # debug(.ctree_test_internal)
          # Passes the influence to LinStatExpCov, which goes to doTest, which gives you an index.  
          # debug(doTest)
  
  # Over-write our ctree_control object even further
  if (!is.null(split_fn)) {
    # TODO: Come back to this after I'm done investigating native
    # Sets the splitter to one of the three non-native functions that we're looking for. 
    # Otherwise, we just use the built-in one that ctree uses. See above for more details. 
    cc$args <- args # How to refit the cSEM models
    cc$indicators <- indicators 
    cc$collector <- collector
    cc$splitfun <- function( # See partykit::ctree_control()$splitfun, partykit:::.split() and partykit:::.ctree_test() for what this function needs to accept and return.
      model,
      trafo,
      data,
      subset,
      weights,
      whichvar,
      ctrl
    ) {
      browser()
      stopifnot("case weights are not supported by doTrees(). partykit has changed since initial doTrees development, please report to developers." = length(weights) == 0L)
      argmax_split(
        splitter = split_fn,
        collector = collector,
        model = model,        
        trafo = trafo, # Only ctrl$lookahead reads this, and it needs the very trafo that fits the node -- partykit hands the splitfun the same one.
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
  
  # ytrafo is not necessarily called at terminal nodes, so sometimes terminal nodes have no csem objects. Here, we refit just the terminal nodes.
  ## What determines whether it will or will not have a csem object depends on the reason why the tree stopped growing. (No significant covariate -> ytrafo is called. Sample size smaller that minbucket, etc -> ytrafo is not called, tree is stopped, so no csem model fitted) 
  ret <- attach_leaf_fits(ret, ret$data, args, indicators, collector)

  ## Deletes the fitted csem objects in all non-terminal nodes. Can be helpful to save on RAM and potentially IGSCA forests. Commented out for later
  ##   ret <- drop_inner_node_objects(ret)
  
  attr(ret, "igsca_info") <- list(
    n_fail_full = collector$n_fail_full, # Convergence failure on root
    n_fail_node = collector$n_fail_node, # Convergence failure on a node. Similar to n_fail_leaf, but includes inner nodes.
    n_fail_candidate = collector$n_fail_candidate, # Relevant to non-native splitters (FIT, DLi, DGi)
    n_split_scan = collector$n_split_scan, # Number of scans in a selected covariate
    n_fail_split = collector$n_fail_split,  # Failed to split on non-native splitters (FIT, DLi, DGi)
    n_fail_leaf = collector$n_fail_leaf, # Convergence failure on a terminal node, fitted via attach_leaf_fits
    root_criteria = root_criteria(ret) # Solely on the root node: The different permutation statistics, their p-value and the criterion that was used to decide whether a split should occur.
  )
  return(ret)
}