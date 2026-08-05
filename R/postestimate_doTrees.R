#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#' Grows a tree that partitions the rows `.object` was fitted on into subgroups
#' whose model estimates differ. Every node is refit by replaying `.object`'s own
#' [csem()] arguments, so the estimator, modes and convergence settings of the
#' tree are those of the fit it was given -- there is nothing to keep in sync and
#' no estimator hard-coded here.
#'
#' `.influence` selects the family: `"mat"` and `"vec"` run conditional-inference
#' (COIN) variable selection on casewise GSCA residuals, `"FIT"`, `"DLi"` and
#' `"DGi"` run permutation tests on a model-comparison statistic. `.splitter`
#' chooses the cutpoint within the selected covariate, and `"native"` (partykit's
#' own maxstat scan) is available only to the COIN family.
#'
#' @param .object A single-group `cSEMResults` object, as returned by [csem()].
#'   Its data must contain the `.covariates` columns; [csem()] ignores
#'   non-indicator columns, so they can simply ride along in the original call.
#' @param .covariates Character vector of columns of `.object`'s data to
#'   partition on.
#' @param .influence Node statistic and, with it, the algorithm family. One of
#'   "mat", "vec", "FIT", "DLi" or "DGi".
#' @param .splitter Cutpoint rule. One of "native", "FIT", "DLi" or "DGi";
#'   "native" is only available when `.influence` is "mat" or "vec".
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
  .influence = c("mat", "vec", "FIT", "DLi", "DGi"),
  .splitter = c("native", "FIT", "DLi", "DGi"),
  .control = igsca_tree_control()
) {
  # Preparation
  influence <- match.arg(.influence)
  splitter <- match.arg(.splitter)
  control <- .control

  split_fn <- switch(
    splitter,
    native = NULL,
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )

  ## Everything the tree needs -- the rows to partition, the covariate columns,
  ## the model, and the estimator settings for every node refit -- comes from
  ## the fit itself, so a node cannot be estimated differently from its root.
  args <- csem_tree_args(.object)
  data <- as.data.frame(args$.data)
  indicators <- parseModel(args$.model)$indicators
  validate_tree_input(data, indicators, .covariates, influence, splitter, args)

  fml <- paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(.covariates, collapse = " + ")
  ) |>
    stats::as.formula()

  collector <- new_collector()


  if (influence %in% c("mat", "vec")) {
    influence_fn <- switch(influence, mat = influence_mat, vec = influence_vec)
    # Conditional Tree Route --------------
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
        was_root <- !collector$root_seen
        collector$root_seen <- TRUE
        ft <- try_fit(mf[subset, indicators, drop = FALSE], args)
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
    ## partykit's `testtype` conflates two orthogonal choices: where the
    ## p-value comes from ("MonteCarlo" = nresample permutations, otherwise the
    ## asymptotic chi-squared limit of the same conditional null) and whether
    ## it is adjusted for multiplicity. `ctree_control()` accepts a length-2
    ## vector to ask for both, so building it this way reaches all four
    ## combinations of our two arguments inside the documented API -- and
    ## `bonferroni = "Bonferroni" %in% testtype` is then derived correctly by
    ## ctree_control() itself, with no post-hoc surgery on the control object.
    testtype <- c(
      if (isTRUE(control$bonferroni)) "Bonferroni",
      if (control$coin_distribution == "approximate") "MonteCarlo"
    )
    if (is.null(testtype)) {
      testtype <- "Univariate"
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

    if (!is.null(split_fn)) {
      # Sets the splitter to one of the three functions that we're looking for.
      cc$args <- args
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
    # TODO: Do I need to pass an explicit converged function?
    ret <- partykit::ctree(fml, data = data, ytrafo = ytrafo, control = cc)
    class(ret) <- c("igsca_tree", class(ret))
    warn_dead_splitter(collector, splitter)

    ## partykit never calls the trafo at a node it does not try to split, so
    ## leaves arrive with info = NULL and carry no IGSCA fit. Refit them.
    ret <- attach_leaf_fits(ret, ret$data, args, indicators, collector)

    ## MEMORY OPTIMIZATION (opt-in): each inner node's info holds a full
    ## cSEMResults object, so a maxdepth = 3 tree keeps up to 7 of them.
    ## Uncomment the next line to drop them; leaf fits are unaffected.
    ##   ret <- drop_inner_node_objects(ret)
    ## Do NOT instead delete `object = ft$fit` from the ytrafo above: a
    ## mixed-pair splitter reads model$object while the tree is still growing,
    ## so removing it there silently disables every non-native splitter.

    attr(ret, "igsca_info") <- list(
      n_fail_full = collector$n_fail_full,
      n_fail_node = collector$n_fail_node,
      ## 0 on the native path (COIN resampling lives inside libcoin); counts
      ## failed candidate MGA fits when a mixed-pair splitter is in play.
      n_fail_resample = collector$n_fail_resample,
      ## Split-kernel scans: n_fail_split counts the scans that produced no
      ## finite statistic. Both stay 0 when splitter = "native".
      n_split_scan = collector$n_split_scan,
      n_fail_split = collector$n_fail_split,
      ## Failed IGSCA refits at terminal nodes (attach_leaf_fits).
      n_fail_leaf = collector$n_fail_leaf,
      root_criteria = root_criteria(ret)
    )
    return(ret)
  } else if (influence %in% c("FIT", "DLi", "DGi")) {
    # NPT, dGi and dLi split -------------------------------------------------

    selector <- switch(
      influence,
      FIT = select_npt,
      DLi = select_ndt_dli,
      DGi = select_ndt_dgi
    )

    stopifnot(
      "splitter should be any one of 'FIT', 'DLi' or 'DGi' when the influence (selector) function is one of 'FIT', 'DLi' or 'DGi'" = splitter %in%
        c("FIT", "DLi", "DGi")
    )

    d <- partykit::extree_data(
      fml,
      data = data,
      yx = "none",
      nmax = c(yx = Inf, z = Inf)
    ) # no binning

    mf <- model.frame(d)
    trafo <- function(
      subset,
      weights,
      info = NULL,
      estfun = TRUE,
      object = TRUE,
      ...
    ) {
      was_root <- !collector$root_seen
      collector$root_seen <- TRUE
      ft <- try_fit(mf[subset, indicators, drop = FALSE], args)
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
    ctrl$args <- args
    ctrl$indicators <- indicators
    ctrl$collector <- collector

    nodes <- grow_extree(d, trafo, selector, split_fn, ctrl)
    fitted <- data.frame(
      "(fitted)" = partykit::fitted_node(nodes, mf),
      "(weights)" = rep.int(1L, nrow(mf)),
      check.names = FALSE
    )
    ret <- partykit::party(
      nodes,
      data = mf,
      fitted = fitted,
      info = list(model = args$.model, covariates = .covariates)
    )
    ret$terms <- d$terms$all
    class(ret) <- c("igsca_tree", "constparty", class(ret))
    warn_dead_splitter(collector, splitter)

    ## As on the ctree path: leaves are never visited by the trafo, and here the
    ## trafo does not compute objfun at all (it is NA_real_ on inner nodes), so
    ## the leaf refit is the only source of a per-node IGSCA fit and SSR.
    ret <- attach_leaf_fits(ret, mf, args, indicators, collector)

    ## MEMORY OPTIMIZATION (opt-in): see the note on the ctree path above.
    ##   ret <- drop_inner_node_objects(ret)
    ## Do NOT instead delete `object = ft$fit` from the trafo above: the
    ## partition selectors read model$object on every candidate partition
    ## (partition_stat), so removing it there breaks the whole branch.

    attr(ret, "igsca_info") <- list(
      n_fail_full = collector$n_fail_full,
      n_fail_node = collector$n_fail_node,
      n_fail_resample = collector$n_fail_resample,
      n_split_scan = collector$n_split_scan,
      n_fail_split = collector$n_fail_split,
      n_fail_leaf = collector$n_fail_leaf,
      root_criteria = root_criteria(ret)
    )
    return(ret)
  }
}