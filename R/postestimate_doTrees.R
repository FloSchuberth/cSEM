#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#'
#' @inheritParams csem_arguments
#' @inheritParams csem
#' @return cSEM Tree model of class `modelparty` and `party`.
#' @export
#' @importFrom partykit ctree_control ctree
#' @references
#'   \insertAllCited{}
doTrees <- function(
  data,
  model,
  covariates,
  influence = c("mat", "vec", "FIT", "DLi", "DGi"),
  splitter = c("native", "FIT", "DLi", "DGi"),
  control = igsca_tree_control()
) {
  # Preparation
  influence <- match.arg(influence)
  splitter <- match.arg(splitter)

  split_fn <- switch(
    splitter,
    native = NULL,
    FIT = split_max_fitdiff,
    DLi = split_max_dli,
    DGi = split_max_dgi
  )



  indicators <- parseModel(model)$indicators
  fml <- paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(covariates, collapse = " + ")
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
        ft <- try_fit(mf[subset, indicators, drop = FALSE], model)
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
    # TODO: Figure out the different test types and what it has to do with the coin_distribution and bonferroni
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
    # TODO: Return to this bonferroni problem
    cc$bonferroni <- isTRUE(control$bonferroni)

    if (!is.null(split_fn)) {
      # Sets the splitter to one of the three functions that we're looking for.
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
      info = list(model = model, covariates = covariates)
    )
    ret$terms <- d$terms$all
    class(ret) <- c("igsca_tree", "constparty", class(ret))
    warn_dead_splitter(collector, splitter)
    attr(ret, "igsca_info") <- list(
      n_fail_full = collector$n_fail_full,
      n_fail_node = collector$n_fail_node,
      n_fail_resample = collector$n_fail_resample,
      n_split_scan = collector$n_split_scan,
      n_fail_split = collector$n_fail_split,
      root_criteria = root_criteria(ret)
    )
    return(ret)
  }
}