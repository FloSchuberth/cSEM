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
  influence = influence_ssr,
  splitter = NULL,
  control = igsca_tree_control()
) {

  # Preparation
  indicators <- parseModel(model)$indicators
  fml <- paste(
    paste(indicators, collapse = " + "),
    "~",
    paste(covariates, collapse = " + ")
  ) |>
    stats::as.formula()

  collector <- new_collector()

  # Conditional Tree Route --------------
  # TODO: Implement igsca_ctree
  ytrafo <- function(data, weights, control) {
    mf <- model.frame(data)
    function(subset, weights, info = NULL, estfun = TRUE, object = TRUE, ...) {
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
  # TODO: Figure out the different test types and what it has to do with the coin_distribution and bonferroni
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
  # TODO: Return to this bonferroni problem
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
  # TODO: Do I need to pass an explicit converged function?
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
 return(ret)

  # TODO: Implement igsca_tree and add conditional switch 
  

  # class(tree) <- c(class(tree), "cSEMResults")

  # return(tree)
}