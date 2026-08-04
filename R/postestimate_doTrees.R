#' Recursive Partitioning Algorithm for *csem* GSCA Models
#'
#'
#' @inheritParams csem_arguments
#' @inheritParams csem
#' @return cSEM Tree model of class `modelparty` and `party`.
#' @export
#' @importFrom partykit ctree_control
#' @importFrom Formula as.Formula
#' @references
#'   \insertAllCited{}
doTrees <- function(
  .object,
  .splitvars,
  .model,
  .data = .object$Information$Arguments$.data,
  .ctree_control = partykit::ctree_control(),
  .igsca_tree_control = igsca_tree_control(),
  .approach_weights = .object$Information$Arguments$.approach_weights,
  .iter_max = .object$Information$Arguments$.iter_max,
  .tolerance = .object$Information$Arguments$.tolerance,
  .disattenuate = .object$Information$Arguments$.disattenuate,
  .dominant_indicators = .object$Information$Arguments$.dominant_indicators,
  .conv_criterion = .object$Information$Arguments$.conv_criterion
) {
  
  # Warning Checks
  stopifnot(
    'This function only works on single-group models of class "cSEMResults_default" and not cSEMResults_multi' = !inherits(
      .object,
      "cSEMREsults_multi"
    ),
    'This function only supports GSCA models' = .object$Information$Arguments$.approach_weights ==
      "GSCA"
  )

  # Preparation
  ex_formula <- paste(
    paste(.object$Information$Model$indicators, collapse = " + "),
    "~",
    paste(.splitvars, collapse = " + ")
  ) |>
    stats::as.formula()

  partied_dat <- partykit::extree_data(
    formula = ex_formula,
    data = .data,
    yx = "none", # TODO: What is this for?
    nmax = c(yx = Inf, z = Inf) # TODO: What is this for?
  )

  browser()
  # model_frame <- model.frame(partied_dat)

  # TODO: Implement igsca_ctree
  ytrafo <- function(data, weights, control) {
    mf <- model.frame(data)
    function(subset, weights, info = NULL, estfun = TRUE, object = TRUE, ...) {
      # TODO: Add fitting function to cSEM here. 
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

  # TODO: Implement igsca_tree and add conditional switch 
  

  class(tree) <- c(class(tree), "cSEMResults")

  return(tree)
}


#' csem Fitting Function for Interfacing with `partykit::mob()`
#' 
#' @inheritParams doTrees
#' @inheritParams csem_arguments
#' @keywords internal
csem_fit <- function(.object,
                     .model,
                     .approach_weights,
                     .iter_max,
                     .tolerance,
                     .disattenuate,
                     .dominant_indicators,
                     .conv_criterion) {
  # TODO: Rethink the arguments and etc
  # Why an unnamed function within a function?
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = TRUE) {
    
    fitted_model <- csem(
      .model = .model,
      .data = y,
      .approach_weights = .approach_weights,
      .iter_max = .iter_max,
      .tolerance = .tolerance,
      .conv_criterion =  .conv_criterion,
      .disattenuate =  .disattenuate,
      .dominant_indicators = .dominant_indicators
    )
    
## Get output --------------------------------------------------------------

    # Coefficients
    
    out_coef <-tidy(fitted_model)
    res_coef <- out_coef[!(out_coef$op %in% c("Direct_effect", "Indirect_effect", "Total_effect")),"estimate"] 
    names(res_coef) <- out_coef[!(out_coef$op %in% c("Direct_effect", "Indirect_effect", "Total_effect")),"term"]
    
    # Objective Function
    if(.object$Information$Arguments$.approach_weights == "GSCA") {
      res_obj <- calculateIgscaObjectiveFunction(fitted_model)
      
    } else {
      res_obj <- NULL
      
    }
    
          
    return(list(coefficients = res_coef,
         objfun = res_obj, # The minimized objective function
         estfun = if(estfun) {} else NULL,
         object = if(object) fitted_model else NULL))
  }
}

#' Prune a grown tree from doTrees
#'
#' @param .tree Fitted tree
#'
#' @return A pruned tree
prune.cSEMResults <- function(.tree) {
  if (!all(inherits(.tree) %in% c("modelparty", "party", "cSEMResults"))) {
    stop("Please pass a completed tree from cSEM::doTrees().")
  }
  
  # TODO: When pruning with mob(), there's a few approaches.
  # We could make the objective function into FIT itself, then when pruning, we
  # sum the FIT of each model and divide by the number of models. This seems to
  # get the FIT of the block-diagonalized multigroup model for some reason
  # 
  # Alternatively, we could have our own pruning function and routine to get the paired bootstrap t-test.
  # 
  
  return(pruned_tree)
}