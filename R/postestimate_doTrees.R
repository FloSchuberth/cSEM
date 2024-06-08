#' Model-Based Recursive Partitioning Algorithm for *csem* Models
#' 
#' This implementation of Model-Based Recursive Partitioning for *csem* models
#' closely follows the implementations described for the `lavaan` package
#' \insertCite{zeileis2020AchimZeileis}{cSEM} and the `networktree` package
#' \insertCite{jonesNetworkTreesMethod2020}{cSEM} using `partykit::mob()` as
#' described in hothornzeileis2015J.Mach.Learn.Res and
#' hothornetal2006J.Comput.Graph.Stat and zeileisetal2008J.Comput.Graph.Stat.
#' 
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] 
#' @inheritParams csem_arguments 
#' @inheritParams csem
#' @return cSEM Tree model of class `modelparty` and `party`.
#' @export
#' @importFrom partykit mob mob_control
#' @importFrom Formula as.Formula
#' @references
#'   \insertAllCited{}
doTrees <- function(.object,
                    .splitvars,
                    .model,
                    .maxdepth = Inf,
                    .minsize = nrow(tidy.cSEMResults(.object)),
                    .data = .object$Information$Arguments$.data,
                    .approach_weights = .object$Information$Arguments$.approach_weights,
                    .iter_max = .object$Information$Arguments$.iter_max,
                    .tolerance = .object$Information$Arguments$.tolerance,
                    .disattenuate = .object$Information$Arguments$.disattenuate,
                    .dominant_indicators = .object$Information$Arguments$.dominant_indicators,
                    .conv_criterion = .object$Information$Arguments$.conv_criterion) {
  
## Warning Checks ----------------------------------------------------------
  stopifnot(
    'This function only works on single-group models of class "cSEMResults_default" and not cSEMResults_multi' = !inherits(.object, "cSEMREsults_multi"),
    'This function is not supported for non I-GSCA models' = .object$Information$Arguments$.approach_weights == "IGSCA"
  )
  

## Prepare for MOB ---------------------------------------------------------
## Formula as it is done in networktree:::networktree.default(), by Payton J.
## Jones and Achim Zeileis
  mob_form  <- paste(
    paste(.object$Information$Model$indicators, collapse = " + "),
    "~",
    paste(.splitvars, collapse = " + ")
  ) |>
    stats::as.formula()
## 
  tree <- partykit::mob(
    formula = mob_form,
    data = .data,
    fit = csem_fit(.object=.object,
                   .model=.model,
                   .approach_weights=.approach_weights,
                   .iter_max=.iter_max,
                   .tolerance=.tolerance,
                   .disattenuate=.disattenuate,
                   .dominant_indicators=.dominant_indicators,
                   .conv_criterion=.conv_criterion),
    control = partykit::mob_control(
      ytype = "data.frame",
      prune = NULL,
      maxdepth = .maxdepth,
      minsize = .minsize,
      terminal = "object",
      restart = FALSE,
      verbose = TRUE
    )
  )
  
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
    if(.object$Information$Arguments$.approach_weights == "IGSCA") {
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


#' Calculate I-GSCA's Objective Function
#' 
#' Numerator of the fraction of unexplained variance in the data of the FIT statistic for GSCA models.
#' 
#' @param .object 
#'
#' @return Sum of Squares of Unexplained Variance
#' @keywords internal
calculateIgscaObjectiveFunction <- function(.object = NULL) {
  
  # As shown in the GSCA_m publication (Hwang et al., 2017)
  Gamma <- .object$Estimates$Construct_scores
  Psi <- cbind(.object$Information$Data, Gamma)
  # I am fairly confident the transpose of B is what's needed
  # See Gamma[1,] %*% t(...$Path_estimates)
  if (!is.null(.object$Estimates$Path_estimates)) {
    # If there's a structural model
    A <- cbind(.object$Estimates$Loading_estimates,
               t(.object$Estimates$Path_estimates))
  }
  else if (is.null(.object$Estimates$Path_estimates) | (!exists(".object$Estimates$Path_estimates"))) {
    # If no structural model
    A <- cbind(.object$Estimates$Loading_estimates,
               matrix(data = 0,
                      nrow = nrow(.object$Estimates$Loading_estimates),
                      ncol = nrow(.object$Estimates$Loading_estimates))
    )
  }
  
  if (!is.null(.object$Estimates$UniqueComponent)) {
    S <- cbind(.object$Estimates$UniqueComponent, matrix(data = 0, nrow = nrow(Gamma), ncol = ncol(Gamma)))
    
  } else if (is.null(.object$Estimates$UniqueComponent)) {
    # UniqueComponent should be NULL when GSCA and not GSCA_m/I-GSCA is run 
    S <- matrix(data = 0, nrow(Psi), ncol = ncol(Psi))  
    
  }
  
  SS_unexplained_variance <- sum(diag(t(Psi - Gamma %*% A - S) %*% (Psi - Gamma %*% A - S)))
  
  return(SS_unexplained_variance)
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


#' Prototype of doTrees
#'
#' Pending the computation of the gradient, here I use a manual greedy approach
#' towards finding optimal splits.
#' 
#' @inheritParams doTrees
#' @return
#' @export
#'
#' @examples
doTreesBeta <- function(.object, .splitvars, .model,
         .maxdepth = Inf,
         .minsize = nrow(tidy.cSEMResults(.object)),
         .data = .object$Information$Arguments$.data,
         .approach_weights = .object$Information$Arguments$.approach_weights,
         .iter_max = .object$Information$Arguments$.iter_max,
         .tolerance = .object$Information$Arguments$.tolerance,
         .disattenuate = .object$Information$Arguments$.disattenuate,
         .dominant_indicators = .object$Information$Arguments$.dominant_indicators,
         .conv_criterion = .object$Information$Arguments$.conv_criterion) {
  browser()
}