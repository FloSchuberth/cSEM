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
#' Pending the computation of the gradient, here I use a manual multi-group
#' versus single group hypothesis test of the difference in FIT statistic,
#' similar to semtrees.
#' 
#' The current implementation assumes that once you use a split-variable, you
#' cannot re-use it. The final implementation will re-use variables.
#' 
#' `.splitvars` currently only supports binary variables
#' `.maxdepth` doesn't do anything
#' `.minsize` does not currently work
#' 
#' @param .splitvars List of the column names to try splitting on.
#' 
#' @inheritParams doTrees
#' @inheritParams csem_arguments
#' @return List of whether split occured and with what variable, the associated
#'   paired t-test and the resulting model
#' @export
#' @importFrom boot boot
doTreesBeta <- function(.object,
                        .splitvars,
                        .model,
                        # .maxdepth = Inf,
                        .R = 100 #,
                        # .minsize = nrow(tidy.cSEMResults(.object)[!(tidy.cSEMResults(.object)$op %in% c("Indirect_effect", "Direct_effect", "Total_effect")),])
                        ) {
                        
  
## Load Arguments ----------------------------------------------------------

  
  .data <- .object$Information$Arguments$.data
  .approach_weights <- .object$Information$Arguments$.approach_weights
  .conv_criterion <- .object$Information$Arguments$.conv_criterion
  .disattenuate <- .object$Information$Arguments$.disattenuate
  .dominant_indicators <- .object$Information$Arguments$.dominant_indicators
  .iter_max <- .object$Information$Arguments$.iter_max
  .tolerance <- .object$Information$Arguments$.tolerance
  
  # FIXME: I don't think .minsize is working correctly
  # .minsize <- nrow(tidy.cSEMResults(.object)[!(tidy.cSEMResults(.object)$op %in% c("Indirect_effect", "Direct_effect", "Total_effect")),])
  
  # if (.minsize > nrow(.data)) {
  #   return(list("split" = NA,
  #          "p.value" = NA,
  #          "FIT" = NA,
  #          "t.test" = NA,
  #          "fitted.model" = NA,
  #          "daughter.model" = NA))
  # }
  
## Figure out what splits to try and pre-process split variables------------
  
  # TODO: Test for whether every column referred to by .splitvars is binary or not
  # Checks through .splitvars
  #.splitvars 
  
## Evaluate Split ----------------------------------------------------------
  
  bootFIT <- boot::boot(data = .data,
                        statistic = calculateFITForSplit,
                        R = .R,
                        sim = "ordinary",
                        .model = .model,
                        .id = .splitvars,
                        .approach_weights = .approach_weights,
                        .conv_criterion = .conv_criterion,
                        .disattenuate = .disattenuate,
                        .dominant_indicators = .dominant_indicators,
                        .iter_max = .iter_max,
                        .tolerance = .tolerance)
  
  # FIXME: Why are the replications failing so frequently
  # browser()
  # If there's fit failure, remove the row
  bootFIT[["t"]] <- na.omit(bootFIT[["t"]])
### Hypothesis Testing ------------------------------------------------------
  if (length(.splitvars) == 1) {
    paired_bootstrap_t_test <- t.test(
      x = bootFIT[["t"]][, 2], # Multi-group FIT
      y = bootFIT[["t"]][, 1], # Single-Group FIT
      paired =  TRUE,
      var.equal = FALSE,
      alternative = "greater"
    )
    # browser()
    # Return result to user
    if (paired_bootstrap_t_test$p.value <= 0.05 & paired_bootstrap_t_test$statistic > 0) {
      return(
        list(
          "split" = .splitvars,
          "p.value" = paired_bootstrap_t_test$p.value,
          "FIT" = bootFIT$t0[2],
          "t.test" = paired_bootstrap_t_test,
          "fitted.model" = csem(
            .data = .data,
            .id = .splitvars,
            .model = .model,
            .approach_weights = .approach_weights,
            .tolerance = .tolerance,
            .conv_criterion = .conv_criterion
          ),
          "daughter.model" = NA
        )
      )
    } else {
      return(
        list(
          "split" = .splitvars, # An attempted split in this case
          "p.value" = paired_bootstrap_t_test$p.value,
          "FIT" = bootFIT$t0[1],
          "t.test" = paired_bootstrap_t_test,
          "fitted.model" = NA,
          "daughter.model" = NA
        )
      )
    }
  } else if (length(.splitvars) > 1) {
    paired_bootstrap_t_tests <- apply(
      X =bootFIT[["t"]][, 2:ncol(bootFIT[["t"]])],
      MARGIN = 2,
      FUN = function(.x) {
        return(t.test(
          x = .x,
          y = bootFIT[["t"]][, 1],
          paired = TRUE,
          var.equal = FALSE,
          alternative = "greater"
        ))
      }
    )
    
    names(paired_bootstrap_t_tests) <- names(bootFIT[["t0"]])[2:length(bootFIT[["t0"]])]
    
    Ps <- unlist(lapply(
      paired_bootstrap_t_tests,
      FUN = function(.x)
        .x$p.value
    ))
    Ts <- unlist(lapply(
      paired_bootstrap_t_tests,
      FUN = function(.x)
        return(unname(.x$statistic))
    ))
    PsTemp <- Ps
    
    
    for (i in length(Ps)) {
      if (Ts[which.min(Ps)] >= 0 && Ps[which.min(Ps)] <= 0.05) {
        
        fitted.model <-csem(
          .data = .data,
          .id = names(which.min(Ps)),
          .model = .model,
          .approach_weights = .approach_weights,
          .tolerance = .tolerance,
          .conv_criterion = .conv_criterion
        )
        
        daughter.model <- lapply(fitted.model, function(.obj_it)
          return(
            doTreesBeta(
              .object = .obj_it,
              .splitvars = .splitvars[!(.splitvars %in% names(which.min(Ps)))],
              .model = .model,
              .R = .R
            )
          ))
        return(
          list(
            "split" = names(which.min(Ps)),
            "p.value" = Ps[which.min(Ps)],
            "FIT" = bootFIT$t0[names(which.min(Ps))],
            "t.test" = paired_bootstrap_t_tests[[names(which.min(Ps))]],
            "fitted.model" = fitted.model,
            "daughter.model" = daughter.model)
        )
      } else if (length(Ts) > 0 && length(Ps) > 0) { # FIXME: The conditions here are likely incorrect in some subtle way, creating strange results
        Ts <- Ts[!(names(Ts) %in% names(which.min(Ps)))]
        Ps <- Ps[!(names(Ps) %in% names(which.min(Ps)))]
        
      } else {
        return(list(
          "split" = NA,
          "p.value" = PsTemp,
          "FIT" = bootFIT$t0[1],
          "t.test" = paired_bootstrap_t_tests),
          "fitted.model" = NA,
          "daughter.model" = NA
        )
      }
    }
    warning("Split failure")
    return(NA) 
  } else {
    stop(".splitvars argument is incorrect")
  }

    # TODO: If .maxdepth has not been reached and some split was significant, then fit the multigroup model and apply doTreesBeta to each sub-group
    # TODO: If .maxdepth has been reached, then stop and return the stopping point.
  
}


#' Calculate the difference in FIT
#'
#'
#' @param data Data passed by boot
#' @param idx Row-indices passed by boot
#' @inheritParams csem
#' @return Named vector of the delta between the FIT of the multi-group model
#'   versus the FIT of the single-group model--positive values favor
#'   multi-group, negative, single-group; the FIT of the single-group model; and
#'   the FIT of the multi-group model.
#' @keywords internal
#'
calculateFITForSplit <- function(data,
                                 idx,
                                 .model,
                                 .id,
                                 .approach_weights,
                                 .conv_criterion,
                                 .disattenuate,
                                 .dominant_indicators,
                                 .iter_max,
                                 .tolerance) {
  
  # Robust Boot by-pass by Thomas on June 11,2013 https://stackoverflow.com/a/17040580
  ret <- tryCatch({
    sg_mod <- csem(
      .data = data[idx, ],
      .model = .model,
      .approach_weights = .approach_weights,
      .tolerance = .tolerance,
      .conv_criterion = .conv_criterion
    )
    sg_fit <- calculateFIT(sg_mod)
    
    if (length(.id) == 1) {
      mg_mod <- csem(
        .data = data[idx, ],
        .id = .id,
        .model = .model,
        .approach_weights = .approach_weights,
        .tolerance = .tolerance,
        .conv_criterion = .conv_criterion
      )
      
      mg_fit <- bdiagonalizeMultiGroupIgscaEstimates(mg_mod) |>
        calculateFIT()
      
      names(mg_fit) <- .id
      
      return(c("sg_fit" = sg_fit, mg_fit))
      
    } else if (length(.id) > 1) {
      mg_fits <- lapply(.id, function(.id_iter) {
        mg_mod <- csem(
          .data = data[idx, ],
          .id = .id_iter,
          .model = .model,
          .approach_weights = .approach_weights,
          .tolerance = .tolerance,
          .conv_criterion = .conv_criterion
        )
        
        mg_fit <- bdiagonalizeMultiGroupIgscaEstimates(mg_mod) |>
          calculateFIT()
        names(mg_fit) <- .id_iter
        
        return(mg_fit)
      })
      
      return(c("sg_fit" = sg_fit, unlist(mg_fits)))
      
    } else {
      stop("Inappropriate length or type of .splitvars passed")
      
    }
  }, error = function(error) {
    return(c("sg_fit" = NA, "mg_fit" = NA))
  })
  
  return(ret)
  
}
