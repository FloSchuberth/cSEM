#' Internal: Second/Third stage of the two-stage approach for second order constructs
#'
#' Performs the second and third stage for a model containing second order 
#' constructs.
#'
#' @usage calculate2ndStage(
#'  .csem_model          = args_default()$.csem_model,
#'  .first_stage_results = args_default()$.first_stage_results,
#'  .original_arguments  = args_default()$.original_arguments,
#'  .approach_2ndorder   = args_default()$.approach_2ndorder
#'   )
#'
#' @inheritParams csem_arguments
#
#' @return A cSEMResults object.
#'
#' @keywords internal
#'
calculate2ndStage <- function(
  .csem_model          = args_default()$.csem_model,
  .first_stage_results = args_default()$.first_stage_results,
  .original_arguments  = args_default()$.original_arguments,
  .approach_2ndorder   = args_default()$.approach_2ndorder
) {
  
  original_model <- .csem_model
  
  ## Parse second stage model
  model2 <- convertModel(
    .csem_model        = original_model,
    .approach_2ndorder = "2stage",
    .stage             = "second")
  
  ## Get scores for the second stage
  scores         <- .first_stage_results$Estimates$Construct_scores
  
  # Add scores for nonlinear terms if a second order construct is attached
  # to a nonlinear term
  # Note: 23.07.2019: We dont allow for first order nonlinear terms to build/
  #                   measure a second-order construct. Hence the code below 
  #                   is superfluos.
  # if(any(grepl("\\.", colnames(model2$measurement)))) {
  #   
  #   temp_names <- grep("\\.", colnames(model2$measurement), value = TRUE)
  #   temp <- lapply(strsplit(temp_names, "\\."), function(x) {
  #     matrixStats::rowProds(scores[, x, drop = FALSE])
  #   })
  #   
  #   temp <- do.call(cbind, temp)
  #   colnames(temp) <- temp_names
  #   scores <- cbind(scores, temp)
  # }
  
  ## Collect all necessary sets
  vars_2nd <- original_model$vars_2nd 
  vars_attached_to_2nd <- original_model$vars_attached_to_2nd
  vars_not_attached_to_2nd <- original_model$vars_not_attached_to_2nd
  
  ## Only linear terms required for subsetting reliabilities
  vars_not_attached_to_2nd <- unique(unlist(strsplit(vars_not_attached_to_2nd, "\\.")))
  
  ## Remove 2nd order terms if there are any in the interaction terms
  vars_not_attached_to_2nd <- setdiff(vars_not_attached_to_2nd, vars_2nd)
  
  ## Get all reliabilities from first stage 
  rel_all_1step  <- .first_stage_results$Estimates$Reliabilities
  
  # Select reliabilities for those constructs that are not attached
  # to a second order factor
  if(!is.null(vars_not_attached_to_2nd)) {
    rel_not_attached_to_2nd <- rel_all_1step[vars_not_attached_to_2nd]
    names(rel_not_attached_to_2nd) <- paste0(vars_not_attached_to_2nd, "_temp")
  } else {
    rel_not_attached_to_2nd <- NULL
  }

  ## Save arguments of the first stage 
  args <- .first_stage_results$Information$Arguments
  
  ## Estimate second stage
  out2 <- csem(
    .data                        = scores, 
    .model                       = model2,
    .approach_cor_robust         = args$.approach_cor_robust,
    .approach_nl                 = args$.approach_nl,
    .approach_paths              = args$.approach_paths,
    .approach_weights            = args$.approach_weights,
    .conv_criterion              = args$.conv_criterion,
    .disattenuate                = args$.disattenuate,
    .dominant_indicators         = NULL,
    .estimate_structural         = args$.estimate_structural,
    .id                          = NULL,
    .iter_max                    = args$.iter_max,
    .normality                   = args$.normality,
    .PLS_approach_cf             = args$.PLS_approach_cf,
    .PLS_ignore_structural_model = args$.PLS_ignore_structural_model,
    .PLS_modes                   = NULL,
    .PLS_weight_scheme_inner     = args$.PLS_weight_scheme_inner,
    .reliabilities               = rel_not_attached_to_2nd,
    .starting_values             = NULL,
    .tolerance                   = args$.tolerance
    )
  
  ## Correct loadings (this basically rebases loadings)
  out2$Estimates$Loading_estimates <-  t(apply(out2$Estimates$Loading_estimates, 1, function(x) {
    x / sqrt(rel_all_1step[colnames(out2$Information$Model$measurement)])
  }))
  colnames(out2$Estimates$Loading_estimates) <- colnames(out2$Information$Model$measurement)
  
  ### Third stage --------------------------------------------------------------
  # If a 2nd order construct is a composite build by at least one common factor 
  # an additional correction is necessary. Mainly based on:
  # VanRiel (2017)- "Estimating hierarchical constructs using consistent 
  # partial least squares: The case of second-order composites of common factors"
  
  if(any(original_model$construct_type[vars_attached_to_2nd] == "Common factor")) {
    
    # Which second order constructs are composites?
    vars_2nd_composites <- original_model$construct_type[vars_2nd]
    vars_2nd_composites <- names(vars_2nd_composites[
      vars_2nd_composites == "Composite"])
    
    if(length(vars_2nd_composites) != 0) {
      
      rel_2nd_order        <- rep(1, length(vars_2nd_composites))
      names(rel_2nd_order) <- vars_2nd_composites
      
      for(i in vars_2nd_composites) {
        
        col_names   <- colnames(original_model$measurement[
          i, original_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
        w           <- out2$Estimates$Weight_estimates[i, col_names, drop = FALSE]
        Sstar       <- out2$Estimates$Indicator_VCV[col_names, col_names,drop=FALSE]
        diag(Sstar) <- rel_all_1step[col_names]
        
        rel_2nd_order[i] <- c(w %*% Sstar %*% t(w))
      }
      
      rel                         <- out2$Estimates$Reliabilities
      rel[vars_2nd_composites] <- rel_2nd_order 
      
      ## Redo second stage including new reliabilities (= third stage)
      out2 <- csem(
        .data                        = scores, 
        .model                       = model2,
        .approach_cor_robust         = args$.approach_cor_robust,
        .approach_nl                 = args$.approach_nl,
        .approach_paths              = args$.approach_paths,
        .approach_weights            = args$.approach_weights,
        .conv_criterion              = args$.conv_criterion,
        .disattenuate                = args$.disattenuate,
        .dominant_indicators         = NULL,
        .estimate_structural         = args$.estimate_structural,
        .id                          = NULL,
        .iter_max                    = args$.iter_max,
        .normality                   = args$.normality,
        .PLS_approach_cf             = args$.PLS_approach_cf,
        .PLS_ignore_structural_model = args$.PLS_ignore_structural_model,
        .PLS_modes                   = NULL,
        .PLS_weight_scheme_inner     = args$.PLS_weight_scheme_inner,
        .reliabilities               = rel,
        .starting_values             = NULL,
        .tolerance                   = args$.tolerance
      )
      
      for(i in vars_2nd_composites) {
        
        col_names <- colnames(original_model$measurement[
          i, original_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
        
        ## Compute consistent weights
        lambda <- out2$Estimates$Loading_estimates[i, col_names]
        q      <- lambda/sqrt(rel_all_1step[col_names])
        S      <- .first_stage_results$Estimates$Construct_VCV[col_names, col_names]
        v      <- solve(S) %*% q
        
        ## Standardize weights
        w        <- c(v / sqrt(c(t(v)%*%S%*%v)))
        names(w) <- col_names
        
        ## Compute consistent loadings
        lambda        <- c(S%*%t(t(w)))
        names(lambda) <- col_names
        
        ## Collect and replace
        out2$Estimates$Loading_estimates[i, col_names] <- lambda
        out2$Estimates$Weight_estimates[i, col_names]  <- w
      }
      
      ## Correct single indicator loadings
      if(!is.null(vars_not_attached_to_2nd)) {
        for(i in paste0(vars_not_attached_to_2nd, "_temp")) {
          out2$Estimates$Loading_estimates[i, ] <- 
            out2$Estimates$Loading_estimates[i, ] / 
            sqrt(rel_all_1step[colnames(out2$Information$Model$measurement)])
        } # END correct single indicator loadings 
      } # END if
    } # END if length(vars_2nd_composites) != 0
  } # END third stage
  
  out <- list("First_stage" = .first_stage_results, "Second_stage" = out2)
  
  ## Append original arguments needed as they are required by e.g. testOMF.
  # Since 
  .original_arguments[[".model"]] <- original_model
  out$Second_stage$Information$Arguments_original <- .original_arguments
  
  ## Add second order approach to $Information
  out$First_stage$Information$Approach_2ndorder  <- .approach_2ndorder
  out$Second_stage$Information$Approach_2ndorder <- .approach_2ndorder
  
  class(out) <- c("cSEMResults", "cSEMResults_2ndorder")
  
  return(out)
}