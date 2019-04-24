#' Second order constructs
#'
#' Performs the second and third stage for a model containing second order 
#' constructs.
#'
#' @usage calculate2ndOrder(.csem_model, .first_stage_results)
#'
#' @inheritParams csem_arguments
#
#' @return A cSEMResults object.
#'
#' @keywords internal
#'
calculate2ndOrder <- function(
  .csem_model          = args_default()$.csem_model,
  .first_stage_results = args_default()$.first_stage_results
) {

  original_model <- .csem_model

  ## Parse second stage model
  model2 <- convertModel(
    .csem_model        = original_model,
    .approach_2ndorder = "3stage",
    .stage             = "second")
  
  ## Get scores for the second stage
  scores         <- .first_stage_results$Estimates$Construct_scores
  
  # Add scores for nonlinear terms if a second order construct is attached
  # to a nonlinear term
  if(any(grepl("\\.", colnames(model2$measurement)))) {
    
    temp_names <- grep("\\.", colnames(model2$measurement), value = TRUE)
    temp <- lapply(strsplit(temp_names, "\\."), function(x) {
      matrixStats::rowProds(scores[, x, drop = FALSE])
    })
    
    temp <- do.call(cbind, temp)
    colnames(temp) <- temp_names
    scores <- cbind(scores, temp)
  }
  
  ## Collect all necessary sets
  # All linear constructs of the original model
  c_linear_original     <- rownames(original_model$structural)
  # All constructs used in the first step (= all first order constructs)
  c_linear_1step        <- names(original_model$construct_order[
    original_model$construct_order == "First order"])
  # All second order constructs
  c_2nd_order           <- setdiff(c_linear_original, c_linear_1step)
  # All indicators of the original model (including linear and nonlinear 
  # constructs that form/measure a second order construct)
  i_original            <- colnames(original_model$measurement)
  i_linear_original     <- intersect(c_linear_original, i_original)
  i_nonlinear_original  <- grep("\\.", i_original, value = TRUE) 
  # Linear constructs not attatched to second order constructs
  c_not_attached_to_2nd <- setdiff(c_linear_1step, i_linear_original)
  # Linear constructs attached to second order constructs
  c_attached_to_2nd     <- intersect(c_linear_1step, i_linear_original)
  
  ## Get all reliabilities from first stage 
  rel_all_1step  <- .first_stage_results$Estimates$Construct_reliabilities
  
  # Select reliabilities for those constructs that are not attached
  # to a second order factor
  rel_not_attached_to_2nd <- rel_all_1step[c_not_attached_to_2nd]
  names(rel_not_attached_to_2nd) <- paste0(c_not_attached_to_2nd, "_temp")
  
  ## Perform second stage
  out2 <- csem(.data          = scores, 
               .model         = model2, 
               .reliabilities = rel_not_attached_to_2nd)
  
  ## Correct loadings (this basically rebases loadings)
  out2$Estimates$Loading_estimates <-  t(apply(out2$Estimates$Loading_estimates, 1, function(x) {
    x / sqrt(rel_all_1step[colnames(out2$Information$Model$measurement)])
  }))
  colnames(out2$Estimates$Loading_estimates) <- colnames(out2$Information$Model$measurement)
  
  ### Third stage --------------------------------------------------------------
  # If a 2nd order construct is a composite build by at least on common factor 
  # an additional correction is necessary. Mainly based on:
  # VanRiel (2017)- "Estimating hierarchical constructs using consistent 
  # partial least squares: The case of second-order composites of common factors"
  
  if(any(original_model$construct_type[c_attached_to_2nd] == "Common factor")) {
    
    # Which second order constructs are composites?
    c_2nd_order_composites <- original_model$construct_type[c_2nd_order]
    c_2nd_order_composites <- names(c_2nd_order_composites[
      c_2nd_order_composites == "Composite"])
    
    if(length(c_2nd_order_composites) != 0) {
      
      rel_2nd_order        <- rep(1, length(c_2nd_order_composites))
      names(rel_2nd_order) <- c_2nd_order_composites
      
      for(i in c_2nd_order_composites) {
        
        col_names   <- colnames(original_model$measurement[
          i, original_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
        w           <- out2$Estimates$Weight_estimates[i, col_names, drop = FALSE]
        Sstar       <- out2$Estimates$Indicator_VCV[col_names, col_names]
        diag(Sstar) <- rel_all_1step[col_names]
        
        rel_2nd_order[i] <- c(w %*% Sstar %*% t(w))
      }
      
      rel                         <- out2$Estimates$Construct_reliabilities
      rel[c_2nd_order_composites] <- rel_2nd_order 
      
      ## Redo second stage including new reliabilities (= third stage)
      out2 <- csem(.data          = scores, 
                   .model         = model2, 
                   .reliabilities = rel)
      
      for(i in c_2nd_order_composites) {
        
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
      
      ## Correct sinlge indicator loadings
      for(i in paste0(c_not_attached_to_2nd, "_temp")) {
        out2$Estimates$Loading_estimates[i, ] <- 
          out2$Estimates$Loading_estimates[i, ] / 
          sqrt(rel_all_1step[colnames(out2$Information$Model$measurement)])
      } # END correct single indicator loadings
    } # END if length(c_2nd_order_composites) != 0
  } # END third stage

  return(out2)
}