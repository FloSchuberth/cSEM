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
  
  # Add scores for nonlinear terms of the if a second order construct is attached
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
  # ###  Combine Results for final output ----------------------------------------
  # out1 <- .first_stage_results
  # n <- rownames(original_model$structural)
  # 
  # # n <- c(setdiff(names_constructs, rownames(model_ordered)), rownames(model_ordered))
  # # m <- order(which(model_measurement[n, ] == 1, arr.ind = TRUE)[, "row"])
  # # structural_ordered <- model_structural[n, c(n, setdiff(colnames(model_ordered), n))]
  # 
  # ## Path coefficients
  # n_structural <- gsub("_temp", "", rownames(out2$Estimates$Path_estimates))
  # dimnames(out2$Estimates$Path_estimates) <- list(n_structural, n_structural)
  # 
  # ## Loadings
  # tmp <- out2$Estimates$Loading_estimates[c_2nd_order, c_attached_to_2nd, drop = FALSE]
  # 
  # l1 <- nrow(out1$Estimates$Loading_estimates)
  # l2 <- ncol(out1$Estimates$Loading_estimates)
  # l3 <- nrow(tmp)
  # l4 <- ncol(tmp)
  # 
  # # Combine
  # Lambda <- rbind(cbind(out1$Estimates$Loading_estimates, 
  #                       matrix(0, l1, l4, dimnames = list(NULL, colnames(tmp)))), 
  #                 cbind(matrix(0, l3, l2), tmp))
  # # Reorder
  # m <- order(which(Lambda[n, ] != 0, arr.ind = TRUE)[, "row"])
  # Lambda <- Lambda[n, m]
  # 
  # ## Weights
  # tmp <- out2$Estimates$Weight_estimates[c_2nd_order, c_attached_to_2nd, drop = FALSE]
  # 
  # W <- rbind(cbind(out1$Estimates$Weight_estimates, 
  #                  matrix(0, l1, l4, dimnames = list(NULL, colnames(tmp)))), 
  #            cbind(matrix(0, l3, l2), tmp))
  # W <- W[n, m]
  # 
  # ## Inner weight estimates
  # dimnames(out2$Estimates$Inner_weight_estimates) <- list(n_structural, n_structural)
  # 
  # ## Construct scores
  # colnames(out2$Estimates$Construct_scores) <- n_structural
  # 
  # ## Indicator VCV
  # # not sure yet
  # 
  # ## Proxy VCV
  # dimnames(out2$Estimates$Proxy_VCV) <- list(n_structural, n_structural)
  # 
  # ## Construct VCV
  # dimnames(out2$Estimates$Construct_VCV) <- list(n_structural, n_structural)
  # 
  # ## Cross loadings
  # # not sure yet
  # 
  # ## Construct reliabilities
  # names(out2$Estimates$Construct_reliabilities) <- n_structural
  # 
  # 
  # out <- list(
  #   "Estimates"   = list(
  #     "Path_estimates"         = out2$Estimates$Path_estimates,
  #     "Loading_estimates"      = Lambda,
  #     "Weight_estimates"       = W,
  #     "Construct_scores"       = out2$Estimates$Construct_scores,
  #     "Indicator_VCV"          = cor(out1$Information$Data),
  #     "Proxy_VCV"              = out2$Estimates$Proxy_VCV,
  #     "Construct_VCV"          = out2$Estimates$Construct_VCV,
  #     "Cross_loadings"         = NULL,
  #     "Construct_reliabilities"= out2$Estimates$Construct_reliabilities,
  #     "Correction_factors"     = NULL,
  #     "R2"                     = out2$Estimates$R2
  #   ),
  #   "Information" = list(
  #     "Data"          = out1$Information$Data,
  #     "Model"         = original_model,
  #     "Arguments"     = out1$Information$Arguments,
  #     "Weight_info"   = list(
  #       "Modes"                  = NULL,
  #       "Number_iterations"      = NULL,
  #       "Inner_weight_estimates" = out2$Estimates$Inner_weight_estimates,
  #       "Convergence_status"     = out2$Information$Weight_info$Convergence_status
  #     )
  #   )
  # )
  return(out2)
}