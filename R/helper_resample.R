#' Helper for resamplecSEMResults.cSEMResults
#' @noRd
#' 
selectAndVectorize <- function(.object) {
  
  ## Select relevant statistics/parameters/quantities
  x1 <- list()
  summary_temp <- summarize(.object)
  
  if(inherits(.object,"cSEMResults_2ndorder")) {
    
    est1_temp <- summary_temp$First_stage$Estimates
    est2_temp <- summary_temp$Second_stage$Estimates
    
    # Path estimates
    x1[["Path_estimates"]] <- est2_temp$Path_estimates$Estimate
    names(x1[["Path_estimates"]]) <- est2_temp$Path_estimates$Name
    
    # Loading estimates
    x1[["Loading_estimates"]] <- c(est1_temp$Loading_estimates$Estimate,
                                   est2_temp$Loading_estimates$Estimate)
    names(x1[["Loading_estimates"]]) <- c(est1_temp$Loading_estimates$Name,
                                          est2_temp$Loading_estimates$Name)
    
    # Weight estimates
    x1[["Weight_estimates"]] <- c(est1_temp$Weight_estimates$Estimate,
                                  est2_temp$Weight_estimates$Estimate)
    names(x1[["Weight_estimates"]]) <- c(est1_temp$Weight_estimates$Name,
                                         est2_temp$Weight_estimates$Name)
    
    if(nrow(est1_temp$Residual_correlation) != 0) {
      # Residual correlation
      x1[["Residual_correlation"]] <- est1_temp$Residual_correlation$Estimate
      names(x1[["Residual_correlation"]]) <- est1_temp$Residual_correlation$Name 
    }
    
    if(nrow(est1_temp$Indicator_correlation) != 0) {
      # Residual correlation
      x1[["Indicator_correlation"]] <- est1_temp$Indicator_correlation$Estimate
      names(x1[["Indicator_correlation"]]) <- est1_temp$Indicator_correlation$Name 
    }
    
    if(nrow(est2_temp$Exo_construct_correlation) != 0) {
      # Residual correlation
      x1[["Exo_construct_correlation"]] <- est2_temp$Exo_construct_correlation$Estimate
      names(x1[["Exo_construct_correlation"]]) <- est2_temp$Exo_construct_correlation$Name 
    }
    
    if(.object$Second_stage$Information$Model$model_type == "Linear") {
      if(any(names(est2_temp$Effect_estimates) == "Indirect_effect")) {
        # Indirect effects
        x1[["Indirect_effect"]] <- est2_temp$Effect_estimates$Indirect_effect$Estimate
        names(x1[["Indirect_effect"]]) <- est2_temp$Effect_estimates$Indirect_effect$Name
      }
      if(any(names(est2_temp$Effect_estimates) == "Total_effect")) {
        # Total effect
        x1[["Total_effect"]] <- est2_temp$Effect_estimates$Total_effect$Estimate
        names(x1[["Total_effect"]]) <- est2_temp$Effect_estimates$Total_effect$Name         
      }
    }
  } else {
    
    est_temp <- summary_temp$Estimates
    
    # Path estimates
    x1[["Path_estimates"]] <- est_temp$Path_estimates$Estimate
    names(x1[["Path_estimates"]]) <- est_temp$Path_estimates$Name
    
    # Loading estimates
    x1[["Loading_estimates"]] <- est_temp$Loading_estimates$Estimate
    names(x1[["Loading_estimates"]]) <- est_temp$Loading_estimates$Name
    
    # Weight estimates
    x1[["Weight_estimates"]] <- est_temp$Weight_estimates$Estimate
    names(x1[["Weight_estimates"]]) <- est_temp$Weight_estimates$Name
    
    if(nrow(est_temp$Residual_correlation) != 0) {
      # Residual correlation
      x1[["Residual_correlation"]] <- est_temp$Residual_correlation$Estimate
      names(x1[["Residual_correlation"]]) <- est_temp$Residual_correlation$Name 
    }
    
    if(nrow(est_temp$Indicator_correlation) != 0) {
      # Residual correlation
      x1[["Indicator_correlation"]] <- est_temp$Indicator_correlation$Estimate
      names(x1[["Indicator_correlation"]]) <- est_temp$Indicator_correlation$Name 
    }
    
    if(nrow(est_temp$Exo_construct_correlation) != 0) {
      # Residual correlation
      x1[["Exo_construct_correlation"]] <- est_temp$Exo_construct_correlation$Estimate
      names(x1[["Exo_construct_correlation"]]) <- est_temp$Exo_construct_correlation$Name 
    }
    
    if(.object$Information$Model$model_type == "Linear") {
      if(any(names(est_temp$Effect_estimates) == "Indirect_effect")){
        # Indirect effects
        x1[["Indirect_effect"]] <- est_temp$Effect_estimates$Indirect_effect$Estimate
        names(x1[["Indirect_effect"]]) <- est_temp$Effect_estimates$Indirect_effect$Name 
      }
      
      if(any(names(est_temp$Effect_estimates) == "Total_effect")) {
      # Total effect
      x1[["Total_effect"]] <- est_temp$Effect_estimates$Total_effect$Estimate
      names(x1[["Total_effect"]]) <- est_temp$Effect_estimates$Total_effect$Name
      }
    }
  }
  
  # Return 
  x1
}

#' Helper for resamplecSEMResults.cSEMResults
#' @noRd
#' 
applyUserFuns <- function(.object, .user_funs, ...) {
  if(!(is.function(.user_funs) | is.list(.user_funs))) {
    stop2("The following error occured in the `resamplecSEMResults()` function:\n",
          "The arguments to `.user_funs` must be a single function or a list",
          " of functions")
  } 
  
  if(is.function(.user_funs)) {
    if(names(formals(.user_funs))[1] != ".object") {
      stop2("The following error occured in the `resamplecSEMResults()` function:\n",
            "The first argument of the function(s) provided to `.user_funs` must",
            " always be `.object`.")
    } else {
      # Delete unused ... arguments
      args <- list(...)[intersect(names(list(...)), names(formals(.user_funs)))]
      args$.object <- .object

      # Compute functions
      fun_out <- do.call(.user_funs, args)
      
      ## Check that the output is a vector or a matrix. Currently, no other
      ## format is allowed 
      if(is.vector(fun_out) & !is.list(fun_out)) {

        list("User_fun" = fun_out)
        
      } else if(is.matrix(fun_out)) {
        
        list("User_fun" = c(fun_out))
        
      } else {
        stop2(
          "The following error occured in the `resamplecSEMResults()` function:\n",
          "The output of the function(s) provided to `.user_funs` must",
          " always be an atomic vector (not a list) or matrix.")
      }
    }
  } else {
    x <- lapply(.user_funs, function(f) {
      # Delete unused ... arguments
      args <- list(...)[intersect(names(list(...)), names(formals(f)))]
      args$.object <- .object
      
      # Compute functions
      fun_out <- do.call(f, args)
      
      ## Check that the output is a vector or a matrix. Currently, no other
      ## format is allowed 
      if(is.vector(fun_out) & !is.list(fun_out)) { 
        
        fun_out
        
      } else if(is.matrix(fun_out))  {
        c(fun_out)
        
      } else {
        stop2(
          "The following error occured in the `resamplecSEMResults()` function:\n",
          "The output of the function(s) provided to `.user_funs` must",
          " always be an atomic vector (not a list) or matrix.")
      }
    })
    
    if(is.null(names(x))) {
      names(x) <- paste0("User_fun", 1:length(x))
    }
    # Return result
    x
  }
}

#' Helper for resamplecSEMResults.cSEMResults
#' @noRd
#' 
reverseSign <- function(.Est_temp, .summary_original) {
  
  x1 <- list()
  summary_temp <- summarize(.Est_temp)
  
  if(inherits(.Est_temp, "cSEMResults_2ndorder")) {
    est1 <- .summary_original$First_stage$Estimates
    est2 <- .summary_original$Second_stage$Estimates
    
    est1_temp <- summary_temp$First_stage$Estimates
    est2_temp <- summary_temp$Second_stage$Estimates
    
    ## Which signs differ:
    sign_diff_path <- sign(est2_temp$Path_estimates$Estimate) != 
      sign(est2$Path_estimates$Estimate)
    
    sign_diff_loadings <- sign(c(est1_temp$Loading_estimates$Estimate,
                                 est2_temp$Loading_estimates$Estimate)) !=
      sign(c(est1$Loading_estimates$Estimate,
             est2$Loading_estimates$Estimate))
    sign_diff_weights <- sign(c(est1_temp$Weight_estimates$Estimate,
                                est2_temp$Weight_estimates$Estimate)) !=
      sign(c(est1$Weight_estimates$Estimate,
             est2$Weight_estimates$Estimate))
    sign_diff_total_effect <- sign(est2_temp$Effect_estimates$Total_effect$Estimate) != 
      sign(est2$Effect_estimates$Total_effect$Estimate)
    
    sign_diff_indirect_effect <- sign(est2_temp$Effect_estimates$Indirect_effect$Estimate) != 
      sign(est2$Effect_estimates$Indirect_effect$Estimate)
    
  } else {
    est <- .summary_original$Estimates
    est_temp <- summary_temp$Estimates
    ## Which signs differ:
    sign_diff_path <- sign(est_temp$Path_estimates$Estimate) != 
      sign(est$Path_estimates$Estimate)
    
    sign_diff_loadings <- sign(est_temp$Loading_estimates$Estimate) !=
      sign(est$Loading_estimates$Estimate)
    sign_diff_weights  <- sign(est_temp$Weight_estimates$Estimate) !=
      sign(est$Weight_estimates$Estimate)
    
    sign_diff_total_effect <- sign(est_temp$Effect_estimates$Total_effect$Estimate) != 
      sign(est$Effect_estimates$Total_effect$Estimate)
    
    sign_diff_indirect_effect <- sign(est_temp$Effect_estimates$Indirect_effect$Estimate) != 
      sign(est$Effect_estimates$Indirect_effect$Estimate)
  }
  
  ## Multiply the coefficients for which the sign differs by (-1)
  
  # Path estimates
  x1[["Path_estimates"]][sign_diff_path] <- 
    x1[["Path_estimates"]][sign_diff_path] * (-1)
  
  # Loading estimates
  x1[["Loading_estimates"]][sign_diff_loadings] <- 
    x1[["Loading_estimates"]][sign_diff_loadings] * (-1)
  
  # Weight estimates
  x1[["Weight_estimates"]][sign_diff_weights] <-
    x1[["Weight_estimates"]][sign_diff_weights] * (-1)
  
  # Total effect estimates
  x1[["Total_effect"]][sign_diff_weights] <-
    x1[["Total_effect"]][sign_diff_weights] * (-1)
  
  # Total effect estimates
  x1[["Indirect_effect"]][sign_diff_weights] <-
    x1[["Indirect_effect"]][sign_diff_weights] * (-1)
  
  ## Return
  x1
}
