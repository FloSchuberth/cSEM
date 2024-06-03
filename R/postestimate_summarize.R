#' Summarize model
#' 
#' \lifecycle{stable}
#'
#' The summary is mainly focused on estimated parameters. For quality criteria 
#' such as the average variance extracted (AVE), reliability estimates, 
#' effect size estimates etc., use [assess()].
#' 
#' If `.object` contains resamples, standard errors, t-values and p-values
#' (assuming estimates are standard normally distributed) are printed as well.
#' By default the percentile confidence interval is given as well. For other
#' confidence intervals use the `.ci` argument. See [infer()] for possible choices
#' and a description.
#'
#' @return An object of class `cSEMSummarize`. A `cSEMSummarize` object has
#'   the same structure as the [cSEMResults] object with a couple differences:
#'   \enumerate{
#'   \item{Elements `$Path_estimates`, `$Loadings_estimates`, `$Weight_estimates`,
#'     `$Weight_estimates`, and  `$Residual_correlation` are standardized data frames instead of matrices.}
#'   \item{Data frames `$Effect_estimates`, `$Indicator_correlation`, and
#'     `$Exo_construct_correlation` are added to `$Estimates`.}
#'   } 
#'   The data frame format is usually much more convenient if users intend to
#'   present the results in e.g., a paper or a presentation.
#'   
#' @usage summarize(
#'  .object = NULL, 
#'  .alpha  = 0.05,
#'  .ci     = NULL,
#'  ...
#'  )
#'
#' @inheritParams csem_arguments
#' @param ... Further arguments to `summarize()`. Currently ignored.
#'
#' @seealso [csem], [assess()], [cSEMResults], [exportToExcel()]
#'
#' @example inst/examples/example_summarize.R
#' @export
#'
summarize <- function(
  .object                = NULL, 
  .alpha                 = 0.05,
  .ci                    = NULL,
  ...
  ) {
  UseMethod("summarize")
}

#' @export

summarize.cSEMResults_default <- function(
  .object                = NULL, 
  .alpha                 = 0.05,
  .ci                    = NULL,
  ...
  ) {
  
  x1  <- .object$Estimates
  x2  <- .object$Information

  ## Add estimation Status
  .object$Information$Estimation_status <- unclass(verify(.object))
  
  ### Structure output =========================================================
  ## Path estimates ------------------------------------------------------------
  if(!all(x2$Model$structural == 0)) {
    # Get construct type for relevant variables
    type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$structural))
    
    # Build names
    temp <- outer(rownames(x1$Path_estimates), colnames(x1$Path_estimates), 
                  FUN = function(x, y) paste(x, y, sep = " ~ "))
    
    path_estimates <- data.frame(
      "Name"           = t(temp)[t(x2$Model$structural) != 0],
      "Construct_type" = type,
      "Estimate"       = t(x1$Path_estimates)[t(x2$Model$structural) != 0 ],
      "Std_err"        = NA,
      "t_stat"         = NA,
      "p_value"        = NA,
      stringsAsFactors = FALSE)
    
    # Delete rownames
    rownames(path_estimates) <- NULL 
  } else {
    path_estimates <- NULL
  }
  
  ## Loading estimates ---------------------------------------------------------
  # Get construct type for relevant variables
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$measurement))
  
  # Build names
  temp <- rep(rownames(x1$Loading_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, " =~ ", colnames(x1$Loading_estimates))
  
  loading_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type,
    "Estimate"       = x1$Loading_estimates[x2$Model$measurement != 0 ],
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  # Delete rownames
  rownames(loading_estimates) <- NULL
  
  ## Weight estimates ----------------------------------------------------------
  temp <- rep(rownames(x1$Weight_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, " <~ ", colnames(x1$Weight_estimates))
  
  weight_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type, 
    "Estimate"       = x1$Weight_estimates[x2$Model$measurement != 0 ],
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  # Delete rownames
  rownames(weight_estimates) <- NULL
  
  ## Inner weight estimates ----------------------------------------------------
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    i <- rownames(x1$Inner_weight_estimates)
    if(all(x2$Model$structural == 0)) {
      D <- x2$Model$cor_specified[i, i , drop = FALSE]
    } else {
      D <- x2$Model$structural[i, i , drop = FALSE] + t(x2$Model$structural[i, i, drop = FALSE]) 
    }
    
    # Note: June 2019
    if(any(D == 2)) { # non recursive model
      # Set elements back to 1 
      D[D == 2] <- 1 
    }
    
    temp <- outer(i, i, FUN = function(x, y) paste(x, y, sep = " -- "))
    type <- rep(x2$Model$construct_type, times = colSums(D))
    
    inner_weight_estimates <- data.frame(
      "Name"           = t(temp)[t(D) != 0],
      "Construct_type" = type, 
      "Estimate"       = t(x1$Inner_weight_estimates)[t(D) != 0 ], 
      stringsAsFactors = FALSE)
    
    # Delete rownames
    rownames(inner_weight_estimates) <- NULL
  }
  
  ## Residual correlation ------------------------------------------------------
  # Set up empty matrix to select the relevant residual correlations 
  cor_selector <- x2$Model$error_cor

  # Variances are not of interest and each correlation should appear only once
  cor_selector[upper.tri(cor_selector, diag = TRUE)] <- 0
  
  # Build names
  temp <- outer(rownames(x1$Residual_correlation), colnames(x1$Residual_correlation), 
                FUN = function(x, y) paste(x, y, sep = " ~~ "))
  
  ## For composites we need the indicator correlations not the residuals
  if(sum(cor_selector) != 0) {
    residual_correlation <- data.frame(
      "Name"           = t(temp)[cor_selector != 0],
      "Estimate"       = x1$Residual_correlation[cor_selector != 0],
      "Std_err"        = NA,
      "t_stat"         = NA,
      "p_value"        = NA,
      stringsAsFactors = FALSE)
    
    # Delete rownames
    rownames(residual_correlation) <- NULL 
  } else {
    residual_correlation <- data.frame(NULL)
  }
  

  ## D2 ----------------------------------------------------------------------
  # Build names
  if (!is.null(x1$D2)) {
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$measurement))
  temp <- paste0(names(x1$D2), " ~~ ", names(x1$D2))
  
  D2 <- data.frame(
    "Name"           = temp,
    "Construct_type" = type,
    "Estimate"       = x1$D2,
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  # Delete rownames
  rownames(D2) <- NULL
  } else if (is.null(x1$D2)) {
    D2 <- x1$D2
  }
  
  
  ## Indicator correlation ------------------------------------------------------
  # Set up empty matrix to select the relevant residual correlations 
  cor_selector[] <- 0
  
  ## Modify and fill Lambda
  for(i in names(x2$Model$construct_type)) {
    indicators <- colnames(x2$Model$measurement[i, x2$Model$measurement[i, ] != 0, drop = FALSE])
    
    # For common factors only the measurment errors specified by the user
    # are returned. Since they are already contained in cor_selector nothing needs
    # to be done
    if(x2$Model$construct_type[i] == "Composite") {
      cor_selector[indicators, indicators] <- 1
    }
  }
  
  # Variances are not of interest and each correlation should appear only once
  cor_selector[upper.tri(cor_selector, diag = TRUE)] <- 0
  
  # Build names
  temp <- outer(rownames(x1$Indicator_VCV), colnames(x1$Indicator_VCV), 
                FUN = function(x, y) paste(x, y, sep = " ~~ "))
  
  ## For composites we need the indicator correlations not the residuals
  if(sum(cor_selector) != 0) {
    indicator_correlation <- data.frame(
      "Name"           = t(temp)[cor_selector != 0],
      "Estimate"       = x1$Indicator_VCV[cor_selector != 0],
      "Std_err"        = NA,
      "t_stat"         = NA,
      "p_value"        = NA,
      stringsAsFactors = FALSE)
    
    # Delete rownames
    rownames(indicator_correlation) <- NULL 
  } else {
    indicator_correlation <- data.frame(NULL)
  }
  
  ## (Exogenous) construct correlations ----------------------------------------
  # Set up empty matrix to select the relevant residual correlations 
  exo_cors <- as.matrix(x1$Construct_VCV[x2$Model$cons_exo, x2$Model$cons_exo])
  
  # Build names
  temp <- outer(rownames(x1$Construct_VCV[x2$Model$cons_exo, x2$Model$cons_exo]), 
                colnames(x1$Construct_VCV[x2$Model$cons_exo, x2$Model$cons_exo]), 
                FUN = function(x, y) paste(x, y, sep = " ~~ "))
  
  ## For composites we need the indicator correlations not the residuals
  if(nrow(exo_cors) > 1) {
    exo_construct_correlation <- data.frame(
      "Name"           = t(temp)[lower.tri(exo_cors, diag = FALSE)],
      "Estimate"       = exo_cors[lower.tri(exo_cors, diag = FALSE)],
      "Std_err"        = NA,
      "t_stat"         = NA,
      "p_value"        = NA,
      stringsAsFactors = FALSE)
    
    # Delete rownames
    rownames(exo_construct_correlation) <- NULL 
  } else {
    exo_construct_correlation <- data.frame(NULL)
  }

  
  ## Construct scores ----------------------------------------------------------
  construct_scores <- as.data.frame(x1$Construct_scores)
  
  ## Effects -------------------------------------------------------------------
  if(x2$Model$model_type == "Linear" & !all(x2$Model$structural == 0)) {
    ## Effects are currently only for linear models. See issue #252 on github.
    effects <- calculateEffects(.object, .output_type = "data.frame")
    # Add effects to .object
    .object$Estimates$Effect_estimates <- effects
  }
  
  ## Inference =================================================================
  
  if(inherits(.object, "cSEMResults_resampled")) {
    ## Check arguments
    match.arg(.ci, args_default(.choices = TRUE)$.ci, several.ok = TRUE)
    
    ## The 95% percentile CI is returned by default (BCa is actually better 
    ##  but requires a second resampling run and should therefore not be the 
    ##  default).
    if(is.null(.ci)) .ci <- "CI_percentile"
    
    quantity <- c("sd", .ci)
    
    infer_out <- infer(
      .object          = .object,
      .alpha           = .alpha,
      .quantity        = quantity,
      ...
    )
    # Note (19.05.2021): infer() also computes standard deviations and confidence
    # intervals for constant values. Because of floating point impressions
    # its possible that the sd is not exactly zero in some instances. This
    # needs to be kept in mind for the print method of summarize()
    
    
    ## Path estimates ----------------------------------------------------------
    if(!all(x2$Model$structural == 0)) {
      path_estimates <- addInfer(infer_out$Path_estimates, path_estimates, .ci)
    }
    
    ## Loading estimates -------------------------------------------------------
    loading_estimates <- addInfer(infer_out$Loading_estimates, loading_estimates, .ci)
    
    ## Weight estimates --------------------------------------------------------
    weight_estimates <- addInfer(infer_out$Weight_estimates, weight_estimates, .ci)
    
    ## Residual correlation ----------------------------------------------------
    if(!is.null(infer_out$Residual_correlation)) {
      residual_correlation <- addInfer(infer_out$Residual_correlation, 
                                       residual_correlation, .ci)
    }
    
    ## Exogenous construct correlation -----------------------------------------
    if(!is.null(infer_out$Exo_construct_correlation)) {
      exo_construct_correlation <- addInfer(infer_out$Exo_construct_correlation, 
                                            exo_construct_correlation, .ci)
    }
    
    ## Indicator correlation ---------------------------------------------------
    if(!is.null(infer_out$Indicator_correlation)) {
      indicator_correlation <- addInfer(infer_out$Indicator_correlation, 
                                        indicator_correlation, .ci)
    }
    
    ## Effects -----------------------------------------------------------------
    if(x2$Model$model_type == "Linear" & !all(x2$Model$structural == 0)) {
      which_effects <- intersect(names(effects), c("Total_effect", "Indirect_effect"))
      if(length(which_effects) > 0) {
        effects <- lapply(which_effects, function(nx) {
          temp   <- infer_out[[nx]]
          x <- effects[[nx]]
          t_temp <- x$Estimate / temp$sd
          
          x["Std_err"] <- temp$sd
          x["t_stat"]  <- t_temp
          x["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
          
          if(!is.null(.ci)) {
            ## Add CI's
            # Column names
            ci_colnames <- paste0(rep(names(temp[.ci]), sapply(temp[.ci], function(x) nrow(x))), ".",
                                  unlist(lapply(temp[.ci], rownames)))
            
            # Add cis to data frame and set names
            x <- cbind(x, t(do.call(rbind, temp[.ci])))
            rownames(x) <- NULL
            colnames(x)[(length(colnames(x)) - 
                           (length(ci_colnames) - 1)):length(colnames(x))] <- ci_colnames
          }
          ## Return
          x
        })
        names(effects) <- which_effects
        
        if(any(which_effects == "Total_effect")) {
          .object$Estimates$Effect_estimates$Total_effect <- effects[["Total_effect"]] # Total effect!
        }
        if(any(which_effects == "Indirect_effect")) {
          .object$Estimates$Effect_estimates$Indirect_effect <- effects[["Indirect_effect"]] # Indirect effect 
        }
      }          
    }
  }
  
  ### Modify relevant .object elements =========================================
  .object$Estimates$Path_estimates    <- path_estimates
  .object$Estimates$Loading_estimates <- loading_estimates
  .object$Estimates$Weight_estimates  <- weight_estimates
  .object$Estimates$D2                <- D2
  .object$Estimates$Residual_correlation <- residual_correlation
  .object$Estimates$Indicator_correlation <- indicator_correlation
  .object$Estimates$Exo_construct_correlation <- exo_construct_correlation
  
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    .object$Estimates$Inner_weight_estimates <- inner_weight_estimates
  }
  .object$Estimates$Construct_scores <- construct_scores
  
  ## Set class for printing and return
  
  class(.object) <- if(inherits(.object, "cSEMResults_resampled")) {
    c("cSEMSummarize", "cSEMSummarize_default", "cSEMSummarize_resampled")
  } else {
    c("cSEMSummarize", "cSEMSummarize_default")
  }
  
  return(.object)
}

#' @export

summarize.cSEMResults_multi <- function(
  .object                = NULL, 
  .alpha                 = 0.05,
  .ci                    = NULL,
  ...
) {
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, summarize.cSEMResults_2ndorder, 
           .alpha, .ci, ...)
  } else {
    lapply(.object, summarize.cSEMResults_default, 
           .alpha, .ci, ...)
  }
}

#' @export

summarize.cSEMResults_2ndorder <- function(
  .object                = NULL, 
  .alpha                 = 0.05,
  .ci                    = NULL,
  ...
) {
  
  ## Run summarize for each stage
  x <- lapply(.object, summarize.cSEMResults_default)

  x11 <- x$First_stage$Estimates
  x12 <- x$First_stage$Information
  
  x21 <- x$Second_stage$Estimates
  x22 <- x$Second_stage$Information
  
  ### Modify second stage estimates ============================================
  ## Path estimates 
  # Only second order path model estimates are relevant. Delete the "_temp" 
  # suffix.
  x21$Path_estimates$Name <- gsub("_temp", "", x21$Path_estimates$Name)
  
  ## Loading estimates 
  # Loadings for first stage (=all loadings except those of the 2nd order constructs)
  # Loadings for second stage (all loadings for 2nd order constructs)
  i <- x21$Loading_estimates$Name[!grepl("_temp", x21$Loading_estimates$Name)]
  x21$Loading_estimates <- x21$Loading_estimates[x21$Loading_estimates$Name %in% i, , drop = FALSE]
  
  ## Weights estimates 
  # Weights for first stage (=all weights except those of the 2nd order constructs)
  # Weights for second stage (all weights for 2nd order constructs)
  i <- x21$Weight_estimates$Name[!grepl("_temp", x21$Weight_estimates$Name)]
  x21$Weight_estimates <- x21$Weight_estimates[x21$Weight_estimates$Name %in% i, , drop = FALSE]
  
  ## Exogenous construct correlation matrix
  x21$Exo_construct_correlation$Name <- gsub("_temp", "", x21$Exo_construct_correlation$Name)
  
  ## Delete _temp suffix wherever it appears ----------------------------------
  ## Inner weight estimates
  x21$Inner_weight_estimates$Name <- gsub("_temp", "", x21$Inner_weight_estimates$Name)
  
  ## Construct scores
  colnames(x21$Construct_scores) <- gsub("_temp", "", colnames(x21$Construct_scores))
  
  ## Proxy VCV
  colnames(x21$Proxy_VCV) <- gsub("_temp", "", colnames(x21$Proxy_VCV))
  rownames(x21$Proxy_VCV) <- gsub("_temp", "", rownames(x21$Proxy_VCV))
  
  ## Construct VCV
  colnames(x21$Construct_VCV) <- gsub("_temp", "", colnames(x21$Construct_VCV))
  rownames(x21$Construct_VCV) <- gsub("_temp", "", rownames(x21$Construct_VCV))
  
  ## Reliabilities
  names(x21$Reliabilities) <- gsub("_temp", "", names(x21$Reliabilities))
  
  ## R2
  names(x21$R2) <- gsub("_temp", "", names(x21$R2))
  names(x21$R2adj) <- gsub("_temp", "", names(x21$R2adj))
  
  ## Vif
  names(x21$VIF) <- gsub("_temp", "", names(x21$VIF))
  lapply(x21$VIF, function(x) {
    names(x) <- gsub("_temp", "", names(x)) 
    x
  })
  
  ## Model
  rownames(x22$Model$structural) <- gsub("_temp", "", rownames(x22$Model$structural))
  colnames(x22$Model$structural) <- gsub("_temp", "", colnames(x22$Model$structural))
  
  rownames(x22$Model$measurement) <- gsub("_temp", "", rownames(x22$Model$measurement))
  colnames(x22$Model$measurement) <- gsub("_temp", "", colnames(x22$Model$measurement))
  
  names(x22$Model$construct_type) <- gsub("_temp", "", names(x22$Model$construct_type))
  names(x22$Model$construct_order) <- gsub("_temp", "", names(x22$Model$construct_order))
  ## Effects -------------------------------------------------------------------
  if(x22$Model$model_type == "Linear") {
    effects <- calculateEffects(.object, .output_type = "data.frame")
    # Add effects to second stage
    x21$Effect_estimates <- effects 
  }
  
  ## Inference =================================================================
  if(inherits(.object, "cSEMResults_resampled")) {
    
    ## Check arguments
    match.arg(.ci, args_default(.choices = TRUE)$.ci, several.ok = TRUE)
    
    ## The 95% percentile CI is returned by default (BCa is acutally better 
    ##  but requires a second resampling run and should therefore not be the 
    ##  default).
    if(is.null(.ci)) .ci <- "CI_percentile"
    
    quantity <- c("sd", .ci)
    
    infer_out <- infer(
      .object          = .object,
      .alpha           = .alpha,
      .quantity        = quantity,
      ...
    )
    
    ### Path estimates ---------------------------------------------------------
    temp   <- infer_out$Path_estimates
    t_temp <- x21$Path_estimates$Estimate / temp$sd
    
    x21$Path_estimates["Std_err"] <- temp$sd
    x21$Path_estimates["t_stat"]  <- t_temp
    x21$Path_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    
    if(!is.null(.ci)) {
      ## Add CI's
      # Column names
      ci_colnames <- paste0(rep(names(temp[.ci]), sapply(temp[.ci], function(x) nrow(x))), ".",
                            unlist(lapply(temp[.ci], rownames)))
      
      # Add cis to data frame and set names
      x21$Path_estimates <- cbind(x21$Path_estimates, t(do.call(rbind, temp[.ci])))
      rownames(x21$Path_estimates) <- NULL
      colnames(x21$Path_estimates)[(length(colnames(x21$Path_estimates)) - 
                                      (length(ci_colnames) - 1)):length(colnames(x21$Path_estimates))] <- ci_colnames
    }
    
    ### Loadings ----------------------
    L <- lapply(list(x11$Loading_estimates, x21$Loading_estimates), function(y) {
      
      temp   <- infer_out$Loading_estimates
      t_temp <- y$Estimate / temp$sd[y$Name]
      
      y["Std_err"] <- temp$sd[y$Name]
      y["t_stat"]  <- t_temp
      y["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
      
      if(!is.null(.ci)) {
        ## Add CI's
        # Add cis to data frame and set names
        ci_temp <- t(do.call(rbind, temp[.ci]))[y$Name, , drop =FALSE]
        y <- cbind(y, ci_temp)
        rownames(y) <- NULL
        colnames(y)[(length(colnames(y)) - (length(ci_colnames) - 1)):length(colnames(y))] <- ci_colnames
      }
      ## Return y
      y
    }) 
    
    ### Weights --------------------
    W <- lapply(list(x11$Weight_estimates, x21$Weight_estimates), function(y) {
      
      temp   <- infer_out$Weight_estimates
      t_temp <- y$Estimate / temp$sd[y$Name]
      
      y["Std_err"] <- temp$sd[y$Name]
      y["t_stat"]  <- t_temp
      y["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
      
      if(!is.null(.ci)) {
        ## Add CI's
        # Add cis to data frame and set names
        ci_temp <- t(do.call(rbind, temp[.ci]))[y$Name, , drop =FALSE]
        y <- cbind(y, ci_temp)
        rownames(y) <- NULL
        colnames(y)[(length(colnames(y)) - (length(ci_colnames) - 1)):length(colnames(y))] <- ci_colnames
      }
      ## Return y
      y
    }) 
    
    ## Update
    x11$Loading_estimates <- L[[1]]
    x21$Loading_estimates <- L[[2]]
    x11$Weight_estimates  <- W[[1]]
    x21$Weight_estimates  <- W[[2]]
    
    ### Exo construct correlation -----------------------------------------------
    temp   <- infer_out$Exo_construct_correlation
    if(!is.null(temp)) {
      t_temp <- x21$Exo_construct_correlation$Estimate / temp$sd
      
      x21$Exo_construct_correlation["Std_err"] <- temp$sd
      x21$Exo_construct_correlation["t_stat"]  <- t_temp
      x21$Exo_construct_correlation["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
      
      if(!is.null(.ci)) {
        ## Add CI's
        # Add cis to data frame and set names
        x21$Exo_construct_correlation <- cbind(x21$Exo_construct_correlation, t(do.call(rbind, temp[.ci])))
        rownames(x21$Exo_construct_correlation) <- NULL
        colnames(x21$Exo_construct_correlation)[(length(colnames(x21$Exo_construct_correlation)) - 
                                              (length(ci_colnames) - 1)):length(colnames(x21$Exo_construct_correlation))] <- ci_colnames
      }
    }
    
    ### Residual correlation ---------------------------------------------------
    temp   <- infer_out$Residual_correlation
    if(!is.null(temp)) {
      t_temp <- x11$Residual_correlation$Estimate / temp$sd
      
      x11$Residual_correlation["Std_err"] <- temp$sd
      x11$Residual_correlation["t_stat"]  <- t_temp
      x11$Residual_correlation["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
      
      if(!is.null(.ci)) {
        ## Add CI's
        # Add cis to data frame and set names
        x11$Residual_correlation <- cbind(x11$Residual_correlation, t(do.call(rbind, temp[.ci])))
        rownames(x11$Residual_correlationn) <- NULL
        colnames(x11$Residual_correlation)[(length(colnames(x11$Residual_correlation)) - 
                                          (length(ci_colnames) - 1)):length(colnames(x11$Residual_correlation))] <- ci_colnames
      }
    }
    
    ## Indicator correlation ---------------------------------------------------
    temp   <- infer_out$Indicator_correlation
    if(!is.null(temp)) {
      t_temp <- x11$Indicator_correlation$Estimate / temp$sd
      
      x11$Indicator_correlation["Std_err"] <- temp$sd
      x11$Indicator_correlation["t_stat"]  <- t_temp
      x11$Indicator_correlation["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
      
      if(!is.null(.ci)) {
        ## Add CI's
        # Add cis to data frame and set names
        x11$Indicator_correlation <- cbind(x11$Indicator_correlation, t(do.call(rbind, temp[.ci])))
        rownames(x11$Indicator_correlationn) <- NULL
        colnames(x11$Indicator_correlation)[(length(colnames(x11$Indicator_correlation)) - 
                                              (length(ci_colnames) - 1)):length(colnames(x11$Indicator_correlation))] <- ci_colnames
      }
    }
    
    ## Effects -----------------------------------------------------------------
    if(x22$Model$model_type == "Linear") {
      which_effects <- intersect(names(effects), c("Total_effect", "Indirect_effect"))
      if(length(which_effects) > 0) {
        effects <- lapply(which_effects, function(nx) {
          temp   <- infer_out[[nx]]
          x <- effects[[nx]]
          t_temp <- x$Estimate / temp$sd
          
          x["Std_err"] <- temp$sd
          x["t_stat"]  <- t_temp
          x["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
          
          if(!is.null(.ci)) {
            ## Add CI's
            # Column names
            ci_colnames <- paste0(rep(names(temp[.ci]), sapply(temp[.ci], function(x) nrow(x))), ".",
                                  unlist(lapply(temp[.ci], rownames)))
            
            # Add cis to data frame and set names
            x <- cbind(x, t(do.call(rbind, temp[.ci])))
            rownames(x) <- NULL
            colnames(x)[(length(colnames(x)) - 
                           (length(ci_colnames) - 1)):length(colnames(x))] <- ci_colnames
          }
          ## Return
          x
        })
        names(effects) <- which_effects
        
        if(any(which_effects == "Total_effect")) {
          x21$Effect_estimates$Total_effect <- effects[["Total_effect"]] # Total effect!
        }
        if(any(which_effects == "Indirect_effect")) {
          x21$Effect_estimates$Indirect_effect <- effects[["Indirect_effect"]] # Indirect effect 
        }
      }
    }
  } # END resampled
  
  out <- list("First_stage"  = list("Estimates" = x11, "Information" = x12), 
              "Second_stage" = list("Estimates" = x21, "Information" = x22))
  ## Set class for list elements
  class(out$First_stage) <- class(out$Second_stage) <- c("cSEMSummarize", "cSEMSummarize_default")
  
  ## Set class for printing and return
  class(out) <- if(inherits(.object, "cSEMResults_resampled")) {
    c("cSEMSummarize", "cSEMSummarize_2ndorder", "cSEMSummarize_resampled")
  } else {
    c("cSEMSummarize", "cSEMSummarize_2ndorder")
  }
  
  return(out)
}
