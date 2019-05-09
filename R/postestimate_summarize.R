#' Summarize model
#'
#' Summarize the model (TODO). 
#' 
#' Summary (TODO)
#'
#' @usage summarize(
#'  .object = NULL, 
#'  .alpha  = args_default()$.alpha,
#'  .ci     = NULL,
#'  ...
#'  )
#'
#' @inheritParams csem_arguments
#' @param ... Further arguments to `summarize()`. Currently ignored.
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
summarize <- function(
  .object                = NULL, 
  .alpha                 = args_default()$.alpha,
  .ci                    = NULL,
  ...
  ) {
  UseMethod("summarize")
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_default <- function(
  .object                = NULL, 
  .alpha                 = args_default()$.alpha,
  .ci                    = NULL,
  ...
  ) {
  
  x1  <- .object$Estimates
  x2  <- .object$Information

  ### Structure output =========================================================
  ## Path estimates ------------------------------------------------------------
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
  
  ## Inner weight estimates ----------------------------------------------------
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    i <- rownames(x1$Inner_weight_estimates)
    D <- x2$Model$structural[i, i , drop = FALSE] + t(x2$Model$structural[i, i, drop = FALSE])
    
    temp <- outer(i, i, FUN = function(x, y) paste(x, y, sep = " -- "))
    type <- rep(x2$Model$construct_type, times = colSums(D))
    
    inner_weight_estimates <- data.frame(
      "Name"           = t(temp)[t(D) != 0],
      "Construct_type" = type, 
      "Estimate"       = t(x1$Inner_weight_estimates)[t(D) != 0 ], 
      stringsAsFactors = FALSE)
  }

  ## Construct scores ----------------------------------------------------------
  construct_scores <- as.data.frame(x1$Construct_scores)

  ## Inference =================================================================
  
  if(inherits(.object, "cSEMResults_resampled")) {
    ## Check arguments
    match.arg(.ci, args_default(.choices = TRUE)$.ci, several.ok = TRUE)
    
    quantity <- c("sd", .ci)
    
    infer_out <- infer(
      .object = .object,
      .alpha           = .alpha,
      .quantity        = quantity,
      ...
    )
    
    ## Path estimates ----------------------------------------------------------
    temp   <- infer_out$Path_estimates
    t_temp <- t(x1$Path_estimates)[t(x2$Model$structural) != 0 ] / temp$sd
    
    path_estimates["Std_err"] <- temp$sd
    path_estimates["t_stat"]  <- t_temp
    path_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    
    if(!is.null(.ci)) {
      ## Add CI's
      # Column names
      ci_colnames <- paste0(rep(names(temp[.ci]), sapply(temp[.ci], function(x) nrow(x))), ".",
             unlist(lapply(temp[.ci], rownames)))
      
      # Add cis to data frame and set names
      path_estimates <- cbind(path_estimates, t(do.call(rbind, temp[.ci])))
      rownames(path_estimates) <- NULL
      colnames(path_estimates)[(length(colnames(path_estimates)) - 
                                  (length(ci_colnames) - 1)):length(colnames(path_estimates))] <- ci_colnames
    }
    
    ## Loading estimates -------------------------------------------------------
    temp   <- infer_out$Loading_estimates
    t_temp <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0 ] / temp$sd
    
    loading_estimates["Std_err"] <- temp$sd
    loading_estimates["t_stat"]  <- t_temp
    loading_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    
    if(!is.null(.ci)) {
      ## Add CI's
      # Add cis to data frame and set names
      loading_estimates <- cbind(loading_estimates, t(do.call(rbind, temp[.ci])))
      rownames(loading_estimates) <- NULL
      colnames(loading_estimates)[(length(colnames(loading_estimates)) - 
                                  (length(ci_colnames) - 1)):length(colnames(loading_estimates))] <- ci_colnames
    }
    
    ## Weight estimates --------------------------------------------------------
    temp   <- infer_out$Weight_estimates
    t_temp <- t(x1$Weight_estimates)[t(x2$Model$measurement) != 0 ] / temp$sd
    
    weight_estimates["Std_err"] <- temp$sd
    weight_estimates["t_stat"]  <- t_temp
    weight_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    
    if(!is.null(.ci)) {
      ## Add CI's
      # Add cis to data frame and set names
      weight_estimates <- cbind(weight_estimates, t(do.call(rbind, temp[.ci])))
      rownames(weight_estimates) <- NULL
      colnames(weight_estimates)[(length(colnames(weight_estimates)) - 
                                  (length(ci_colnames) - 1)):length(colnames(weight_estimates))] <- ci_colnames
    }
  }
  
  ### Modify relevant .object elements =========================================
  .object$Estimates$Path_estimates    <- path_estimates
  .object$Estimates$Loading_estimates <- loading_estimates
  .object$Estimates$Weight_estimates  <- weight_estimates
  if(x2$Arguments$.approach_weights == "PLS-PM") {
  .object$Estimates$Inner_weight_estimates <- inner_weight_estimates
  }
  .object$Estimates$Construct_scores <- construct_scores
  
  
  ## Set class for printing and return
  
  class(.object) <- if(inherits(.object, "cSEMResults_resampled")) {
    c("cSEMSummarize_default", "cSEMSummarize_resampled")
  } else {
    "cSEMSummarize_default"
  }
  
  return(.object)
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_multi <- function(
  .object                = NULL, 
  .alpha                 = args_default()$.alpha,
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

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_2ndorder <- function(
  .object                = NULL, 
  .alpha                 = args_default()$.alpha,
  .ci                    = NULL,
  ...
) {
  
  ## Run summarize for each stage
  x <- lapply(.object, summarize.cSEMResults_default)
  
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
  # Wieghts for first stage (=all weights except those of the 2nd order constructs)
  # Weights for second stage (all weights for 2nd order constructs)
  i <- x21$Weight_estimates$Name[!grepl("_temp", x21$Weight_estimates$Name)]
  x21$Weight_estimates <- x21$Weight_estimates[x21$Weight_estimates$Name %in% i, , drop = FALSE]
  
  ## Inner weight estimates
  # Rename. Delete the "_temp" suffix.
  # x21$Inner_weight_estimates$Name <- gsub("_temp", "", x21$Inner_weight_estimates$Name)
  
  ## Construct scores
  # colnames(x21$Construct_scores) <- gsub("_temp", "", colnames(x21$Construct_scores))
  
  if(any(class(.object) == "cSEMResults_resampled")) {
    
    ## Check arguments
    match.arg(.ci, args_default(.choices = TRUE)$.ci, several.ok = TRUE)
    
    quantity <- c("sd", .ci)
    
    infer_out <- infer(
      .object = .object,
      .alpha           = .alpha,
      .quantity        = quantity,
      ...
    )
    
    ## Path coefficients
    
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
    # 
    # ## Loadings
    # temp   <- lapply(infer_out$Loading_estimates, function(x) x[]
    # t_temp <- x21$Loading_estimates$Estimate / temp$sd
    # 
    # x21$Loading_estimates["Std_err"] <- temp$sd
    # x21$Loading_estimates["t_stat"]  <- t_temp
    # x21$Loading_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    # 
    # if(!is.null(.ci)) {
    #   ## Add CI's
    #   # Add cis to data frame and set names
    #   x21$Loading_estimates <- cbind(x21$Loading_estimates, t(do.call(rbind, temp[.ci])))
    #   rownames(x21$Loading_estimates) <- NULL
    #   colnames(x21$Loading_estimates)[(length(colnames(x21$Loading_estimates)) - 
    #                                  (length(ci_colnames) - 1)):length(colnames(x21$Loading_estimates))] <- ci_colnames
    # }
    # 
    # ## Weight estimates
    # temp   <- infer_out$Weight_estimates
    # t_temp <- x21$Weight_estimates / temp$sd
    # 
    # x21$Weight_estimates["Std_err"] <- temp$sd
    # x21$Weight_estimates["t_stat"]  <- t_temp
    # x21$Weight_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
    # 
    # if(!is.null(.ci)) {
    #   ## Add CI's
    #   # Add cis to data frame and set names
    #   x21$Weight_estimates <- cbind(x21$Weight_estimates, t(do.call(rbind, temp[.ci])))
    #   rownames(x21$Weight_estimates) <- NULL
    #   colnames(x21$Weight_estimates)[(length(colnames(x21$Weight_estimates)) - 
    #                                 (length(ci_colnames) - 1)):length(colnames(x21$Weight_estimates))] <- ci_colnames
    # }
  }
  
  ## Set class for printing and return
  out <- list("First_stage" = x$First_stage, 
              "Second_stage" = list("Estimates" = x21, "Information" = x22))
  
  class(out) <- "cSEMSummarize_2ndorder"
  
  if(inherits(.object) == "cSEMResults_resampled") {
    class(.object) <- c(class(.object), "cSEMSummarize_resampled" )
  }
  
  return(out)
}
