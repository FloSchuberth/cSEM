#' Calculate direct, indirect and total effect
#' Helper for generic function summarize()
#' @noRd
#' 
calculateEffects <- function(.object, .output_type = c("data.frame", "matrix")) {
  
  output_type <- match.arg(.output_type)
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    m <- .object$Second_stage$Information$Model
    
    ## Matrix of direct effects (second stage):
    direct <- .object$Second_stage$Estimates$Path_estimates
    
    ## Rename
    colnames(direct) <- gsub("_temp", "", colnames(direct))
    rownames(direct) <- gsub("_temp", "", rownames(direct))
    
  } else {
    m <- .object$Information$Model
    
    ## Matrix of direct effects:
    direct <- .object$Estimates$Path_estimates
  }

  ## Endogenous (lhs) variables
  vars_endo <- rownames(m$structural)[rowSums(m$structural) != 0]
  
  ## Matrix of total total effects: B = direct
  # Note: eta = B x eta + zeta
  #       (I - B)*eta = zeta
  
  B_star <- diag(nrow(direct)) - direct
  total <- solve(B_star) - diag(nrow(direct))
  
  ## Matrix of indirect effects:
  indirect <- total - direct
  
  # MAtrix containing the variance accounted for (VAR)
  VAR <- indirect/total
  VAR[which(is.nan(VAR), arr.ind = TRUE)] <- 0
  
  out <- list(
    "Direct_effect"          = direct, 
    "Indirect_effect"        = indirect, 
    "Total_effect"           = total,
    "Variance_accounted_for" = VAR
  )
  
  ## Convert to data frame
  if(output_type == "data.frame") {
    
    # Lookup table for the names
    temp <- outer(rownames(direct), colnames(direct), 
                  FUN = function(x, y) paste(x, y, sep = " ~ "))
    
    lout <- lapply(out, function(x) {
      # Get construct type for relevant variables
      # Note: round() in the formula below is necessary as calculateEffect() can
      #       sometimes produce values that are not exactly 0 due to floating point
      #       imprecisions. To round to 10 digits is rather arbitrary. Maybe
      #       there is a better way to check if a number is "acutally zero".
      type <- rep(m$construct_type, times = rowSums(round(x, 10) != 0))
      
      # If there are no indirect effects the matrix "indirect" is the zero matrix
      # return an empty data frame in this case
      if(all(round(x, 10) == 0)) {
        data.frame(
          "Name"           = character(),
          "Construct_type" = character(),
          "Estimate"       = double(),
          "Std_err"        = double(),
          "t_stat"         = double(),
          "p_value"        = double(),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      } else {
        data.frame(
          "Name"           = t(temp)[round(t(x), 10) != 0 ],
          "Construct_type" = type,
          "Estimate"       = t(x)[round(t(x), 10) != 0 ],
          "Std_err"        = NA,
          "t_stat"         = NA,
          "p_value"        = NA,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }

    })
    
    out[["Direct_effect"]]   <- lout[[1]]
    out[["Indirect_effect"]] <- lout[[2]] 
    out[["Total_effect"]]    <- lout[[3]]
    out[["Variance_accounted_for"]] <- lout[[4]]
  }
  
  ## Return
  out
}