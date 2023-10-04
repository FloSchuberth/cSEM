#' Internal: Calculate direct, indirect and total effect
#' 
#' The direct effects are equal to the estimated coefficients. The total effect 
#' equals (I-B)^{-1}Gamma. The indirect effect equals the difference between
#' the total effect and the indirect effect. In addition, the variance accounted
#' for (VAF) is calculated. The VAF is defined as the ratio of a variables
#' indirect effect to its total effect. Helper for generic functions [summarize()] and [assess()].
#' 
#' @usage calculateEffects(
#'  .object       = NULL,
#'  .output_type  = c("data.frame", "matrix")
#' )
#'
#' @return A matrix or a data frame of effects.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess()], [summarize()] [cSEMResults]
#'
#' @keywords internal

calculateEffects <- function(.object = NULL, .output_type = c("data.frame", "matrix")) {
  
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
  
  # check whether system is non-recursive (Bollen, 1987):
  
  # k: number of variables in the structural model
  k=nrow(m$structural)
  temp=m$structural
  for (i in 1:k) {
    temp <- temp %*% m$structural
  }
 
  nonrec=all(temp == 0)
  
   
  # # Create indicator matrix for total and indirect effects; if model without reciprocal relationship
  # Otherwise it does not work
  
  if(nonrec){
  Btemp <- diag(nrow(m$structural)) - m$structural

  totaltemp <- solve(Btemp) - diag(nrow(m$structural))
  totalInd = matrix(0,nrow=nrow(totaltemp), ncol=ncol(totaltemp))
  totalInd[totaltemp!=0] <- 1

  indirecttemp = totaltemp - m$structural
  indirectInd <- matrix(0,nrow=nrow(indirecttemp), ncol=ncol(indirecttemp))
  indirectInd[indirecttemp!=0] <- 1
  }
  # Matrix containing the variance accounted for (VAF)
  VAF <- indirect/total
  VAF[which(is.nan(VAF), arr.ind = TRUE)] <- 0
  
  out <- list(
    "Direct_effect"          = direct, 
    "Indirect_effect"        = indirect, 
    "Total_effect"           = total,
    "Variance_accounted_for" = VAF
  )
  
  ## Convert to data frame
  if(output_type == "data.frame") {
    
    # Lookup table for the names
    temp <- outer(rownames(direct), colnames(direct), 
                  FUN = function(x, y) paste(x, y, sep = " ~ "))
    
    
    lout <- lapply(1:length(out),function(x){
      
        # Get construct type for relevant variables
        # Note: round() in the formula below is necessary as calculateEffects() can
        #       sometimes produce values that are not exactly 0 due to floating point
        #       imprecision. To round to 10 digits is rather arbitrary. Maybe
        #       there is a better way to check if a number is "actually zero".
        # type <- rep(m$construct_type, times = rowSums(round(x, 10) != 0))

        # Note (07.08.2019): Rounding may be confusing. I was using the repeated
        #       indicators approach for the model
        #       c4 ~ eta1
        #       eta2 ~ eta1 + c4
        #       and expecting there to be an effect of eta1 on c4 (knowingly that
        #       it should be zero). Round killed it, but it should be there as it
        #       is an effect that is falsely estimated to zero.
        # type <- rep(m$construct_type, times = rowSums(x != 0))
        
        # Note (04.10.2023): To address the rounding problem, the binary indicator matrices are added in case of non-recursive models
        # before rounding. In case of recursive models, it is just rounded. 

      # Choose binary indicator matrix
        IndMat = if(nonrec){
          if(names(out[x]) %in% c("Direct_effect")){
          m$structural
        } else if(names(out[x]) %in% c("Indirect_effect")){
          indirectInd
        } else if(names(out[x]) %in% c("Total_effect")){
          totalInd
        }else{
          0
        }} else {
          0
        }
      
        type <- rep(m$construct_type, times = rowSums(round(out[[x]] + IndMat, 10) != 0))

      if(all(round(out[[x]]+ IndMat, 10) == 0)) {
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
          "Name"           = t(temp)[round(t(out[[x]] + IndMat), 10) != 0 ],
          "Construct_type" = type,
          "Estimate"       = t(out[[x]])[round(t(out[[x]] + IndMat), 10) != 0 ],
          # "Name"           = t(temp)[t(x) != 0 ],
          # "Construct_type" = type,
          # "Estimate"       = t(x)[t(x) != 0 ],
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
  # Delete matrices/data frame that are empty
  if(output_type == "matrix") {
    sumout <- function(x) sum(x) != 0
  } else {
    sumout <- function(x) nrow(x) > 0
  }
  out <- Filter(sumout, out)
  out
}

#' Helper for summarize
#' @noRd

addInfer <- function(.what = NULL, .estimates = NULL, .ci = NULL) {
  temp <- .what
  t_temp <- .estimates$Estimate / temp$sd
  
  .estimates["Std_err"] <- temp$sd
  .estimates["t_stat"]  <- t_temp
  .estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
  
  if(!is.null(.ci)) {
    ## Add CI's
    # Column names
    ci_colnames <- paste0(rep(names(temp[.ci]), sapply(temp[.ci], function(x) nrow(x))), ".",
                          unlist(lapply(temp[.ci], rownames)))
    
    # Add cis to data frame and set names
    .estimates <- cbind(.estimates, t(do.call(rbind, temp[.ci])))
    rownames(.estimates) <- NULL
    colnames(.estimates)[(length(colnames(.estimates)) - 
                                (length(ci_colnames) - 1)):length(colnames(.estimates))] <- ci_colnames
  }
  .estimates
}