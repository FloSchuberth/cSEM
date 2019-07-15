################################################################################
#
#   Purpose: central file to be sourced to several test-xxx.R files
#
################################################################################
## Function to compare path and/or loading and/or weight estimates from a cSEM object
## to a vector of population parameters
comparecSEM <- function(.object, .what, .pop_parameters) {
  # .object: cSEM object
  # .what: what to compare
  # .pop_parameters: a vector of population values
  
  x <- cSEM::summarize(.object)
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    x1 <- x$First_stage$Estimates
    x2 <- x$Second_stage$Estimates
  } else {
    x1 <- NULL
    x2 <- x$Estimates
  }
  
  if(.what == "Path_estimates") {
    est <- x2$Path_estimates[, c("Name", "Estimate")]
    
  } else if(.what == "Loading_estimates") {
    est <- rbind(x1$Loading_estimates[, c("Name", "Estimate")],
                 x2$Loading_estimates[, c("Name", "Estimate")])
    
  } else if(.what == "Weight_estimates") {
    ## Compare only weights for composites, since only those weights are
    ## specified when creating the DGP
    x1$Weight_estimates
    est <- rbind(
      x1$Weight_estimates[x1$Weight_estimates$Construct_type == "Composite", 
                          c("Name", "Estimate")],
      x2$Weight_estimates[x2$Weight_estimates$Construct_type == "Composite", 
                          c("Name", "Estimate")])
    
  } else {
    stop("Error") 
  }
  
  data.frame(est, 
             "Pop_value" = unname(.pop_parameters),
             "Pop_value_name" = names(.pop_parameters),
             stringsAsFactors = FALSE)
}