# Empirical Bayes correction for loadings suggested by Dijkstra (2018) 
#
# Last modified: 07.2019

## Notes
# 1. Function is a post-hoc function. Written to be applied to a cSEMResults
#    object that produced inadmissble results in a sense that the loadings
#    are larger than 1 and (although only partly implemented) the construct
#    correlation matrix is semi-positve definite.
# 2. Not clear how to tackle inadmissible results due to not semi-positive 
#    construct correlation matrix. not semi-positve definite does not mean
#    that correlations are larger/smaller than 1 /-1.
#    Idea by Flo: do a singluar value decomposition and correct the eigen-
#                 values, reassemble the matrix and continue. Whether this 
#                 is statistically sound and desierable at all is an open question.
# 2. Not clear how to handle reliabilities. 
#    - if loadings are corrected: reliabilities change --> construct correlations
#      change --> path coefficients change. Not clear if reliabilities should be
#      adjusted as well or whether we simply correct the loadings and nothing else
#    - if construct correlations are corrected, path coefficients change but for 
#      nonlinear models we would need the reliabilities.

correctInadmissibles <- function(
  .object = NULL, 
  .method = c("median","mean")
) {
  
  .object <- res1
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    
    stop("Not implemented yet")
    
  } else {
    x2 <- .object$Estimates
    sd <- cSEM:::SdResample(.object$Estimates$Estimates_resample$Estimates1,
                            .resample_method = "bootstrap") 
  }
  
  if(verify(.object)["2"]) {
    ### Loadings larger than 1 in absolute value
    L    <- x2$Loading_estimates[.object$Information$Model$measurement != 0]
    L_sd <- sd$Loading_estimates
    
    for(i in which(abs(L) > 1)) {
      A <- pnorm((-1 - L[i]) / L_sd[i])
      B <- pnorm((1 - L[i]) / L_sd[i]) 
      
      ## Replace VCV elements
      if(.method == "median"){
        ## Overwrite the old loadings
        L[i] <- L[i] + L_sd[i] * qnorm(mean(c(A , B)))
      }
      
      if(.method == "mean"){
        ## Overwrite the old loadings
        L[i] <- L[i] + L_sd[i] / (B - A) * 
          (dnorm((-1 - L[i]) / L_sd[i]) - 
             dnorm((1 - L[i]) / L_sd[i]))
      }
    }
    # Overwrite original construct VCV
    .object$Estimates$Loading_estimates[.object$Information$Model$measurement != 0] <- L
    
    # Recompute reliabilities based on the new loadings
    L  <- .object$Estimates$Loading_estimates
    W  <- .object$Estimates$Weight_estimates
    Q2 <- .object$Estimates$Reliabilities
    
    for(j in rownames(L)) {
      Q2[j] <- c(W[j, ] %*% L[j, ])^2
    }
    
    .object$Estimates$Reliabilities <- Q2
  }
  
  if(verify(.object)["3"]) {
    ### Construct VCV not positive semi-definite (because construct correlations 
    ### are larger than 1)
    VCV <- x2$Construct_VCV
    VCV_sd <- matrix(sd$User_fun, nrow = nrow(VCV), ncol = ncol(VCV),
                     dimnames = list(rownames(VCV), colnames(VCV)))
    
    for(j in rownames(VCV)) {
      for(i in colnames(VCV)) {
        if(abs(VCV[j, i]) > 1) {
          A <- pnorm((-1 - VCV[j, i]) / VCV_sd[j, i])
          B <- pnorm((1 - VCV[j, i]) / VCV_sd[j, i]) 
          
          ## Replace VCV elements
          if(.method == "median"){
            ## Overwrite the old loadings
            VCV[j, i] <- VCV[j, i] + VCV_sd[j, i] * qnorm(mean(c(A , B)))
          }
          
          if(.method == "mean"){
            ## Overwrite the old loadings
            VCV[j, i] <- VCV[j, i] + VCV_sd[j, i] / (B - A) * 
              (dnorm((-1 - VCV[j, i]) / VCV_sd[j, i]) - 
                 dnorm((1 - VCV[j, i]) / VCV_sd[j, i]))
          }
        }
      }
    }
    # Overwrite original construct VCV
    .object$Estimates$Construct_VCV <- VCV
    
    # Reestimate path model
    P <- .object$Estimates$Construct_VCV
    
    path <- cSEM:::estimatePath(
      .approach_nl    = .object$Information$Arguments$.approach_nl,
      .approach_paths = .object$Information$Arguments$.approach_paths,
      .csem_model     = .object$Information$Arguments$.model,
      .H              = .object$Estimates$Construct_scores,
      .normality      = .object$Information$Arguments$.normality,
      .P              = .object$Estimates$Construct_VCV,
      .Q              = .object$Estimates$Reliabilities # only for nonlinear models
    )
    
    .object$Estimates$Path_estimates <- path
  }
  
  ## Return
  .object
}