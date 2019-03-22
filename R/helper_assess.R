#' Internal: AVE
#'
#' Calculate the average variance extracted (AVE). See (TODO) for details.
#'
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateAVE(
#'  .object              = NULL,
#'  .only_common_factors = TRUE,
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

calculateAVE <- function(
  .object              = NULL,
  .only_common_factors = args_default()$.only_common_factors
  ){

  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  c_types <- .object$Information$Model$construct_type
  c_names <- names(c_types)
  
  ## Extract loadings
  Lambda  <- .object$Estimates$Loading_estimates
  x <- c_names
  
  ## Calculate AVE (extracted variance / total variance)
  # Note: Within the cSEM package indicators are always standardized (e.g, have
  #       a variance of 1). Therefore the term (sum(lambda^2) + sum(1 - lambda^2) is
  #       simply equal to the number of indicator attached to a construct j.
  #       Hence AVE is simply the average over all lambda^2_k of construct j.
  #       Since for x_k standardized --> lambda^2_k := (indicator) reliability
  #       the AVE in cSEM is simply the average indicator reliability
  AVEs <- sapply(c_names, function(x){
    lambda <- c(Lambda[x, Lambda[x,] != 0])
    ave    <- sum(lambda^2) / (sum(lambda^2) + sum(1 - lambda^2))
    ave
  })
  
  names(AVEs) <- c_names
  
  # #By default AVE's for constructs modeled as composites are not returned
  if(.only_common_factors){
    co_names <- names(c_types[c_types == "Composite"])
    AVEs     <- AVEs[setdiff(c_names, co_names)]
  }
  
  # Return vector of AVEs
  return(AVEs) 
}



#' Internal: Reliability
#'
#' Compute several reliability measures. See the vignette for details.
#' 
#' Reliability is the consistency of measurement. Practically, reliability 
#' is empirically derived based on the classical true score measurement theory 
#' and defined as the ratio of the true score variance relative to the proxy 
#' (test score) variance, where the later is a weighted linear combinations of 
#' the indicators (i.e. a proxy or stand-in for the true score). 
#' 
#' The current literatue provides numerous reliability measures, based on 
#' (seemingly) different formulae with no systematic naming conventions. 
#' This package follows \insertCite{Cho2016;textual}{cSEM}'s attempt for a 
#' systematic approach to reliability by adopting his proposed unified naming 
#' convention. See \insertCite{Cho2016;textual}{cSEM} for details.
#' 
#' Since reliability is defined with respect to a classical true score measurement
#' model only constructs modeled as common factors are considered by default.
#' For constructs modeled as composites reliability may be estimated by setting
#' `.only_common_factors = FALSE`, however, it is unclear how to
#' interpret reliability in this case.
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#'
#' @references 
#' 
#' \insertAllCited{}
#' @keywords internal
#' @name reliability
NULL

#' @describeIn reliability Calculate the congeneric reliability, also known as
#'                         composite reliability or rho_A.
calculateRhoC <- function(
  .object              = NULL,
  .only_common_factors = args_default()$.only_common_factors) {

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  W <- .object$Estimates$Weight_estimates
  rhoC     <- rep(1, times = nrow(W))
  names(rhoC) <- rownames(W)
  
  # Get loadings
  Lambda  <- .object$Estimates$Loading_estimates
  
  # Compute the congeneric reliability by block 
  for(j in rownames(Lambda)) {
    rhoC[j]  <- c(W[j, ] %*% Lambda[j, ])^2 
  }
  
  c_types <- .object$Information$Model$construct_type
  c_names <- names(c_types)
  
  # By default only reliabilities for constructs model common factors are returned
  if(.only_common_factors){
    co_names <- names(c_types[c_types == "Composite"])
    rhoC     <- rhoC[setdiff(c_names, co_names)]
  }

  # Return named vector of reliabilities
  return(rhoC)
}

#' @describeIn reliability Calculate the tau-equivalent reliability, also known
#'                           as Cronbach alpha or coefficient alpha. Since
#'                           indicators are always standarized in `cSEM`, 
#'                           tau-equivalent reliability is identical to the
#'                           parallel reliability.
#' 
calculateRhoT <- function(
  .object              = NULL,
  .only_common_factors = args_default()$.only_common_factors) {
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  W <- .object$Estimates$Weight_estimates
  rhoT     <- rep(1, times = nrow(W))
  names(rhoT) <- rownames(W)
  
  # Get correlation matrix
  S  <- .object$Estimates$Indicator_VCV
  
  # Calculate tau-equivalent reliability by block
  for(j in rownames(W)) {
    indicators_names <- colnames(W[j, W[j,] != 0, drop = FALSE])
    S_jj             <- S[indicators_names, indicators_names]
    rho_bar_j <- mean(S_jj[upper.tri(S_jj)])
    rhoT[j]  <- rho_bar_j * sum(W[j, ])^2  
  }
  
  c_types <- .object$Information$Model$construct_type
  c_names <- names(c_types)
  
  # By default only reliabilities for constructs modeled as common factors are returned
  if(.only_common_factors){
    co_names <- names(c_types[c_types == "Composite"])
    rhoT     <- rhoT[setdiff(c_names, co_names)]
  }
  
  # Return named vector of reliabilities
  return(rhoT)
}



#' Internal: HTMT
#'
#' Compute the heterotrait-monotrait ratio of correlations 
#' \insertCite{Henseler2015}{cSEM}.
#'
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateHTMT(
#'  .object              = NULL,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateHTMT <- function(
  .object              = NULL,
  .only_common_factors = args_default()$.only_common_factors
){
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Get relevant quantities
  m <- .object$Information$Model
  
  if(isTRUE(.only_common_factors)) {
    cf_names <- names(m$construct_type[m$construct_type == "Common factor"])
    
    ## Stop if there are no common factors
    if(length(cf_names) < 2) {
      stop2("Computation of the HTMT requires at least two common factors, ",
            "unless `.only_common_factors = FALSE`.")
    }
  } else {
    cf_names <- names(m$construct_type)
  }
  
  cf_measurement <- m$measurement[cf_names, colSums(m$measurement[cf_names, ]) != 0, drop = FALSE]
  
  ## HTMT can only be calculated for constructs with more than one indicator
  x <- rowSums(cf_measurement) > 1
  cf_measurement <- cf_measurement[x, colSums(cf_measurement[x, ]) != 0, drop = FALSE]
  
  ## At least two multi-indicator constructs required
  if(length(which(x)) < 2) {
    stop2("Computation of the HTMT requires at least two multi indicator constructs.")
  }
  
  i_names <- colnames(cf_measurement)
  S       <- .object$Estimates$Indicator_VCV[i_names, i_names]
  
  ## Average correlation of the indicators of a block 
  avrg_cor <- cf_measurement %*% (S - diag(diag(S))) %*% t(cf_measurement) /
    cf_measurement %*% (1 - diag(nrow(S))) %*% t(cf_measurement)
  
  ## Compute HTMT
  out <- avrg_cor*lower.tri(avrg_cor) / sqrt(diag(avrg_cor) %o% diag(avrg_cor))
  
  # Return
  return(out)
}


#' Internal: Calculate difference between S and Sigma_hat
#'
#' Calculates the difference between the empirical 
#' and the model-implied indicator variance-covariance matrix 
#' using different distance measures. See vignette for assess.
#'
#' The functions are only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @inheritParams csem_arguments
#'
#' @keywords internal
#' @name distance_measures
NULL

#' @describeIn distance_measures The standardized root means squared residual (SRMR).

calculateSRMR <- function(.object = NULL) {
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object)

  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 

#' @describeIn distance_measures The geodesic distance (dG).

calculateDG <- function(.object = NULL) {

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object) 
  
  # Not sure if logarithm naturalis is used or logarithm with base 10. 
  Eigen            <- eigen(solve(S) %*% Sigma_hat)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   
  
  ## Calculate distance
  0.5 * sum(logEigenvaluessq)
}

#' @describeIn distance_measures The squared euclidian distance

calculateDL <- function(.object = NULL) {
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object) 
  
  ## Calculate distance
  0.5 * sum((S - Sigma_hat)[lower.tri(S, diag = FALSE)]^2)
}

#' @describeIn  distance_measures The distance measure used by FIML

calculateDML <- function(.object = NULL){

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  n         <- nrow(.object$Information$Data)
  S         <- .object$Estimates$Indicator_VCV
  p         <- dim(S)[1]
  Sigma_hat <- fit(.object, .saturated = FALSE, .type_vcv = 'indicator')
  
  (n - 1)*(log(det(Sigma_hat)) 
              + sum(diag(S %*% solve(Sigma_hat))) 
              - log(det(S)) - p)
}



#' Internal: RMS
#'
#' Compute the RMS.
#'
#' @usage RMS(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateRMS <- function(.object) {
  print("not yet implemented")
}



#' Internal: Calculate effect size
#'
#' Calculate the effect size
#'
#' @usage calculateEffectSize(.object = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateEffectSize <- function(.object = NULL) {
  
  # Get relevant quantities
  H <- .object$Estimates$Construct_scores
  Q <- sqrt(.object$Estimates$Construct_reliabilities)
  P <- .object$Estimates$Construct_VCV
  csem_model  <- .object$Information$Model
  normality   <- .object$Information$Arguments$.normality
  approach_nl <- .object$Information$Arguments$.approach_nl
  
  s <- csem_model$structural
  
  vars_endo <- rownames(s)[rowSums(s) != 0]
  outer_out <- lapply(vars_endo, function(x) {
    
    # get colnames
    indep_vars <- colnames(s[x , s[x, ] != 0, drop = FALSE])
    
    inner_out <- lapply(indep_vars, function(i) {
      # update csem_model
      model_temp <- csem_model
      model_temp$structural[x, i] <- 0 
      
      out <- estimatePathOLS(
        .H = H,
        .Q = Q,
        .P = P,
        .csem_model = model_temp,
        .normality = normality,
        .approach_nl = approach_nl
      )
      
      # calculate 
      r2_excluded <- out$R2[x]
      r2_included <- .object$Estimates$R2[x]
      
      effect_size <- (r2_included - r2_excluded)/(1 - r2_included)
      # list("r2_ex" = r2_excluded, "r2_in" = r2_included, "eff_size" = effect_size)
      list("effect_size" = effect_size)
    })
    names(inner_out) <- indep_vars
    inner_out
  })
  names(outer_out) <- vars_endo
  outer_out
}

