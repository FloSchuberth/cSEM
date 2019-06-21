#' Internal: AVE
#'
#' Calculate the average variance extracted (AVE) as proposed by 
#' \insertCite{Fornell1981;textual}{cSEM}. For details see the
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#ave}{cSEM website} 
#'
#' The AVE is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret AVE results in the context of a 
#' composite model. It is possible to report the AVE for concepts
#' modeled as composites by setting `.only_common_factors = FALSE`, 
#' however, we explicitly warn to interpret results with caution, 
#' as they may not even have a conceptual meaning.
#' 
#' The function is only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @usage calculateAVE(
#'  .object              = NULL,
#'  .only_common_factors = TRUE
#' )
#'
#' @return A named vector of numeric values of length equal to the number of constructs
#'   in the model.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

calculateAVE <- function(
  .object              = NULL,
  .only_common_factors = TRUE
  ){

  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Extract loadings
  Lambda  <- .object$Estimates$Loading_estimates
  c_names <- rownames(Lambda)
  
  ## Calculate AVE (extracted variance / total variance)
  # Note: Within the cSEM package indicators are always standardized (i.e, have
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
  
  # By default AVE's for constructs modeled as composites are not returned
  if(.only_common_factors){
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    AVEs      <- AVEs[names_cf]
  }
  
  # Return vector of AVEs
  return(AVEs) 
}


#' Internal: Degrees of freedom
#' 
#' Calculates the degrees of freedom from a [cSEMResults] object.
#' 
#' @usage calculateDF(.object = NULL)
#'
#' @return A single numeric value.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [cSEMResults]
#' @keywords internal

calculateDf <- function(.object) {
  
  x1 <- .object$Estimates
  x2 <- .object$Information
  
  ## Df for composite models ---------------------------------------------------
  if(all(x2$Model$construct_type == "Composite")) {
    if(inherits(.object, "cSEMResults_default")) {
      # Number of non-redundant off-diagonal elements of the indicator covariance 
      # matrix
      nS <- sum(lower.tri(x1$Indicator_VCV)) # same as dim_nS *(dim_nS - 1) / 2
      
      # Number of free correlations among the composites
      n_cor_composites <- sum(lower.tri(x1$Construct_VCV))
      
      # Number of free covariances between the composites and indicators not forming
      # a composite
      # Note 21.06.2019: I think this is not nessary yet, as indicators are always attached 
      #                  to a composite in cSEM
      
      # Number pf free non-redundant off-diagonal element of each intra-block
      # covariance matrix
      n_intra_block <- c()
      for(i in 1:nrow(x2$Model$measurement)) {
        dim <- sum((x2$Model$measurement[i, ] == 1))
        n_intra_block[i] <- dim * (dim - 1) / 2
      }
      
      # Number of weights
      n_weights <- ncol(x2$Model$measurement)
      
      # Number of blocks
      n_blocks  <- nrow(x2$Model$measurement)
      
      ## Calculate Dfs
      df <- nS - n_cor_composites - sum(n_intra_block) - n_weights + n_blocks
    } else {
      stop("not yet implemented")
    }
  } else {
    # DF for common factor model -----------------------------------------------
    
    stop2("Not yet implement for common factor models.")
  }
  
  df
}



#' Internal: GoF
#'
#' Calculate the Goodness of Fit (GoF) proposed by \insertCite{Tenenhaus2004;textual}{cSEM}. 
#' Note that, contrary to what the name suggests, the GoF is **not** a 
#' measure of model fit in a Chi-square fit test sense. See e.g. \insertCite{Henseler2012a;textual}{cSEM}
#' for a discussion.
#'
#' The GoF is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret the GoF in the context of a 
#' composite model. Hence, only constructs modeled as common factors are 
#' considered. It is possible to force computation of the GoF including constructs
#' modeled as composites as well by setting `.only_common_factors = FALSE`,
#' however, we explicitly discourage to do so as the result may not even 
#' have a conceptual meaning.
#' 
#' The function is only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @usage calculateGoF(
#'  .object              = NULL,
#'  .only_common_factors = TRUE
#' )
#'
#' @return A single numeric value.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

calculateGoF <- function(
  .object              = NULL,
  .only_common_factors = TRUE
){
  
  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Get relevant quantities
  Lambda    <- .object$Estimates$Loading_estimates
  R2        <- .object$Estimates$R2
  
  # The GoF is defined as the sqrt of the mean of the R^2s of the structural model 
  # times the variance in the indicators that is explained by the construct (lambda^2).
  # For the latter, only constructs modeled as common factors are considered
  # as they explain their indicators in contrast to a composite where 
  # indicators acutally build the construct.
  
  if(.only_common_factors) {
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    Lambda    <- Lambda[names_cf, ]
  }
  
  # Select only non-zero loadings
  L <- Lambda[Lambda != 0]
  

  gof <- sqrt(mean(L^2) * mean(R2))
  
  return(gof)
}



#' Internal: Reliability
#'
#' Compute several reliability estimates. See the 
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#reliability}{Reliability}
#' section of the \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website}
#' for details.
#' 
#' Since reliability is defined with respect to a classical true score measurement
#' model only constructs modeled as common factors are considered by default.
#' For constructs modeled as composites reliability may be estimated by setting
#' `.only_common_factors = FALSE`, however, it is unclear how to
#' interpret reliability in this case.
#' 
#' Reliability is are traditionally based on a test score (proxy) based on unit weights.
#' To compute congeneric and tau-equivalent reliability based on a score that 
#' uses the weights of the weight approach used to obtain `.object` use `.weighted = TRUE` 
#' instead.
#' 
#' For the the tau-equivalent reliability (rhoT or Cronbach alpha) a closed-form 
#' confidence interval may be computed \insertCite{Trinchera2018}{cSEM} by setting
#' `.closed_form_ci = TRUE` (default is `FALSE`). If `.alpha` is a vector
#' several CI's are returned.
#' 
#' The function is only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#'
#' @return For `calculateRhoC()` a named numeric vector containing the reliability
#'   estimates.
#'   For `calculateRhoT()` a `data.frame` with as many rows as there are
#'   constructs modeled as common factors in the model (unless 
#'   `.only_common_factors = FALSE` in which case the number of rows equals the
#'   total number of constructs in the model). The first column contains the name of the construct.
#'   The second column the reliability estimate.
#'   If `.closed_form_ci = TRUE` the remaining columns contain lower and upper bounds
#'   for the (1 - `.alpha`) confidence interval(s).
#'   
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @seealso [assess], [cSEMResults]
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
  .only_common_factors = TRUE,
  .weighted            = args_default()$.weighted
  ) {

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Get relevant objects
  if(.weighted) {
    W <- .object$Estimates$Weight_estimates
  } else {
    W <- .object$Information$Model$measurement
    W <- scaleWeights(.object$Estimates$Indicator_VCV, W)
  }
  
  rhoC <- rep(1, times = nrow(W))
  names(rhoC) <- rownames(W)
  
  # Get loadings
  Lambda  <- .object$Estimates$Loading_estimates
  
  # Compute the congeneric reliability by block 
  for(j in rownames(Lambda)) {
    rhoC[j]  <- c(W[j, ] %*% Lambda[j, ])^2 
  }
  
  # By default only reliabilities for constructs model common factors are returned
  if(.only_common_factors){
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    rhoC      <- rhoC[names_cf]
  }

  # Return named vector of reliabilities
  return(rhoC)
}

#' @describeIn reliability Calculate the tau-equivalent reliability, also known
#'                         as Cronbach's alpha or coefficient alpha. Since
#'                         indicators are always standarized in `cSEM`, 
#'                         tau-equivalent reliability is identical to the
#'                         parallel reliability.
#' 
calculateRhoT <- function(
  .object              = NULL,
  .alpha               = args_default()$.alpha,
  .closed_form_ci      = args_default()$.closed_form_ci,
  .only_common_factors = TRUE,
  .weighted            = args_default()$.weighted,
  ...
  ) {
  
  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Get relevant objects
  S <- .object$Estimates$Indicator_VCV
  
  if(.weighted) {
    W <- .object$Estimates$Weight_estimates
  } else {
    W <- .object$Information$Model$measurement
    W <- scaleWeights(S, W)
  }
  
  ## Calculate tau-equivalent reliability by block/construct
  out <- data.frame()
  for(j in rownames(W)) {
    indicator_names <- colnames(W[j, W[j,] != 0, drop = FALSE])
    S_jj            <- S[indicator_names, indicator_names]
    rho_bar_j <- mean(S_jj[upper.tri(S_jj)])
    rhoT      <- rho_bar_j * sum(W[j, ])^2  
    
    ## Add to data.frame
    out[j, "Construct"] <- j
    out[j, "Estimate"]  <- rhoT
    
    ## Calculate confidence interval for rhoT
    if(.closed_form_ci) {
      # Calculation of the CIs are based on Trinchera et al. (2018).
      # The code for the calculation of the CIs is addpated from the paper.

      ## If CI's are computed:
      # Calculation of the CIs are based on Trinchera et al. (2018).
      # In the paper, the CI was proposed and studied using Monte Carlo simulation
      # assuming scores are build as sum scores. Therefore a warning is
      # given if a weightning scheme other than "unit" is used, since they have
      # not been formally studied yet.
      
      if(.object$Information$Arguments$.approach_weights != "unit" & .weighted) {
        warning2("Calculation of the confidence interval (CI) for rhoT (Cronbach alpha)",
                 " was proposed and studied assuming sum scores.\n",
                 "Statistical properties of the CI for rhoT based on a weighted composite",
                 " obtained using weight approach ",
                 paste0("`", .object$Information$Arguments$.approach_weights, "`"),
                 " have not been formally examined yet.")
      }
      
      
      K <- length(indicator_names)
      X <- .object$Information$Data[, indicator_names, drop=FALSE]
      N <- nrow(X)
      H <- .object$Estimates$Construct_scores[, j]
      
      # In cSEM X and H (the scores) are always standardized. As a consequence, 
      # many of the quantities in Trichera et al. (2018)'s proposed 
      # formula for the standard error become zero or one and therefore drop.
      # Quantities that are zero (using their notation):
      # - bar(X)_p; bar(X)_q;  bar(S)_X
      # Quantities that are one:
      # - var(H) and var(X) and therefore sigma^2_Sx, sigma^4_Sx, sigma^6, 
      #   sigma^8 and sigma_Xp
      
      ## Standard error formula for cSEM (for original formula see the paper)
      #
      #   theta_hat = (K / (K - 1))^2 * *(A + 2B + C)
      #
      #   where 
      #   A := sum_k sum_l (-1 + 1/(N-1) sum_i (X^2_ik * X^2_il))
      #   B := K * sum_k (-1 + 1/(N-1) sum_i (H^2_i * X^2_ik))
      #   C := K^2 * (1 + 1/(N-1) sum_i H^4_i)
      # Note also: in the code provided in the paper N was used instead of N-1
      
      A <- sum(t(X^2) %*% X^2/(N-1) - 1)
      B <- K * sum(H^2 %*% X^2/(N-1) - 1)
      C <- K^2 * (1/(N-1)*sum(H^4) - 1)

      ## Calculate the variance and the se
      var_rhoT <- (K^2/(K - 1)^2) * (A - 2*B + C)
      se_rhoT  <- sqrt(var_rhoT/N)
      
      ## Order the alphas provided, compute CI and add to data frame
      .alpha <- .alpha[order(.alpha)]
      for(i in .alpha) { 
        z_value <- qnorm((1 - i/2), mean = 0, sd = 1)
        up  <- rhoT + z_value * se_rhoT
        low <- rhoT - z_value * se_rhoT
        name_L <- sprintf("%.6g%%L", 100* (1 - i))
        name_U <- sprintf("%.6g%%U", 100* (1 - i))
        
        out[j, name_L] <- low
        out[j, name_U] <- up
      }
    }
  }
  
  # By default only reliabilities for constructs modeled as common factors are returned
  if(.only_common_factors){
    con_types <- .object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    out       <- out[names_cf, ]
  }
  
  # Return named vector of reliabilities
  return(out)
}



#' Internal: HTMT
#'
#' Compute the heterotrait-monotrait ratio of correlations (HTMT) based on 
#' \insertCite{Henseler2015;textual}{cSEM}.
#'
#' The HTMT is used to assess discriminant validity.
#' 
#' The function is only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @usage calculateHTMT(
#'  .object              = NULL,
#'  .only_common_factors = TRUE
#' )
#'
#' @inheritParams csem_arguments
#'
#' @return A lower tringular matrix of HTMT values.
#' 
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateHTMT <- function(
  .object              = NULL,
  .only_common_factors = TRUE
){
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
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
  if(nrow(cf_measurement) < 2) {
    stop2("Computation of the HTMT requires at least two multi-indicator common factors.")
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
#' Calculate the difference between the empirical (S) 
#' and the model-implied indicator variance-covariance matrix (Sigma_hat)
#' using different distance measures. See vignette assess.
#'
#' The functions are only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @return A single numeric value.
#' 
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @keywords internal
#' @name distance_measures
NULL

#' @describeIn distance_measures The geodesic distance (dG).

calculateDG <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv,
  ...
  ){

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv) 
  
  # Not sure if logarithm naturalis is used or logarithm with base 10. 
  Eigen            <- eigen(solve(S) %*% Sigma_hat)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   
  
  ## Calculate distance
  0.5 * sum(logEigenvaluessq)
}

#' @describeIn distance_measures The squared Euclidian distance

calculateDL <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv,
  ...
  ){
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)
  
  ## Calculate distance
  0.5 * sum((S - Sigma_hat)[lower.tri(S, diag = FALSE)]^2)
}

#' @describeIn  distance_measures The distance measure used by FIML

calculateDML <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv,
  ...
  ){

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  n         <- nrow(.object$Information$Data)
  S         <- .object$Estimates$Indicator_VCV
  p         <- dim(S)[1]
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)
  
  # This is the distance function. The test statistic is T_ML = (n-1) or n * DML! 
  sum(diag(S %*% solve(Sigma_hat))) - log(det(S%*%solve(Sigma_hat))) - p
  
}

#' Internal: Fit measures
#' 
#' Calculate common fit measures
#' 
#' The functions are only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @return A single numeric value.
#' 
#' @inheritParams csem_arguments
#'
#' @keywords internal
#' @name fit_measures 

#' @describeIn fit_measures The goodness of fit index (GFI).

calculateGFI <- function(.object) {
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object)
  
  1 - matrixcalc::matrix.trace(t(S - Sigma_hat) %*% (S - Sigma_hat)) / 
    matrixcalc::matrix.trace(t(S) %*% S)
}

#' @describeIn fit_measures The normed fit index (NFI).

calculateNFI <- function(.object) {
  
  n <- nrow(.object$Information$Data)
  S <- .object$Estimates$Indicator_VCV
  p <- dim(S)[1]
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  (F0 - FT) / F0
}

#' @describeIn fit_measures The root mean square error of approximation (RMSEA).

calculateRMSEA <- function(.object) {
  
  n  <- nrow(.object$Information$Data)
  df <- calculateDf(.object)
  
  F0 <- max(calculateDML(.object) - calculateDf(.object)/(n - 1), 0)
  
  sqrt(F0 / df) # RMSEA
}

#' @describeIn fit_measures The standardized root mean square residual (SRMR).

calculateSRMR <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv,
  ...
) {
  
  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)
  
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 

#' @describeIn fit_measures The RMS theta (SRMR).

calculateRMSTheta <- function(.object) {
  S      <- .object$Estimates$Indicator_VCV
  W      <- .object$Estimates$Weight_estimates
  Lambda <- .object$Estimates$Loading_estimates
  P      <- .object$Estimates$Construct_VCV
  
  Theta <- S - S %*% t(W) %*% Lambda - t(S %*% t(W) %*% Lambda) + t(Lambda) %*% P %*% Lambda
  
  ## For compsites, within block indicator correlations should be excluded as 
  ## they are allowed to freely covary.
  
  comp <- which(.object$Information$Model$construct_type == "Composite")
  
  for(i in comp) {
    indi <- which(.object$Information$Model$measurement[i, ] == 1)
    Theta[indi, indi] <- NA
  }
  
  sqrt(mean(Theta[lower.tri(Theta1)]^2, na.rm = TRUE))
}

#' Internal: Calculate effect size
#'
#' Calculate the effect size.
#'
#' @usage calculateEffectSize(.object = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @return A matrix with as many rows as there are structural equations. The 
#'   number of columns is equal to the total number of right-hand side variables
#'   of these equations.
#' 
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateEffectSize <- function(.object = NULL) {
  
  ## Get relevant quantities
  approach_nl    <- .object$Information$Arguments$.approach_nl
  approach_paths <- .object$Information$Arguments$.approach_paths
  csem_model     <- .object$Information$Model
  H <- .object$Estimates$Construct_scores
  normality <- .object$Information$Arguments$.normality
  P <- .object$Estimates$Construct_VCV
  Q <- sqrt(.object$Estimates$Reliabilities)
  
  s <- csem_model$structural
  
  vars_endo <- rownames(s)[rowSums(s) != 0]
  ## Start loop
  outer_out <- lapply(vars_endo, function(x) {
    
    # Get names of the independent variables
    indep_vars <- colnames(s[x , s[x, ] != 0, drop = FALSE])
    
    inner_out <- lapply(indep_vars, function(i) {
      # Update csem_model
      model_temp <- csem_model
      model_temp$structural[x, i] <- 0 
      
      out <- estimatePath(
        .approach_nl    = approach_nl,
        .approach_paths = approach_paths,
        .csem_model     = model_temp,
        .H              = H,
        .normality      = normality,
        .P              = P,
        .Q              = Q
      )
      
      ## Calculate effect size
      # For equations with only one independet variable the R2_excluded is 0
      r2_excluded <- ifelse(is.na(out$R2[x]), 0, out$R2[x])
      names(r2_excluded) <- x
      r2_included <- .object$Estimates$R2[x]
      
      effect_size <- unname((r2_included - r2_excluded)/(1 - r2_included))

      # list("r2_ex" = r2_excluded, "r2_in" = r2_included, "eff_size" = effect_size)
    })
    inner_out <- unlist(inner_out)
    names(inner_out) <- indep_vars
    inner_out
  })
  names(outer_out) <- vars_endo
  
  ## Make output a matrix
  # Note: this is necessary for calculateEffectSize to work
  #       when supplied to the .user_funs argument. Currently, .user_funs functions 
  #       need to return a vector or a matrix. I may change that in the future.
  ss <- s[vars_endo, , drop = FALSE]
  tm <- t(ss)
  tm[which(tm == 1)] <- unlist(outer_out)
  # Return
  t(tm)
}



#' Internal: Calculate VIF values for PLS Mode B
#'
#' Calculate the variance inflation factor (VIF) for indicators of constructs whose
#' outer PLS estimation scheme is Mode B.
#'
#' Weights obtained using Mode B often suffer from multicollinearity. VIF values
#' are commonly used to assess the severity of multicollinearity. 
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateVIFModeB(.object = NULL)
#'
#' @return A named list of vectors containing the VIF values. Each list name
#'   is the name of a construct whose weights were obtained using Mode B. 
#'   The vectors contain the VIF values obtained from a regression of each 
#'   explanatory variable of a given construct on the remaining explanatory 
#'   variables of that construct.
#'   
#'   If the weightning approach was not `"PLS-PM"` or non of the constructs had Mode B,
#'   the function silently returns `NA`.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

calculateVIFModeB <- function(.object = NULL) {
  
  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## Get the modes
  modes  <- .object$Information$Weight_info$Modes
  modesB <- modes[modes == "modeB"] 
  
  # Only compute if 1. PLS-PM is the weight approach, 2. at least
  # one construct has outer weightning scheme Mode B 
  if (.object$Information$Arguments$.approach_weights == "PLS-PM" &&
      length(modesB) > 0) {
    
    m <- .object$Information$Model$measurement
    X <- .object$Information$Data
    
    VIF <- lapply(names(modesB), function(j) {
      indicator_names <- colnames(m[j, m[j,] != 0, drop = FALSE])
      
      # If j is a single indicator construct set VIF values to NA otherwise
      # continue with the computation.
      if(length(indicator_names) > 1) {

        VIF_j <- c()
        for(k in indicator_names) {
          
          y  <- X[, k]
          Xk <- X[, setdiff(indicator_names, k)]
          
          beta <- solve(t(Xk) %*% Xk) %*% t(Xk) %*% y
          R2   <- cor(Xk %*% beta, y)
          VIF_j[k] <- 1 / (1 - R2) 
        }
        
        VIF_j
        
      } else {
        VIF_j <- NA
      }
    })
    names(VIF) <- names(modesB)
  } else {
    VIF <- NA
  }
  
  return(VIF)
}



#' Internal: Calculate a redundancy analysis
#'
#' Calculate/do a redundancy analysis (RA) as proposed by Hair (2014) with
#' reference to Chin (1998).
#'
#' According to Hair (2014), redundancy analysis (RA)
#' is the process of regressing the scores of a reflectivly measured construct
#' on the scores of a formatively measured construct in order to gain empirical
#' evidence for convergent validity of a formatively measured construct. 
#' RA is therefore confined to PLS, specifically PLS with at least one construct
#' whose mode is Mode B. This is the case if the construct is modeled as a 
#' composite or if the construct was explicitly given Mode B.
#' Hence RA is only conducated if `.approach_weights = "PLS-PM"` and if at least
#' one construct's mode is Mode B.
#'  
#' The principal idea of RA is to take two different measures of the 
#' same construct and regress the scores obtained for each measure on each
#' other. If they are similar they are likely to measure the same "thing"
#' which is then taken as evidence that both measurement models actually
#' measure what they are supposed to measure (validity). 
#' 
#' There are several issues with the termionlogy and the reasoning behind this logic.
#' RA is therefore only implemented since reviewers are likely to demand
#' its compuation, however, its actual application for validitiy assessment 
#' is discouraged.
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' 
#' @usage calculateRA(.object = NULL)
#'
#' @return A named numeric vector of correlations obtained by the RA. If 
#'   the weightning approach is not `"PLS-PM"` or non of the PLS outer modes
#'   was mode B, the function silently returns `NA`.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

calculateRA <- function(.object = NULL) {
  
  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  if(.object$Information$Arguments$.approach_weights == "PLS-PM") {
    
    modes <- .object$Information$Weight_info$Modes
    modesB <- modes[modes == "modeB"]
    
    if(length(modesB) > 0) {
      
      args <- .object$Information$Arguments
      
      # Functions resampleData() and testMICOm require the .id argument.
      # It is therefore present in the Arguments list although the 
      # data set in the Arguments list does not contain the id column anymore.
      # Therefore .id needs to be set to NULL
      args[[".id"]] <- NULL
      new_modes <- as.list(modes) 
      
      beta <- c()
      for(j in names(modesB)) {
        new_modes_j <- new_modes
        new_modes_j[j] <- "modeA"
        
        args[[".PLS_modes"]] <- new_modes_j
        
        res_reflective <- do.call(csem, args)
        
        Y <- res_reflective$Estimates$Construct_scores
        X <- .object$Estimates$Construct_scores
        
        beta[j] <- c(solve(t(X[,j]) %*% X[,j]) %*% t(X[,j]) %*% Y[,j])
      }
    } else {
      beta <- NA
    }
  } else {
    beta <- NA
  }
  
  return(beta)
}