#' Internal: AVE
#'
#' Calculate the average variance extracted (AVE) as proposed by 
#' \insertCite{Fornell1981;textual}{cSEM}. For details see the
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#ave}{cSEM website} 
#'
#' The AVE is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret AVE results in the context of a 
#' composite model. It is possible to report the AVE for composites by
#'  setting `.only_common_factors = FALSE`, 
#' however, result should be interpreted with caution 
#' as they may not have a conceptual meaning.
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
  
  # By default AVE's for composites are not returned
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
#' Calculates the degrees of freedom for a given model from a [cSEMResults] object.
#' 
#' @usage calculateDf(
#'   .object     = NULL,
#'   .null_model = args_default()$.null_model,
#'   ...
#'   )
#' 
#' @return A single numeric value.
#'   
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @seealso [assess], [cSEMResults]
#' @keywords internal

calculateDf <- function(
  .object     = NULL, 
  .null_model = args_default()$.null_model,
  ...
  ) {
  
  if(inherits(.object, "cSEMResults_default")) { 
    x1 <- .object$Estimates
    x2 <- .object$Information$Model
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Degrees of freedom for 2nd order composites not implemented yet.")
  } else {
    
    df <- lapply(.object, calculateDf)
    
    ## Return
    names(df) <- names(.object)
    return(df)
  }
  
  ## Number of non-redundant off-diagonal elements of the indicator covariance 
  ## matrix (S)
  vS <- sum(lower.tri(x1$Indicator_VCV)) # same as dim_nS *(dim_nS - 1) / 2
  
  ## Caculate the degrees of freedom of a null model (Sigma_hat = I)
  if(.null_model) {
    # The degrees of freedom of the null model are identical to 
    # the number of non-redundant elements of S (since there is nothing to
    # estimate. Everything is set to 0 a priori)
    return(vS)
  }
  
  ## Number of correlations between exogenous constructs
  temp <- sum(rowSums(x2$structural) == 0)
  n_exo <- temp * (temp - 1) / 2
  
  ## Number of structural parameters
  n_structural <- sum(x2$structural)
  
  if(.object$Information$Arguments$.approach_weights %in% c("bartlett", "regression", "unit")) {
    # If one of these estimators is used degrees of freedom are counted 
    # the same way as one would count in CB-SEM: 
    # df = non_redundant_elements_of_S - 
    #      - (structural_parameter + cor_between_exos) 
    #      - loadings
    #      - assumed_measurement_error_cors (usually 0)
    
    ## Number of measurement errors assumed to correlate
    n_error <- sum(x2$error_cor == 1) / 2
    
    n_loadings <- ncol(x2$measurement)
    
    df_total <- vS - (n_exo + n_structural) - n_loadings - n_error
  } else {
    ## Construct names
    names_constructs <- rownames(x2$structural)
    k <- c()
    for(j in names_constructs) {
      # DF for 
      ## Number of free covariances between the composites and indicators not forming
      ## a composite
      # Not relevant, as indicators are always attached to a composite in cSEM
      
      ## Number of free non-redundant off-diagonal element of each intra-block
      #+ covariance matrix
      temp    <- sum((x2$measurement[j, ] == 1))
      n_intra <- temp * (temp - 1) / 2
      
      ## Number of weights minus 1 (since weights are choosen s.t. Var(eta) = 1)
      n_weights <- sum(x2$measurement[j, ]) - 1
      
      ## Calculate Dfs
      k[j] <- n_intra + n_weights
    }
    
    df_total <- vS - (n_exo + n_structural) - sum(k)
  }
  
  # return degrees of freedom
  df_total
}



#' Internal: GoF
#'
#' Calculate the Goodness of Fit (GoF) proposed by \insertCite{Tenenhaus2004;textual}{cSEM}. 
#' Note that, contrary to what the name suggests, the GoF is **not** a 
#' measure of model fit in the sense of SEM. See e.g. \insertCite{Henseler2012a;textual}{cSEM}
#' for a discussion.
#'

# The GoF is inherently tied to the common factor model. It is therefore 
# unclear how to meaningfully interpret the GoF in the context of a 
# composite model. Hence, only constructs modeled as common factors are 
# considered. It is possible to force computation of the GoF including constructs
# modeled as composites as well by setting `.only_common_factors = FALSE`,
# however, we explicitly discourage to do so as the result may not even 
# have a conceptual meaning.
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
#' model only concepts modeled as common factors are considered by default.
#' For concepts modeled as composites reliability may be estimated by setting
#' `.only_common_factors = FALSE`, however, it is unclear how to
#' interpret reliability in this case.
#' 
#' Reliability is are traditionally based on a test score (proxy) based on unit weights.
#' To compute congeneric and tau-equivalent reliability based on a score that 
#' uses the weights of the weight approach used to obtain `.object` use `.weighted = TRUE` 
#' instead.
#' 
#' For the the tau-equivalent reliability (rhoT or Cronbach's alpha) a closed-form 
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
      # given if a weighting scheme other than "unit" is used, since they have
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
#' Similarly to the Cronbach's alpha, the HTMT is a consistent estimator for the
#' construct correlations in case of tau equivalent measurement models. 
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
    
    ## Return NA if there are not at least 2 common factors
    if(length(cf_names) < 2) {
      warning2("Computation of the HTMT requires at least two common factors, ",
               "unless `.only_common_factors = FALSE`. NA is returned.")
      return(NA)
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
#' using different distance measures.
#'
#' The functions are only applicable to objects inheriting class `cSEMResults_default`.
#' or `cSEMResults_2ndorder`. For objects of class `cSEMResults_multi` 
#' use `lapply(.object, calculateXX())`.
#' 
#' The geodesic and the squared Euclidian distance may also be 
#' computed for any two matrices A and B by supplying A and B directly via the 
#' `.matrix1` and `.matrix2` arguments. If A and B are supplied `.object` is ignored.
#' 
#' @return A single numeric value giving the distance between two matrices.
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
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = args_default()$.saturated,
  ...
  ){

  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    # Only applicable to objects of class cSEMResults_default and cSEMResults_2ndorder
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateDG() function:\n",
        "`.object` must be of class `cSEMResults_default` or `cSEMResults_2ndorder`.")
    }
    
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')  
  }
  
  # Not sure if logarithm naturalis is used or logarithm with base 10. 
  Eigen            <- eigen(solve(S) %*% Sigma_hat)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   
  
  ## Calculate distance
  0.5 * sum(logEigenvaluessq)
}

#' @describeIn distance_measures The squared Euclidian distance

calculateDL <- function(
  .object    = NULL, 
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = args_default()$.saturated,
  ...
  ){
  
  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    # Only applicable to objects of class cSEMResults_default and cSEMResults_2ndorder
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateDL() function:\n",
        "`.object` must be of class `cSEMResults_default` or `cSEMResults_2ndorder`.")
    }
    
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')
  }
  
  ## Calculate distance
  0.5 * sum((S - Sigma_hat)[lower.tri(S, diag = FALSE)]^2)
}

#' @describeIn  distance_measures The distance measure used by FIML

calculateDML <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  ...
  ){
  
  # Only applicable to objects of class cSEMResults_default and cSEMResults_2ndorder
  if(inherits(.object, "cSEMResults_default")) {
    S <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    S <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateDML() function:\n",
      "`.object` must be of class `cSEMResults_default` or `cSEMResults_2ndorder`.")
  }

  p         <- dim(S)[1]
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')
  
  # This is the distance function. The test statistic is T_ML = (n-1) or n * DML! 
  sum(diag(S %*% solve(Sigma_hat))) - log(det(S%*%solve(Sigma_hat))) - p
}

#' Internal: Fit measures
#' 
#' Calculate common fit measures.
#' 
#' All functions, except for `calculateSRMR()` are only applicable to 
#' objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#' 
#' @return A single numeric value.
#' 
#' @inheritParams csem_arguments
#'
#' @keywords internal
#' @name fit_measures 
NULL

#' @describeIn fit_measures The comparative fit index (CFI).

calculateCFI <- function(.object) {
  
  n    <- nrow(.object$Information$Data)
  S    <- .object$Estimates$Indicator_VCV
  p    <- dim(S)[1]
  df_T <- calculateDf(.object)
  df_0 <- calculateDf(.object, .null_model = TRUE)
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- max((n-1)*calculateDML(.object) - calculateDf(.object), 0)
  F0 <- max((n-1)*F0 - df_0 , (n-1)*calculateDML(.object) - calculateDf(.object), 0)
  
  1 - FT/F0
}

#' @describeIn fit_measures The goodness of fit index (GFI).

calculateGFI <- function(.object) {
  
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object)
  
  1 - matrixcalc::matrix.trace(t(S - Sigma_hat) %*% (S - Sigma_hat)) / 
    matrixcalc::matrix.trace(t(S) %*% S)
}

#' @describeIn fit_measures The incremental fit index (IFI).

calculateIFI <- function(.object) {
  
  n  <- nrow(.object$Information$Data)
  S  <- .object$Estimates$Indicator_VCV
  p  <- dim(S)[1]
  df <- calculateDf(.object)
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  ((n-1)*F0 - (n-1)*FT) / ((n-1)*F0 - df)
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

#' @describeIn fit_measures The non-normed fit index (NNFI).

calculateNNFI <- function(.object) {
  
  n    <- nrow(.object$Information$Data)
  S    <- .object$Estimates$Indicator_VCV
  p    <- dim(S)[1]
  df_T <- calculateDf(.object)
  df_0 <- calculateDf(.object, .null_model = TRUE)
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  (F0/df_0 - FT/df_T) / (F0/df_0 - 1/(n-1))
}

#' @describeIn fit_measures The root mean square error of approximation (RMSEA).

calculateRMSEA <- function(.object) {
  
  n  <- nrow(.object$Information$Data)
  df <- calculateDf(.object)
  
  F0 <- max(calculateDML(.object) - calculateDf(.object)/(n - 1), 0)
  
  sqrt(F0 / df) # RMSEA
}

#' @describeIn fit_measures The root mean squared residual covariance matrix of the outer model residuals (RMS theta).

calculateRMSTheta <- function(
  .object, 
  .model_implied = args_default()$.model_implied
  ) {
  S      <- .object$Estimates$Indicator_VCV
  W      <- .object$Estimates$Weight_estimates
  Lambda <- .object$Estimates$Loading_estimates
  P      <- .object$Estimates$Construct_VCV
  
  
  if(.model_implied) {
    Theta <- S - S %*% t(W) %*% Lambda - t(S %*% t(W) %*% Lambda) + t(Lambda) %*%  fit(.object, .type_vcv = "construct") %*% Lambda
  } else {
    Theta <- S - S %*% t(W) %*% Lambda - t(S %*% t(W) %*% Lambda) + t(Lambda) %*% P %*% Lambda
  }
  
  # Check how its done in smartpls
  
  ## For compsites, within block indicator correlations should be excluded as 
  ## they are allowed to freely covary.
  
  comp <- which(.object$Information$Model$construct_type == "Composite")
  
  for(i in comp) {
    indi <- which(.object$Information$Model$measurement[i, ] == 1)
    Theta[indi, indi] <- NA
  }
  
  sqrt(mean(Theta[lower.tri(Theta)]^2, na.rm = TRUE))
}

#' @describeIn fit_measures The standardized root mean square residual (SRMR).

calculateSRMR <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  ...
) {
  
  # Only applicable to objects of class cSEMResults_default and cSEMResults_2ndorder
  if(inherits(.object, "cSEMResults_default")) {
    S <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    S <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateSRMR() function:\n",
      "`.object` must be of class `cSEMResults_default` or `cSEMResults_2ndorder`.")
  }
  
  # The SRMR as calculated by us is always based on the the difference 
  # between correlation matrices.
  Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')
  
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 




#' Internal: Calculate effect size
#'
#' Calculate the effect size for regression analysis \insertCite{Cohen1992}{cSEM}.
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
  approach_nl      <- .object$Information$Arguments$.approach_nl
  approach_paths   <- .object$Information$Arguments$.approach_paths
  approach_weights <- .object$Information$Arguments$.approach_weights
  csem_model       <- .object$Information$Model
  H         <- .object$Estimates$Construct_scores
  normality <- .object$Information$Arguments$.normality
  P         <- .object$Estimates$Construct_VCV
  Q         <- sqrt(.object$Estimates$Reliabilities)
  
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
        .approach_nl      = approach_nl,
        .approach_paths   = approach_paths,
        .approach_se      = "none",
        .approach_weights = approach_weights,
        .csem_model       = model_temp,
        .H                = H,
        .normality        = normality,
        .P                = P,
        .Q                = Q
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



#' Internal: Calculate variance inflation factors (VIF) for weights obtained by PLS Mode B
#'
#' Calculate the variance inflation factor (VIF) for weights obtained by PLS-PM's Mode B.
#'
#' Weight estimates obtained by Mode B can suffer from multicollinearity. VIF values
#' are commonly used to assess the severity of multicollinearity. 
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateVIFModeB(.object = NULL)
#'
#' @return A named list of vectors containing the VIF values. Each list name
#'   is the name of a construct whose weights were obtained by Mode B. 
#'   The vectors contain the VIF values obtained from a regression of each 
#'   explanatory variable of a given construct on the remaining explanatory 
#'   variables of that construct.
#'   
#'   If the weighting approach is not `"PLS-PM"` or for none of the constructs Mode B is used,
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
  # one construct has outer weighting scheme Mode B 
  if (.object$Information$Arguments$.approach_weights == "PLS-PM" &&
      length(modesB) > 0) {
    
    m <- .object$Information$Model$measurement
    X <- .object$Information$Data_with_RI
    
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



#' Internal: Do a redundancy analysis
#'
#' Calculate/do a redundancy analysis (RA) as proposed by \insertCite{Hair2016;textual}{cSEM}
#' with reference to \insertCite{Chin1998;textual}{cSEM}.
#'
#' RA is therefore confined to PLS-PM, specifically PLS-PM with at least one construct
#' whose weight are obtained by Mode B. This is the case if the construct is modeled as a 
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
#' @usage doRedundancyAnalysis(.object = NULL)
#'
#' @return A named numeric vector of correlations. If 
#'   the weighting approach is not `"PLS-PM"` or non of the PLS outer modes
#'   was mode B, the function silently returns `NA`.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @keywords internal

doRedundancyAnalysis <- function(.object = NULL) {
  
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

#' Do a floodlight analysis
#'
#' Calculate the the effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a (continous) moderator variable (x) 
#' to perform a floodlight analysis \insertCite{Spiller2013}{cSEM}. Moreover, 
#' the Johnson-Neyman points are calculated, i.e. the value(s) of x for which 
#' lower or upper boundary of the confidence interval of the effect
#' estimate of z for a given x switch signs. 
#' 
#' @usage doFloodlightAnalysis(
#'  .object        = NULL,
#'  .alpha         = args_default()$.alpha,
#'  .y             = args_default()$.y, 
#'  .x             = args_default()$.x,
#'  .z             = args_default()$.z,
#'  .n_spotlights  = args_default()$.n_spotlights
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
doFloodlightAnalysis <- function(
  .object        = NULL,
  .alpha         = args_default()$.alpha,
  .y             = args_default()$.y, 
  .x             = args_default()$.x,
  .z             = args_default()$.z,
  .n_spotlights  = args_default()$.n_spotlights
){
  
  ### Errors and warnings ------------------------------------------------------
  ## Check whether .object is of class cSEMResults_resampled; if not perform
  ## standard resampling (bootstrap with .R = 499 reps)
  if(!inherits(.object, "cSEMResults_resampled")) {
    if(inherits(.object, "cSEMResults_default")) {
      args <- .object$Information$Arguments
    } else {
      args <- .object$Second_stage$Information$Arguments_original
    }
    args[".resample_method"] <- "bootstrap"
    .object <- do.call(csem, args)
  }
  
  ##  Select relevant quantities
  if(inherits(.object, "cSEMResults_default")) {
    m   <- .object$Information$Model
    est <- .object$Estimates$Estimates_resample$Estimates1$Path_estimates
    H   <- .object$Estimates$Construct_scores
  } else {
    m   <- .object$Second_stage$Information$Model
    est <- .object$Second_stage$Information$Resamples$Estimates$Estimates1$Path_estimates
    H   <- .object$Second_stage$Estimates$Construct_scores
  }
  
  # Character string containing the names of the dependent and independent variables
  dep_vars   <- rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE])
  indep_vars <- colnames(m$structural[, colSums(m$structural[.y, ,drop = FALSE]) !=0 , drop = FALSE])
  
  ## Check if model is nonlinear.
  if(m$model_type != "Nonlinear"){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "The structural model must be nonlinear (i.e. contain interactions).")
  }
  
  ## Works only for one significance level, i.e., no vector of significances is allowed
  if(length(.alpha) != 1){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "Currently only a single significance level (.alpha) is allowed.")
  }
  
  ## Check whether dependent, independent, and moderator variable are provided
  if(is.null(.y) | is.null(.x) | is.null(.z)){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "All variables (.y, .x, and .z) must be supplied.")
  }
  
  # Check if the name of the dependent variable is valid
  if(!(.y %in% dep_vars)){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "The dependent variable supplied to `.y` is not a dependent variable in the original model.")
  }
  
  # Check if the name of the moderator (.x) and the dependent variable (.z) supplied are used 
  # in the original model
  if(!all(c(.x, .z) %in% colnames(m$structural))){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "Independent and/or moderator variable are not part of the original model.")
  }
  
  ### Calculation --------------------------------------------------------------
  # Possible names of the interaction terms
  possible_names <- c(paste(.x, .z, sep = '.'), paste(.z, .x, sep= '.'))
  
  # Name of the interaction term
  xz <- possible_names[possible_names %in% indep_vars]
  
  if(length(xz) != 1){
    stop2(
      "The defined interaction term does not exist in the model or is not ",
      " part of the equation of the dependent variable.")
  }
  
  
  # Effect names
  beta_z  <- paste(.y, .x, sep = ' ~ ')
  beta_xz <- paste(.y, .z, sep = ' ~ ')
  
  ## Compute spotlights (= effects of independent (z) on dependent (y) for given 
  ## values of moderator (x))
  steps <- seq(min(H[, .x]), max(H[, .x]), length.out = .n_spotlights)
  
  dataplot_temp <- lapply(steps, function(step){
    ## Note:
    #   
    #   y = a + bz + cx + dxz
    #   The marginal effect delta y = b + dx evaluated at x0 is identical to the
    #   b' of the regression that uses x' = x - x0, since
    #   y = a' + b'z + c'x' + d'x'z
    #     = (a - cx0) + (b - dx0)z + cx + dxz
    #   ---> b' = b - dx0
    #
    # Resamples of the effect of z on y at different levels of x 
    effect_resampled <- est$Resampled[ , beta_z] + est$Resampled[, beta_xz] * step 
    
    # Value of the originally estimated effect of z on y at different levels of x
    effect_original <- est$Original[beta_z] + est$Original[beta_xz] * step
    
    # Compute empirical quantile based on resamples
    bounds <- quantile(effect_resampled , c(.alpha/2, 1 - .alpha/2))
    
    # Return output
    c(effect_original, step, bounds[1],  bounds[2])
  })
  
  out <- do.call(rbind, dataplot_temp)
  colnames(out) <- c('direct_effect', 'value_z', 'lb', 'ub')
  
  # Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  pos_ub <- which(diff(sign(out[, 'ub'])) != 0)
  pos_lb <- which(diff(sign(out[, 'lb'])) != 0) 
  
  # Prepare and return output
  out <- list(
    "out"                   = out, 
    "Johnson_Neyman_points" = list(
      JNlb = c(x = out[,'value_z'][pos_lb], 
               y = out[,'lb'][pos_lb]),
      JNub = c(x = out[,'value_z'][pos_ub], 
               y = out[,'ub'][pos_ub]
      )
    ),
    "Information" = list(
      alpha       = .alpha, 
      dependent   = .y,
      independent = .z,
      moderator   = .x
    )
  )
  
  class(out) = "cSEMFloodlight"
  return(out)
} 

#' `cSEMFloodlight` method for `plot()`
#'
#' Plot the direct effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a moderator variable (x), including
#' the confidence interval and the Johnson-Neyman points. 
#' 
#' @param x An R object of class `cSEMFloodlight` resulting from a call to [doFloodlightAnalysis].
#' @param ... ignored.
#' 
plot.cSEMFloodlight <- function(x, ...) {
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  plot1 <- ggplot2::ggplot(as.data.frame(x$out), ggplot2::aes(x = x$out[, 'value_z'], 
                                                              y = x$out[, 'direct_effect'])) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = x$out[,'lb'], ymax = x$out[, 'ub']), alpha = 0.2) +
    ggplot2::labs(
      "x" = paste('Level of ', x$Information['moderator']), 
      "y" = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent'])) +
    ggplot2::theme_bw() +
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  # Add Johnson-Neyman points, if they exist in the considered range
  if(length(x$Johnson_Neyman_points$JNlb) == 2){
    JN <- x$Johnson_Neyman_points$JNlb
    plot1 <- plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  if(length(x$Johnson_Neyman_points$JNub) == 2){
    JN = x$Johnson_Neyman_points$JNub
    plot1 = plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  # Plot
  plot1
}