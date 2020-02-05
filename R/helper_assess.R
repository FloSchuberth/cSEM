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
#' @seealso [assess()], [cSEMResults]
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
#' Calculate the degrees of freedom for a given model from a [cSEMResults] object.
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
#' @seealso [assess()], [cSEMResults]
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
  
  ## Number of correlations between (exogenous) constructs
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
#' @seealso [assess()], [cSEMResults]
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
#' Reliability is traditionally based on a test score (proxy) based on unit weights.
#' To compute congeneric and tau-equivalent reliability based on a score that 
#' uses the weights of the weight approach used to obtain `.object` use `.weighted = TRUE` 
#' instead.
#' 
#' For the tau-equivalent reliability ("`rho_T`" or "`cronbachs_alpha`") a closed-form 
#' confidence interval may be computed \insertCite{Trinchera2018}{cSEM} by setting
#' `.closed_form_ci = TRUE` (default is `FALSE`). If `.alpha` is a vector
#' several CI's are returned.
#' 
#' The function is only applicable to objects inheriting class `cSEMResults_default`.
#' For objects of class `cSEMResults_multi` and `cSEMResults_2ndorder` use [assess()].
#'
#' @return For `calculateRhoC()` and `calculateRhoT()` (if `.output_type = "vector"`) 
#'   a named numeric vector containing the reliability estimates.
#'   If `.output_type = "data.frame"` `calculateRhoT()` returns a `data.frame` with as many rows as there are
#'   constructs modeled as common factors in the model (unless 
#'   `.only_common_factors = FALSE` in which case the number of rows equals the
#'   total number of constructs in the model). The first column contains the name of the construct.
#'   The second column the reliability estimate.
#'   If `.closed_form_ci = TRUE` the remaining columns contain lower and upper bounds
#'   for the (1 - `.alpha`) confidence interval(s).
#'   
#' @inheritParams csem_arguments
#' @param .output_type Character string. The type of output. One of "vector" or
#'   "data.frame". Defaults to "vector".
#' @param .model_implied Logical. Should weights be scaled using the model-implied
#'   indicator correlation matrix? Defaults to `TRUE`.
#' @param ... Ignored.
#'
#' @seealso [assess()], [cSEMResults]
#'
#' @references 
#' 
#' \insertAllCited{}
#' 
#' @keywords internal
#' @name reliability
NULL

#' @describeIn reliability Calculate the congeneric reliability
calculateRhoC <- function(
  .object              = NULL,
  .model_implied       = TRUE,
  .only_common_factors = TRUE,
  .weighted            = args_default()$.weighted
  ) {

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to objects of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  # Note (21.01.2020): The actual formula for congeneric reliability assuming
  #   1. constructs modeled as common factors
  #   2. scores build using unit weights
  #   is given by (notation from cSEM website):
  #
  #      rhoC = Var(eta_bar)/ Var(eta_hat_k) = (sum lambda_k)^2 / (sum lambda_k)^2 + Var(epsilon_bar)
  #
  # Reliability in general is given by
  # 
  #  rhoC = Var(eta_bar)/ Var(eta_hat_k) = (w'lambda)^2 / w'S w
  # We can write the formula es
  # In cSEM weights are chosen such that Var(eta_hat_k) is 1. We always use
  
  ## Get relevant objects
  if(.weighted) {
    W <- .object$Estimates$Weight_estimates
  } else {
    W <- .object$Information$Model$measurement
  }
  
  if(.model_implied) {
    ## Redefine the construct types: all composites to common factor
    ## Reason: to compute Joereskogs rho we need the model-implied indicator
    ## correlation matrix as if the model was a common factor model.
    c_type_original <- .object$Information$Model$construct_type 
    c_type <- c_type_original
    c_type[c_type == "Composite"] <- "Common factor" 
    # Replace
    .object$Information$Model$construct_type <- c_type
    
    indicator_vcv <- fit(.object)

    # Reset
    .object$Information$Model$construct_type <- c_type_original
  } else {
    indicator_vcv <- .object$Estimates$Indicator_VCV
  }
  
  W <- scaleWeights(indicator_vcv, W)
  
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

#' @describeIn reliability Calculate the tau-equivalent reliability
#' 
calculateRhoT <- function(
  .object              = NULL,
  .alpha               = args_default()$.alpha,
  .closed_form_ci      = args_default()$.closed_form_ci,
  .only_common_factors = TRUE,
  .output_type         = c("vector", "data.frame"),
  .weighted            = args_default()$.weighted,
  ...
  ) {
  .output_type <- match.arg(.output_type)
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
  out_df  <- data.frame()
  out_vec <- c()
  
  for(j in rownames(W)) {
    indicator_names <- colnames(W[j, W[j,] != 0, drop = FALSE])
    S_jj            <- S[indicator_names, indicator_names]
    rho_bar_j <- mean(S_jj[upper.tri(S_jj)])
    rhoT      <- rho_bar_j * sum(W[j, ])^2  
    
    # Add to vector 
    out_vec[j] <- rhoT
    
    ## Add to data.frame
    out_df[j, "Construct"] <- j
    out_df[j, "Estimate"]  <- rhoT
    
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
      
      ## If .closed_form_ci == TRUE then .output_type must be "data.frame"
      if(.output_type != "data.frame") {
        stop2(
          "The following error occured in the `calculateRhoT()` function:\n",
          "Output type must be 'data.frame' if `.closed_form_ci = TRUE`"
        )
      }
      
      K <- length(indicator_names)
      X <- .object$Information$Data[, indicator_names, drop = FALSE]
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
        
        out_df[j, name_L] <- low
        out_df[j, name_U] <- up
      }
    }
  }
  
  # By default only reliabilities for constructs modeled as common factors are returned
  if(.only_common_factors){
    con_types <- .object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    out_vec   <- out_vec[names_cf]
    out_df    <- out_df[names_cf, ]
  }
  
  # Return reliabilities
  if(.output_type == "vector") {
    return(out_vec) 
  } else {
    return(out_df)
  }
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
#' @seealso [assess()], [csem], [cSEMResults]
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
  0.5 * sum((S - Sigma_hat)^2)
}

#' @describeIn  distance_measures The distance measure used by FIML

calculateDML <- function(
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
        "The following error occured in the calculateDML() function:\n",
        "`.object` must be of class `cSEMResults_default` or `cSEMResults_2ndorder`.")
    }
    
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')  
  }

  p <- dim(S)[1]
  
  # This is the distance function. The test statistic is T_ML = (n-1) or n * DML! 
  sum(diag(S %*% solve(Sigma_hat))) - log(det(S%*%solve(Sigma_hat))) - p
}

#' Internal: Fit measures
#' 
#' Calculate fit measures.
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

#' @describeIn fit_measures The ChiSquare statistic.

calculateChiSquare <- function(.object) {
  
  n    <- nrow(.object$Information$Data)
  F0   <- calculateDML(.object)
  
  (n - 1) * F0
}

#' @describeIn fit_measures The ChiSquare statistic divided by its degrees of freedom.

calculateChiSquareDf <- function(.object) {
  
  n    <- nrow(.object$Information$Data)
  F0   <- calculateDML(.object)
  
  ((n - 1) * F0) / calculateDf(.object)
}

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
  
  S <- .object$Estimates$Indicator_VCV
  p <- dim(S)[1]
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  (F0 - FT) / F0
}

#' @describeIn fit_measures The non-normed fit index (NNFI). Also called the Tucker-Lewis index (TLI).

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
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = args_default()$.saturated,
  ...
) {

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
    
    # The SRMR as calculated by us is always based on the the difference 
    # between correlation matrices.
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator') 
  }
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 




#' Internal: Calculate Cohens f^2
#'
#' Calculate the effect size for regression analysis \insertCite{Cohen1992}{cSEM}
#' known as Cohen's f^2
#'
#' @usage calculatef2(.object = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @return A matrix with as many rows as there are structural equations. The 
#'   number of columns is equal to the total number of right-hand side variables
#'   of these equations.
#' 
#' @seealso [assess()], [csem], [cSEMResults]
#' @keywords internal

calculatef2 <- function(.object = NULL) {
  
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
  
  ## Loop over each structural equation
  outer_out <- lapply(vars_endo, function(x) {
    
    # Get names of the independent variables of equation x
    indep_vars <- colnames(s[x , s[x, ] != 0, drop = FALSE])
    
    # Compute R_excluded by regressing variable x on all other variables of
    # that equation except the independent variable i
    inner_out <- lapply(indep_vars, function(i) {
      # Update csem_model
      model_temp <- csem_model
      model_temp$structural[x, i] <- 0 
      
      # If there is only one independent variable the R^2_excluded is 0
      # since s_u_dach^2 = s_y^2
      if(sum(model_temp$structural[x, ]) > 0) { 
        out <- estimatePath(
          .approach_nl      = approach_nl,
          .approach_paths   = approach_paths,
          .approach_weights = approach_weights,
          .csem_model       = model_temp,
          .H                = H,
          .normality        = normality,
          .P                = P,
          .Q                = Q
        )
        ## Calculate effect size
        r2_excluded <- out$R2[x]
      } else {
        r2_excluded <- 0
      }
      names(r2_excluded) <- x
      r2_included <- .object$Estimates$R2[x]
      
      effect_size <- unname((r2_included - r2_excluded)/(1 - r2_included))
    })
    inner_out <- unlist(inner_out)
    names(inner_out) <- indep_vars
    inner_out 
  })
  names(outer_out) <- vars_endo
  
  ## Remove 
  
  ## Make output a matrix
  # Note: this is necessary for calculatef2() to work
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
#' @seealso [assess()], [cSEMResults]
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

#' Helper for assess()
#' @noRd
printEffects <- function(.effect, .ci_colnames, .what = "Total effect") {
  l <- max(nchar(.effect[, "Name"]), nchar(.what))
  
  if(length(.ci_colnames) != 0) {
    xx <- regmatches(.ci_colnames, regexpr("\\.", .ci_colnames), invert = TRUE)
    interval_names    <- unique(sapply(xx, `[`, 1))
    sig_level_names   <- unique(gsub("[LU]", "", sapply(xx, `[`, 2)))
    
    cat2("\n  ",  col_align("", width = max(l, nchar(.what)) + 44))
    for(i in interval_names) {
      cat2(col_align(i, width = 20*length(sig_level_names), align = "center"))
    }
  }
  cat2(
    "\n  ", 
    col_align(.what, l + 2), 
    col_align("Estimate", 10, align = "right"), 
    col_align("Std. error", 12, align = "right"),
    col_align("t-stat.", 10, align = "right"), 
    col_align("p-value", 10, align = "right")
  )
  if(length(.ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat2(col_align(i, 20, align = "center"))
    } 
  }
  
  for(i in 1:nrow(.effect)) {
    cat2(
      "\n  ", 
      col_align(.effect[i, "Name"], l + 2), 
      col_align(sprintf("%.4f", .effect[i, "Estimate"]), 10, align = "right"),
      col_align(sprintf("%.4f", .effect[i, "Std_err"]), 12, align = "right"),
      col_align(sprintf("%.4f", .effect[i, "t_stat"]), 10, align = "right"),
      col_align(sprintf("%.4f", .effect[i, "p_value"]), 10, align = "right")
    )
    if(length(.ci_colnames) != 0) {
      for(j in seq(1, length(.ci_colnames), by = 2) + 6) {
        cat2(
          col_align(
            paste0("[", sprintf("%7.4f", .effect[i, j]), ";", 
                   sprintf("%7.4f", .effect[i, j+1]), " ]"), 20, align = "center")
        )
      } 
    }
  }
}