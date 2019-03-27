#' Internal: AVE
#'
#' Calculate the average variance extracted (AVE) as proposed by 
#' \insertCite{Fornell1981;textual}{cSEM}.
#'
#' The AVE is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret AVE results in the context of a 
#' composite model. It is possible to force computation of the AVE for constructs
#' modeled as composites as well by setting `.only_common_factors = FALSE`, 
#' however, we explicitly warn to interpret results with caution, 
#' as they may not even have a conceptual meaning.
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateAVE(
#'  .object              = NULL,
#'  .only_common_factors = TRUE,
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
  
  # By default AVE's for constructs modeled as composites are not returned
  if(.only_common_factors){
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    AVEs      <- AVEs[names_cf]
  }
  
  # Return vector of AVEs
  return(AVEs) 
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
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateGoF(
#'  .object              = NULL,
#'  .only_common_factors = TRUE,
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
#' Compute several reliability measures. See the vignette for details.
#' 
#' Reliability is the consistency of measurement. Practically, reliability 
#' is empirically assessed based on the classical true score measurement theory 
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
#' For the the tau-equivalent reliability (rhoT) a closed-form confidence 
#' interval may be computed \insertCite{Trinchera2018}{cSEM} by setting
#' `.calculate_ci = TRUE` (default is `FALSE`). If `.alpha` is a vector
#' several CI's are returned.
#' 
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#'
#' @return For `calculateRhoC()` a named numeric vector containing the reliability
#'   estimates.
#'   For `calculateRhoT()` a `data.frame` with as many rows as there are
#'   constructs modeled as common factors in the model (unless 
#'   `.only_common_factors = FALSE` in which case the number of rows equals the
#'   total number of constructs in the model). The first column contains the name of the construct.
#'   The second column the reliability estimate.
#'   If `.calculate_ci = TRUE` the remaining columns contain lower and upper bounds
#'   for the (1 - `.alpha`) confidence interval(s).
#'   
#' @inheritParams csem_arguments
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
  .only_common_factors = TRUE) {

  # Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  W <- .object$Estimates$Weight_estimates
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
#'                         as Cronbach alpha or coefficient alpha. Since
#'                         indicators are always standarized in `cSEM`, 
#'                         tau-equivalent reliability is identical to the
#'                         parallel reliability.
#' 
calculateRhoT <- function(
  .object              = NULL,
  .alpha               = args_default()$.alpha,
  .calculate_ci        = args_default()$.calculate_ci,
  .only_common_factors = TRUE) {
  
  ## Only applicable to objects of class cSEMResults_default
  if(!any(class(.object) == "cSEMResults_default")) {
    stop2("`", match.call()[1], "` only applicable to object of",
          " class `cSEMResults_default`. Use `assess()` instead.")
  }
  
  ## If CI's are computed:
  # Calculation of the CIs are based on Trinchera et al. (2018).
  # In the paper, the CI was proposed and studied using Monte Carlo simulation
  # assuming scores are build as sum scores. Therefore a warning is
  # given if a weightning scheme other than "unit" is used, since they have
  # not been formally studied yet.
  
  if(.object$Information$Arguments$.approach_weights != "unit" & .calculate_ci == TRUE) {
    warning2("Calculation of the confidence interval (CI) for rhoT (Cronbach alpha)",
             " was proposed and studied assuming sum scores.\n",
             "Statistical properties of the CI for rhoT based on a weighted composite",
             " obtained using weight approach ", 
             paste0("`", .object$Information$Arguments$.approach_weights, "`"), 
             " have not been formally examined yet.")
  }
  
  ## Get relevant objects
  W <- .object$Estimates$Weight_estimates
  S <- .object$Estimates$Indicator_VCV
  
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
    if(.calculate_ci) {
      # Calculation of the CIs are based on Trinchera et al. (2018).
      # The code for the calculation of the CIs is addpated from the paper.

      K <- length(indicator_names)
      X <- .object$Information$Data[, indicator_names]
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
        z_value <- stats::qnorm((1 - i/2), mean = 0, sd = 1)
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
#' The function is only applicable to objects of class `cSEMResults_default`.
#' For other object classes use [assess()].
#' 
#' @usage calculateHTMT(
#'  .object              = NULL,
#'  .only_common_factors = TRUE
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [assess], [csem], [cSEMResults]
#' @keywords internal

calculateHTMT <- function(
  .object              = NULL,
  .only_common_factors = TRUE
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
#' Calculate the difference between the empirical (S) 
#' and the model-implied indicator variance-covariance matrix (Sigma_hat)
#' using different distance measures. See vignette assess.
#'
#' The distance functions are only applicable to objects of class 
#' `cSEMResults_default`. For other object classes use [assess()].
#' 
#' @return A single numeric value.
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

#' @describeIn distance_measures The squared Euclidian distance

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
#' is strongly discouraged.
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