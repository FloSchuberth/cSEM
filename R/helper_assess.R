#' Model selection criteria
#' 
#' Calculate several information or model selection criteria (MSC) such as the
#' Akaike information criterion (AIC), the Bayesian information criterion (BIC) or
#' the Hannan-Quinn criterion (HQ).
#'
#' By default, all criteria are calculated (`.ms_criterion == "all"`). To compute only
#' a subset of the criteria a vector of criteria may be given.
#' 
#' If `.by_equation == TRUE` (the default), the criteria are computed for each 
#' structural equation of the model separately, as suggested by 
#' \insertCite{Sharma2019;textual}{cSEM} in the context of PLS. The relevant formula can be found in
#'  Table B1 of the appendix of \insertCite{Sharma2019;textual}{cSEM}. 
#'  
#' If `.by_equation == FALSE` the AIC, the BIC and the HQ for whole model 
#' are calculated. All other criteria are currently ignored in this case! 
#' The relevant formulae are (see, e.g., \insertCite{Akaike1974}{cSEM},
#' \insertCite{Schwarz1978;textual}{cSEM}, 
#' \insertCite{Hannan1979;textual}{cSEM}): 
#' 
#' \deqn{AIC = - 2*log(L) + 2*k}
#' \deqn{BIC = - 2*log(L) + k*ln(n)}
#' \deqn{HQ  = - 2*log(L) + 2*k*ln(ln(n))}
#' 
#' where log(L) is the log likelihood function of the multivariate normal
#' distribution of the observable variables, k the (total) number of estimated parameters,
#' and n the sample size.
#' 
#' If `.only_structural == TRUE`, log(L) is based on the structural model only.
#' The argument is ignored if `.by_equation == TRUE`.
#' 
#' @usage calculateModelSelectionCriteria(
#'   .object          = NULL,
#'   .ms_criterion    = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
#'                        "hqc", "mallows_cp"),
#'   .by_equation     = TRUE, 
#'   .only_structural = TRUE 
#'   )
#'
#' @return If `.by_equation == TRUE` a named list of model selection criteria.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess()], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @export

calculateModelSelectionCriteria <- function(
  .object          = NULL, 
  .ms_criterion    = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                       "hqc", "mallows_cp"),
  .by_equation     = TRUE,
  .only_structural = TRUE
  ){
  
  if(inherits(.object, "cSEMResults_multi")) {
    
    out <- lapply(.object, calculateModelSelectionCriteria,
                  .ms_criterion    = .ms_criterion,
                  .by_equation     = .by_equation,
                  .only_structural = .only_structural)
    return(out)
    
  } else if(inherits(.object, "cSEMResults_default")) {
    
    x1 <- .object$Estimates
    x2 <- .object$Information$Model
    S  <- x1$Indicator_VCV
    n <- nrow(.object$Information$Data)
    
    ## Mean square of the saturated model
    args <- .object$Information$Arguments
    args$.model$structural[lower.tri(args$.model$structural)] <- 1
    
    out_saturated <- do.call(csem, args)
    
    MSE <- (1 - out_saturated$Estimates$R2)[names(x1$R2)]
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    
    x1 <- summarize(.object)$Second_stage$Estimates
    x2 <- summarize(.object)$Second_stage$Information$Model
    S  <- .object$First_stage$Estimates$Indicator_VCV
    # Sample size
    n <- nrow(.object$First_stage$Information$Data)
    
    ## Mean square of the saturated model
    args <- .object$Second_stage$Information$Arguments
    args$.model$structural[lower.tri(args$.model$structural)] <- 1
    
    out_saturated <- do.call(csem, args)
    
    MSE <- (1 - out_saturated$Estimates$R2)
    # remove _temp suffix in second stage
    names(MSE) <- gsub("_temp","", names(MSE))
    MSE <- MSE[names(x1$R2)]
    
  } else {
    stop2(
      "The following error occured in the calculateModelSelectionCriteria() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  ## Model-implied indicator covariance matrix
  Sigma_hat <- fit(.object = .object, 
                   .saturated  = FALSE, 
                   .type_vcv = "indicator"
                   )
  ## Number of estimated parameters
  # Note: # estimated parameters = total number of elements of S - df
  k <- sum(lower.tri(S)) - calculateDf(.object)
  ## Number of indicators
  p <- nrow(S)
  
  ## Set up empty list
  out <- list()
  
  if(!.by_equation) {
    if(.only_structural) {
      # Reassign S, Sigma_hat, k and p. Now:
      # S         := construct correlation matrix
      # Sigma_hat := model-implied construct correlation matrix
      # k         := # of estimated parameters of the structural model (either 
      #              correlations or structural model parameters)
      # p         := # of constructs in the structural model
      S <- x1$Construct_VCV
      Sigma_hat <- fit(.object, 
                       .saturated = FALSE, 
                       .type_vcv = 'construct')
      
      # Count construct correlations and/or structural model parameters
        if(!all(x2$structural == 0)) {
          # Free correlations between exogenous constructs
          n_cor_exo <- length(x2$cons_exo) * (length(x2$cons_exo) - 1) / 2 
          
          # Number of structural parameters
          n_structural <- sum(x2$structural)
        } else {
          n_cor_exo    <- nrow(x2$structural) * (nrow(x2$structural) - 1) / 2
          n_structural <- 0
        }
      k <- n_cor_exo + n_structural
      p <- nrow(S)
    }

    ## Log Likelihood
    logL <- - n*p/2*log(2*pi) -n/2*sum(diag(S %*% solve(Sigma_hat))) - n/2*log(det(Sigma_hat))
    
    if(any(.ms_criterion %in% c("all", "aic"))) {
      ## Compute AIC
      out[["AIC"]] <- -2*logL + 2*k 
    }
    if(any(.ms_criterion %in% c("all", "bic"))) {
      ## Compute BIC
      out[["BIC"]] <- -2*logL + k*log(n) 
    }
    if(any(.ms_criterion %in% c("all", "hq"))) {
      ## Compute HQ
      out[["HQ"]] <- -2*logL + 2*k*log(log(n)) 
    }
  } else {
    ### Implementation based on Sharma et al. (2019), p.391, Table B1
    ## 
    SSE <- (1 - x1$R2)*(n - 1)
    
    ## Update k to contain the number of parameters for each equation instead of
    ## all model parameters.
    ## Note: in line with Sharma et. al (2019) the number of paramters per equation
    ##       is computed as number of regressors + 1 (for the constant).
    ##       Its unclear if the constant is really necessary since constructs are
    ##       standardized. Hence, the constant will always be zero.
    ## 
    k <- rowSums(x2$structural)[names(x1$R2)] + 1
   
    if(any(.ms_criterion %in% c("all", "aic"))) {
      out[["AIC"]] <- n*(log(SSE / n) + 2*k/n)
    }
    if(any(.ms_criterion %in% c("all", "aicc"))) {
      out[["AICc"]] <- n*(log(SSE / n) + (n + k)/(n - k -2))
    }
    if(any(.ms_criterion %in% c("all", "aicu"))) {
      out[["AICu"]] <- n*(log(SSE / (n - k)) + 2*k/n)
    }
    if(any(.ms_criterion %in% c("all", "bic"))) {
      out[["BIC"]] <- n*(log(SSE / n) + k*log(n)/n)
    }
    if(any(.ms_criterion %in% c("all", "fpe"))) {
      out[["FPE"]] <- (SSE / (n-k)) * (1 + k/n)
    }
    if(any(.ms_criterion %in% c("all", "gm"))) {
      out[["GM"]] <- (SSE / MSE) + k*log(n)
    }
    if(any(.ms_criterion %in% c("all", "hq"))) {
      out[["HQ"]] <- n*(log(SSE / n) + (2*k*log(log(n)))/n)
    }
    if(any(.ms_criterion %in% c("all", "hqc"))) {
      out[["HQc"]] <- n*(log(SSE / n) + (2*k*log(log(n)))/(n - k - 2))
    }
    if(any(.ms_criterion %in% c("all", "mallows_cp"))) {
      out[["Mallows_Cp"]] <- (SSE / MSE) - (n - 2*k)
    }
  }
  
  return(out)
}

#' Average variance extracted (AVE)
#'
#' Calculate the average variance extracted (AVE) as proposed by 
#' \insertCite{Fornell1981;textual}{cSEM}. For details see the
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#ave}{cSEM website} 
#'
#' The AVE is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret the AVE in the context of a 
#' composite model. It is possible, however, to force computation of the AVE for constructs 
#' modeled as composites by setting `.only_common_factors = FALSE`.
#' 
#' @usage calculateAVE(
#'  .object              = NULL,
#'  .only_common_factors = TRUE
#' )
#'
#' @return A named vector of numeric values (the AVEs). If `.object` is a list 
#'   of `cSEMResults` objects, a list of AVEs is returned.
#'   
#' @inheritParams csem_arguments
#'
#' @seealso [assess()], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @export

calculateAVE <- function(
  .object              = NULL,
  .only_common_factors = TRUE
  ){

  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateAVE, .only_common_factors = .only_common_factors)
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    ## Extract loadings
    Lambda   <- .object$Estimates$Loading_estimates
    c_names  <- rownames(Lambda)
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    # ## Extract construct type
    c_names2 <- .object$Second_stage$Information$Arguments_original$.model$vars_2nd
    
    out <- lapply(.object, calculateAVE, .only_common_factors = .only_common_factors)
    out$Second_stage <- out$Second_stage[c_names2]
    out <- if(is.na(out$Second_stage)) {
      out$First_stage
    } else {
      out <- c(out$First_stage, out$Second_stage)
    }
    return(out)
    
  } else {
    stop2(
      "The following error occured in the calculateAVE() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
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
    ## Extract construct type
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    AVEs      <- AVEs[names_cf]
  }
  
  # Return vector of AVEs
  return(AVEs) 
}


#' Degrees of freedom
#' 
#' Calculate the degrees of freedom for a given model from a [cSEMResults] object.
#' 
#' Although, composite-based estimators always retrieve parameters of the 
#' postulated models via the estimation of a composite model,
#' the computation of the degrees of freedom depends on the postulated model.
#' 
#' See: \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{cSEM website} 
#' for details on how the degrees of freedom are calculated.
#' 
#' To compute the degrees of freedom of the null model use `.null_model = TRUE`.
#' The degrees of freedom of the null model are identical to the number of
#' non-redundant off-diagonal elements of the empirical indicator correlation matrix.
#' This implicitly assumes a null model with model-implied indicator correlation
#' matrix equal to the identity matrix.
#' 
#' @usage calculateDf(
#'   .object     = NULL,
#'   .null_model = FALSE,
#'   ...
#'   )
#' 
#' @return A single numeric value.
#'   
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @seealso [assess()], [cSEMResults]
#' @export

calculateDf <- function(
  .object     = NULL, 
  .null_model = FALSE,
  ...
) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    
    df <- lapply(.object, calculateDf)
    
    ## Return
    names(df) <- names(.object)
    return(df)
    
  } else if(inherits(.object, "cSEMResults_default")) { 
    
    x1 <- .object$Estimates
    x2 <- .object$Information$Model
    
    if(x2$model_type == "Nonlinear") {
      warning("Computation of the degrees of freedom for models containing nonlinear",
              " terms is experimental. Computation may not be correct.")
    }
    
    # Number of non-redundant off-diagonal elements of the indicator covariance 
    # matrix (S)
    vS <- sum(lower.tri(x1$Indicator_VCV)) # same as dim_nS *(dim_nS - 1) / 2
    
    # Caculate the degrees of freedom of a null model (Sigma_hat = I)
    if(.null_model) {
      # The degrees of freedom of the null model are identical to 
      # the number of non-redundant elements of S (since there is nothing to
      # estimate. Everything, except for the main diagonal element, is set to 0 a 
      # priori).
      return(vS)
    }
    
    # Construct correlations and structural model parameters
    if(!all(x2$structural == 0)) {
      # Free correlations between exogenous constructs
      n_cor_exo <- length(x2$cons_exo) * (length(x2$cons_exo) - 1) / 2 
      
      # Correlations between endogenous constructs
      names_cor_endo <- intersect(x2$cons_endo, rownames(x2$cor_specified))
      if(length(names_cor_endo) != 0) {
        n_cor_endo <- sum(x2$cor_specified[names_cor_endo, names_cor_endo, drop = FALSE]
                          [lower.tri(x2$cor_specified[names_cor_endo, names_cor_endo, drop = FALSE])])
      } else {
        n_cor_endo <- 0
      }
      
      # Number of structural parameters
      n_structural <- sum(x2$structural)
      
    } else {
      n_cor_exo    <- nrow(x2$structural) * (nrow(x2$structural) - 1) / 2
      n_cor_endo   <- 0
      n_structural <- 0
    }
    
    ## Construct names
    names_constructs <- rownames(x2$structural)
    k <- c()
    for(j in names_constructs) {
      if(x2$construct_type[j] == "Composite") {
        
        ## Number of weights minus 1 (since weights are choosen s.t. Var(eta) = 1)
        n_weights <- sum(x2$measurement[j, ]) - 1 
        
        ## Number of free non-redundant off-diagonal element of each intra-block
        #+ covariance matrix
        temp    <- sum(x2$measurement[j, ] == 1)
        n_intra <- temp * (temp - 1) / 2
        
        ## Calculate dfs per construct
        k[j] <- n_weights + n_intra 
        
      } else {
        
        # Number of loadings (since either 1 loading is fixed or the variance 
        # of the construct is fixed)
        n_loadings <- sum(x2$measurement[j, ] == 1)
        
        # Number of measurement errors assumed to correlate
        j_names <- names(which(x2$measurement[j, ] == 1))
        n_error <- sum(x2$error_cor[j_names, j_names, drop = FALSE] == 1) / 2
        
        # Calculate dfs per construct
        k[j] <- n_loadings + n_error
      }
    }
    
    df_total <- vS - (n_cor_exo + n_cor_endo + n_structural) - sum(k)
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    
    x11 <- .object$First_stage$Estimates
    x12 <- .object$First_stage$Information$Model
    
    x21 <- .object$Second_stage$Estimates
    x22 <- .object$Second_stage$Information$Model
  
    if(x22$model_type == "Nonlinear") {
      warning("Computation of the degrees of freedom for models containing nonlinear",
              " terms is experimental. Computation may not be correct.")
    }
    
    # Number of non-redundant off-diagonal elements of the indicator covariance 
    # matrix (S) of the first stage
    vS <- sum(lower.tri(x11$Indicator_VCV)) # same as dim_nS *(dim_nS - 1) / 2
    
    # Caculate the degrees of freedom of a null model (Sigma_hat = I)
    if(.null_model) {
      # The degrees of freedom of the null model are identical to 
      # the number of non-redundant elements of S (since there is nothing to
      # estimate. Everything is set to 0 a priori)
      return(vS)
    }
    
    # Construct correlations and structural model parameters
    if(!all(x22$structural == 0)) {
      # Free correlations between exogenous constructs
      n_cor_exo <- length(x22$cons_exo) * (length(x22$cons_exo) - 1) / 2 
      
      # Correlations between endogenous constructs
      names_cor_endo <- intersect(x22$cons_endo, rownames(x22$cor_specified))
      if(!is.null(names_cor_endo)) {
        n_cor_endo <- sum(x22$cor_specified[names_cor_endo, names_cor_endo, drop = FALSE]
                          [lower.tri(x22$cor_specified[names_cor_endo, names_cor_endo, drop = FALSE])])
      } else {
        n_cor_endo <- 0
      }
      
      # Number of structural parameters
      n_structural <- sum(x22$structural)
      
    } else {
      n_cor_exo    <- nrow(x22$structural) * (nrow(x22$structural) - 1) / 2
      n_cor_endo   <- 0
      n_structural <- 0
    }
    
    ## Construct names
    construct_order  <- .object$First_stage$Information$Model$construct_order 
    names_constructs <- names(construct_order)
    k <- c()
    for(j in names_constructs) {
      if(construct_order[j] == "Second order") {
        if(x22$construct_type[j] == "Composite") {
          ## Number of weights minus 1 (since weights are choosen s.t. Var(eta) = 1)
          n_weights <- sum(x22$measurement[j, ]) - 1 
          
          ## Number of free non-redundant off-diagonal element of each intra-block
          #+ covariance matrix
          temp    <- sum(x22$measurement[j, ] == 1)
          n_intra <- temp * (temp - 1) / 2
          
          ## Calculate dfs per construct
          k[j] <- n_weights + n_intra 
          
        } else {
          
          # Number of loadings -1 (since either 1 loading is fixed or the variance 
          # of the construct is fixed)
          n_loadings <- sum(x22$measurement[j, ] == 1)
          
          # Number of measurement errors assumed to correlate
          n_error <- sum(x22$error_cor == 1) / 2
          
          # Calculate dfs per construct
          k[j] <- n_loadings + n_error
        }
      } else {
        if(x12$construct_type[j] == "Composite") {
          ## Number of weights minus 1 (since weights are choosen s.t. Var(eta) = 1)
          n_weights <- sum(x12$measurement[j, ]) - 1
          
          ## Number of free non-redundant off-diagonal element of each intra-block
          #+ covariance matrix
          temp    <- sum(x12$measurement[j, ] == 1)
          n_intra <- temp * (temp - 1) / 2
          
          ## Calculate dfs per construct
          k[j] <- n_weights + n_intra 
          
        } else {
          
          # Number of loadings -1 (since either 1 loading is fixed or the variance 
          # of the construct is fixed)
          n_loadings <- sum(x12$measurement[j, ] == 1)
          
          # Number of measurement errors assumed to correlate
          n_error <- sum(x12$error_cor == 1) / 2
          
          # Calculate dfs per construct
          k[j] <- n_loadings + n_error
        }
      }
    }
  } else {
    stop2(
      "The following error occured in the calculateDL() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  ## Compute degrees of freedom
  df_total <- vS - (n_cor_exo + n_cor_endo + n_structural) - sum(k)
  
  ## Return degrees of freedom
  df_total
}



#' Fornell-Larcker criterion
#' 
#' Computes the Fornell-Larcker matrix.
#' 
#' The Fornell-Larcker criterion (FL criterion) is a rule suggested by \insertCite{Fornell1981;textual}{cSEM}
#' to assess discriminant validity. The Fornell-Larcker
#' criterion is a decision rule based on a comparison between the squared
#' construct correlations and the average variance extracted (AVE).
#' 
#' The FL criterion is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret the FL criterion in the context of a 
#' model that contains constructs modeled as composites.
#' 
#' @usage calculateFLCriterion(
#'   .object              = NULL,
#'   .only_common_factors = TRUE,
#'   ...
#'   )
#'
#' @return A matrix with the squared construct correlations on the off-diagonal and 
#' the AVE's on the main diagonal.
#'   
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @seealso [assess()], [cSEMResults]
#'
#' @references 
#' \insertAllCited{}
#' @export

calculateFLCriterion <- function(
  .object              = NULL,
  .only_common_factors = TRUE,
  ...
  ) {
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateFLCriterion,
                  .only_common_factors = .only_common_factors)
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    # Fornell-Larcker
    ## Get relevant objects
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    P         <- .object$Estimates$Construct_VCV
    
    if(.only_common_factors) {
      P <- P[names_cf, names_cf]
    }
    
    if(sum(dim(P)) > 0) {
      FL_matrix <- cov2cor(P)^2
      diag(FL_matrix) <- calculateAVE(.object, 
                                      .only_common_factors = .only_common_factors)
      FL_matrix
    } 
    
  } else { # 2nd order
    stop2(
      "Computation of the Fornell-Larcker criterion",
      " not supported for models containing second-order constructs:\n")
  }
}



#' Goodness of Fit (GoF)
#'
#' Calculate the Goodness of Fit (GoF) proposed by \insertCite{Tenenhaus2004;textual}{cSEM}. 
#' Note that, contrary to what the name suggests, the GoF is **not** a 
#' measure of model fit in the sense of SEM. See e.g. \insertCite{Henseler2012a;textual}{cSEM}
#' for a discussion.
#' 
#' The GoF is inherently tied to the common factor model. It is therefore 
#' unclear how to meaningfully interpret the GoF in the context of a 
#' model that contains constructs modeled as composites.
#' 
#' @usage calculateGoF(
#'  .object              = NULL
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
#' @export

calculateGoF <- function(
  .object              = NULL
){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateGoF)
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    ## Get relevant quantities
    Lambda    <- .object$Estimates$Loading_estimates
    R2        <- .object$Estimates$R2
    
    # Select only non-zero loadings
    L  <- Lambda[Lambda != 0]
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    c_names2  <- .object$Second_stage$Information$Arguments_original$.model$vars_2nd
    
    ## Extract loadings
    Lambda   <- .object$First_stage$Estimates$Loading_estimates
    Lambda2  <- .object$Second_stage$Estimates$Loading_estimates 
    R2       <- .object$Second_stage$Estimates$R2
    
    Lambda2   <- Lambda2[c_names2, ]

    L  <- Lambda[Lambda != 0]
    L2 <- Lambda2[Lambda2 != 0]
    
    L <- c(L, L2)
  } else {
    stop2(
      "The following error occured in the calculateGoF() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  # The GoF is defined as the sqrt of the mean of the R^2s of the structural model 
  # times the variance in the indicators that is explained by the construct (lambda^2).

  gof <- sqrt(mean(L^2) * mean(R2))
  
  return(gof)
}



#' Reliability
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
#' @name reliability
NULL

#' @describeIn reliability Calculate the congeneric reliability
#' @export
calculateRhoC <- function(
  .object              = NULL,
  .model_implied       = TRUE,
  .only_common_factors = TRUE,
  .weighted            = FALSE
  ) {

  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateRhoC, 
                  .model_implied       = .model_implied,
                  .only_common_factors = .only_common_factors,
                  .weighted            = .weighted)
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    # continue
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    # ## Extract construct type
    c_names2 <- .object$Second_stage$Information$Arguments_original$.model$vars_2nd
    
    out <- lapply(.object, calculateRhoC, 
                  .model_implied       = .model_implied,
                  .only_common_factors = .only_common_factors,
                  .weighted            = .weighted)
    
    out$Second_stage <- out$Second_stage[c_names2]
    out <- if(is.na(out$Second_stage)) {
      out$First_stage
    } else {
      c(out$First_stage, out$Second_stage)
    }
    return(out)
    
  } else {
    stop2(
      "The following error occured in the calculateRhoC() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
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
#' @export
calculateRhoT <- function(
  .object              = NULL,
  .alpha               = 0.05,
  .closed_form_ci      = FALSE,
  .only_common_factors = TRUE,
  .output_type         = c("vector", "data.frame"),
  .weighted            = FALSE,
  ...
  ) {

  # Match arguments
  .output_type <- match.arg(.output_type)
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateRhoT,
                  .alpha               = .alpha,
                  .closed_form_ci      = .closed_form_ci,
                  .only_common_factors = .only_common_factors,
                  .output_type         = .output_type,
                  .weighted            = .weighted
                  )
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    # continue
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    # ## Extract construct type
    c_names2 <- .object$Second_stage$Information$Arguments_original$.model$vars_2nd
    
    out <- lapply(.object, calculateRhoT,
                  .alpha               = .alpha,
                  .closed_form_ci      = .closed_form_ci,
                  .only_common_factors = .only_common_factors,
                  .output_type         = .output_type,
                  .weighted            = .weighted
    )
    if(.output_type == "vector") {
      out$Second_stage <- out$Second_stage[c_names2]
      out <- if(is.na(out$Second_stage)) {
        out$First_stage
      } else {
        c(out$First_stage, out$Second_stage)
      }
    } else {
      out$Second_stage <- out$Second_stage[out$Second_stage$Construct %in% c_names2, ]
      out <- if(nrow(out$Second_stage) == 0) {
        out$First_stage
      } else {
        rbind(out$First_stage, out$Second_stage)
      }
    }

    return(out)
    
  } else {
    stop2(
      "The following error occured in the calculateRhoT() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
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
    
    ## Add to data.frame and remove rownames
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
    rownames(out_df) <- NULL
    return(out_df)
  }
}



#' HTMT
#'
#' Computes either the heterotrait-monotrait ratio of correlations (HTMT) based on 
#' \insertCite{Henseler2015;textual}{cSEM} or its advancement HTMT2 \insertCite{Roemerinprint;textual}{cSEM}.
#' While the HTMT is a consistent estimator for the construct correlation in 
#' case of tau-equivalent measurement models, the HTMT2 is a consistent estimator
#' for congeneric measurement models. In general, they are used to assess discriminant validity.
#' 
#' Computation of the HTMT assumes that all intra-block and inter-block 
#' correlations between indicators are either all-positive or all-negative.
#' A warning is given if this is not the case. If all intra-block or inter-block
#' correlations are negative the absolute HTMT values are returned (`.absolute = TRUE`).
#' 
#' To obtain the 1-alpha%-quantile of the bootstrap distribution for each HTMT 
#' value set `.inference = TRUE`. To choose the type of confidence interval to use
#' to compute the 1-alpha%-quantile, use `.ci`. To control the bootstrap process,
#' arguments `.handle_inadmissibles`, `.R` and `.seed` are available. 
#' 
#' Since the HTMT and the HTMT2 both assume a reflective measurement
#' model only concepts modeled as common factors are considered by default.
#' For concepts modeled as composites the HTMT may be computed by setting
#' `.only_common_factors = FALSE`, however, it is unclear how to
#' interpret values in this case.
#' 
#' @usage calculateHTMT(
#'  .object               = NULL,
#'  .type_htmt            = c('htmt','htmt2'),
#'  .absolute             = TRUE,
#'  .alpha                = 0.05,
#'  .ci                   = c("CI_percentile", "CI_standard_z", "CI_standard_t", 
#'                            "CI_basic", "CI_bc", "CI_bca", "CI_t_interval"),
#'  .handle_inadmissibles = c("drop", "ignore", "replace"),
#'  .inference            = FALSE,
#'  .only_common_factors  = TRUE,
#'  .R                    = 499,
#'  .seed                 = NULL,
#'  ...
#' )
#'
#' @inheritParams csem_arguments
#' @param .alpha A numeric value giving the significance level. 
#'   Defaults to `0.05`.
#' @param .ci A character strings naming the type of confidence interval to use 
#'   to compute the 1-alpha% quantile of the bootstrap HTMT values. For possible 
#'   choices see [infer()]. Ignored
#'   if `.inference = FALSE`. Defaults to "*CI_percentile*".
#' @param ... Ignored.
#'
#' @return A lower tringular matrix of HTMT values. If `.inference = TRUE`
#'   the upper tringular part is the 1-.alpha%-quantile of the HTMT's bootstrap
#'   distribution.
#' 
#' @seealso [assess()], [csem], [cSEMResults]
#' 
#' @references 
#' 
#' \insertAllCited{}
#' 
#' @export

calculateHTMT <- function(
  .object               = NULL,
  .type_htmt            = c('htmt','htmt2'),
  .absolute             = TRUE,
  .alpha                = 0.05,
  .ci                   = c("CI_percentile", "CI_standard_z", "CI_standard_t", 
                            "CI_basic", "CI_bc", "CI_bca", "CI_t_interval"),
  .handle_inadmissibles = c("drop", "ignore", "replace"),
  .inference            = FALSE,
  .only_common_factors  = TRUE,
  .R                    = 499,
  .seed                 = NULL,
  ...
){
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .ci                   <- match.arg(.ci) # allow only one CI
  .type_htmt            <- match.arg(.type_htmt)

  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateHTMT,
                  .type_htmt     = .type_htmt,
                  .absolute = .absolute,
                  .alpha                = .alpha,
                  .handle_inadmissibles = .handle_inadmissibles,
                  .inference            = .inference,
                  .only_common_factors  = .only_common_factors,
                  .R                    = .R,
                  .seed                 = .seed
                  )
    return(out)
  } else if(inherits(.object, "cSEMResults_default")) {
    ## Get relevant quantities
    m <- .object$Information$Model
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    ## Get relevant quantities
    m <- .object$Second_stage$Information$Arguments_original$.model
  } else {
    stop2(
      "The following error occured in the calculateHTMT() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  if(.only_common_factors) {
    cf_names <- names(m$construct_type[m$construct_type == "Common factor"])
    
    ## Return NA if there are not at least 2 common factors
    if(length(cf_names) < 2) {
      return(NULL)
    }
  } else {
    cf_names <- names(m$construct_type)
  }
  
  cf_measurement <- m$measurement[cf_names, colSums(m$measurement[cf_names, ]) != 0, drop = FALSE]
  
  if(inherits(.object, "cSEMResults_default")) {
    
    i_names <- colnames(cf_measurement)
    S       <- .object$Estimates$Indicator_VCV[i_names, i_names]
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2(
      "The following error occured in the calculateHTMT() function:\n",
      "The HTMT is not (yet) implemented for models containing second-order constructs."
    )
  }
  
  ## Warning if S contains negative and positive correlations within a block
  S_signs    <- cf_measurement %*% (sign(S) - diag(nrow(S))) %*% t(cf_measurement)
  S_elements <- cf_measurement %*% (1 - diag(nrow(S))) %*% t(cf_measurement)
  
  if(any(abs(S_signs) != S_elements)) {
    warning(
      "The following warning occured in the calculateHTMT() function:\n",
      "Intra-block and inter-block correlations between indicators", 
      " must be either all-positive or all-negative.", call. = FALSE)
  }
  
  # Extract names of the indicators belnging to one block
  ind_blocks<-lapply(rownames(cf_measurement),function(x){
    colnames(cf_measurement)[cf_measurement[x,]==1]
  })
  names(ind_blocks)<-rownames(cf_measurement)
  
  # Create pairs of indicator blocks 
  block_pairs <- utils::combn(ind_blocks, 2, simplify = FALSE)
  
  # extract monotrait and heterotrait correlations
  correlations <- lapply(block_pairs, function(x){
    monocortemp<-S[x[[1]],x[[1]],drop=FALSE]
    # Set monotrait-heteromethod cor to one for single-indicator constructs 
    if(nrow(monocortemp)==1){
     monocor1=1 
    }else{
      monocor1<-monocortemp[lower.tri(monocortemp)]
    }
    monocortemp<-S[x[[2]],x[[2]],drop=FALSE]
    # Set monotrait-heteromethod cor to one for single-indicator constructs 
    if(nrow(monocortemp)==1){
      monocor2=1
    }else{
      monocor2<-monocortemp[lower.tri(monocortemp)]
    }
    # take always the absolute value of the heterotrait-heteromethod correlations
    hetcor<-abs(c(S[x[[1]],x[[2]]]))
    
    # return correlations as list of vectors containing correlations
    list(monocor1,monocor2,hetcor)
  })
  
  # check if sign of all heterotrait-heteromethod correlations is negative; 
  # if this i the case -1 else 1
  sign_identification=sapply(block_pairs,function(x){
    S_signs_two_blocks=S_signs[names(x)[1],names(x)[2]]
    S_elements_two_blocks=S_elements[names(x)[1],names(x)[2]]
    
    if(S_signs_two_blocks==S_elements_two_blocks){
      1
    }else{
      -1
    }
  })
  
  
  if(.type_htmt=='htmt2'){
    # calculate geometric mean of the correlations
    avg_cor<-lapply(correlations, function(x){
      sapply(x, function(y){
        prod(y)^(1/length(y))
      })
    })
  }
  
  if(.type_htmt=='htmt'){
    # calculate the arithmetic mean of the correlations
    avg_cor<-lapply(correlations, function(x){
      sapply(x, function(y){
        mean(y)
      })
    })
  }
  
  # Compute HTMT
  htmts <- sapply(avg_cor,function(x){
    x[3]/sqrt(x[1]*x[2])
  })
  
  # 
  htmts <- htmts * sign_identification 
  
  # Sort HTMT values in matrix
  out<-matrix(0,
         nrow=length(names(ind_blocks)),
         ncol=length(names(ind_blocks)),
         dimnames=list(names(ind_blocks),names(ind_blocks)))
  out[lower.tri(out)]<-htmts
  
  
  # ## Average correlation of the indicators within and across blocks
  # ## The eta_i - eta_i (main diagonal) element is the monotrait-heteromethod correlation
  # ## The eta_i - eta_j (off-diagonal) element is the heterotrait-heteromethod correlation
  # avrg_cor <- cf_measurement %*% (S - diag(diag(S))) %*% t(cf_measurement) / S_elements
  # 
  # # Single-indicator constructs monotrait-heteromethod correlation is set to 1
  # x <- rowSums(cf_measurement) == 1
  # 
  # if(sum(x) == 1) {
  #   avrg_cor[x, x] <- 1
  # } else if(sum(x) > 1) {
  #   diag(avrg_cor[x, x]) <- 1 
  # } # else: dont do anything
  # 
  # ## Compute HTMT
  # # HTMT_ij = Average heterotrait-heteromethod correlation between i and j divided by 
  # # the geometric means of the average monotrait-heteromethod correlation of 
  # # eta_i with the average monotrait-heteromethod correlation of construct eta_j
  # # (can be negative if some indicators are negatively correlated)
  # tryCatch({sqrt(diag(avrg_cor) %o% diag(avrg_cor))},
  #          warning = function(w) {
  #            warning(
  #              "The following warning occured in the calculateHTMT() function:\n",
  #              "The geometric mean of the average monotrait-heteromethod",
  #              " correlation of at least one construct with",
  #              " the average monotrait-heteromethod correlation of the",
  #              " other constructs is negative. NaNs produced.",
  #              call. = FALSE)
  #          }
  # )
  # out <- avrg_cor*lower.tri(avrg_cor) / suppressWarnings(sqrt(diag(avrg_cor) %o% diag(avrg_cor))) 
  
  if(.absolute) {
    out <- abs(out)
  }
  
  if(.inference) {
    # Bootstrap if necessary
    out_resample <- resamplecSEMResults(
      .object, 
      .user_funs = list("HTMT" = calculateHTMT), 
      .type_htmt = .type_htmt,
      .absolute = .absolute,
      .handle_inadmissibles = .handle_inadmissibles,
      .inference = FALSE,
      .only_common_factors = .only_common_factors,
      .force = TRUE, # to force computation even if .object already contains resamples
      .R = .R,
      .seed = .seed
    )
    
    # Compute quantile
    if(length(.alpha) == 1) {
      out_infer <- infer(out_resample, .alpha = .alpha*2, .quantity = .ci)
      quants <- out_infer$HTMT[[1]][2, ] 
    } else {
      stop2(
        "The following error occured in the calculateHTMT() function:\n",
        "Only a single numeric probability accepted. You provided:", paste(.alpha, sep = ", "))
    }
    
    ## Reassemble matrix
    htmt_quantiles <- out
    htmt_quantiles[] <- quants
    
    htmt_inference <- out + t(htmt_quantiles)
    
    # Return
    diag(htmt_inference) <- 1
    return(htmt_inference)
  }
  # Return
  diag(out) <- 1
  out
}


#' Calculate difference between S and Sigma_hat
#'
#' Calculate the difference between the empirical (S) 
#' and the model-implied indicator variance-covariance matrix (Sigma_hat)
#' using different distance measures.
#' 
#' The distances may also be computed for any two matrices A and B by supplying 
#' A and B directly via the `.matrix1` and `.matrix2` arguments. 
#' If A and B are supplied `.object` is ignored.
#' 
#' @return A single numeric value giving the distance between two matrices.
#' 
#' @inheritParams csem_arguments
#' @param ... Ignored.
#'
#' @name distance_measures
NULL

#' @describeIn distance_measures The geodesic distance (dG).
#' @export

calculateDG <- function(
  .object    = NULL, 
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = FALSE,
  ...
){
  
  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    if(inherits(.object, "cSEMResults_multi")) {
      out <- lapply(.object, calculateDG, .saturated = .saturated)
      return(out)
    }
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateDG() function:\n",
        "`.object` must be of class `cSEMResults`."
      )
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
#' @export

calculateDL <- function(
  .object    = NULL, 
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = FALSE,
  ...
){
  
  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    if(inherits(.object, "cSEMResults_multi")) {
      out <- lapply(.object, calculateDL, .saturated = .saturated)
      return(out)
    }
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateDL() function:\n",
        "`.object` must be of class `cSEMResults`."
      )
    }
    
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')
  }
  
  ## Calculate distance
  0.5 * sum((S - Sigma_hat)^2)
}

#' @describeIn  distance_measures The distance measure (fit function) used by ML
#' @export

calculateDML <- function(
  .object    = NULL, 
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = FALSE,
  ...
){
  
  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    if(inherits(.object, "cSEMResults_multi")) {
      out <- lapply(.object, calculateDML, .saturated = .saturated)
      return(out)
    }
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateDML() function:\n",
        "`.object` must be of class `cSEMResults`."
      )
    }
    
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator')  
  }
  
  p <- dim(S)[1]
  
  # This is the distance function. The test statistic is T_ML = (n-1) or n * DML! 
  sum(diag(S %*% solve(Sigma_hat))) - log(det(S%*%solve(Sigma_hat))) - p
}

#' Model fit measures
#' 
#' Calculate fit measures.
#' 
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#fit_indices}{Fit indices}
#' section of the \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website}
#' for details on the implementation.
#' 
#' @return A single numeric value.
#' 
#' @inheritParams csem_arguments
#' @param ... Ignored.
#' 
#' @name fit_measures 
NULL

#' @describeIn fit_measures The chi square statistic.
#' @export

calculateChiSquare <- function(.object) {
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateChiSquare)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
  } else {
    stop2(
      "The following error occured in the calculateChiSquare() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  F0   <- calculateDML(.object)
  
  (n - 1) * F0
}

#' @describeIn fit_measures The ChiSquare statistic divided by its degrees of freedom.
#' @export

calculateChiSquareDf <- function(.object) {
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateChiSquareDf)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
  } else {
    stop2(
      "The following error occured in the calculateChiSquareDf() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  F0   <- calculateDML(.object)
  
  ((n - 1) * F0) / calculateDf(.object)
}

#' @describeIn fit_measures The comparative fit index (CFI).
#' @export

calculateCFI <- function(.object) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateCFI)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
    S    <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
    S    <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateCFI() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
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
#' @export

calculateGFI <- function(.object, .type_gfi = c("ML", "GLS", "ULS"), ...) {
  .type_gfi <- match.arg(.type_gfi)

  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateGFI, .type_gfi = .type_gfi)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    S    <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    S    <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateGFI() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  Sigma_hat <- fit(.object)
  
  if(.type_gfi == "ML") {
    # If ML; See Mulaik (1989, p. 345)
    1 - matrixcalc::matrix.trace(t(solve(Sigma_hat) %*% S - diag(nrow(S))) %*% 
                                   (solve(Sigma_hat) %*% S - diag(nrow(S)))) / 
      matrixcalc::matrix.trace(t(solve(Sigma_hat) %*% S) %*% (solve(Sigma_hat) %*% S))
  } else if(.type_gfi == "GLS") {
    # If GLS; See Tanaka & Huba (1985, p. 199; Eq. 16)
    1 - matrixcalc::matrix.trace(t(diag(nrow(S)) - Sigma_hat %*% solve(S)) %*% 
                                   (diag(nrow(S)) - Sigma_hat %*% solve(S))) / nrow(S)
  } else if(.type_gfi == "ULS") {
    # If ULS; See Mulaik (1989, p. 345)
    1 - matrixcalc::matrix.trace(t(S - Sigma_hat) %*% (S - Sigma_hat)) / 
      matrixcalc::matrix.trace(t(S) %*% S)
  }
}

#' @describeIn fit_measures The Hoelter index alias Hoelter's (critical) N (CN).
#' @export

calculateCN <- function(.object, .alpha = 0.05, ...) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateCN)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    N    <- nrow(.object$Information$Data)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    N    <- nrow(.object$First_stage$Information$Data)
  } else {
    stop2(
      "The following error occured in the calculateHoeltersN() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  chi_square      <- calculateChiSquare(.object)
  df              <- calculateDf(.object)
  chi_square_crit <- qchisq(1 - .alpha, df)
  # z <- qnorm(1 - .alpha/2ap)                              
  
  # (z + sqrt(2*df - 1))^2 / (2*chi_square/(N-1)) + 1 
  # Formula given in Hoelters (1983), p.331 and Hu & Bentler (1998), p.428; the
  # formula is less exact than the one below. See Bollen & Liang (1988) - Some properties
  # of Hoelter's CN; especially footnote 2

  chi_square_crit / (chi_square/(N-1)) + 1 # formula given on David Kennys website and Bollen & Liang (1988)
  # Motivation: chi_crit = (CN - 1)*F = (CN -1)*Chi_square/(N-1)
}

#' @describeIn fit_measures The incremental fit index (IFI).
#' @export

calculateIFI <- function(.object) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateIFI)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
    S    <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
    S    <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateIFI() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  p  <- dim(S)[1]
  df <- calculateDf(.object)
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  ((n-1)*F0 - (n-1)*FT) / ((n-1)*F0 - df)
}

#' @describeIn fit_measures The normed fit index (NFI).
#' @export

calculateNFI <- function(.object) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateNFI)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    S    <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    S    <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateNFI() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  p <- dim(S)[1]
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  
  (F0 - FT) / F0
}

#' @describeIn fit_measures The non-normed fit index (NNFI). Also called the Tucker-Lewis index (TLI).
#' @export

calculateNNFI <- function(.object) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateNNFI)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
    S    <- .object$Estimates$Indicator_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
    S    <- .object$First_stage$Estimates$Indicator_VCV
  } else {
    stop2(
      "The following error occured in the calculateNNFI() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  p    <- dim(S)[1]
  df_T <- calculateDf(.object)
  df_0 <- calculateDf(.object, .null_model = TRUE)
  
  F0 <- log(det(diag(nrow(S)))) + 
    sum(diag(S %*% solve(diag(nrow(S))))) - log(det(S)) - p
  
  FT <- calculateDML(.object)
  # Note: "If the index is greater than one, it is set at one" 
  # Source: http://www.davidakenny.net/cm/fit.htm
  
  min((F0/df_0 - FT/df_T) / (F0/df_0 - 1/(n-1)), 1)
}

#' @describeIn fit_measures The root mean square error of approximation (RMSEA).
#' @export

calculateRMSEA <- function(.object) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateRMSEA)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    n    <- nrow(.object$Information$Data)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    n    <- nrow(.object$First_stage$Information$Data)
  } else {
    stop2(
      "The following error occured in the calculateRMSEA() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  
  df <- calculateDf(.object)
  
  F0 <- max(calculateDML(.object) - calculateDf(.object)/(n - 1), 0)
  
  sqrt(F0 / df) # RMSEA
}

#' @describeIn fit_measures The root mean squared residual covariance matrix of the outer model residuals (RMS theta).
#' @export

calculateRMSTheta <- function(
  .object
) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateRMSTheta)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    S      <- .object$Estimates$Indicator_VCV
    W      <- .object$Estimates$Weight_estimates
    Lambda <- .object$Estimates$Loading_estimates
    P      <- .object$Estimates$Construct_VCV
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Not yet implemented.")
  } else {
    stop2(
      "The following error occured in the calculateRMSTheta() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
  }
  # This is the "non-model-implied" error correlation matrix 
  Theta <- S - S %*% t(W) %*% Lambda - t(S %*% t(W) %*% Lambda) + t(Lambda) %*% P %*% Lambda
  
  ## For compsites, within block indicator correlations should be excluded as 
  ## they are allowed to freely covary.
  if(inherits(.object, "cSEMResults_default")) {
    comp <- which(.object$Information$Model$construct_type == "Composite")
    
    for(i in comp) {
      indi <- which(.object$Information$Model$measurement[i, ] == 1)
      Theta[indi, indi] <- NA
    }
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2("Not yet implemented.")
  }
  
  sqrt(mean(Theta[lower.tri(Theta)]^2, na.rm = TRUE))
}

#' @describeIn fit_measures The standardized root mean square residual (SRMR).
#' @export

calculateSRMR <- function(
  .object    = NULL, 
  .matrix1   = NULL,
  .matrix2   = NULL,
  .saturated = FALSE,
  ...
) {
  
  if(!is.null(.matrix1) & !is.null(.matrix2)) {
    S         <- .matrix1
    Sigma_hat <- .matrix2
  } else {
    if(inherits(.object, "cSEMResults_multi")) {
      out <- lapply(.object, calculateSRMR, .saturated = .saturated)
      return(out)
    }
    if(inherits(.object, "cSEMResults_default")) {
      S <- .object$Estimates$Indicator_VCV
    } else if(inherits(.object, "cSEMResults_2ndorder")) {
      S <- .object$First_stage$Estimates$Indicator_VCV
    } else {
      stop2(
        "The following error occured in the calculateSRMR() function:\n",
        "`.object` must be of class `cSEMResults`."
      )
    }
    
    # The SRMR as calculated by us is always based on the the difference 
    # between correlation matrices.
    Sigma_hat <- fit(.object, .saturated = .saturated, .type_vcv = 'indicator') 
  }
  
  # Perhaps in the future we allow to estimate unstandardized coefficients
  C_diff    <- cov2cor(S) -  cov2cor(Sigma_hat)
  
  sqrt(sum(C_diff[lower.tri(C_diff, diag = T)]^2) / sum(lower.tri(C_diff, diag = T)))
} 




#' Calculate Cohens f^2
#'
#' Calculate the effect size for regression analysis \insertCite{Cohen1992}{cSEM}
#' known as Cohen's f^2. 
#'
#' @usage calculatef2(.object = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @return A matrix with as many rows as there are structural equations. The 
#'   number of columns is equal to the total number of right-hand side variables
#'   of these equations.
#' 
#' @references 
#' \insertAllCited{}
#' @seealso [assess()], [csem], [cSEMResults]
#' @export

calculatef2 <- function(.object = NULL) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculatef2)
    return(out)
    
  } else if(inherits(.object, "cSEMResults_default")) {
    info <- .object
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    info <-  .object$Second_stage
    
  } else {
    stop2(
      "The following error occured in the calculatef2() function:\n",
      "`.object` must be a `cSEMResults` object."
    )
  }
  
  ## Get relevant quantities
  approach_nl      <- info$Information$Arguments$.approach_nl
  approach_paths   <- info$Information$Arguments$.approach_paths
  approach_weights <- info$Information$Arguments$.approach_weights
  csem_model       <- info$Information$Model
  H         <- info$Estimates$Construct_scores
  normality <- info$Information$Arguments$.normality
  P         <- info$Estimates$Construct_VCV
  Q         <- sqrt(info$Estimates$Reliabilities)
  
  s <- csem_model$structural
  
  ## The R2 and the VIF for the 2SLS approach are not implemented yet, hence,
  ## the f2 statistic cannot be calculated 
  if(approach_paths != "OLS") {
    stop2(
      "The following error occured in the calculatef2() function:\n",
      "Calculation of the effect size (f2) only implemented for .approach_path = 'OLS'."
    )
  }
  
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
      
      # Obtain the R^2 included
      r2_included <- info$Estimates$R2[x]
      
      effect_size <- unname((r2_included - r2_excluded)/(1 - r2_included))
    })
    inner_out <- unlist(inner_out)
    names(inner_out) <- indep_vars
    inner_out 
  })
  names(outer_out) <- vars_endo
  
  ## Make output a matrix
  # Note: this is necessary for calculatef2() to work
  #       when supplied to the .user_funs argument. Currently, .user_funs functions 
  #       need to return a vector or a matrix. I may change that in the future.
  ss <- s[vars_endo, , drop = FALSE]
  tm <- t(ss)
  tm[which(tm == 1)] <- unlist(outer_out)
  
  # Remove "_temp" suffix if it appears
  rownames(tm) <- gsub("_temp", "", rownames(tm))
  colnames(tm) <- gsub("_temp", "", colnames(tm))
  
  # Return
  t(tm)
}



#' Calculate variance inflation factors (VIF) for weights obtained by PLS Mode B
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
#' @export

calculateVIFModeB <- function(.object = NULL) {
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, calculateVIFModeB)
    return(out)
  }
  if(inherits(.object, "cSEMResults_default")) {
    # continue
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    out <- lapply(.object, calculateVIFModeB)
    return(out)
  } else {
    stop2(
      "The following error occured in the calculateVIFModeB() function:\n",
      "`.object` must be of class `cSEMResults`."
    )
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
          R2   <- cor(Xk %*% beta, y)^2
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
  
  ## Make output a matrix
  # Note: this is necessary for calculateVIFModeB() to work
  #       when supplied to the .user_funs argument of resamplecSEMResults(). 
  #       Currently, the .user_funs functions 
  #       need to return a vector or a matrix.
  
  if(!anyNA(VIF)) {
    mm <- m[names(modesB), colSums(m[names(modesB), , drop = FALSE]) != 0 , drop = FALSE]
    tm <- t(mm)
    tm[which(tm == 1)] <- unlist(VIF)
    
    # Remove "_temp" suffix if it appears
    rownames(tm) <- gsub("_temp", "", rownames(tm))
    colnames(tm) <- gsub("_temp", "", colnames(tm))
    
    # Return
    VIF <- t(tm)
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
