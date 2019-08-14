### 
# Date: 14.08.2019
# This is the fit function when i was trying to make it work with the repeated
# indicators approach. Decided to delete the appraoch. unless we decided to 
# put it back in, this will not be required anymore.

fit.cSEMResults_default <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv
) {
  
  ### For maintenance: ---------------------------------------------------------
  ## Cons_exo  := (J_exo x 1) vector of exogenous constructs names.
  ## Cons_endo := (J_endo x 1) vector of endogenous constructs names.
  ## S         := (K x K) Empirical indicator VCV matrix: V(X).
  ## B         := (J_endo x J_endo) matrix of (estimated) path coefficients 
  ##              from endogenous to endogenous constructs. (zero if there is no path)
  ## Gamma     := (J_endo x J_exo) matrix of (estimated) path coefficients from
  ##              exogenous to endogenous constructs.
  ## Lambda    := (J X K) matrix of factor (dissatenuated if requested) 
  ##              and/or composite loadings.
  ## Phi       := (J_exo x J_exo) empirical construct correlation matrix 
  ##              between exogenous constructs (attenuated if requested).
  ## I         := (J_endo x J_endo) identity matrix.
  ## Theta     := (K x K) diagonal matrix of measurement model error variances.
  ## Psi       := (J_endo x J_endo) diagonal matrix of structural model error 
  ##              variances (zetas).
  ## Corr_exo_endo := (J_exo x J_endo) model-implied correlation matrix between 
  ##                  exogenous and endogenous constructs.
  ## Corr_endo     := (J_endo x J_endo)  model-implied correlation matrix between
  ##                  endogenous constructs.
  
  ### Preparation ==============================================================
  ## Check if linear
  if(.object$Information$Model$model_type != "Linear"){
    stop2("`fit()` currently not applicable to nonlinear models.")
  }
  
  mod       <- .object$Information$Model
  S         <- .object$Estimates$Indicator
  Lambda    <- .object$Estimates$Loading_estimates
  
  # Prune S, Lambda and Theta if there are repeated indicators.
  # Its important to have the if clause here because fit.cSEMResults_default
  # is also called on the second stage of the 2stage/mixed approach. There the
  # rows for "cons_2nd" must not be deleted!
  # if(.object$Information$Approach_2ndorder %in% c("RI_original", "RI_extended")) {
  
  # # Which constructs are second orders (NULL for models containing no 2nd orders)
  # cons_2nd      <- .object$Information$Model_original$vars_2nd
  # # Which constructs are constructs attachted to a second order (atto2nd)
  # cons_atto2nd  <- .object$Information$Model_original$vars_attached_to_2nd
  # # Which constructs are constructs not attached to a second order (natto2nd)
  # cons_natto2nd <- .object$Information$Model_original$vars_not_attached_to_2nd
  # 
  # # Select only columns/rows that are not repeated indicators (if there are no
  # # repeated indicators this will simply select all columns of S and Lambda)
  # # Also delete rows in Lambda that are second orders
  # selector <- !grepl("_2nd_", colnames(S))
  # 
  # S        <- S[selector, selector]
  # Lambda   <- Lambda[setdiff(rownames(Lambda), cons_2nd), selector]
  
  ## The model-implied construct VCV of a model containing second order
  ## constructs is:
  #
  #     V(eta) = (         V_atto2nd                    lambda_2nd %*% V_2nd_to_natto2nd')
  #              ((lambda_2nd %*% V_2nd_to_natto2nd)'           V_natto2nd               )
  #
  #  2nd      := second order 
  #  atto2nd  := attached to second order construct
  #  natto2nd := not attachted to second order construct
  
  #  V_atto2nd  := VCV of cons_atto2nd
  #  lambda_2nd := correlation ("loadings") between cons_atto2n and cons_2nd
  #  V_2nd_to_natto2nd := correlation between cons_2nd and cons_natto2nd
  #  V_natto2nd := VCV of cons_natto2nd
  
  # ## Extract lambda_2nd:
  # lambda_2nd <- .object$Estimates$Construct_VCV[cons_atto2nd, cons_2nd, drop = FALSE]
  # 
  # ## Compute V_atto2nd
  # V_atto2nd  <- .object$Estimates$Construct_VCV[cons_atto2nd, cons_atto2nd, drop = FALSE] 
  # # fallunterscheidung composite common factor; measurement errors
  # 
  # ## Construct names required to estimate V_2nd_to_natto2nd and V_natto2nd
  # m         <- mod$structural
  # Cons_endo <- setdiff(rownames(m)[rowSums(m) != 0], cons_atto2nd)
  # Cons_exo  <- setdiff(colnames(m), c(Cons_endo, cons_atto2nd))
  # } else {
  
  m         <- mod$structural
  Cons_endo <- rownames(m)[rowSums(m) != 0]
  Cons_exo  <- setdiff(colnames(m), Cons_endo)
  Theta     <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(Theta) <- dimnames(S)
  # }
  
  ## Check if recursive, otherwise return a warning
  if(any(m[Cons_endo, Cons_endo] + t(m[Cons_endo, Cons_endo]) == 2)){
    warning2("`fit()` currently not applicable to non-recursive models.",
             " The model-implied indicator covariance matrix is likely to be wrong.")
  }
  
  ### VCV of the constructs ====================================================
  if(.saturated) {
    # If a saturated model is assumed the structural model is ignored in
    # the calculation of the construct VCV (i.e. a full graph is estimated). 
    # Hence: V(eta) = WSW' (ppssibly disattenuated)
    vcv_construct <- .object$Estimates$Construct_VCV
    
  } else {
    
    B      <- .object$Estimates$Path_estimates[Cons_endo, Cons_endo, drop = FALSE]
    Gamma  <- .object$Estimates$Path_estimates[Cons_endo, Cons_exo, drop = FALSE]
    Phi    <- .object$Estimates$Construct_VCV[Cons_exo, Cons_exo, drop = FALSE]
    I      <- diag(length(Cons_endo))
    
    ## Calculate variance of the zetas
    # Note: this is not yet fully correct, athough it does not currently affect 
    # the results. This may have to be fixed in the future to avoid potential 
    # problems that might arise in setups we have not considered yet.
    vec_zeta <- 1 - rowSums(.object$Estimates$Path_estimates * 
                              .object$Estimates$Construct_VCV)
    names(vec_zeta) <- rownames(.object$Estimates$Construct_VCV)
    
    vcv_zeta <- matrix(0, nrow = nrow(I), ncol = ncol(I))
    diag(vcv_zeta) <- vec_zeta[Cons_endo]
    
    ## Correlations between exogenous and endogenous constructs
    Corr_exo_endo <- Phi %*% t(Gamma) %*% t(solve(I-B))
    ## Correlations between endogenous constructs 
    Cor_endo <- solve(I-B) %*% (Gamma %*% Phi %*% t(Gamma) + vcv_zeta) %*% t(solve(I-B))
    diag(Cor_endo) <- 1
    
    vcv_construct <- rbind(
      cbind(Phi, Corr_exo_endo),
      cbind(t(Corr_exo_endo), Cor_endo)
    ) 
    ## Make symmetric
    vcv_construct[lower.tri(vcv_construct)] <- t(vcv_construct)[lower.tri(vcv_construct)]
  }
  
  ## If repeated indicators appraoch: assemble construct vcv matrix
  # if(.object$Information$Approach_2ndorder %in% c("RI_original", "RI_extended")) { 
  #   # V_2nd_to_natto2nd := correlation between cons_2nd and cons_natto2nd
  #   V_2nd_to_natto2nd <- vcv_construct[cons_natto2nd, cons_2nd]
  #   # V_natto2nd := VCV of cons_natto2nd
  #   V_natto2nd <- vcv_construct[cons_natto2nd, cons_natto2nd]
  #   
  #   vcv_construct <- rbind(
  #     cbind(V_atto2nd, lambda_2nd %*% t(V_2nd_to_natto2nd)),
  #     cbind(t(lambda_2nd %*% t(V_2nd_to_natto2nd)), V_natto2nd)
  #   )
  # 
  #   ## Reoder to match Lambda
  #   selector2 <- setdiff(colnames(.object$Estimates$Construct_VCV), cons_2nd)
  #   vcv_construct <- vcv_construct[selector2, selector2]
  # 
  #   Theta           <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  #   dimnames(Theta) <- dimnames(S)
  # }
  
  ## If only the fitted construct VCV is needed, return it now
  if(.type_vcv == "construct") {
    return(vcv_construct)
  }
  
  ## Calculate model-implied VCV of the indicators
  vcv_ind <- t(Lambda) %*% vcv_construct %*% Lambda
  
  Sigma <- vcv_ind + Theta
  
  ## Make symmetric
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  ## Replace indicators connected to a composite by their correponding elements of S.
  composites <- names(mod$construct_type[mod$construct_type == "Composite"])
  index  <- t(mod$measurement[composites, , drop = FALSE]) %*% mod$measurement[composites, , drop = FALSE]
  
  # if(.object$Information$Approach_2ndorder %in% c("RI_original", "RI_extended")) { 
  #   index <- index[selector, selector]
  #   mod$error_cor <- mod$error_cor[selector, selector]
  # }
  
  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  # Replace indicators whose measurement errors are allowed to be correlated by s_ij
  Sigma[mod$error_cor == 1] = S[mod$error_cor == 1]
  
  return(Sigma)
}