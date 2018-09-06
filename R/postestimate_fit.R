#' Model-implied indicator covariance matrix
#'
#' Calculate the model-implied indicator variance-covariance (VCV) matrix. 
#' Currently only the model-implied VCV for linear model is implemented.
#' 
#' Notation is taken from \insertCite{Bollen1989;textual}{cSEM}.
#' By default the model-implied VCV matrix is based on the structural model. 
#' If `.saturated = TRUE` the structural model is ignored (i.e. a full graph is estimated). 
#' Hence: V(eta) = WSW' (possibly disattenuated).
#'
#' @usage fit(
#'   .object    = args_default()$.object, 
#'   .saturated = args_default()$.saturated
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
fit <- function(
  .object    = args_default()$.object, 
  .saturated = args_default()$.saturated
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
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  ## Check if linear
  if(.object$Information$Model$model_type != "Linear"){
    stop("`fit()` currently not applicable to nonlinear models.",
         call. = FALSE)
  }
  
  ## Collect matrices
  S      <- .object$Estimates$Indicator_VCV
  Lambda <- .object$Estimates$Loading_estimates
  Theta  <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(Theta) <- dimnames(S)
  
  if(.saturated) {
    # If a saturated model is assumed the structural model is ignored in
    # the calculation of the construct VCV (i.e. a full graph is estimated). 
    # Hence: V(eta) = WSW' (ppssibly disattenuated)
    vcv_construct <- .object$Estimates$Construct_VCV
    
  } else {
    
    Cons_exo  <- .object$Information$Model$vars_exo
    Cons_endo <- .object$Information$Model$vars_endo
    
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
    
    vcv_construct <- rbind(cbind(Phi, Corr_exo_endo),
                           cbind(t(Corr_exo_endo), Cor_endo)) 
  }
  
  ## Calculate model-implied VCV of the indicators
  vcv_ind <- t(Lambda) %*% vcv_construct %*% Lambda
  
  Sigma <- vcv_ind + Theta
  
  ## Make symmetric
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  ## Replace indicators connected to a composite by their correponding elements of S.
  
  mod <- .object$Information$Model
  composites <- names(mod$construct_type[mod$construct_type == "Composite"])
  index <- t(mod$measurement[composites, , drop = FALSE]) %*% mod$measurement[composites, , drop = FALSE]
  
  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  # Replace indicators whose measurement errors are allowed to be correlated by s_ij
  Sigma[.object$Information$Model$error_cor == 1] = S[.object$Information$Model$error_cor == 1]
  
  return(Sigma)
}
