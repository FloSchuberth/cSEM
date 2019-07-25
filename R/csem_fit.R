#' Model-implied indicator and construct variance-covariance matrix
#'
#' Calculate the model-implied indicator or construct variance-covariance (VCV) 
#' matrix. Currently only the model-implied VCV for recursive linear models 
#' is implemented (including models containing second order constructs).
#' 
#' Notation is taken from \insertCite{Bollen1989;textual}{cSEM}.
#' If `.saturated = TRUE` the model-implied variance-covariance matrix is calculated 
#' for a saturated structural model (i.e., the VCV of the constructs is replaced 
#' by their correlation matrix). Hence: V(eta) = WSW' (possibly disattenuated).
#'
#' @usage fit(
#'   .object    = NULL, 
#'   .saturated = args_default()$.saturated,
#'   .type_vcv  = args_default()$.type_vcv
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return Either a (K x K) matrix or a (J x J) matrix depending on the `*type_vcv*`.
#' 
#' @examples 
#' \dontrun{
#' res <- csem(.data = dataset, .model = model)
#' 
#' fit(.object = res, .saturated = FALSE, .type_vcv = "indicator")
#' }
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export

fit <- function(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv
  ) {
  UseMethod("fit")
}

#' @describeIn fit (TODO)
#' @export

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
  
  ## Collect relevant 
  m         <- .object$Information$Model$structural
  Cons_endo <- rownames(m)[rowSums(m) != 0]
  Cons_exo  <- setdiff(colnames(m), Cons_endo)
  S      <- .object$Estimates$Indicator_VCV
  Lambda <- .object$Estimates$Loading_estimates
  Theta  <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(Theta) <- dimnames(S)
  
  ## Check if recursive, otherwise return a warning
  if(any(m[Cons_endo, Cons_endo] + t(m[Cons_endo, Cons_endo]) == 2)){
    warning2("`fit()` currently not applicable to non-recursive models.",
             " The model-implied indicator covariance matrix is likely to be wrong.")
  }
  
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
    
    vcv_construct <- rbind(cbind(Phi, Corr_exo_endo),
                           cbind(t(Corr_exo_endo), Cor_endo)) 
    ## Make symmetric
    vcv_construct[lower.tri(vcv_construct)] <- t(vcv_construct)[lower.tri(vcv_construct)]
  }
  
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
  
  mod <- .object$Information$Model
  composites <- names(mod$construct_type[mod$construct_type == "Composite"])
  index <- t(mod$measurement[composites, , drop = FALSE]) %*% mod$measurement[composites, , drop = FALSE]
  
  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  # Replace indicators whose measurement errors are allowed to be correlated by s_ij
  Sigma[.object$Information$Model$error_cor == 1] = S[.object$Information$Model$error_cor == 1]
  
  return(Sigma)
}

#' @describeIn fit (TODO)
#' @export

fit.cSEMResults_multi <- function(
  .object    = NULL,
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv
  ) {
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, fit.cSEMResults_2ndorder, 
           .saturated = .saturated,
           .type_vcv  = .type_vcv)
  } else {
    lapply(.object, fit.cSEMResults_default, 
           .saturated = .saturated,
           .type_vcv  = .type_vcv)
  }
  
}

#' @describeIn fit (TODO)
#' @export

fit.cSEMResults_2ndorder <- function(
  .object    = NULL,
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv
  ) {
  
  ## Get relevant quantities
  S             <- .object$First_stage$Estimates$Indicator_VCV
  vcv_construct <- fit.cSEMResults_default(.object$Second_stage, 
                                           .saturated = .saturated,
                                           .type_vcv  = .type_vcv)
  Lambda        <- .object$First_stage$Estimates$Loading_estimates
  Theta         <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  
  # Reorder dimnames to match the order of Lambda and ensure symmetrie
  vcv_construct <- vcv_construct[rownames(Lambda), rownames(Lambda)]
  
  ## If only the fitted construct VCV is needed, return it now
  if(.type_vcv == "construct") {
    return(vcv_construct)
  }
  
  # Compute VCV and ensure symmetrie
  Sigma <- t(Lambda) %*% vcv_construct %*% Lambda + Theta
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  # Replace composite blocks by corresponding elements of S
  m          <- .object$First_stage$Information$Model
  composites <- names(m$construct_type[m$construct_type == "Composite"])
  index      <- t(m$measurement[composites, , drop = FALSE]) %*% m$measurement[composites, , drop = FALSE]
  
  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  # Replace indicators whose measurement errors are allowed to be correlated by s_ij
  Sigma[.object$Information$Model$error_cor == 1] = S[.object$Information$Model$error_cor == 1]
  Sigma
  
  return(Sigma)
}
  