#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()] prints basic 
#' information on a [cSEMResults] object including potentially
#' useful postestimation functions to call.
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#' 
print.cSEMResults <- function(.object) {
  
  cat(cli::rule(line = "bar2"), "\n")
  cat(cli::rule(center = "Overview"), "\n\n")
  # if(length(.object) > 1) {
  #   cat("Estimation status:\n\n\t")
  #       for(i in names(.object)) {
  #         cat(crayon::col_align(i, 10), ": ", 
  #             ifelse(.object[[i]]$Information$Weight_info$Convergence_status == TRUE, "successful\n\t", "not successful\n\t"), sep = "")
  #       }
  # } else {
  cat("Estimation was successful. The result is a named list of class " %+% bold("cSEMResults") %+%" \n",
      "containing the following list elements (:\n\n\t",
      "- ", crayon::green("Estimates\n\t"),
      "- ", crayon::green("Information\n\n"), sep = "")
  cat("To get an overview or help type: \n\n\t",
      "- ", crayon::magenta("?"), crayon::cyan("cSEMResults"),"\n\t",
      "- ", crayon::magenta("str"), "(", crayon::cyan("<listname>"), ")\n\t",
      "- ", crayon::magenta("summary"), "(", crayon::cyan("<listname>"), ") or \n\t",
      "- ", crayon::magenta("listviewer"), crayon::yellow("::"), crayon::magenta("jsondedit"),
      "(", crayon::cyan("<listname>"), ", ", crayon::red("mode"), " = ", crayon::cyan("'view'"), ")\n\n", sep = "")
  cat("If you wish to access the list elements directly type e.g. \n\t",
      "- ", crayon::cyan("<listname>"), crayon::yellow("$"), crayon::green("Estimates"), "\n", sep = "")
  cat(cli::rule(line = "bar2"), "\n")
  # }
}

#' `cSEMResultssummary` method for `print()`
#'
#' The [cSEMResultssummary] method for the generic function [print()]. 
#' 
#' Prints a summary of the results obtained from runinng [csem], [cca], or
#' [foreman].
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
print.cSEMResultssummary <- function(x, ...) {
  cat(cli::rule(line = "bar2"), "\n",
      cli::rule(center = "General Information"), "\n\n\t", sep = "")
  
  cat(crayon::col_align("Number of Observations", 25), "= ", x$Number_of_observations, "\n\t",
      crayon::col_align("Weight estimator", 25), "= ", x$Weight_estimator, "\n\t", sep = "")
  if(x$Weight_estimator == "PLS") {
    cat(crayon::col_align("Inner Weighting Scheme ", 25), "= ", x$PLS_weight_scheme_inner, "\n\t", sep = "")
  }
  cat(
    crayon::col_align("Path estimator", 25), "= ", x$Path_estimator, "\n\t",
    crayon::col_align("Convergence Status", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Overall model fit", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Degrees of Freedom", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Computation Time", 25), "= ", c("not yet implemented"), "\n\n\t",
    sep = "")
  
  cat("Construct Types:\n\t","----------------","\n\t", sep = "")
  
  for(i in seq_along(x$Construct_types$Name)) {
    cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$Construct_types$Type[i],"\n\t", sep = "")
  }
  cat("\n")
  
  if(x$Weight_estimator == "PLS") {
    cat("\tPLS Modes:\n\t","----------------","\n\t", sep = "")
    
    for(i in seq_along(x$Construct_type$Name)) {
      cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$PLS_modes[i],"\n\t", sep = "")
    }
    cat("\n")
  }
  
  cat(cli::rule(center = "Estimates"), "\n\n", sep = "")
  
  cat("Estimated Path Coefficients:\n============================\n", sep = "")
  print(x$Path_estimates, row.names = FALSE)
  
  
  cat("\nEstimated Loadings:\n===================\n", sep = "")
  print(x$Loading_estimates, row.names = FALSE)
  
  cat("\nEstimated Weights:\n==================\n", sep = "")
  print(x$Weight_estimates, row.names = FALSE)
  
  if(x$Weight_estimator == "PLS") {
    cat("\nEstimated Correction Factors:\n=============================\n", sep = "")
    print(x$Correction_factors)
  }
  
  cat("\n\n", cli::rule(center = "Other output"), "\n\n\t", sep = "")
  
  cat("<not yet implemented>")
  
  cat("\n\n", cli::rule(center = "Fit Indices"), "\n\n\t", sep = "")
  
  cat(crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\n",
      sep = "")
  
  cat(cli::rule(line = "bar2"))
}

#' `cSEMResults` method for `summary()`
#'
#' The [cSEMResults] method for the generic function [summary()]. 
#' 
#' Computes a summary of the results obtained from runinng [csem], [cca], or
#' [foreman].
#'
#' @usage summary(.object, .what = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#' 
#' @inherit csem_resultssummary return
#'
#' @export
#'
summary.cSEMResults <- function(.object, .what = NULL) {
  
  ## Structure loadings output
  temp <- x$Estimates$Loading_estimates
  names_loadings <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                           " =~ ", colnames(temp))
  loading_estimates <- data.frame("Loading" = names_loadings,
                                  "Estimate" = unlist(t(temp)[t(temp) != 0 ]),
                                  stringsAsFactors = FALSE)
  
  ## Structure weights output
  temp <- x$Estimates$Weight_estimates
  names_weights <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                          " -- ", colnames(temp))
  weight_estimates <- data.frame("Weights" = names_weights,
                                 "Estimate" = unlist(t(temp)[t(temp) != 0 ]),
                                 stringsAsFactors = FALSE)
  
  ## Create summary list
  summary_out <- list(
    "Construct_types"        = x$Meta_information$Construct_types,
    "Correction_factors"     = x$Estimates$Correction_factors,
    "Loading_estimates"      = loading_estimates,
    "Names_endogenous_var"   = x$Meta_information$modellist$vars_endo,
    "Number_of_observations" = x$Meta_information$Number_of_observations,
    "Path_estimates"         = x$Estimates$Path_estimates,
    "Path_estimator"         = x$Meta_information$Path_approach,
    "PLS_modes"              = x$Meta_information$PLS_Modes,
    "PLS_weight_scheme_inner"= x$Meta_information$PLS_Inner_Weightning_scheme,
    "Weight_estimates"       = weight_estimates,
    "Weight_estimator"       = x$Meta_information$Weight_approach
  )
  
  ## Set class for printing and return
  class(summary_out) <- "cSEMResultssummary"
  return(summary_out)
}

#' Model implied indicator variance-covariance matrix
#'
#' Compute the model implied indicator variance-covariance (VCV) matrix. Usually 
#' this is called \eqn{\hat\Sigma}. Currently only the model implied VCV for 
#' linear model for constructs modeled as common factors or composites are 
#' implemented. An error is given if the model is not linear.
#'
#' @usage fitted(.object, test = NULL)
#'
#' @inheritParams csem_arguments
#' @param test not used
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
fitted.cSEMResults <- function(.object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  # Function to compute the model implied VCV matrix of the indicators
  # The implementation is based on an approach proposed by
  # Bentler, P. M., & Weeks, D. G. (1980) - Linear Structural Equations with Latent Variables 
  
  ### For maintenance: ---------------------------------------------------------
  ## S      := (K x K) Empirical indicator VCV matrix: V(X)  
  ## P      := (J x J) Empirical construct covariance/correlation matrix (attenuated) 
  ## B      := (J x J) Matrix of (estimated) path coefficients (zero if there is no path)
  ## Lambda := (J X K) Matrix of factor and/or composite loadings
  ## Lambda_cross := (J x K) Empirical indicator-construct CV matrix (attenuated)
  ##
  ## ---- Bentler & Weeks notation
  ##
  ## p      := Number of observed dependent variables (usually all indicators)
  ## q      := Number of observed independent variables (0 in a complete latent-variable model)
  ## r      := Number of observed variables, r = p + q
  ## m      := Number of dependent variables ("endogenous" constructs + indicators)
  ## n      := Number of independent variables ("exogenous" constructs + zetas + deltas)
  ## s      := Total number of variables, s = m + n
  ## beta   := (m x m) Coefficient matrix of dep. variables on dep. variables.
  ## gamma  := (m x n) Coefficient matrix of indep. variables on dep. variables.
  ## Phi    := (n x n) VCV of independent variables.
  ## Beta   := (s x s) Relationship super-matrix of dependent variables on
  ##                   dependent variables. Path from independent and dependent 
  ##                   on independent are zero.
  ## Gamma  := (s x n) Relationship super-matrix of independent on dependent 
  ##                   and independent variables.
  ## G      := (r x s) Selector matrix selecting only observable variables.
  ### --------------------------------------------------------------------------
  ### Preparation ==============================================================
  ## Check if linear
  if(.object$Information$Model$model_type != "Linear") {
    stop("Model is nonlinear. Currently the model-implied indicator covariance",
         " matrix can only be computed for linear models.", call. = FALSE)
  }
  
  Cons_exo <- .object$Information$Model$vars_exo
  Cons_endo <- .object$Information$Model$vars_endo
  
  ## Get relevant matrices
  S <- .object$Estimates$Indicator_VCV
  Phi <- .object$Estimates$Construct_VCV[Cons_exo,Cons_exo,drop=FALSE]
  B <- .object$Estimates$Path_estimates[Cons_endo,Cons_endo,drop=FALSE] 
  Gamma = .object$Estimates$Path_estimates[Cons_endo,Cons_exo,drop=FALSE]
  Lambda <- .object$Estimates$Cross_loadings*.object$Information$Model$measurement
  I = diag(length(Cons_endo))
  
  
  vcv_delta <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(vcv_delta)=dimnames(S)
  
  # calculate the vcv of zeta
  #  Currently, the calculation of the VCV of the zetas is not 100%, this should be adjusted, 
  # Calculate the variance of zeta equation by equation with replacement
  vec_zeta=1-rowSums(.object$Estimates$Path_estimates*
                       .object$Estimates$Construct_VCV)
  names(vec_zeta)=rownames(.object$Estimates$Construct_VCV)
  
  vcv_zeta=matrix(0,nrow=nrow(I),ncol=ncol(I))
  diag(vcv_zeta)=vec_zeta[Cons_endo]
  
  # Correlation among the exogenous constructs
  Corr_exo=Phi
  # Correlations  between exogenous and endogenous constructs
  Corr_exo_endo=Phi %*% t(Gamma) %*%t(solve(I-B))
  # Correlations among endogenous cosntructs 
  Cor_endo=solve(I-B)%*%(Gamma%*%Phi%*%t(Gamma)+vcv_zeta)%*%t(solve(I-B))
  diag(Cor_endo)=1
  
  VCV_construct=rbind(cbind(Corr_exo,Corr_exo_endo),
                      cbind(t(Corr_exo_endo),Cor_endo))
  
  # calculate model-implied VCV of the indicators
  VCV_ind=t(Lambda)%*%VCV_construct%*%Lambda
  
  Sigma=VCV_ind+vcv_delta
  
  # Make symmetric
  Sigma[lower.tri(Sigma)] = t(Sigma)[lower.tri(Sigma)]
  
  # ### Replace indicators connected to a composite by S
  
  mod <- .object$Information$Model
  composites <- names(mod$construct_type[mod$construct_type == "Composite"])
  index <- t(mod$measurement[composites, , drop = FALSE]) %*% mod$measurement[composites, , drop = FALSE]
  
  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  # Replace indicators which measurement errors are allowed to be correlated by s_ij
  
  return(Sigma)
}



#' Effects
#'
#' Compute direct, indirect and total effects.
#'
#' @usage effects(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export
#'
effects.cSEMResults <- function(.object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  ## Endogenous (lhs) variables
  vars_endo <- .object$Information$Model$vars_endo
  
  ## Matrix of direct effects:
  direct <- .object$Estimates$Path_estimates
  
  ## Matrix of total total effects: B = direct
  # Note: eta = B x eta + zeta
  #       (I - B)*eta = zeta
  
  B_star <- diag(nrow(direct)) - direct
  total <- solve(B_star) - diag(nrow(direct))
  
  ## Matrix of indirect effects:
  indirect <- total - direct
  
  list(direct = direct[vars_endo, , drop = FALSE], 
       indirect = indirect[vars_endo, , drop = FALSE], 
       total = total[vars_endo, , drop = FALSE])
}

#' `cSEMResults` method to check the validity of the results
#'
#' Check the validity of the estimated quantities.
#' 
#' @usage status(object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @return A NULL vector or a vector of status codes pointing to specific problems. 
#' Problems are:
#'  \itemize{
#' \item 1 Algorithm has not not converged
#' \item 2: at least one absolute standardized loading estimate is larger than 1,
#' which implies either a negative veriance of the measurement error or
#' a correlation larger than 1
#' \item 3: construct VCV is not positive semi-definit
#' \item 4 model-implied indicators VCV is not positive semi-definit
#' }
#' 
#' @export
#'

status <- function(.object){
  UseMethod("status")
}

#' @describeIn status Some text
#' @export
status.cSEMResults <- function(.object){
  
  # NULL: Everything is fine
  # 1 Algorithm has not not converged
  # 2: at least one absolute standardized loading estimate is larger than 1,
  # which implies either a negative veriance of the measurement error or
  # a correlation larger than 1
  # 3: construct VCV is not positive semi-definit
  # 4: model-implied indicators VCV is not positive semi-definit
  # 5: at least one construct reliability is larger than 1
  
  if(.object$Information$Model$model_type != "Linear"){
    stop("Currently, the status function only works for linear models.",
         call. = FALSE)}
  
  stat <- c("1" = FALSE, "2" = FALSE, "3" = FALSE, "4" = FALSE, "5" = FALSE)
  
  # if(.object$Information$Arguments$.approach_weights == "PLS") {
  
  if(!(is.null(.object$Information$Weight_info$Convergence_status) || .object$Information$Weight_info$Convergence_status)) {
    stat["1"] <- TRUE
  }
  
  if(max(abs(.object$Estimates$Cross_loadings)) > 1) {
    stat["2"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(.object$Estimates$Construct_VCV)) {
    stat["3"] <- TRUE
  }
  
  if(!matrixcalc::is.positive.semi.definite(fitted(.object))) {
    stat["4"] <- TRUE
  }
  
  if(max(.object$Estimates$Construct_reliabilities)>1) {
    stat["5"] <- TRUE
  }
  
  
  # If if no problem occured, it seems that the estimation is fine.
  if(sum(stat) == 0) {
    stat <- NULL
  }
  
  # } else {
  #   stop("Only applicable if PLS is used.", call. = FALSE)
  # }
  
  return(stat) 
}
