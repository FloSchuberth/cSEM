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

  cat(rule(line = "bar2"), "\n")
  cat(rule(center = "Overview"), "\n\n")

  if(!all(names(.object) %in% c("Estimates", "Information")) ) {
    cat("Estimation status by group/data set:\n\n\t")
        for(i in names(.object)) {
          cat(col_align(cyan(i), 15), ": ",
              ifelse(.object[[i]]$Information$Weight_info$Convergence_status == TRUE, 
                     green("successful\n\t"), red("not successful\n\t")), sep = "")
        }
  } else {
    cat("Estimation was ", ifelse(.object$Information$Weight_info$Convergence_status == TRUE, 
                                 green("successful.\n"), red("not successful.\n")), sep = "") 
  }
    cat("\nThe result is a list of class " %+% bold("cSEMResults") %+%" with list elements:\n\n\t",
        "- ", green("Estimates\n\t"),
        "- ", green("Information\n\n"), sep = "")
    cat("To get an overview or help type:\n\n\t",
        "- ", magenta("?"), cyan("cSEMResults"),"\n\t",
        "- ", magenta("str"), "(", cyan("<object-name>"), ")\n\t",
        "- ", magenta("listviewer"), yellow("::"), magenta("jsondedit"),
        "(", cyan("<object-name>"), ", ", red("mode"), " = ", cyan("'view'"), ")\n\n", sep = "")
    cat("If you wish to access the list elements directly type e.g. \n\n\t",
        "- ", cyan("<object-name>"), yellow("$"), green("Estimates"), "\n\n", sep = "")
    cat("Frequently used available postestimation commands:\n\n\t",
        "- ", magenta("summary"), "(", cyan("<object-name>"), ")\n\t",
        "- ", magenta("test"), "(", cyan("<object-name>"), ")\n\t"  ,
        "- ", magenta("status"), "(", cyan("<object-name>"), ")\n", sep = "")
    cat(rule(line = "bar2"), "\n")
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
#' Notation anis taken from Bollen (1989), page 319 f.
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
  if(.object$Information$Model$model_type != "Linear") {
    stop("Model is nonlinear. Currently the model-implied indicator covariance",
         " matrix can only be computed for linear models.", call. = FALSE)
  }
  
  Cons_exo  <- .object$Information$Model$vars_exo
  Cons_endo <- .object$Information$Model$vars_endo
  
  ## Get relevant matrices
  S      <- .object$Estimates$Indicator_VCV
  B      <- .object$Estimates$Path_estimates[Cons_endo, Cons_endo, drop = FALSE]
  Gamma  <- .object$Estimates$Path_estimates[Cons_endo, Cons_exo, drop = FALSE]
  Lambda <- .object$Estimates$Loading_estimates
  Phi    <- .object$Estimates$Construct_VCV[Cons_exo, Cons_exo, drop = FALSE]
  I      <- diag(length(Cons_endo))
  Theta  <- diag(diag(S) - diag(t(Lambda) %*% Lambda))
  dimnames(Theta) <- dimnames(S)
  
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
  
  VCV_construct <- rbind(cbind(PH, Corr_exo_endo),
                         cbind(t(Corr_exo_endo), Cor_endo))
  
  ## Calculate model-implied VCV of the indicators
  VCV_ind <- t(Lambda) %*% VCV_construct %*% Lambda
  
  Sigma <- VCV_ind + Theta
  
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
  # 4 model-implied indicators VCV is not positive semi-definit
  
  if(.object$Information$Model$model_type != "Linear"){
    stop("Currently, the status function only works for linear models.",
         call. = FALSE)}

  stat <- c("1" = FALSE, "2" = FALSE, "3" = FALSE, "4" = FALSE)
  
  # if(.object$Information$Arguments$.approach_weights == "PLS") {
    
    if(!.object$Information$Weight_info$Convergence_status) {
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
    # If if no problem occured, it seems that the estimation is fine.
    if(sum(stat) == 0) {
      stat <- NULL
    }
    
  # } else {
  #   stop("Only applicable if PLS is used.", call. = FALSE)
  # }
  
  return(stat) 
}
