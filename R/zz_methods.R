#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()]. 
#' 
#' Prints basic information on the [cSEMResults] object and suggested next 
#' steps to the console.
#'
#' @usage print(object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
print.cSEMResults <- function(object) {

  cat(cli::rule(line = "bar2"), "\n")
  cat(cli::rule(center = "Overview"), "\n\n")
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
}

#' `cSEMResultssummary` method for `print()`
#'
#' The [cSEMResultssummary] method for the generic function [print()]. 
#' 
#' Prints a summary of the results obtained from runinng [csem], [cca], or
#' [workhorse].
#'
#' @usage print(object)
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
#' [workhorse].
#'
#' @usage summary(object, .what = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#' 
#' @inherit csem_resultssummary return
#'
#' @export
#'
summary.cSEMResults <- function(object, .what = NULL) {

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
#' @usage fitted(object, test = NULL)
#'
#' @inheritParams csem_arguments
#' @param test not used
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
fitted.cSEMResults <- function(object) {
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
  if(object$Information$Model$model_type != "Linear") {
    stop("Model is nonlinear. Currently the model-implied indicator covariance",
         " matrix can only be computed for linear models.", call. = FALSE)
  }
    
  ## Get relevant matrices
  
  S <- object$Estimates$Indicator_VCV
  P <- object$Estimates$Construct_VCV
  B <- object$Estimates$Path_estimates 
  Lambda <- object$Estimates$Loading_estimates 
  Lambda_cross <- object$Estimates$Cross_loadings
  
  a1 <- object$Information$Model$vars_exo
  a2 <- object$Information$Model$vars_endo
  
  ## Compute variances of the errors (diagonal matrices Var(delta) and Var(zeta)).
  
  vcv_delta <- diag(diag(S) - diag(t(Lambda) %*% P %*% Lambda))
  Delta <- vcv_delta
  diag(Delta) <- 1

  tmp <- diag(1 - diag(solve(diag(length(a2)) - B[a2, a2]) %*% B[a2, a1]%*% 
    P[a1, a1] %*% t(B[a2, a1]) %*% t(solve(diag(length(a2)) - B[a2, a2]))))
  rownames(tmp) <- colnames(tmp) <- a2
  
  # vcv_zeta <- diag(1 - diag(B%*% P %*% t(B)))
  # vcv_zeta[1:length(a), 1:length(a)] <- 0
  vcv_zeta <- matrix(0, nrow(B), ncol(B), dimnames = list(rownames(B), colnames(B)))
  vcv_zeta[a2, a2] <- tmp

  Zeta <- vcv_zeta
  Zeta[Zeta != 0] <- 1
  
  ## Get names of all variables (dependent and independent)
  names <- c(colnames(S), rownames(Lambda), paste0("del", 1:nrow(vcv_delta)),
             paste0("zeta", 1:nrow(vcv_zeta)))
  
  ## Matrix of structural relationships
  # Note: Here, "structural" refers to realtionsips between variables in general
  #       not just constructs. The matrix is (s x s), where s is the total number
  #       of variables.
  
  A1 <- rbind(
    cbind(S * 0, t(Lambda), Delta, matrix(0, nrow(S), ncol(Zeta))),
    cbind(Lambda * 0, B, matrix(0, nrow(B), ncol(Delta)), Zeta)
  )

  A1 <- rbind(A1,
    matrix(0, nrow(Delta), ncol(A1)),
    matrix(0, nrow(Zeta), ncol(A1))
  )
  rownames(A1) <- colnames(A1) <- names   
  ## Set rows of variables who are only related to an error to 0 .
  # A1[rowSums(A1) == 1, ] <- 0  ## Set rows of variables who are only related to an error to 0 .
  
  ## Matrix of variances and covariances between all variables 
  # Note: is necessary for the computation of Phi. The matrix is also (s x s)
  
  A2 <- rbind(
    cbind(S, t(Lambda_cross), matrix(0, nrow(S), sum(ncol(vcv_delta), ncol(vcv_zeta)))),
    cbind(Lambda_cross, P, matrix(0, nrow(P), sum(ncol(vcv_delta), ncol(vcv_zeta)))),
    cbind(matrix(0, nrow(vcv_delta), sum(ncol(S) + nrow(Lambda_cross))), 
          vcv_delta,
          matrix(0, nrow(vcv_delta), ncol(vcv_zeta))),
    cbind(matrix(0, nrow(vcv_zeta), 
                 sum(ncol(S) + nrow(Lambda_cross) + ncol(vcv_delta))), vcv_zeta) 
  )
  rownames(A2) <- colnames(A2) <- names 

  ## Distinguish and get dimensions --------------------------------------------
  indep_vars <- rownames((A1[rowSums(A1) == 0, ]))
  dep_vars   <- setdiff(rownames(A1), indep_vars)
  observed   <- colnames(S)
  # unobserved <- setdiff(rownames(A1), observed)
  
  m <- length(dep_vars)
  n <- length(indep_vars)
  s <- m + n
  p <- length(intersect(dep_vars, observed)) 
  q <- length(intersect(indep_vars, observed))
  r <- p + q

  ### Define the required matrices ---------------------------------------------
  
  Beta <- matrix(0, s, s, dimnames = list(names, names))  # (s x s)
  Beta[dep_vars, dep_vars] <- A1[dep_vars, dep_vars] # beta0 (m x m)
  
  Gamma <- matrix(0, s, n, dimnames = list(names, indep_vars))  # (s x n)
  Gamma[dep_vars, ]   <- A1[dep_vars, indep_vars] # gamma (m x n)
  Gamma[indep_vars, ] <- diag(n) # I (n x n)
  
  G <- matrix(0, r, s, dimnames = list(observed, names))  # (r x s)
  G[observed, observed] <- diag(r) # (r x s)
  
  Phi   <- A2[indep_vars, indep_vars] # (n x n)
  I     <- diag(s)          # (s x s)
  
  ### Calculate the model implied covariance matrix ============================

  Sigma <- G %*% solve(I - Beta) %*% Gamma %*% Phi %*% t(Gamma) %*% t(solve(I - Beta))%*% t(G)
  rownames(Sigma) <- colnames(Sigma) <- rownames(S)

  # ### Replace indicators connected to a composite by S

  mod <- object$Information$Model
  composites <- names(mod$construct_type[mod$construct_type == "Composite"])
  index <- t(mod$measurement[composites, , drop = FALSE]) %*% mod$measurement[composites, , drop = FALSE]

  Sigma[which(index == 1)] <- S[which(index == 1)]
  
  diag(Sigma) <- 1
  
  
  # make symmetric; it may happen that the matrix is symmetric due to maschine calculations
  Sigma[lower.tri(Sigma)] = t(Sigma)[lower.tri(Sigma)]
  
  return(Sigma)
}

#' Effects
#'
#' Compute direct, indirect and total effects.
#'
#' @usage effects(object)
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
effects.cSEMResults <- function(object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  ## Endogenous (lhs) variables
  vars_endo <- object$Information$Model$vars_endo
  
  ## Matrix of direct effects:
  direct <- object$Estimates$Path_estimates

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

status <- function(object){
  UseMethod("status")
}

status.cSEMResults <- function(object){
  
  # 0: Everything is fine
  # 1 Algorithm has not not converged
  # 2: at least one absolute standardized loading estimate is larger than 1,
  # which implies either a negative veriance of the measurement error or
  # a correlation larger than 1
  # 3: construct VCV is not positive semi-definit
  # 4 model-implied indicators VCV is not positive semi-definit

  stat <- c("1" = FALSE, "2" = FALSE, "3" = FALSE, "4" = FALSE)
  
  if(object$Information$Weight_approach == "PLS") {
    
    if(!object$Information$Convergence_status) {
      stat["1"] <- TRUE
    }
    
    if(max(abs(object$Estimates$Cross_loadings)) > 1) {
      stat["2"] <- TRUE
    }
    
    if(!matrixcalc::is.positive.semi.definite(object$Estimates$Construct_VCV)) {
      stat["3"] <- TRUE
    }
    
    if(!matrixcalc::is.positive.semi.definite(fitted(object))) {
      stat["4"] <- TRUE
    }
    # If if no problem occured, it seems that the estimation is fine.
    if(sum(stat) == 0) {
      stat <- NULL
    }
    
  } else {
    stop("Only applicable if PLS is used.", call. = FALSE)
  }
  
  return(stat) 
}
