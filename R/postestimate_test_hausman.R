#' Regression-based Hausman test
#' 
#' Calculates the regression-based Hausman test to be used to compare 
#' OLS to 2SLS estimates or 2SLS to 3SLS estimates.
#' 
#' @usage testHausman(
#'  .object               = NULL,
#'  .alpha                = 0.05,
#'  .eval_plan            = c("sequential", "multiprocess"),
#'  .handle_inadmissibles = c("drop", "ignore", "replace"),
#'  .R                    = 499,
#'  .resample_method      = c("bootstrap", "jackknife"),
#'  .seed                 = NULL
#'  )
#' 
#' @inheritParams csem_arguments
#' 
#' @seealso [csem()], [cSEMResults]
#'   
#' @export

testHausman <- function(
  .object               = NULL,
  .alpha                = 0.05,
  .eval_plan            = c("sequential", "multiprocess"),
  .handle_inadmissibles = c("drop", "ignore", "replace"),
  .R                    = 499,
  .resample_method      = c("bootstrap", "jackknife"),
  .seed                 = NULL
) {
  
  .eval_plan <- match.arg(.eval_plan)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .resample_method      <- match.arg(.resample_method)
  
  if(!inherits(.object, "cSEMResults_default")) {
    stop2("Hausman test currently only available for linear models.")
  }
  
  # Define function to bootstrap
  fun_residuals <- function(.object) {
    
    ## Extract relevant quantities
    m  <- .object$Information$Model
    P  <- .object$Estimates$Construct_VCV
    
    # Dependent (LHS) variables of the structural equations
    dep_vars  <- rownames(m$structural)[rowSums(m$structural) != 0]
    
    res <- lapply(dep_vars, function(y) {
      # Which of the variables in dep_vars have instruments specified, i.e.
      # have endogenous variables on the RHS. By default: FALSE.
      endo_in_RHS <- FALSE
      
      if(!is.null(m$instruments)) {
        endo_in_RHS <- y %in% names(m$instruments)
      }
      
      # All independent variables (X) of the structural equation of construct y
      # (including the endogenous RHS variables)
      names_X <-  colnames(m$structural[y, m$structural[y, ] != 0, drop = FALSE])
      
      ## Only for equations with endogenous variables on the RHS do we need to compute
      ##  the Hausman test.
      if(endo_in_RHS & (.object$Information$Arguments$.approach_paths == "2SLS")) {
        ## First stage
        # Note: Technically, we only need to regress the P endogenous variables 
        #       (y1) on the L instruments and the K exogenous independent variables 
        #       (which must be part of Z).
        #       Therefore: y1 (N x P) and Z (N x (L + K)) and
        #       beta_1st = (Z'Z)^-1*(Z'y1)
        #       
        #       beta_1st would be ((L + K) x P), i.e. each columns represents the 
        #       first stage estimates for a regression of the instruments on
        #       on the p'th endogenous variable.
        names_endo <- rownames(m$instruments[[y]]) # names of the endogenous variables (y1's)
        names_Z <- colnames(m$instruments[[y]]) # names of the instruments 
        
        # Assuming that P (the construct correlation matrix) also contains 
        # the instruments (ensured if only internal instruments are allowed)
        
        beta_1st <- solve(P[names_Z, names_Z, drop = FALSE], 
                          P[names_Z, names_endo, drop = FALSE])
        
        ## Second stage
        # Note: we have to use (construct) correlations here since using 
        #       the plain score would yield inconsistent estimates for concepts
        #       modeled as common factors. Hence we need to "translate" the
        #       regression based second stage estimation 
        #                  y = X * beta_1 + v_hat*beta_2 + u
        #       and the estimator
        #                  beta_2nd = (W'W)^-1W'y
        #       into correlations.
        #
        #   y1_hat = Z*beta_1st
        #   v_hat  = y1 - y1_hat = y1 - Z*beta_1st
        #
        #   y = X * beta_1 + v_hat*beta_2 + u
        #
        #   Define W = [X; v_hat]
        #
        #   E(W'W ) = E[X'X ; X'v_hat
        #               v_hat'X ; v_hat'v_hat]
        #   E(W'y)  = E[X'y ; v_hat'y]'
        
        ww11 <- P[names_X, names_X, drop = FALSE]
        ww12 <- P[names_X, names_endo, drop = FALSE] - P[names_X, names_Z, drop = FALSE] %*% beta_1st
        ww21 <- t(ww12)
        ww22 <- P[names_endo, names_endo, drop = FALSE] - 
          P[names_endo, names_Z, drop = FALSE] %*% beta_1st -
          t(beta_1st) %*% P[names_Z, names_endo, drop = FALSE] +
          t(beta_1st) %*% P[names_Z, names_Z, drop = FALSE] %*% beta_1st
        
        WW <- rbind(cbind(ww11, ww12), cbind(ww21, ww22))
        
        wy1 <- P[names_X, y, drop = FALSE]
        wy2 <- P[names_endo, y, drop = FALSE] - t(beta_1st) %*% P[names_Z, y, drop = FALSE]
        
        Wy <- rbind(wy1, wy2)
        
        ## Estimate the second stage estimation
        
        beta_2nd <- solve(WW, Wy)
        rownames(beta_2nd) <- c(names_X, paste0("Resid_", names_endo))
        colnames(beta_2nd) <- y
        
        return(beta_2nd)
        
      } else {
        NA
      }
    })
    
    names(res) <- dep_vars
    res        <- Filter(Negate(anyNA), res) 
    
    # Vectorize "res" to be able to use resamplecSEMResults (which requires
    # output to be a vector or a matrix)
    # 1. Assign a unique name to each element
    res <- lapply(res, function(x) {
      rownames(x) <- paste0(colnames(x), "_", rownames(x))
      x
    })
    # 2. Combine, vectorize and assign names
    res2 <- do.call(rbind, res)
    res <- c(res2)
    names(res) <- rownames(res2)
    res
  }
  
  ## Resample to get standard errors
  out <- resamplecSEMResults(
    .object               = .object,
    .resample_method      = .resample_method,
    .R                    = .R,
    .handle_inadmissibles = .handle_inadmissibles,
    .user_funs            = fun_residuals,
    .eval_plan            = .eval_plan,
    .seed                 = .seed
  )
  
  ## Get relevant quantities
  out_infer <- infer(out, .alpha = .alpha)
  beta      <- out$Estimates$Estimates_resample$Estimates1$User_fun$Original
  se        <- out_infer$User_fun$sd
  t         <- beta/se
  p_normal  <- 2*pnorm(abs(t), mean = 0, sd = 1, lower.tail = FALSE) 
  ci_percentile <- out_infer$User_fun$CI_percentile 
  
  out_data_frame <- data.frame(
    "Name"       = names(beta),
    "Estimate"   = beta,
    "Std. error" = se,
    "t-stat."    = t,
    "p-value"    = p_normal,
    "Ci_perc. L" = ci_percentile[1, ],
    "Ci_perc. H" = ci_percentile[2, ],
    stringsAsFactors = FALSE
  )
  rownames(out_data_frame) <- NULL
  return(out_data_frame)
}
