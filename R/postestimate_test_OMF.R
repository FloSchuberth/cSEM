#' Test for overall model fit
#'
#' Test the overall model fit.
#' 
#' Description: TODO.
#' Test is based on \insertCite{Dijkstra2015;textual}{cSEM}.
#' Data is transformed according to \insertCite{Beran1985;textual}{cSEM}.
#' After 10000 iterations it stops automatically
#' 
#' @usage testOMF(
#'  .object                = args_default()$.object, 
#'  .alpha                 = args_default()$.alpha, 
#'  .handle_inadmissibles  = args_default()$.handle_inadmissibles, 
#'  .runs                  = args_default()$.runs, 
#'  .saturated             = args_default()$.saturated,
#'  .verbose               = args_default()$.verbose
#' )
#' 
#' @inheritParams  csem_arguments
#' 
#' @inherit csem_test return
#' 
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export

testOMF <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .runs                  = args_default()$.runs,
  .saturated             = args_default()$.saturated,
  .verbose               = args_default()$.verbose
) {
  # Implementation is based on:
  # Dijkstra & Henseler (2015) - Consistent Paritial Least Squares Path Modeling
  
  if(.verbose) {
    cat(rule(center = "Test for overall model fit based on Dijkstra & Henseler (2015)",
             line = "bar3"), "\n\n")
  }
  
  UseMethod("testOMF")
}
  
#' @describeIn testOMF (TODO)
#' @export

testOMF.cSEMResults_default <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .runs                  = args_default()$.runs,
  .saturated             = args_default()$.saturated,
  .verbose               = args_default()$.verbose
){
  ## Check arguments
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
  
  ## Check if initial results are inadmissible
  if(sum(verify(.object)) != 0) {
    stop("Initial estimation results are inadmissible.\n", 
         "See `verify(.object)` for details.",
         call. = FALSE)
  }
  
  ## Extract required information 
  X         <- .object$Information$Data
  S         <- .object$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object,
                   .saturated = .saturated,
                   .type_vcv  = "indicator")
  
  ## Calculate test statistic
  teststat <- c(
    "dG"   = dG(S, Sigma_hat),
    "SRMR" = SRMR(S, Sigma_hat),
    "dL"   = dL(S, Sigma_hat)
  )
  
  ## Transform dataset, see Beran & Srivastava (1985)
  S_half   <- solve(expm::sqrtm(S))
  Sig_half <- expm::sqrtm(Sigma_hat)
  
  X_trans           <- X %*% S_half %*% Sig_half
  colnames(X_trans) <- colnames(X)
  
  ## Collect arguments
  arguments <- .object$Information$Arguments
  
  # Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .runs, style = 3)
  }
  
  ## Calculate reference distribution
  ref_dist         <- list()
  n_inadmissibles  <- 0
  i <- 0
  repeat{
    # Counter
    i <- i + 1
    
    # Draw dataset
    X_temp <- X_trans[sample(1:nrow(X), replace = TRUE), ]
    
    # Replace the old dataset by the new one
    arguments[[".data"]] <- X_temp
    
    # Estimate model
    Est_temp <- do.call(foreman, arguments) # its important to use foreman here 
    # instead of csem() to allow for lapply(x, testOMF.cSEMResults_default) when x 
    # is of class cSEMResults_2ndorder.
    
    # Check status
    status_code <- sum(verify(Est_temp))
    
    # Distinguish depending on how inadmissibles should be handled
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
      # Compute if status is ok or .handle inadmissibles = "ignore" AND the status is 
      # not ok
      S_temp         <- Est_temp$Estimates$Indicator_VCV
      Sigma_hat_temp <- fit(Est_temp,
                            .saturated = .saturated,
                            .type_vcv  = "indicator")
      
      ref_dist[[i]] <- c(
        "dG"   = dG(S_temp, Sigma_hat_temp),
        "SRMR" = SRMR(S_temp, Sigma_hat_temp),
        "dL"   = dL(S_temp, Sigma_hat_temp)
      ) 
      
    } else if(status_code != 0 & .handle_inadmissibles == "drop") {
      # Set list element to zero if status is not okay and .handle_inadmissibles == "drop"
      ref_dist[[i]] <- NULL
      
    } else {# status is not ok and .handle_inadmissibles == "replace"
      # Reset counter and raise number of inadmissibles by 1
      i <- i - 1
      n_inadmissibles <- n_inadmissibles + 1
    }
    
    # Break repeat loop if .runs results have been created.
    if(length(ref_dist) == .runs) {
      break
    } else if(i + n_inadmissibles == 10000) { 
      ## Stop if 10000 runs did not result in insufficient admissible results
      stop("Not enough admissible result.", call. = FALSE)
    }
    if(.verbose){
      setTxtProgressBar(pb, i)
    }
    
  } # END repeat 
  
  # close progress bar
  if(.verbose){
    close(pb)
  }
  
  # Combine
  ref_dist_matrix <- do.call(cbind, ref_dist) # note: this drops the NULL elements,
  # so only the admissibles remain.
  ## Compute critical values (Result is a (2 x p) matrix, where n is the number
  ## of quantiles that have been computed (1 by default)
  .alpha <- .alpha[order(.alpha)]
  critical_values <- matrixStats::rowQuantiles(ref_dist_matrix, 
                                               probs =  1-.alpha, drop = FALSE)
  
  ## Compare critical value and teststatistic
  decision <- teststat < critical_values # a logical (2 x p) matrix with each column
                                         # representing the decision for one
                                         # significance level. TRUE = no evidence 
                                         # against the H0 --> not reject
                                         # FALSE --> reject
  
  # Return output
  out <- list(
    "Test_statistic"     = teststat,
    "Critical_value"     = critical_values, 
    "Decision"           = decision, 
    "Information"        = list(
      "Number_admissibles" = ncol(ref_dist_matrix),
      "Total_runs"         = i + n_inadmissibles 
    )
  )
  
  class(out) <- "cSEMTestOMF"
  return(out)
}

#' @describeIn testOMF (TODO)
#' @export

testOMF.cSEMResults_multi <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .runs                  = args_default()$.runs,
  .saturated             = args_default()$.saturated,
  .verbose               = args_default()$.verbose
){
  lapply(.object, testOMF.cSEMResults_default,
         .alpha                = .alpha,
         .handle_inadmissibles = .handle_inadmissibles,
         .runs                 = .runs,
         .saturated            = .saturated,
         .verbose              = .verbose
         )
}

#' @describeIn testOMF (TODO)
#' @export

testOMF.cSEMResults_2ndorder <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .runs                  = args_default()$.runs,
  .saturated             = args_default()$.saturated,
  .verbose               = args_default()$.verbose
){
  ## Check arguments
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)

  ## Check if initial results are inadmissible
  if(sum(unlist(verify(.object))) != 0) {
    stop("Initial estimation results are inadmissible.\n", 
         "See `verify(.object)` for details.",
         call. = FALSE)
  }
  
  ## Extract required information 
  x1 <- .object$First_stage
  x2 <- .object$Second_stage
  X         <- x1$Information$Data
  S         <- x1$Estimates$Indicator_VCV
  Sigma_hat <- fit(.object,
                   .saturated = .saturated,
                   .type_vcv  = "indicator")
  
  ## Calculate test statistic
  teststat <- c(
    "dG"   = dG(S, Sigma_hat),
    "SRMR" = SRMR( S, Sigma_hat),
    "dL"   = dL(S, Sigma_hat)
  )
  
  ## Transform dataset, see Beran & Srivastava (1985)
  S_half   <- solve(expm::sqrtm(S))
  Sig_half <- expm::sqrtm(Sigma_hat)
  
  X_trans           <- X %*% S_half %*% Sig_half
  colnames(X_trans) <- colnames(X)
  
  ## Collect arguments
  arguments <- x2$Information$Arguments_original
  
  # Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .runs, style = 3)
  }
  
  ## Calculate reference distribution
  ref_dist         <- list()
  n_inadmissibles  <- 0
  i <- 0
  repeat{
    # Counter
    i <- i + 1
    
    # Draw dataset
    X_temp <- X_trans[sample(1:nrow(X), replace = TRUE), ]
    
    # Replace the old dataset by the new one
    arguments[[".data"]] <- X_temp
    
    # Estimate model
    Est_temp <- do.call(csem, arguments)               
    
    # Check status (Note: output of verify for second orders is a list)
    status_code <- sum(unlist(verify(Est_temp)))
    
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
      # Compute if status is ok or .handle inadmissibles = "ignore" AND the status is 
      # not ok
      S_temp         <- Est_temp$First_stage$Estimates$Indicator_VCV
      Sigma_hat_temp <- fit(Est_temp,
                            .saturated = .saturated,
                            .type_vcv  = "indicator")
      
      ref_dist[[i]] <- c(
        "dG"   = dG(S_temp, Sigma_hat_temp),
        "SRMR" = SRMR(S_temp, Sigma_hat_temp),
        "dL"   = dL(S_temp, Sigma_hat_temp)
      )  
      
    } else if(status_code != 0 & .handle_inadmissibles == "drop") {
      # Set list element to zero if status is not okay and .handle_inadmissibles == "drop"
      ref_dist[[i]] <- NULL
      
    } else {# status is not ok and .handle_inadmissibles == "replace"
      # Reset counter and raise number of inadmissibles by 1
      i <- i - 1
      n_inadmissibles <- n_inadmissibles + 1
    }
    
    # Break repeat loop if .runs results have been created.
    if(length(ref_dist) == .runs) {
      break
    } else if(i + n_inadmissibles == 10000) { 
      ## Stop if 10000 runs did not result in insufficient admissible results
      stop("Not enough admissible result.", call. = FALSE)
    }
    if(.verbose){
      setTxtProgressBar(pb, i)
    }
    
  } # END repeat 
  
  # close progress bar
  if(.verbose){
    close(pb)
  }

  # Combine
  ref_dist_matrix <- do.call(cbind, ref_dist) # note: this drops the NULL elements,
  # so only the admissibles remain.
  ## Compute critical values (Result is a (2 x p) matrix, where n is the number
  ## of quantiles that have been computed (1 by default)
  .alpha <- .alpha[order(.alpha)]
  critical_values <- matrixStats::rowQuantiles(ref_dist_matrix, 
                                               probs =  1-.alpha, drop = FALSE)
  
  ## Compare critical value and teststatistic
  decision <- teststat < critical_values # a logical (2 x p) matrix with each column
                                         # representing the decision for one
                                         # significance level. TRUE = no evidence 
                                         # against the H0 --> not reject
                                         # FALSE --> reject
  
  # Return output
  out <- list(
    "Test_statistic"     = teststat,
    "Critical_value"     = critical_values, 
    "Decision"           = decision, 
    "Information"        = list(
      "Number_admissibles" = ncol(ref_dist_matrix),
      "Total_runs"         = i + n_inadmissibles 
    )
  ) 
  
  class(out) <- "cSEMTestOMF"
  return(out)
  
}