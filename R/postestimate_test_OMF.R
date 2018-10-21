#' Test for overall model fit
#'
#' Test the overall model fit.
#' 
#' Description: TODO.
#' Test is based on \insertCite{Dijkstra2015;textual}{cSEM}.
#' Data is transformed according to \insertCite{Beran1985;textual}{cSEM}.
#' 
#' @usage testOMF(
#'  .object              = args_default()$.object, 
#'  .alpha               = args_default()$.alpha, 
#'  .drop_inadmissibles  = args_default()$.drop_inadmissibles, 
#'  .runs                = args_default()$.runs, 
#'  .saturated           = args_default()$.saturated,
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
  .object              = args_default()$.object,
  .alpha               = args_default()$.alpha,
  .drop_inadmissibles  = args_default()$.drop_inadmissibles,
  .runs                = args_default()$.runs,
  .saturated           = args_default()$.saturated
  ){
  
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  ## Check if single
  if(attr(.object, "single") == FALSE) {
    stop("`testOMF()` not applicable to multiple groups or data sets.\n",
         "Use `lapply(.object, testOMF)` instead.",
         call. = FALSE)
  }
  
  ## Extract required information 
  X         <- .object$Information$Data
  S         <- .object$Estimates$Indicator_VCV
  # nObs      <- nrow(X)
  Sigma_hat <- fit(.object,
                   .saturated = .saturated,
                   .type_vcv  = 'indicator')
  
  
  ## Calculate test statistic
  teststat <- c(
    "dG"   = dG(.matrix1 = S, .matrix2 = Sigma_hat),
    "SRMR" = SRMR(.object, .saturated = .saturated),
    "dL"   = dL(.matrix1 = S, .matrix2 = Sigma_hat)
  )
  
  ## Transform dataset, see Beran & Srivastava (1985)
  S_half   <- solve(expm::sqrtm(S))
  Sig_half <- expm::sqrtm(Sigma_hat)
  
  X_trans           <- X %*% S_half %*% Sig_half
  colnames(X_trans) <- colnames(X)
  
  ## Collect arguments
  arguments <- .object$Information$Arguments
  
  ## Calculate reference distribution
  ref_dist <- lapply(replicate(.runs, X_trans[sample(1:nrow(X), replace = TRUE), ], simplify = FALSE), function(x) {
    
    # Replace data
    arguments[[".data"]] <- x
    
    # Run estimation
    Est_temp <- do.call(csem, arguments)
    
    # Check if admissible
    status_code <- verify(Est_temp)
    
    # if it is controlled for inadmissible
    if(.drop_inadmissibles){
      if(sum(status_code) == 0){
        
        c("dG"   = dG(Est_temp$Estimates$Indicator_VCV, fit(Est_temp, 
                                                            .saturated = .saturated,
                                                            .type_vcv= 'indicator')),
          "SRMR" = SRMR(Est_temp, .saturated = .saturated),
          "dL"   = dL(Est_temp$Estimates$Indicator_VCV, fit(Est_temp,
                                                            .saturated = .saturated,
                                                            .type_vcv='indicator'))
        ) 
      } else {
        NULL
      }
    } else { 
      c("dG"   = dG(Est_temp$Estimates$Indicator_VCV, fit(Est_temp,
                                                          .saturated = .saturated,
                                                          .type_vcv='indicator')),
        "SRMR" = SRMR(Est_temp, .saturated = .saturated),
        "dL"   = dL(Est_temp$Estimates$Indicator_VCV, fit(Est_temp,
                                                          .saturated = .saturated,
                                                          .type_vcv='indicator'))
      ) 
    }
  }) # lapply
  
  ref_dist_matrix <- do.call(cbind, ref_dist)
  critical_value  <- matrix(apply(ref_dist_matrix, 1, quantile, 1-.alpha), 
                            ncol =length(teststat),
                            dimnames = list(paste(.alpha*100, sep = "","%"), 
                                            names(teststat))
                            )
  
  if (length(.alpha) > 1) {
    decision <- t(apply(critical_value, 1, function(x) {teststat < x}))
  }
    
  if (length(.alpha) == 1) {
    decision <- teststat < critical_value
  }
  
  out <- list(
    "Test_statistic"     = teststat,
    "Critical_value"     = critical_value, 
    "Decision"           = decision, 
    "Number_admissibles" = ncol(ref_dist_matrix)
    ) 
  
  class(out) <- "cSEMTestOMF"
  return(out)
}