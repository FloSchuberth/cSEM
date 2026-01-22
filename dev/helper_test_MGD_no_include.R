#' Decision fot test of multigroup difference is based on the critical values instead
#' of p-values 
#' 
#' @usage decide_by_critical=function(
#' .test_MGD = NULL,
#' .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
#' .alpha = args_default()$.alpha)
#' 
#' @inheritParams csem_arguments
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [cSEMResults]
#'
#' @export
decide_by_critical=function(
  .test_MGD = NULL,
  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
  .alpha = args_default()$.alpha){
  
  # Check whether .test_MGD is of class "cSEMTestMGD"
  if(class(.test_MGD) != "cSEMTestMGD"){
    stop2("Please provide an object of class cSEMTestMGD to the function decide_by_crtical")
  }
  
  .alpha <- .alpha[order(.alpha)]
  
  ### Approach suggested by Klesel et al. (2019) -------------------------------
  if("Klesel" %in% names(.test_MGD)){
    # Extract test statistic
    teststat_Klesel <- .test_MGD$Klesel$Test_statistic
    
    # Extract reference distribution
    ref_dist_matrix_Klesel <- .test_MGD$Information$Permutation_values$Klesel
    
    # Compute critical values (Result is a (2 x p) matrix, where n is the number
    # of quantiles that have been computed (1 by default)
    critical_values_Klesel <- matrixStats::rowQuantiles(ref_dist_matrix_Klesel, 
                                                        probs =  1-.alpha, drop = FALSE)
    
    ## Compare critical value and test statistic
    decision_Klesel <- teststat_Klesel < critical_values_Klesel 
    # a logical (2 x p) matrix with each column
    # representing the decision for one significance level. 
    # TRUE = no evidence against the H0 --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
  }
  
  ### Approach suggested by Chin & Dibbern (2010) ------------------------------
  if("Chin" %in% names(.test_MGD)){
    # Extract test statistic
    teststat_Chin <- .test_MGD$Chin$Test_statistic
    
    #Extract reference distribution
    ref_dist_matrices_Chin <- .test_MGD$Information$Permutation_values$Chin
    
    # Calculation of adjusted alphas:
    # Number of comparisons equals the number of parameter times the number of group comparisons
    alpha_Chin <- lapply(.approach_alpha_adjust, function(x){
      cSEM:::adjustAlpha(
        .alpha = .alpha,
        .approach_alpha_adjust = x,
        .nr_comparison = sum(sapply(teststat_Chin,length)))
    })
    
    names(alpha_Chin) <- .approach_alpha_adjust 
    
    # Compute the critical values for all alphas
    critical_values_Chin <- lapply(alpha_Chin, function(alpha_list) {
      res <- lapply(alpha_list, function(alpha) {
        probs_Chin <- c(alpha/2,1-alpha/2)
        temp<- lapply(ref_dist_matrices_Chin, function(x){
          matrixStats::rowQuantiles(x, probs = probs_Chin, drop = FALSE)
        })
      })
      names(res) <- paste0(alpha_list*100, '%')
      return(res)
    })
    
    names(critical_values_Chin) <- .approach_alpha_adjust 
    
    decision_Chin <- lapply(critical_values_Chin, function(critical_list){
      lapply(critical_list,function(critical){# goes over the different significance levels
        temp=mapply(function(stat,crit){# goes over the different group comparisons
          crit[,1]< stat & stat < crit[,2]
          
        },stat=teststat_Chin,crit=critical,SIMPLIFY = FALSE)
        return(temp)
      })
    })
    
    # Overall decision, i.e., was any of the test belonging to one significance levl rejected
    decision_overall_Chin = lapply(decision_Chin, function(decision_Chin_list){
      lapply(decision_Chin_list,function(x){
        all(unlist(x))
      })
    })
  }
  # Approach suggested by Sarstedt et al. (2011)-------------------------------------------
  if("Sarstedt" %in% names(.test_MGD)){
    
    # Extract test statistic
    teststat_Sarstedt <- .test_MGD$Sarstedt$Test_statistic 
    
    #Extract reference distribution 
    ref_dist_matrix_Sarstedt <- .test_MGD$Information$Permutation_values$Sarstedt
    
    # Calculation of adjusted alphas: 
    # Number of the comparisons equals the number of parameters that are compared 
    alpha_Sarstedt <- lapply(.approach_alpha_adjust, function(x){
      cSEM:::adjustAlpha(
        .alpha = .alpha,
        .approach_alpha_adjust = x,
        .nr_comparison = nrow(ref_dist_matrix_Sarstedt))
    })
    
    names(alpha_Sarstedt) <- .approach_alpha_adjust 
    
    critical_values_Sarstedt <- lapply(alpha_Sarstedt, function(alpha_list) {
      matrixStats::rowQuantiles(ref_dist_matrix_Sarstedt, 
                                probs =  1-alpha_list, drop = FALSE)
    })
    
    names(critical_values_Sarstedt) <- .approach_alpha_adjust 
    
    # Compare critical value and test statistic    
    decision_Sarstedt <- lapply(critical_values_Sarstedt, function(critical){
      teststat_Sarstedt < critical
    })
    # FALSE: Reject
    # TRUE: don't reject
    
    names(decision_Sarstedt) <- .approach_alpha_adjust
    
    # Overall decision, i.e., was any of the test belonging to one significance levl rejected
    decision_overall_Sarstedt = lapply(decision_Sarstedt, function(x){
      all(unlist(x))
    })
    
  }
  
  ### Return output ------------------------------------------------------------
  out <- list()
  if("Klesel" %in% names(.test_MGD)){
    out[["Klesel"]] <- list(
      "Test_statistic"     = teststat_Klesel,
      "Critical_value"     = critical_values_Klesel, 
      "Decision"           = decision_Klesel) 
  }
  if("Chin" %in% names(.test_MGD)){
    out[["Chin"]] <- list(
      "Test_statistic"     = teststat_Chin,
      "Critical_value"     = critical_values_Chin, 
      "Decision"           = decision_Chin,
      "Decision_overall"   = decision_overall_Chin,
      "Alpha_adjusted"     = alpha_Chin
    )
  }
  # "Information"        = list(
  #   "Number_admissibles"    = ncol(ref_dist_matrix_Klesel),
  #   "Total_runs"            = counter + n_inadmissibles,
  #   "Group_names"           = names(.object),
  #   "Number_of_observations"= sapply(X_all_list, nrow),
  #   "Permutation_values"      = list(
  #     "Klesel" = ref_dist_Klesel,
  #     "Chin"   = ref_dist_matrices_Chin),
  #   "Approach" = .approach_mgd,
  #   "Seed"     = .seed,
  #   "Alpha"    = .alpha
  # )
  
  
  if("Sarstedt" %in% names(.test_MGD)){
    out[["Sarstedt"]] <- list(
      "Test_statistic"   = teststat_Sarstedt,
      "Critical_value"     = critical_values_Sarstedt, 
      "Decision"           = decision_Sarstedt,
      "Decision_overall"   = decision_overall_Sarstedt,
      "Alpha_adjusted"     = alpha_Sarstedt
    )
  }  
  
  return(out)
  
}