#' Tests for multiple groups. 
#'
#' 
#' The geodesic distance (dG) and the Euclidean distance (dE) is used
#' to assess group differences. Permutation is used to generate a 
#' reference distribtuion. 
#' 
#' @usage testMGD(
#'  .object               = args_default()$.object,
#'  .approach_mgd          = args_default()$.approach_mgd,
#'  .comparison            = args_default()$comparison,
#'  .alpha                = args_default()$.alpha,
#'  .alpha_adjustment     = args_default()$.alpha_adjustment,
#'  .handle_inadmissibles = args_default()$.handle_inadmissibles,
#'  .R                    = args_default()$.R,
#'  .saturated            = args_default()$.saturated,
#'  .type_vcv             = args_default()$.type_vcv,
#'  .verbose              = args_default()$.verbose
#'  ) 
#' 
#' @inheritParams csem_arguments
#' 
#' @inherit csem_test return
#'
#' @seealso [cSEMResults]
#'
#' @examples
#' \dontrun{
#' require(cSEM)
#' data(satisfaction)
#'
#' model <- "
#' # Structural model
#' QUAL ~ EXPE
#' EXPE ~ IMAG
#' SAT  ~ IMAG + EXPE + QUAL + VAL
#' LOY  ~ IMAG + SAT
#' VAL  ~ EXPE + QUAL
#'
#' # Measurement model
#'
#' EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
#' IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
#' LOY  =~ loy1  + loy2  + loy3  + loy4
#' QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
#' SAT  <~ sat1  + sat2  + sat3  + sat4
#' VAL  <~ val1  + val2  + val3  + val4
#' "
#' 
#' listData <- list(satisfaction[-3,], satisfaction[-5, ], satisfaction[-10, ])
#' out.cSEM <- csem(listData, model) 
#'
#' testMGD(.object = out.cSEM, .R = 20, .type_vcv= 'construct')
#' }
#'
#' @export

testMGD <- function(
  .object                = args_default()$.object,
  .approach_mgd          = args_default()$.approach_mgd,
  .comparison            = args_default()$comparison,
  .alpha                 = args_default()$.alpha,
  .alpha_adjustment     = args_default()$.alpha_adjustment,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){
  
  if(.verbose) {  
    if(.approach_mgd == "Klesel"){
        # Implementation is based on:
        # Klesel et al. (2019) - (TODO) name
      cat(rule(center = "Test for multigroup differences based on Klesel et al. (2019)",
               line = "bar3"), "\n\n")
    
    }else if(.approach_mgd == "Chin"){
      # Implementation is based on:
      # Chin & Dibbern (2010) - (TODO) name
      cat(rule(center = "Test for multigroup differences based on Chin & Dibbern (2010)",
               line = "bar3"), "\n\n")
    }else if(.approach_mgd == "Sarstedt"){
      # Implementation is based on:
      # Sarstedt et al. (2011) - (TODO) name
      cat(rule(center = "Test for multigroup differences based on Sarstedt et al. (2011)",
               line = "bar3"), "\n\n")
    }
  }
  UseMethod("testMGD")
  
}

#' @describeIn testMGD (TODO)
#' @export

testMGD.cSEMResults_default <- function(
  .object                = args_default()$.object,
  .approach_mgd          = args_default()$.approach_mgd,
  .comparison            = args_default()$comparison,
  .alpha                 = args_default()$.alpha,
  .alpha_adjustment     = args_default()$.alpha_adjustment,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
) { 
  stop("At least two groups required.", call. = FALSE)
}

#' @describeIn testMGD (TODO)
#' @export

testMGD.cSEMResults_multi <- function(
  .object                = args_default()$.object,
  .approach_mgd          = args_default()$.approach_mgd,
  .comparison            = args_default()$comparison,
  .alpha                 = args_default()$.alpha,
  .alpha_adjustment     = args_default()$.alpha_adjustment,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){
  if(inherits(.object, "cSEMResults_2ndorder")) {
    
    out <- lapply(.object, testMICOM.cSEMResults_2ndorder)
    
  } else {
    
    ### Checks and errors ========================================================
    match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
    match.arg(.type_vcv, args_default(.choices = TRUE)$.type_vcv)
    
    # Check if any of the group estimates are inadmissible
    if(sum(unlist(verify(.object))) != 0) {
      stop("Initial estimation results for at least one group are inadmissible.\n", 
           "See `verify(.object)` for details.",  call. = FALSE)
    }
    
    # Check if data for different groups is identical
    if(.verbose) {
      if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                          function(x){ identical(x[[1]], x[[2]])})){
        warning("At least two groups are identical.", call. = FALSE)
      } 
    }
    
    

    ### Calculation of the test statistics======================================
    

    ## Get the model-implied VCV
    fit <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)
    
    ## Compute the test statistics 
    teststat <- list(
      #  Approach suggested by Klesel et al. (2019)
      Klesel=c(
      "dG" = cSEM:::calculateDistance(.matrices = fit, .distance = "geodesic"),
      "dL" = cSEM:::calculateDistance(.matrices = fit, .distance = "squared_euclidian")
    ),
    # Approach suggested by Chin & Dibbern (2010)
    Chin=cSEM:::parameter_difference(.object=.object,.comparison = .comparison))
      
      


    ## Start Permutation
    # Put data of each groups in a list and combine
    X_all_list  <- lapply(.object, function(x) x$Information$Data)
    X_all       <- do.call(rbind, X_all_list)
    
    # Collect initial arguments (from the first object, but could be any other)
    arguments <- .object[[1]]$Information$Arguments
    
    # Create a vector "id" to be used to randomly select groups (permutate) and
    # set id as an argument in order to identify the groups.
    id <- rep(1:length(X_all_list), sapply(X_all_list, nrow))
    arguments[[".id"]] <- "id"
    
    # Start progress bar if required
    if(.verbose){
      pb <- txtProgressBar(min = 0, max = .R, style = 3)
    }
    
    ## Calculate reference distribution
    ref_dist        <- list()
    n_inadmissibles  <- 0
    counter <- 0
    repeat{
      # Counter
      counter <- counter + 1
      
      # Permutate data
      X_temp <- cbind(X_all, id = sample(id))
      
      # Replace the old dataset by the new permutated dataset
      arguments[[".data"]] <- X_temp
      
      # Estimate model
      Est_temp <- do.call(csem, arguments)   
      
      # Check status
      status_code <- sum(unlist(verify(Est_temp)))
      
      # Distinguish depending on how inadmissibles should be handled
      if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
        # Compute if status is ok or .handle inadmissibles = "ignore" AND the status is 
        # not ok

        # Calculate test statistic for permutation sample

        fit_temp <- fit(Est_temp, .saturated = .saturated, .type_vcv = .type_vcv)
        
        ref_dist[[counter]] <- list(Klesel=c(
          "dG" = cSEM:::calculateDistance(.matrices = fit_temp, .distance = "geodesic"),
          "dL" = cSEM:::calculateDistance(.matrices = fit_temp, .distance = "squared_euclidian")
        ),Chin=cSEM:::parameter_difference(.object=Est_temp,.comparison = .comparison))

      # Approach suggested by Chin & Dibbern
          # Calculate differences based on the estimation based on the permutation sample
          # ref_dist[[counter]] <- parameter_difference(.object=Est_temp,.comparison = .comparison)
          
          
      } else if(status_code != 0 & .handle_inadmissibles == "drop") {
        # Set list element to zero if status is not okay and .handle_inadmissibles == "drop"
        ref_dist[[counter]] <- NA
        
      } else {# status is not ok and .handle_inadmissibles == "replace"
        # Reset counter and raise number of inadmissibles by 1
        counter <- counter - 1
        n_inadmissibles <- n_inadmissibles + 1
      }
      
      # Break repeat loop if .R results have been created.
      if(length(ref_dist) == .R) {
        break
      } else if(counter + n_inadmissibles == 10000) { 
        ## Stop if 10000 runs did not result in insufficient admissible results
        stop("Not enough admissible result.", call. = FALSE)
      }
      
      if(.verbose){
        setTxtProgressBar(pb, counter)
      }
      
    } # END repeat 
    
    # close progress bar
    if(.verbose){
      close(pb)
    }
    
    # Delete potential NA's
    ref_dist1 <- Filter(Negate(anyNA), ref_dist)
    
    
    # Appraoch suggest by Klesel et al. (2019)
    ref_dist_Klesel = lapply(ref_dist1,function(x){
      x$Klesel
    })
    
    # Combine
    ref_dist_matrix_Klesel <- do.call(cbind, ref_dist_Klesel)
    
    
    ## Compute critical values (Result is a (2 x p) matrix, where n is the number
    ## of quantiles that have been computed (1 by default)
    .alpha <- .alpha[order(.alpha)]
    critical_values_Klesel <- matrixStats::rowQuantiles(ref_dist_matrix_Klesel, 
                                                 probs =  1-.alpha, drop = FALSE)
    teststat_Klesel <- teststat$Klesel
    ## Compare critical value and teststatistic
    decision_Klesel <- teststat_Klesel < critical_values_Klesel # a logical (2 x p) matrix with each column
    # representing the decision for one
    # significance level. TRUE = no evidence 
    # against the H0 --> not reject
    # FALSE --> reject
    
    # Approach suggested by Chin & Dibbern (2010)
    teststat_Chin = lapply(teststat$Chin, function(x){
      x[!is.na(x)]
    })
    
    ref_dist_Chin = lapply(ref_dist1,function(x){
      x$Chin
    })
    
    ref_dist_Chin_temp=purrr::transpose(ref_dist_Chin)
    
    # Outcome of the following needs to be matrix
    ref_dist_matrices_Chin =lapply(ref_dist_Chin_temp,function(x){
      temp = do.call(cbind,x)
      temp_ind = stats::complete.cases(temp)
      temp[temp_ind,,drop=FALSE]
      }
      )
    
    # Here we can correct alpha for multiple test
    alpha_Chin=.alpha
    

    critical_values_Chin = lapply(alpha_Chin,function(alpha){
      probs_Chin = c(alpha/2,1-alpha/2)
      temp=lapply(ref_dist_matrices_Chin, function(x){
        matrixStats::rowQuantiles(x, probs = probs_Chin, drop = FALSE)
        })
    })
    
    names(critical_values_Chin) = paste0(alpha_Chin*100,'%')

    decision_Chin <- mapply(function(teststat,critical){
      temp=lapply(critical, function(x){
        x[,1]< teststat & teststat < x[,2]

        })
    },teststat = teststat_Chin, critical = critical_values_Chin,SIMPLIFY = FALSE)
    
    names(decision_Chin) = paste0(alpha_Chin*100,'%')
    
    # Return output
    out <- list(
      "Klesel"=list(
      "Test_statistic"     = teststat_Klesel,
      "Critical_value"     = critical_values_Klesel, 
      "Decision"           = decision_Klesel), 

      "Chin" = list(
        "Test_statistic"     = teststat_Chin,
        "Critical_value"     = critical_values_Chin, 
        "Decision"           = decision_Chin,
        "Alpha adjusted"     =  NULL),
    "Information"        = list(
      "Number_admissibles"    = ncol(ref_dist_matrix_Klesel),
      "Total_runs"            = counter + n_inadmissibles,
      "Group_names"           = names(.object),
      "Number_of_observations"= sapply(X_all_list, nrow),
      "Bootstrap_values"      = list(
        "Klesel" = ref_dist_Klesel,
        "Chin"   = ref_dist_matrices_Chin),
      "Approach" = .approach_mgd
      )
    )
  }
  
  class(out) <- "cSEMTestMGD"
  return(out)
}

#' @describeIn testMGD (TODO)
#' @export

testMGD.cSEMResults_2ndorder <- function(
  .object                = args_default()$.object,
  .approach_mgd          = args_default()$.approach_mgd,
  .comparison            = args_default()$comparison,
  .alpha                 = args_default()$.alpha,
  .alpha_adjustment     = args_default()$.alpha_adjustment,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){
  stop2("Currently, second-order models are not supported by `testMGD()`.")
}
