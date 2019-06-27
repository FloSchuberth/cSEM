#' Estimates the distance between multiple groups. 
#'
#' 
#' The geodesic distance (dG) and the Euclidean distance (dE) is used
#' to assess group differences. Permutation is used to generate a 
#' reference distribtuion. 
#' 
#' @usage testMGD(
#'  .object               = args_default()$.object,
#'  .approach_mgd          = args_default()$.approach_mgd,
#'  .alpha                = args_default()$.alpha,
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
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){
  
  # Implementation is based on:
  # Klesel & Schuberth (forthcoming) - (TODO) name
  
  if(.verbose) {
    cat(rule(center = "Test for multigroup differences based on Klesel et al. (forthcoming)",
             line = "bar3"), "\n\n")
  }
  UseMethod("testMGD")
  
}

#' @describeIn testMGD (TODO)
#' @export

testMGD.cSEMResults_default <- function(
  .object                = args_default()$.object,
  .approach_mgd          = args_default()$.approach_mgd,
  .alpha                 = args_default()$.alpha,
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
  .alpha                 = args_default()$.alpha,
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
    
    ### Calculation===============================================================
    ## Get the fitted values
    fit <- fit(.object, .saturated = .saturated, .type_vcv = .type_vcv)
    
    ## Compute the test statistics
    teststat <- c(
      "dG" = calculateDistance(.matrices = fit, .distance = "geodesic"),
      "dL" = calculateDistance(.matrices = fit, .distance = "squared_euclidian")
    )
    
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
    ref_dist         <- list()
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
        fit_temp <- fit(Est_temp, .saturated = .saturated, .type_vcv = .type_vcv)
        ref_dist[[counter]] <- c(
          "dG" = calculateDistance(.matrices = fit_temp, .distance = "geodesic"),
          "dL" = calculateDistance(.matrices = fit_temp, .distance = "squared_euclidian")
        )
        
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
    
    # Combine
    ref_dist_matrix <- do.call(cbind, ref_dist1)
    
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
        "Number_admissibles"    = ncol(ref_dist_matrix),
        "Total_runs"            = counter + n_inadmissibles,
        "Group_names"           = names(.object),
        "Number_of_observations"= sapply(X_all_list, nrow),
        "Bootstrap_values"      = ref_dist
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
  .alpha                 = args_default()$.alpha,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R                     = args_default()$.R,
  .saturated             = args_default()$.saturated,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){
  stop2("Currently, second-order models are not supported by `testMGD()`.")
}
