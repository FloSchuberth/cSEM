#' Test measurement invariance of composites
#'
#' This functions performs the test for measurement invariance of composites
#' proposed by Henseler et al. (2016).
#'
#' The test is only meaningful for composite models.
#'
#' If more than two groups are to be compared issues related to multiple testing
#' should be taken into account.
#'
#' The number of permutation runs is defaults to `args_default()$.runs` for performance reasons.
#' According to Henseler et al. (2016) the number of permutations should be at least 5000 for
#' assessment to be reliable.
#'
#' @usage testMICOM(
#'  .object               = args_default()$.object,
#'  .alpha                = args_default()$.alpha,
#'  .handle_inadmissibles = args_default()$.handle_inadmissibles,
#'  .runs                 = args_default()$.runs,
#'  .verbose              = args_default()$.verbose
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @return An object of class "testMICOM"
#'
#' An object of class "testMICOM" is a named list containing the following list element:
#' \describe{
#' \item{**Step2**}{A list containing the results of the test for compositional invariance (Step 2).}
#' \item{**Step3**}{A list containing the results of the test for mean and variance equality (Step 3).}
#' \item{**Meta_information**}{A list of additional information on the test.}
#' }
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export
#'

testMICOM <- function(
  .object               = args_default()$.object,
  .alpha                = args_default()$.alpha,
  .handle_inadmissibles = args_default()$.handle_inadmissibles,
  .runs                 = args_default()$.runs,
  .verbose              = args_default()$.verbose
) {
  # Implementation is based on:
  # Henseler et al. (2016) - Testing measurement invariance of composites using
  #                          partial least squares
  
  if(.verbose) {
    cat(rule(center = "Test for measurement invariance based on Henseler et al. (2016)",
             line = "bar3"), "\n\n")
  }
  UseMethod("testMICOM")

}

#' @describeIn testMICOM (TODO)
#' @export

testMICOM.cSEMResults_default <- function(.object = args_default()$.object) {
  stop("At least 2 groups required for the MICOM test.", call. = FALSE)
}

#' @describeIn testMICOM (TODO)
#' @export

testMICOM.cSEMResults_multi <- function(
  .object               = args_default()$.object,
  .alpha                = args_default()$.alpha,
  .handle_inadmissibles = args_default()$.handle_inadmissibles,
  .runs                 = args_default()$.runs,
  .verbose              = args_default()$.verbose
) {
  ### Checks and errors ========================================================
  if(sum(unlist(verify(.object))) != 0) {
    stop("Initial estimation results for at least one group are inadmissible.\n", 
         "See `verify(.object)` for details.",  call. = FALSE)
  }
  
  if(.verbose & all(.object[[1]]$Information$Model$construct_type == "Common factor")) {
    warning("\tAll constructs are modelled as common factors.\n",
            "\tTest results are only meaningful for composite models!",
            call. = FALSE)
  }

  if(.verbose & length(.object) > 2) {
    warning("\tComparing more than two groups inflates the familywise error rate.\n",
            "\tInterpret statistical significance with caution.\n",
            "\n\tFuture versions of the package will likely include\n",
            "\tappropriate correction options.",
            call. = FALSE)
  }
  
  ### Preparation ==============================================================
  ## Get data (pooled, potentially unstandardized data)
  X <- .object$Data_1$Information$Data_pooled 
  
  # Collect initial arguments (from the first object, but could be any other)
  arguments <- .object[[1]]$Information$Arguments
  
  # Create a vector "id" to be used to randomly select groups (permutate) and
  # set id as an argument in order to identify the groups.
  X_list <- lapply(.object, function(x) x$Information$Data)
  id <- rep(1:length(X_list), sapply(X_list, nrow))
  arguments[[".id"]] <- "id"
  
  ### Step 1 - Configural invariance ===========================================

  # Has to be assessed by the user prior to using the testMICOM function. See
  # the original paper for details.

  ### Step 2 - Compositional invariance ========================================
  # Procedure: see page 414 and 415 of the paper
  ## Compute proxies/scores using the pooled data 
  # (it does not matter if the data is scaled or unscaled as this does not 
  # affect the correlation)

  H <- lapply(.object, function(x) X %*% t(x$Estimates$Weight_estimates))

  ## Compute the correlation of the scores for all group combinations
  # Get the scores for all group combinations
  H_combn <- utils::combn(H, 2, simplify = FALSE)

  # Compute the correlation c for each group combination
  c <- lapply(H_combn, function(x) diag(cor(x[[1]], x[[2]])))
  
  # Set the names for each group combination
  names(c) <- utils::combn(names(.object), 2, FUN = paste0, 
                           collapse = "_", simplify = FALSE)
 
## Permutation ---------------------------------------------------------------
  # Start progress bar
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
    
    # Permutate data
    X_temp <- cbind(X, id = sample(id))
    
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
      
      ## Compute weights for each group and use these to compute proxies/scores using
      # the pooled data (= the original combined data). Note that these
      # scores are unstandardized, however since we consider the correlation 
      # it does not matter whether we consider standardized or unstandardized 
      # proxies
      H_temp <- lapply(Est_temp, function(x) X %*% t(x$Estimates$Weight_estimates))
      
      ## Compute the correlation of the scores for all group combinations
      # Get the scores for all group combinations
      H_combn_temp <- utils::combn(H_temp, 2, simplify = FALSE)
      
      # Compute the correlation c for each group combination
      c_temp <- lapply(H_combn_temp, function(x) diag(cor(x[[1]], x[[2]])))
      
      # Set the names for each group combination
      names(c_temp) <- utils::combn(names(H_temp), 2, FUN = paste0, 
                                    collapse = "_", simplify = FALSE)
      
      ref_dist[[i]] <- c_temp
      
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
  
  ## Bring data to form and compute quantiles
  # Delete potential NULL entries
  ref_dist <- Filter(Negate(is.null), ref_dist)
  # Bind 
  temp <- do.call(rbind, lapply(ref_dist, function(x) do.call(rbind, x)))
  temp <- split(as.data.frame(temp), rownames(temp))
  
  critical_values_step2 <- lapply(lapply(temp, as.matrix), matrixStats::colQuantiles, 
                      probs = 1-.alpha, drop = FALSE)
  
  ## Decision
  decision <- mapply(function(x, y) x < y,
                     x = c,
                     y = critical_values_step2,
                     SIMPLIFY = FALSE)
  
  # if(.verbose) {
  #   cat("\nCompositional invariance test finished.\n",
  #       "Proceding to test equality of composite mean values and variances.\n", sep = "")
  # }
  
  ### Step 3 - Equal mean values and variances==================================
  # Update arguments
  arguments[[".data"]] <- X
  arguments[[".id"]]   <- NULL

  
  # Estimate model using pooled data set,
  Est <- do.call(csem, arguments)
  
  # Procedure below:
  # 1. Create list of original group id's + .runs permutated id's
  # 2. Extract construct scores, attach ids and convert to data.frame (for split)
  # 3. Split construct scores by its id column
  # --- Now the structure is a list of length 1 + .runs each containing a list
  #     of the same length as there are numbers of groups (often just 2).
  # 4. Convert group data set to a matrix and delete id column
  # 5. Compute the mean and the variance for each replication (+ the original data of course)
  #    and data set for each group.
  # 6. Transpose list
  # --- Now we have a list of length 1 + .runs containing two list elements
  #     "Mean" and "Var" which in turn contain as many list elements as there are
  #     groups
  meanvar <- c(list(id), replicate(.runs, sample(id), simplify = FALSE)) %>% 
    lapply(function(x) as.data.frame(cbind(Est$Estimates$Construct_scores, id = x))) %>% 
    lapply(function(x) split(x, f = x$id)) %>% 
    lapply(function(x) lapply(x, function(y) as.matrix(y[, -ncol(y), drop = FALSE]))) %>% 
    lapply(function(x) lapply(x, function(y) list(
      "Mean" = colMeans(y), 
      "Var"  = {tt <- matrixStats::colVars(y); names(tt) <- colnames(y); tt})
    )) %>% 
    lapply(function(x) purrr::transpose(x))
  
  # Compute the difference of the means for each possible group combination 
  # and attach a name
  m <- meanvar %>% 
    lapply(function(x) {
      tt <- utils::combn(x[["Mean"]], m = 2, 
                   FUN = function(y) y[[1]] - y[[2]], 
                   simplify = FALSE)
      names(tt) <- utils::combn(names(.object), m = 2, FUN = paste0, 
                         collapse = "_", simplify = FALSE)
      tt
    })
  
  # Compute the log difference of the variances for each possible group combination
  # and attach a name
  v <- meanvar %>% 
    lapply(function(x) {
      tt <- utils::combn(x[["Var"]], m = 2, 
                         FUN = function(y) log(y[[1]]) - log(y[[2]]), 
                         simplify = FALSE)
      names(tt) <- utils::combn(names(.object), m = 2, FUN = paste0, 
                                collapse = "_", simplify = FALSE)
      tt
    })
  
  # Combine in a list and bind (log) differences for each replication in a list
  mv <- list("Mean" = m, "Var"= v) %>% 
    lapply(function(x) lapply(x, function(y) do.call(rbind, y))) %>% 
    lapply(function(x) cbind(as.data.frame(do.call(rbind, x), row.names = FALSE), 
                             id = rownames(do.call(rbind, x)))) %>% 
    lapply(function(x) split(x, f = x$id)) %>% 
    lapply(function(x) lapply(x, function(y) as.matrix(y[, -ncol(y), drop = FALSE])))
  
  # Extract the original group differences. This is the first row of each group
  # dataset
  mv_o <- lapply(mv, function(x) lapply(x, function(y) y[1, ]))

  # Compute quantiles/critical values
  critical_values_step3 <- lapply(mv, function(x) lapply(x, function(y) y[-1, ])) %>% 
    lapply(function(x) lapply(x, function(y) matrixStats::colQuantiles(y, probs =  1-.alpha, drop = FALSE)))

  ## Compare critical value and teststatistic
  # For Mean
  decision_m <- mapply(function(x, y) x < y,
                    x = mv_o[[1]],
                    y = critical_values_step3[[1]],
                    SIMPLIFY = FALSE)
  # For Var
  decision_v <- mapply(function(x, y) x < y,
                       x = mv_o[[2]],
                       y = critical_values_step3[[2]],
                       SIMPLIFY = FALSE)

  # Return output
  out <- list(
    "step2" = list(
      "Test_statistic"     = c,
      "Critical_value"     = critical_values_step2, 
      "Decision"           = decision 
    ),
    "step3" = list(
      "Mean" = list(
        "Test_statistic"     = mv_o$Mean,
        "Critical_value"     = critical_values_step3$Mean, 
        "Decision"           = decision_m
      ),
      "Var" = list(
        "Test_statistic"     = mv_o$Var,
        "Critical_value"     = critical_values_step3$Var, 
        "Decision"           = decision_v
      )
    ),
    "Information" = list(
      "Number_admissibles"    = length(ref_dist),
      "Total_runs"            = i + n_inadmissibles,
      "Group_names"           = names(.object),
      "Number_of_observations"= sapply(X_list, nrow)
    ) 
  )
  
  class(out) <- "cSEMTestMICOM"
  return(out)
}

#' @describeIn testMICOM (TODO)
#' @export

testMICOM.cSEMResults_2ndorder <- function(
  .object               = args_default()$.object,
  .alpha                = args_default()$.alpha,
  .handle_inadmissibles = args_default()$.handle_inadmissibles,
  .runs                 = args_default()$.runs,
  .verbose              = args_default()$.verbose
) {
  stop("Not yet implemented", call. = FALSE)
}
