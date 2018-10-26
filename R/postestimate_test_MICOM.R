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
    warning("\tComparing multiple groups inflates the familywise error rate.\n",
            "\tInterpret statistical significance with caution.\n",
            "\n\tFuture versions of the package will likely include\n",
            "\tappropriate correction options.",
            call. = FALSE)
  }
  
  ### Preparation ==============================================================
  # Put data of each groups in a list and combine
  X_all_list  <- lapply(.object, function(x) x$Information$Data)
  X_all       <- do.call(rbind, X_all_list)
  
  # Collect initial arguments (from the first object, but could be any other)
  arguments <- .object[[1]]$Information$Arguments
  
  # Create a vector "id" to be used to randomly select groups (permutate) and
  # set id as an argument in order to identify the groups.
  id <- rep(1:length(X_all_list), sapply(X_all_list, nrow))
  arguments[[".id"]] <- "id"
  
  ### Step 1 - Configural invariance ===========================================

  # Has to be assessed by the user prior to using the testMICOM function. See
  # the original paper for details.

  ### Step 2 - Compositional invariance ========================================
  # Procedure: see page 414 and 415 of the paper
  ## Compute weights for each group and use these to compute proxies/scores using
  # the pooled data

  H <- lapply(.object, function(x) X_all %*% t(x$Estimates$Weight_estimates))
  
  ## Compute the correlation of the scores for all group combinations
  # Get the scores for all group combinations
  H_combn <- utils::combn(H, 2, simplify = FALSE)

  # Compute the correlation c for each group combination
  c <- lapply(H_combn, function(x) diag(cor(x[[1]], x[[2]])))

  # Set the names for each group combination
  names(c) <- utils::combn(names(H), 2, FUN = paste0, collapse = "_", simplify = FALSE)

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
      
      ## Compute weights for each group and use these to compute proxies/scores using
      # the pooled data (= the original combined data)
      H_temp <- lapply(Est_temp, function(x) X_all %*% t(x$Estimates$Weight_estimates))
      
      ## Compute the correlation of the scores for all group combinations
      # Get the scores for all group combinations
      H_combn_temp <- utils::combn(H_temp, 2, simplify = FALSE)
      
      # Compute the correlation c for each group combination
      c_temp <- lapply(H_combn_temp, function(x) diag(cor(x[[1]], x[[2]])))
      
      # Set the names for each group combination
      names(c_temp) <- utils::combn(names(H_temp), 2, FUN = paste0, collapse = "_", simplify = FALSE)
      
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
  
  step2_out <- lapply(lapply(temp, as.matrix), matrixStats::colQuantiles, 
                      probs = 1-.alpha, drop = FALSE)
  step2_out <- mapply(function(x, y) cbind("c" = y, x),
                      x = step2_out,
                      y = c,
                      SIMPLIFY = FALSE)
  
  # if(.verbose) {
  #   cat("\nCompositional invariance test finished.\n",
  #       "Proceding to test equality of composite mean values and variances.\n", sep = "")
  # }
  
  ### Step 3 - Equal mean values and variances==================================
# 
#   H <- workhorse(
#     .data                    = data_ordered[, -which(names(data_ordered) == .group_var)],
#     .model                   = csem_model,
#     .approach_cf             = "dist_euclid",
#     .approach_nl             = NULL,
#     .approach_paths          = NULL,
#     .approach_weights        = approach_weights,
#     .disattenuate            = FALSE,
#     .estimate_structural     = FALSE,
#     .ignore_structural_model = FALSE,
#     .iter_max                = 100,
#     .normality               = TRUE,
#     .PLS_mode                = NULL,
#     .PLS_weight_scheme_inner = "centroid",
#     .tolerance               = 1e-06
#   )$Estimates$Proxies
# 
#   # H <- workhorse(.data  = data_ordered[, -which(names(data_ordered) == .group_var)],
#   #                .model = csem_model,
#   #                .approach_weights    = approach_weights,
#   #                .PLS_weight_scheme_inner = PLS_weight_scheme_inner,
#   #                .estimate_structural = FALSE,
#   #                ...)$Estimates$Proxies
# 
#   temp <- split(x = as.data.frame(H), rep(unique(data_ordered[, .group_var]),
#                                           times = table(data_ordered[, .group_var])))
#   temp_m <- lapply(temp, colMeans)
#   temp_v <- lapply(temp, function(x) {
#     a <- matrixStats::colVars(as.matrix(x))
#     names(a) <- colnames(x)
#     a
#   })
# 
#   # Make the combinations of list elements
#   temp_m <- combn(temp_m, 2, simplify = FALSE)
#   temp_v <- combn(temp_v, 2, simplify = FALSE)
# 
#   # Compute difference for each group combination
#   m <- lapply(temp_m, function(x) x[[1]] - x[[2]])
#   v <- lapply(temp_v, function(x) x[[1]] - x[[2]])
# 
#   # Get the names
#   names(m) <- names(v) <- combn(names(temp), 2, FUN = paste0, collapse = "_", simplify = FALSE)
# 
#   ### Permutation --------------------------------------------------------------
#   if(.verbose) {
#     cat("\nProgress:\n", sep = "")
#   }
# 
#   perm2 <- lapply(1:.runs, function(x) {
# 
#     if(.verbose & x %in% seq(0, .runs, by = 10)) {
#       cat(x, " / ", .runs, "\n", sep = "")
#     }
# 
#     H <- workhorse(
#       .data                    = data_no_group_var[sample(x = nrow(data_no_group_var)), ],
#       .model                   = csem_model,
#       .approach_cf             = "dist_euclid",
#       .approach_nl             = NULL,
#       .approach_paths          = NULL,
#       .approach_weights        = approach_weights,
#       .disattenuate            = FALSE,
#       .estimate_structural     = FALSE,
#       .ignore_structural_model = FALSE,
#       .iter_max                = 100,
#       .normality               = TRUE,
#       .PLS_mode                = NULL,
#       .PLS_weight_scheme_inner = "centroid",
#       .tolerance               = 1e-06
#     )$Estimates$Proxies
# 
#     temp <- split(x = as.data.frame(H), rep(unique(data_ordered[, .group_var]),
#                                             times = table(data_ordered[, .group_var])))
#     temp_m <- lapply(temp, colMeans)
#     temp_v <- lapply(temp, function(x) {
#       a <- matrixStats::colVars(as.matrix(x))
#       names(a) <- colnames(x)
#       a
#     })
# 
#     # Make the combinations of list elements
#     temp_m <- combn(temp_m, 2, simplify = FALSE)
#     temp_v <- combn(temp_v, 2, simplify = FALSE)
# 
#     # Compute difference for each group combination
#     m <- lapply(temp_m, function(x) x[[1]] - x[[2]])
#     v <- lapply(temp_v, function(x) log(x[[1]]) - log(x[[2]]))
# 
#     # Get the names
#     names(m) <- names(v) <- combn(names(temp), 2, FUN = paste0, collapse = "_", simplify = FALSE)
#     list("m" = m, "v" = v)
#   })
# 
#   temp_mv <- lapply(perm2, function(x) lapply(x, function(y) do.call(rbind, y)))
#   temp_mv <- lapply(purrr::transpose(temp_mv), function(x) do.call(rbind, x))
# 
#   temp_m <- split(as.data.frame(temp_mv[["m"]]), names(m))
#   temp_v <- split(as.data.frame(temp_mv[["v"]]), names(v))
# 
#   step3_out_m <- lapply(temp_m, function(x) as.data.frame(t(apply(x, 2, quantile, probs = c(.alpha/2, 1 - .alpha/2)))))
#   step3_out_m <- mapply(function(x, y) cbind("Diff_mean" = y, x),
#                       x = step3_out_m,
#                       y = m,
#                       SIMPLIFY = FALSE)
# 
#   step3_out_v <- lapply(temp_v, function(x) as.data.frame(t(apply(x, 2, quantile, probs = c(.alpha/2, 1 - .alpha/2)))))
#   step3_out_v <- mapply(function(x, y) cbind("Diff_log_var" = y, x),
#                         x = step3_out_v,
#                         y = v,
#                         SIMPLIFY = FALSE)
## Results ==================================================================
  out <- list("Step2"            = step2_out
              # ,
              # "Step3"            = list("Mean_diff" = step3_out_m,
              #                           "Var_diff"  = step3_out_v),
              # "Meta_information" = list("Number_of_observations" = no_obs,
              #                           "Number_of_Groups"       = nrow(no_obs) - 1,
              #                           "Grouping_variable"      = .group_var)
              )
  ## Set class
  class(out) <- "testMICOM"
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