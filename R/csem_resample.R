#' Resample data 
#'
#' Resample from a given cSEMResults object or data.frame using common
#' resampling methods.
#'
#' @usage resampleData(
#'  .object        = args_default()$.object,
#'  .data          = NULL,
#'  .method        = c("bootstrap", "jackknife", "permutation", "cross-validation"),
#'  .R             = 499,
#'  .id            = NULL,
#'  .cv_folds      = 10      
#' )
#'
#' @inheritParams csem_arguments
#' @param .data A `data.frame` or a `matrix` containing the raw data. Possible
#'   data column types or classes are: logical, numeric (double or integer), factor 
#'   (ordered and unordered) or a mix of several types. Data may include
#'   one character column whose name is given by `.id`. This column is assumed
#'   to contain group identifiers used to split the data.
#'   If `.data` is provided `.object` is ignored. Defaults to `NULL`.
#' @param .method Character string. The resampling method to use. One of: 
#'  "*bootstrap*", "*jackknife*", "*permutation*", or "*cross-validation*". 
#'  Defaults to "*bootstrap*".
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
resampleData <- function(
  .object   = args_default()$.object,
  .data     = args_default()$.data,
  .method   = args_default()$.method,
  .R        = args_default()$.R,
  .id       = args_default()$.id,
  .cv_folds = args_default()$.cv_folds
) {
  match.arg(.method, c("bootstrap", "jackknife", "permutation", "cross-validation"))

  ## Get data set
  if(is.null(.data)) {
    data <- as.data.frame(.object$Informatio$Data)
    id   <- .object$Information$Arguments$.id
  } else {
    
    ## Checks
    if(!any(class(.data) %in% c("data.frame", "matrix"))) {
      stop("Data must be provided as a `matrix`, a `data.frame`.", 
           " .data has class: ", 
           class(.data), call. = FALSE)
    }
    if(is.matrix(.data)) {
      data <- as.data.frame(.data)
    }
    
    ## Set id 
    id   <- .id
    
    if(!is.null(id)) {
      # extract id column
      id_col <- data[, .id] # a vector
      # data without id column
      data   <- data[, -which(colnames(data) == .id)] 
    }
  }
  
  ## Choose resampling method
  out <- switch (.method,
    "jackknife"   = {
      lapply(1:nrow(data), function(x) data[-x, ])
    },
    "bootstrap"   = {
      lapply(1:.R, function(x) {
        data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]})
    },
    "permutation" = {
      if(is.null(id)) {
        lapply(1:.R, function(x) {
          data[sample(1:nrow(data), size = nrow(data), replace = FALSE), ]
        })
      } else {
        lapply(1:.R, function(x) cbind(data, sample(id)))
      }
    },
    "cross-validation" = {
      if(is.null(.cv_folds)) {
        # LOOCV (Leave one out CV).
        lapply(1:nrow(data), function(x) data[x, ])
        
      } else {
        # k-fold cross-validation (=draw k samples of equal size.).
        # Note the last sample may contain less observations if equal sized
        # samples are not possible
        suppressWarnings(
          split(as.data.frame(data), rep(1:.cv_folds, 
                                         each = ceiling(nrow(data)/.cv_folds)))
        )
      }
    } # END cross-validation
  ) # END switch
  ## Return samples
  out
}

#' Resample cSEMResults
#'
#' Resample (TODO)
#'
#' @usage resampleResults(
#'  .object        = args_default()$.object,
#'  .data          = NULL,
#'  .method        = c("bootstrap", "jackknife", "permutation", "cross-validation"),
#'  .R         = 499,
#'  .id            = NULL,
#'  .cv_folds      = 10      
#' )
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
resamplecSEMResults <- function(
  .object                = args_default()$.object, 
  .method                = args_default()$.method,
  .R                     = args_default()$.R,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .verbose               = args_default()$.verbose,
  .user_funs             = args_default()$.user_funs,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles2 = args_default()$.handle_inadmissibles2,
  .method2               = args_default()$.method2
  ) {
  
  ## Check arguments
  match.arg(.method, args_default(.choices = TRUE)$.method)
  match.arg(.method2, args_default(.choices = TRUE)$.method2)
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
  match.arg(.handle_inadmissibles2, args_default(.choices = TRUE)$.handle_inadmissibles2)
  
  ## Has the object to use the data to resample from produced admissible results?
  if(sum(verify(.object)) != 0) {
    warning("Estimation based on the data in `.object` has produced",
            " inadmissible results.\n", 
            "This may be a sign that something is wrong.",
            call. = FALSE, immediate. = TRUE)
  }
  
  ## Check for the minimum number of necessary resamples
  if(.R < 3 | .R2 < 3) {
    stop2("The following error occured in the `resamplecSEMResults()` functions:\n",
         "At least 3 resamples required.")
  }
  
  ### Process original data ----------------------------------------------------
  Est_original  <- .object$Estimates
  Info_original <- .object$Information
  
  ## Apply user defined function if specified
  user_funs <- if(!is.null(.user_funs)) {
    if(is.function(.user_funs)) {
      c("User_fun" = .user_funs(.object))
    } else {
      x <- lapply(.user_funs, function(f) c(f(.object)))
      if(is.null(names(x))) {
        names(x) <- paste0("User_fun", 1:length(x))
      }
      x
    }
  }
  
  ## Build names for path coef, loadings, and weights
  # Path estimates
  names_p <- outer(rownames(Est_original$Path_estimates), 
                   colnames(Est_original$Path_estimates), 
                   FUN = function(x, y) paste(x, y, sep = " ~ "))
  # Loadings
  names_l <- rep(rownames(Est_original$Loading_estimates), 
                 times = rowSums(Info_original$Model$measurement))
  names_l <- paste0(names_l, " =~ ", colnames(Est_original$Loading_estimates))
  
  # Weights
  names_w <- rep(rownames(Est_original$Weight_estimates), 
                 times = rowSums(Info_original$Model$measurement))
  names_w <- paste0(names_w, " <~ ", colnames(Est_original$Weight_estimates))
  
  ## Vectorize quantities that are not already vectorized 
  ## and set those that we dont need to NULL.
  Est_original$Path_estimates    <- t(Est_original$Path_estimates)[t(Info_original$Model$structural) != 0]
  names(Est_original$Path_estimates) <- t(names_p)[t(Info_original$Model$structural) != 0]
  Est_original$Loading_estimates <- t(Est_original$Loading_estimates)[t(Info_original$Model$measurement) != 0]
  names(Est_original$Loading_estimates) <- names_l
  Est_original$Weight_estimates  <- t(Est_original$Weight_estimates)[t(Info_original$Model$measurement) != 0]
  names(Est_original$Weight_estimates) <- names_w
  Est_original$Inner_weight_estimates <- NULL
  Est_original$Construct_scores       <- NULL
  Est_original$Indicator_VCV          <- NULL
  Est_original$Proxy_VCV              <- NULL
  Est_original$Construct_VCV          <- NULL
  Est_original$Cross_loadings         <- NULL
  Est_original$Construct_reliabilities<- NULL
  Est_original$Correction_factors     <- NULL
  
  ## Add output of the user functions to Est_original
  if(!is.null(.user_funs)) {
    Est_original <- c(Est_original, user_funs)
  }
  
  ## Get argument list
  args <- .object$Information$Arguments
  
  ### Alternative
  if(.method == "jackknife") {
    
    resample_jack <- resampleData(.object, .method = "jackknife")

    ## Start progress bar if required
    if(.verbose){
      pb <- txtProgressBar(min = 0, max = length(resample_jack), style = 3)
    }
    
    ## Substitute 'replace' by 'drop' if necessary
    ## .handle_inadbmissibles  = "replace" cannot be combined with "jackknife"
    if(.handle_inadmissibles == "replace") {
      warning("`'replace'` is set to `'.drop'` for jackknife resampling. ",
              call. = FALSE, immediate. = TRUE)
      .handle_inadmissibles <- "drop"
    }
    
    Est_ls <- lapply(seq_along(resample_jack), function(i) {
      
      # Replace the old dataset by a resampled data set (resampleData always returns
      # a list so for just one draw we need to pick the first list element)
      args[[".data"]] <- resample_jack[[i]]
      
      # Estimate model
      Est_temp <- do.call(csem, args)   
      
      # Check status
      status_code <- sum(unlist(verify(Est_temp)))
      
      # Distinguish depending on how inadmissibles should be handled
      if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
        
        ## Apply user defined function if specified
        user_funs <- if(!is.null(.user_funs)) {
          if(is.function(.user_funs)) {
            c("User_fun" = .user_funs(Est_temp))
          } else {
            x <- lapply(.user_funs, function(f) c(f(Est_temp)))
            if(is.null(names(x))) {
              names(x) <- paste0("User_fun", 1:length(x))
            }
            x
          }
        }
        
        ## Process
        x1  <- Est_temp$Estimates
        x2  <- Est_temp$Information
        
        x1$Path_estimates    <- t(x1$Path_estimates)[t(x2$Model$structural) != 0]
        x1$Loading_estimates <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0]
        x1$Weight_estimates  <- t(x1$Weight_estimates)[t(x2$Model$measurement) != 0]
        x1$Inner_weight_estimates  <- NULL
        x1$Construct_scores        <- NULL
        x1$Indicator_VCV           <- NULL
        x1$Proxy_VCV               <- NULL
        x1$Construct_VCV           <- NULL
        x1$Cross_loadings          <- NULL
        x1$Construct_reliabilities <- NULL
        x1$Correction_factors      <- NULL
        
        if(!is.null(.user_funs)) {
          x1 <- c(x1, user_funs)
        }
        
        ## If resampling from a bootstrap sample is required (e.g. for the 
        ## bootstraped t-interval CI) the original object needs to be returned. 
        if(.method2 != "none") {
          # Sometimes both jackknife and bootstrap resamples are required
          if(.method2 == "both") {
            Est_resamples2 <- lapply(list("bootstrap", "jackknife"), function(x) {
              resamplecSEMResults(
                .object               = Est_temp,
                .R                    = .R2,
                .handle_inadmissibles = .handle_inadmissibles2,
                .method               = x,
                .verbose              = FALSE,
                .user_funs            = .user_funs,
              )
            })
            x1 <- list("Estimates1" = x1, 
                       "Estimates2" = Est_resamples2[[1]],
                       "Estimates3" = Est_resamples2[[2]])
          } else {
            Est_resamples2 <- resamplecSEMResults(
              .object               = Est_temp,
              .R                    = .R2,
              .handle_inadmissibles = .handle_inadmissibles2,
              .method               = .method2,
              .verbose              = FALSE,
              .user_funs            = .user_funs,
            )
            x1 <- list("Estimates1" = x1, "Estimates2" = Est_resamples2)
          }
        } # END if .R2
      } else if(status_code != 0 & .handle_inadmissibles == "drop") {
        # Retrun NA if status is not okay and .handle_inadmissibles == "drop"
        x1 <- NA
      } 
      
      if(.verbose){
        setTxtProgressBar(pb, i)
      }
      
      ## Return
      x1
    }) # END lapply(seq_along(resample_jack), ...)
    # END if .method = .jackknife
  } else { # BEGIN if "bootstrap"
    
    ## Start progress bar if required
    if(.verbose){
      pb <- txtProgressBar(min = 0, max = .R, style = 3)
    }
    
    Est_ls          <- list()
    n_inadmissibles <- 0
    counter         <- 0
    repeat{
      # Counter
      counter <- counter + 1
      
      # Replace the old dataset by a resampled data set (resampleData always returns
      # a list so for just one draw we need to pick the first list element)
      args[[".data"]] <- resampleData(.object, .method = .method, .R = 1)[[1]]
      
      # Estimate model
      Est_temp <- do.call(csem, args)   
      
      # Check status
      status_code <- sum(unlist(verify(Est_temp)))
      
      # Distinguish depending on how inadmissibles should be handled
      if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
        
        ## Apply user defined function if specified
        user_funs <- if(!is.null(.user_funs)) {
          if(is.function(.user_funs)) {
            c("User_fun" = .user_funs(Est_temp))
          } else {
            x <- lapply(.user_funs, function(f) c(f(Est_temp)))
            if(is.null(names(x))) {
              names(x) <- paste0("User_fun", 1:length(x))
            }
            x
          }
        }
        
        ## Process
        x1  <- Est_temp$Estimates
        x2  <- Est_temp$Information
        
        x1$Path_estimates    <- t(x1$Path_estimates)[t(x2$Model$structural) != 0]
        x1$Loading_estimates <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0]
        x1$Weight_estimates  <- t(x1$Weight_estimates)[t(x2$Model$measurement) != 0]
        x1$Inner_weight_estimates  <- NULL
        x1$Construct_scores        <- NULL
        x1$Indicator_VCV           <- NULL
        x1$Proxy_VCV               <- NULL
        x1$Construct_VCV           <- NULL
        x1$Cross_loadings          <- NULL
        x1$Construct_reliabilities <- NULL
        x1$Correction_factors      <- NULL
        
        if(!is.null(.user_funs)) {
          x1 <- c(x1, user_funs)
        }
        
        ## If resampling from a bootstrap sample is required (e.g. for the 
        ## bootstraped t-interval CI) the original object needs to be returned. 
        if(.method2 != "none") {
          # Sometimes both jackknife and bootstrap resamples are required
          if(.method2 == "both") {
            Est_resamples2 <- lapply(list("bootstrap", "jackknife"), function(x) {
              resamplecSEMResults(
                .object               = Est_temp,
                .R                    = .R2,
                .handle_inadmissibles = .handle_inadmissibles2,
                .method               = x,
                .verbose              = FALSE,
                .user_funs            = .user_funs,
              )
            })
            x1 <- list("Estimates1" = x1, 
                       "Estimates2" = Est_resamples2[[1]],
                       "Estimates3" = Est_resamples2[[2]])
          } else {
            Est_resamples2 <- resamplecSEMResults(
              .object               = Est_temp,
              .R                    = .R2,
              .handle_inadmissibles = .handle_inadmissibles2,
              .method               = .method2,
              .verbose              = FALSE,
              .user_funs            = .user_funs,
            )
            x1 <- list("Estimates1" = x1, "Estimates2" = Est_resamples2)
          }
        }
        
        Est_ls[[counter]] <- x1
        
      } else if(status_code != 0 & .handle_inadmissibles == "drop") {
        # Set list element to NA if status is not okay and .handle_inadmissibles == "drop"
        Est_ls[[counter]] <- NA
        
      } else {# status is not ok and .handle_inadmissibles == "replace"
        # Reset counter and raise number of inadmissibles by 1
        counter <- counter - 1
        n_inadmissibles <- n_inadmissibles + 1
      }
      
      # Break repeat loop if .R results have been created.
      if(length(Est_ls) == .R) {
        break
      } else if(counter + n_inadmissibles == 10000) { 
        ## Stop if 10000 runs did not result in insufficient admissible results
        stop("Not enough admissible result.", call. = FALSE)
      }
      
      if(.verbose){
        setTxtProgressBar(pb, counter)
      }
      
    } # END repeat 
  } # END if "bootstrap"

  # Close progress bar
  if(.verbose){
    close(pb)
  }
  
  ## Process data --------------------------------------------------------------
  # Delete potential NA's
  out <- Filter(Negate(anyNA), Est_ls)
  
  # Check if at least 3 admissible results were obtained
  n_admissibles <- length(out)
  if(n_admissibles < 3) {
    stop("The following error occured in the `resamplecSEMResults()` functions:\n",
         "Less than 2 admissible results produced.", 
         " Consider setting .handle_inadmissibles(2) == 'replace' instead.",
         call. = FALSE)
  }
  # Turn list "inside out" and bind bootstrap samples to matrix 
  # columns are variables
  # rows are bootstrap runs
  out <- purrr::transpose(out) 
  
  if(.method2 != "none") {
    if(.method2 == "both") {
      out_3 <- out$Estimates3
    }
    out_2 <- out$Estimates2
    out   <- purrr::transpose(out$Estimates1)
  }

  out <- out %>% 
    lapply(function(x) do.call(rbind, x))
 
  # Add estimated quantities based on the the original sample/data
  out <- mapply(function(l1, l2) list("Resampled" = l1, "Original" = l2),
                l1 = out,
                l2 = Est_original,
                SIMPLIFY = FALSE)
  
  ## Return --------------------------------------------------------------------
  # The desired output is: a list of two:
  #  1. Resamples := a list of two
  #      1.1 Estimates1 := the .R resamples of the "outer" resampling run.
  #      1.2 Estimates2 := the .R2 resamples of the "inner" run or NA
  #  2. Information := a list of three
  #      2.1 Method  := the method used to obtain the "outer" resamples
  #      2.2 Method2 := the method used to obtain the "inner" resamples or NA
  #      2.3 Number_of_admissibles  := the number of admissibles of the "outer" run
  #      2.4 Number_of_observations := the number of observations.
  #
  # Since resamplecSEMResults is called recursivly within its body it is 
  # difficult to produce the desired output. There is now straightforward way 
  # to define if a call is an recursive call or not
  # My workaround for now is to simply check if the argument provided for ".object"
  # is "Est_temp" since this is the argument for the recurive call. This
  # is of course not general but a dirty workaround. I guess it is save enough
  # though.
  is_recursive_call <- eval.parent(as.list(match.call()))$.object == "Est_temp"
  
  ## Return
  if(is_recursive_call) {
    out
  } else {
    if(.method2 != "none") {
      if(.method2 == "both") {
        out <- list("Estimates1" = out, "Estimates2" = out_2, "Estimates3" = out_3)
        method2a <- "bootstrap"
        method2b <- "jackknife"
      } else {
        out <- list("Estimates1" = out, "Estimates2" = out_2, "Estimates3" = NA)
        method2a <- .method2
        method2b <- NA
      }
    } else {
      out <- list("Estimates1" = out, "Estimates2" = NA, "Estimates3" = NA)
      method2a <- method2b <- NA
    }
    
    out <- list(
      "Resamples" = out,
      "Information" = list(
        "Method"                  = .method,
        "Method2a"                = method2a,
        "Method2b"                = method2b,
        "Number_of_observations"  = nrow(.object$Information$Data),
        "Number_admissibles"      = n_admissibles
      )
    )
    
    ## Set class
    class(out) <- "cSEMResults_resampled"
    return(out)
  }
}

#' Internal: Inference
#'
#' Compute quantities for inference (TODO)
#'
#' @usage infer(
#'  .object        = args_default()$.object
#' )
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

infer <- function(
  .csem_resample   = args_default()$.csem_resample,
  .alpha           = args_default()$.alpha,
  .bias_corrected  = args_default()$.bias_corrected,
  .statistic       = args_default()$.statistic
) {
  
  resamples1 <- .csem_resample$Resamples$Estimates1
  resamples2 <- .csem_resample$Resamples$Estimates2 # jackknife or bootstrap
  resamples3 <- .csem_resample$Resamples$Estimates3 # always jackknife
  info       <- .csem_resample$Information
  
  ## Compute quantiles/critical values -----------------------------------------
  probs <- c()
  for(i in seq_along(.alpha[order(.alpha)])) { 
    probs <- c(probs, round(.alpha[i]/2, 4), round(1 - .alpha[i]/2, 4)) 
  }

  ## Compute statistics and quantities
  l <- list(
    "Mean" = MeanResample(resamples1),
    "Var"  = VarResample(.x      = resamples1, 
                         .method = info$Method, 
                         .n      = info$Number_of_observations
                         ),
    "Sd"   = SdResample(.x = resamples1, 
                        .method = info$Method, 
                        .n = info$Number_of_observations
                        ),
    "Bias" = BiasResample(.x = resamples1, 
                          .method = info$Method, 
                          .n = info$Number_of_observations
                          ),
    "CI_standard_z" = StandardCIResample(
      .x              = resamples1, 
      .bias_corrected = .bias_corrected,
      .df             = NULL,
      .dist           = "z",
      .method         = info$Method, 
      .n              = info$Number_of_observations,
      .probs          = probs
      ),
    "CI_standard_t" = StandardCIResample(
      .x              = resamples1, 
      .bias_corrected = .bias_corrected,
      .df             = "type1",
      .dist           = "t",
      .method         = info$Method, 
      .n              = info$Number_of_observations,
      .probs          = probs
    ),
    "CI_percentile" = PercentilCIResample(.x = resamples1, .probs = probs),
    "CI_t_invertall"= if(anyNA(resamples2)) {
      NA
    } else {
      TStatCIResample(
        .x1      = resamples1, 
        .x2      = resamples2, 
        .bias_corrected = .bias_corrected,
        .method  = info$Method, 
        .method2 = info$Method2a, 
        .n       = info$Number_of_observations, 
        .probs   = probs
        ) 
    }
  )
  if(anyNA(resamples2)) {
    l <- l[-length(l)]
  } else {
    l
  }
  return(purrr::transpose(l))
}

#' @describeIn infer The bootstraped mean
MeanResample <- function(.x) {
  
  lapply(.x, function(x) {
    out        <- colMeans(x$Resampled)
    names(out) <- names(x$Original)
    out
  })
}
#' @describeIn infer (TODO)
VarResample <- function(.x, .method, .n) {

  lapply(.x, function(x) {
    out        <- matrixStats::colVars(x$Resampled)
    names(out) <- names(x$Original)
    if(.method == "jackknife") {
      out <- out * (.n - 1)^2/.n 
    }
    out
  })
}
#' @describeIn infer (TODO)
SdResample <- function(.x, .method, .n) {
  
  lapply(.x, function(x) {
    out        <- matrixStats::colSds(x$Resampled)
    names(out) <- names(x$Original)
    if(.method == "jackknife") {
      out <- out * (.n - 1)/sqrt(.n)
    }
    out
  })
}
#' @describeIn infer (TODO)
BiasResample <- function(.x, .method, .n) {
  
  lapply(.x, function(x) {
    out        <- colMeans(x$Resampled) - x$Original
    names(out) <- names(x$Original)
    if(.method == "jackknife") {
      out <- out * (.n - 1)
    }
    out
  })
}
#' @describeIn infer (TODO)
StandardCIResample <- function(
  .x, 
  .bias_corrected,
  .df = c("type1", "type2"),
  .dist = c("z", "t"), 
  .method, 
  .n,
  .probs
  ) {
  # Standard CI with bootstraped SE's
  # Apparently critical quantiles can be based on both the t or the 
  # standard normal distribution (z). The former may perform better in
  # small samples but there is no clear consenus on what the degrees of freedom
  # should be.
  # CI: [theta_hat + c(alpha/2)*sd_boot ; theta_hat + c(1 - alpha/2)*sd_boot]
  #   = [theta_hat - a ; theta_hat + a]
  boot_sd <- SdResample(.x, .method = .method, .n = .n)
  
  # w = .x := the list containing the estimated values based on .R resamples
  #           and the original values
  # y      := the estimated standard errors (based on .R resamples)
  # z      := a vector of probabilities
  out <- mapply(function(w, y) {
    lapply(.probs, function(z) {
      
      theta_star <- if(.bias_corrected) {
        2*w$Original - colMeans(w$Resampled) 
        # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
      } else {
        w$Original
      }
      
      if(.dist == "t") {
        df <- switch(.df,
                     "type1" = .n - 1,
                     "type2" = nrow(w) - 1
        )
        theta_star + qt(z, df = df) * y
      } else {
        theta_star + qnorm(z) * y
      }
    })
  }, w = .x, y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      rownames(x) <- sprintf("%.2f%%", .probs*100)
      x
    })
  
  return(out)
}
#' @describeIn infer (TODO)
PercentilCIResample <- function(.x, .probs) {
  # Percentile CI 
  # Take the distribution of the bootstrap distribution F* (the CDF) as an estimator for
  # the true distribution of theta_hat. Use the quantilies of that distribution
  # to estimate the CI for theta
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  lapply(.x, function(x) {
    out <- t(matrixStats::colQuantiles(x$Resample, probs = .probs, drop = FALSE))
    colnames(out) <- names(x$Original)
    out
  })
}
#' @describeIn infer (TODO)
TStatCIResample <- function(
  .x1, 
  .x2, 
  .bias_corrected,
  .method, 
  .method2, 
  .n, 
  .probs
  ) {
  # Bootstraped t statistic
  # Bootstrap the t-statisic (since it is roughly pivotal) and compute
  # CI based on bootstraped t-values and bootstraped/jackknife SE
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  
  boot_sd_star <- lapply(.x2, SdResample, .method = .method2, .n = .n) %>% 
    purrr::transpose(.) %>% 
    lapply(function(y) do.call(rbind, y))
  
  boot_sd <- SdResample(.x1, .method = .method, .n = .n)
  
  boot_t_stat <- mapply(function(x, y) t(t(x$Resampled) - x$Original) / y,
                   x = .x1,
                   y = boot_sd_star)

  out <- mapply(function(i, y) {
    lapply(.probs, function(z) {
      qt_boot <- t(matrixStats::colQuantiles(boot_t_stat[[i]], probs = z, drop = FALSE))
      
      theta_star <- if(.bias_corrected) {
        2*.x1[[i]]$Original - colMeans(.x1[[i]]$Resampled) # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
      } else {
        .x1[[i]]$Original
      }
      ## Confidence interval
      theta_star + qt_boot * y
    })
  }, i = seq_along(.x1), y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      rownames(x) <- sprintf("%.2f%%", .probs*100)
      x
    })

  return(out)
}
# bootBcaCI<- function(x) {
#   # Bias corrected and accelerated CI
#   
#   if(is.null(L)){
#     Call <- x$call
#     Call[[1]] <- as.name("jackknife")
#     Call$R <- NULL
#     Call$seed <- NULL
#     Call$sampler <- NULL
#     Call$block.size <- NULL
#     jackObject <- eval(Call, sys.parent())
#     L <- -jackObject$replicates
#   }
#   a <- apply(L, 2, skewness) / (6 * sqrt(nrow(L)))
#   w <- qnorm(colMeans(x$replicates < rep(x$observed, each = x$R)))
#   probs2 <- IfElse(expand, ExpandProbs(probs, min(x$n)), probs)
#   zalpha <- qnorm(probs2)
#   # For now probs in rows, statistics in columns
#   zalpha <- matrix(zalpha, nrow = length(probs), ncol = x$p)
#   probs3 <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
#   result <- probs3 * NA
#   for(j in 1:x$p) {
#     result[, j] <- Quantile(x$replicates[, j], probs3[, j])
#   }
#   dimnames(result) <- list(.FormatProbs(probs), names(x$observed))
#   t(result)
# }
