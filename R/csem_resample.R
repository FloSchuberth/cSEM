#' Resample data 
#'
#' Resample from a given cSEMResults object or data.frame using common
#' resampling methods.
#'
#' @usage resampleData(
#'  .object        = args_default()$.object,
#'  .data          = NULL,
#'  .method        = c("bootstrap", "jackknife", "permutation", "cross-validation"),
#'  .draws         = 499,
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
  .draws    = args_default()$.draws,
  .id       = args_default()$.id,
  .cv_folds = args_default()$.cv_folds
) {
  
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
      lapply(1:.draws, function(x) data[-x, ])
    },
    "bootstrap"   = {
      lapply(1:.draws, function(x) {
        data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]})
    },
    "permutation" = {
      if(is.null(id)) {
        lapply(1:.draws, function(x) {
          data[sample(1:nrow(data), size = nrow(data), replace = FALSE), ]
        })
      } else {
        lapply(1:.draws, function(x) cbind(data, sample(id)))
      }
    },
    "cross-validation" = {
      if(is.null(.cv_folds)) {
        # LOOCV (Leave one out CV). Technically, this is the same a jackknife
        lapply(1:.draws, function(x) data[-x, ])
        
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
#'  .draws         = 499,
#'  .id            = NULL,
#'  .cv_folds      = 10      
#' )
#'
#' @inheritParams csem_arguments
#' @param .method Character string. The resampling method to use. One of: 
#'  "*bootstrap*" or "*jackknife*". Defaults to "*bootstrap*".
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
  .draws                 = args_default()$.draws,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .method                = args_default()$.method2,
  .verbose               = args_default()$.verbose,
  .user_funs             = args_default()$.user_funs,
  .draws2                = args_default()$.draws2,
  .handle_inadmissibles2 = args_default()$.handle_inadmissibles2,
  .method2               = args_default()$.method2
  ) {
  
  ## Has the object to use the data to resample from produced admissible results?
  if(sum(verify(.object)) != 0) {
    warning("Estimation based on the data in `.object` has produced",
            " inadmissible results.\n", 
            "This may be a sign that something is wrong.",
            call. = FALSE, immediate. = TRUE)
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
  
  ## Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .draws, style = 3)
  }
  
  Est_ls          <- list()
  n_inadmissibles <- 0
  counter         <- 0
  args            <- .object$Information$Arguments
  repeat{
    # Counter
    counter <- counter + 1
    
    # Replace the old dataset by a resampled data set (resampleData always returns
    # a list so for just one draw we need to pick the first list element)
    args[[".data"]] <- resampleData(.object, .method = .method, .draws = 1)[[1]]

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
      ## bootstraped t-interval CI) the original object needs to be returned. This
      ## may be problematic for 
      if(!is.null(.draws2)) {
        Est_resamples2 <- resamplecSEMResults(
          .object               = Est_temp,
          .draws                = .draws2,
          .handle_inadmissibles = .handle_inadmissibles2,
          .method               = .method2,
          .verbose              = FALSE,
          .user_funs            = .user_funs,
        )
        
        x1 <- list("Estimates1" = x1, "Estimates2" = Est_resamples2)
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
    
    # Break repeat loop if .draws results have been created.
    if(length(Est_ls) == .draws) {
      break
    } else if(counter + n_inadmissibles == 10000) { 
      ## Stop if 10000 runs did not result in insufficient admissible results
      stop("Not enough admissible result.", call. = FALSE)
    }
    
    if(.verbose){
      setTxtProgressBar(pb, counter)
    }
    
  } # END repeat 
  # Close progress bar
  if(.verbose){
    close(pb)
  }
  
  ## Process data --------------------------------------------------------------
  # Delete potential NA's
  out <- Filter(Negate(anyNA), Est_ls)
  
  # Check if at least 3 admissible results were obtained
  if(length(out) < 3) {
    stop("The following error occured in the `infer()` functions:\n",
         "Less than 2 admissible results produced.", 
         " Consider setting .handle_inadmissibles = 'replace' instead.",
         call. = FALSE)
  }
  # Turn list "inside out" and bind bootstrap samples to matrix 
  # columns are variables
  # rows are bootstrap runs
  out <- purrr::transpose(out) 
  
  if(!is.null(.draws2)) {
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
  
  out <- if(!is.null(.draws2)) {
    list("Estimates1" = out, "Estimates2" = out_2)
  } else {
    out
  }
  
  ## Set class
  class(out) <- "cSEMResults_resampled"
  return(out)
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
  .object          = NULL,
  .csem_resample   = args_default()$.csem_resample,
  .alpha           = args_default()$.alpha,
  .statistic       = c("all", "mean", "var", "sd", "bias", "CI-standard_z",
                       "CI-standard_t", "CI-percentile", "CI-t-intervall",
                       "CI-Bc", "CI-Bca")
) {
 
  x1 <- .csem_resample$Estimates1
  x2 <- .csem_resample$Estimates2
  
  if(anyNA(x2) & any(.statistic %in% c("all", "CI-t-intervall", "CI-Bca"))) {
    stop(paste0("`", intersect(.statistic, c("all", "CI-t-intervall", "CI-Bca")), 
                "`", collapse = ", "), 
         " requires the resampling distribution for each resample.",
         " Please set `.draws2` in `resamplecSEMResults()`.",
         call. = FALSE)
  }
  
  ## Compute quantiles/critical values -----------------------------------------
  ### Note if jackknife was chosen a correction is needed for the variance 
  # but not sure how exactly
  probs <- c()
  for(i in seq_along(.alpha[order(.alpha)])) { 
    probs <- c(probs, round(.alpha[i]/2, 4), round(1 - .alpha[i]/2, 4)) 
  }
  
  ## Compute statistics and quantities
  l <- list(
    "Mean" = bootMean(x1),
    "Var"  = bootVar(x1),
    "Sd"   = bootSd(x1),
    "Bias" = bootBias(x1),
    "CI-standard_z" = bootStandardCI(x1, .dist = "z"),
    "CI-standard_t" = bootStandardCI(x1, .dist = "t"),
    "CI-percentile" = bootPercentilCI(x1),
    "CI-t-invertall"= bootTstatCI(x1, x2)
  )
  
  return(l)
}

#' @describeIn infer The bootstraped mean
bootMean <- function(x) {
  lapply(x, function(y) {
    out        <- colMeans(y$Resampled)
    names(out) <- names(y$Original)
    out
  })
}
#' @describeIn infer (TODO)
bootVar <- function(x) {
  lapply(x, function(y) {
    out        <- matrixStats::colVars(y$Resampled)
    names(out) <- names(y$Original)
    out
  })
}
#' @describeIn infer (TODO)
bootSd <- function(x) {
  lapply(x, function(y) {
    out        <- matrixStats::colSds(y$Resampled)
    names(out) <- names(y$Original)
    out
  })
}
#' @describeIn infer (TODO)
bootBias <- function(x) {
  lapply(x, function(y) {
    out        <- colMeans(y$Resampled) - y$Original
    names(out) <- names(y$Original)
    out
  })
}
#' @describeIn infer (TODO)
bootStandardCI <- function(x, .dist = c("z", "t")) {
  # Standard CI with bootstraped SE's
  # Apparently critical quantiles can be based on both the t or the 
  # standard normal distribution (z). The former may perform better in
  # small samples.
  # CI: [theta_hat + c(alpha/2)*sd_boot ; theta_hat + c(1 - alpha/2)*sd_boot]
  #   = [theta_hat - a ; theta_hat + a]
  boot_sd <- bootSd(x)
  
  # w = x := the list containing the estimated values based on .draw resamples
  #          and the original values
  # y     := the estimated standard errors (based on .draw resamples)
  # z     := a vector of probabilities
  out <- mapply(function(w, y) {
    lapply(probs, function(z) {
      if(.dist == "t") {
        df <- nrow(.object$Information$Data) - 1
        w$Original + qt(z, df = df) * y
      } else {
        w$Original + qnorm(z) * y
      }
    })
  }, w = x, y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      rownames(x) <- sprintf("%.2f%%", probs*100)
      x
    })
  
  return(out)
}
#' @describeIn infer (TODO)
bootPercentilCI <- function(x) {
  # Percentile CI 
  # Take the distribution of the bootstrap distribution F* (CDF) as an estimator for
  # the true distribution of theta_hat. Use the quantilies of that distribution
  # to estimate the CI for theta
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  lapply(x, function(y) {
    out <- t(matrixStats::colQuantiles(y$Resample, probs = probs, drop = FALSE))
    colnames(out) <- names(y$Original)
    out
  })
}
#' @describeIn infer (TODO)
bootTstatCI <- function(x1, x2) {
  # Bootstraped t statistic
  # Bootstrap the t-statisic (since it is roughly pivotal) and compute
  # CI based on bootstraped t-values and bootstraped SE
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]

  boot_sd_star <- lapply(x2, bootSd) %>% 
    purrr::transpose(.) %>% 
    lapply(function(y) do.call(rbind, y))
  
  boot_sd <- bootSd(x1)
  
  boot_t_stat <- mapply(function(x, y) t(t(x$Resampled) - x$Original) / y,
                   x = x1,
                   y = boot_sd_star)

  out <- mapply(function(i, y) {
    lapply(probs, function(z) {
      qt_boot <- t(matrixStats::colQuantiles(boot_t_stat[[i]], probs = z, drop = FALSE))
      
      ## Confidence interval
      x1[[i]]$Original + qt_boot * y
    })
  }, i = seq_along(x1), y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      rownames(x) <- sprintf("%.2f%%", probs*100)
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