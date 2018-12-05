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
    data <- as.data.frame(.object$Information$Arguments$.data)
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

#' Internal: resample cSEMResults 
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
  .method2               = args_default()$.method2,
  .R                     = args_default()$.R,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .user_funs             = args_default()$.user_funs,
  .future_plan           = c("sequential", "transparent", "multiprocess")
  ) {
  
  ## Set plan on how to resolve futures 
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(.future_plan)
    
  ### Checks, warnings and errors -----------
  match.arg(.method, args_default(.choices = TRUE)$.method)
  match.arg(.method2, args_default(.choices = TRUE)$.method2)
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
  
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
  
  ### Resample and compute -----------------------------------------------------
  
  out <- resamplecSEMResultsCore(
    .object                = .object, 
    .method                = .method,
    .method2               = .method2,
    .R                     = .R,
    .R2                    = .R2,
    .handle_inadmissibles  = .handle_inadmissibles,
    .user_funs             = .user_funs
  )
  
  # Check if at least 3 admissible results were obtained
  n_admissibles <- length(out)
  if(n_admissibles < 3) {
    stop("The following error occured in the `resamplecSEMResults()` functions:\n",
         "Less than 2 admissible results produced.", 
         " Consider setting `.handle_inadmissibles == 'replace'` instead.",
         call. = FALSE)
  }
  
  # Turn list "inside out" and bind bootstrap samples to matrix 
  # columns are variables
  # rows are bootstrap runs
  out <- purrr::transpose(out) 
  
  if(.method2 != "none") {
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
      out <- list("Estimates1" = out, "Estimates2" = out_2)
      
    } else {
      out <- list("Estimates1" = out, "Estimates2" = NA)
    }
    
    out <- list(
      "Resamples" = out,
      "Information" = list(
        "Method"                  = .method,
        "Method2"                 = .method2,
        "Number_of_observations"  = nrow(.object$Information$Data)
      )
    )
    
    ## Set class
    class(out) <- "cSEMResults_resampled"
    return(out)
  }
}

#' Core tasks of the resamplecSEMResults function
#' @noRd
#' 
resamplecSEMResultsCore <- function(
  .object                = args_default()$.object, 
  .method                = args_default()$.method,
  .method2               = args_default()$.method2,
  .R                     = args_default()$.R,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .user_funs             = args_default()$.user_funs,
  .future_plan           = c("sequential", "transparent", "multiprocess")
) {
  
  ## Get arguments
  args <- .object$Information$Arguments
  
  ## Resample jackknife
  if(.method == "jackknife") {
    resample_jack <- resampleData(.object, .method = "jackknife") 
    .R <- length(resample_jack)
  }
  
  Est_ls <- future.apply::future_lapply(1:.R, function(i) {
    # Replace the old dataset by a resampled data set (resampleData always returns
    # a list so for just one draw we need to pick the first list element)
    
    data_temp <- if(.method == "jackknife") {
      resample_jack[[i]]
    } else {
      resampleData(.object, .method = "bootstrap", .R = 1)[[1]]
    }
    
    args[[".data"]] <- data_temp
    
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
      
      ## Resampling from a bootstrap sample is required for the
      ## bootstraped t-interval CI, hence the second run
      if(.method2 != "none") {
        
        Est_resamples2 <- resamplecSEMResults(
          .object               = Est_temp,
          .R                    = .R2,
          .handle_inadmissibles = .handle_inadmissibles,
          .method               = .method2,
          .user_funs            = .user_funs
        )
        x1 <- list("Estimates1" = x1, "Estimates2" = Est_resamples2)
      }
    } else if(status_code != 0 & .handle_inadmissibles != "ignore") {
      # Retrun NA if status is not okay and .handle_inadmissibles == "drop" or
      # "replace"
      x1 <- NA
    }
    ## Return
    x1
  })
  
  ## Process data --------------------------------------------------------------
  # Delete potential NA's
  out <- Filter(Negate(anyNA), Est_ls)
  
  ## Replace inaadmissibles
  if(length(out) != .R && .method == "bootstrap" && 
     .handle_inadmissibles == "replace") {
    
    R_new <- .R - length(out)
    
    while (length(out) < .R) {
      Est_replace <- resamplecSEMResultsCore(
        .object               = .object,
        .R                    = R_new,
        .handle_inadmissibles = .handle_inadmissibles,
        .method               = .method,
        .method2              = .method2,
        .R2                   = .R2,
        .user_funs            = .user_funs,
        .future_plan          = ifelse(R_new < 150, "sequential", .future_plan)
      )
      
      out <- c(out, Est_replace)
    }
  }
  
  out
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
  .object          = args_default()$.object,
  .csem_resample   = args_default()$.csem_resample,
  .alpha           = args_default()$.alpha,
  .bias_corrected  = args_default()$.bias_corrected,
  .statistic       = args_default()$.statistic,
  ...
) {
  
  ## Compute the .csem_resample object if none has been provided
  if(is.null(.csem_resample)) {
    .csem_resample <- resamplecSEMResults(
      .object = ..object,
      ...
    )
  }
  
  resamples1 <- .csem_resample$Resamples$Estimates1 # jackknife or bootstrap
  resamples2 <- .csem_resample$Resamples$Estimates2 # jackknife or bootstrap
  info       <- .csem_resample$Information

  ## Compute quantiles/critical values -----------------------------------------
  probs  <- c()
  .alpha <- .alpha[order(.alpha)]
  for(i in seq_along(.alpha)) { 
    # for every alpha: alpha/2 and 1 - alpha/2
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
    "CI_basic"      = BasicCIResample(
      .x              = resamples1, 
      .bias_corrected = .bias_corrected, 
      .probs          = probs
      ),
    "CI_Bc"  = BcCI(resamples1, probs),
    "CI_Bca" = if(is.null(.object)) {
      NA
    } else {
      BcaCI(.object, resamples1, probs)
    },
    "CI_t_invertall"= if(anyNA(resamples2)) {
      NA
    } else {
      TStatCIResample(
        .x1      = resamples1, 
        .x2      = resamples2, 
        .bias_corrected = .bias_corrected,
        .method  = info$Method, 
        .method2 = info$Method2, 
        .n       = info$Number_of_observations, 
        .probs   = probs
      ) 
    }
  )
  
  if(anyNA(l)) {
    l <- l[-which(is.na(l))]
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
  # Standard CI with bootstrap SE's
  # Critical quantiles can be based on both the t- or the 
  # standard normal distribution (z). The former may perform better in
  # small samples but there is no clear consenus on what the degrees of freedom
  # should be.
  # CI: [theta_hat + c(alpha/2)*boot_sd ; theta_hat + c(1 - alpha/2)*boot_sd]
  #   = [theta_hat - c(1 - alpha/2)*boot_sd ; theta_hat + c(1 - alpha/2)*boot_sd]
  # if c() is based on a symmetric distribution.
  
  ## Compute standard deviation
  boot_sd <- SdResample(.x, .method = .method, .n = .n)
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  ## Compute intervals
  # Notation:
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
      rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      x
    })
  return(out)
}
#' @describeIn infer (TODO)
PercentilCIResample <- function(.x, .probs) {
  # Percentile CI 
  # Take the bootstrap distribution F* (the CDF) as an estimator for
  # the true distribution F. Use the quantiles of the estimated distribution
  # to estimate the CI for theta.
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  lapply(.x, function(x) {
    out <- t(matrixStats::colQuantiles(x$Resampled, probs = .probs, drop = FALSE))
    colnames(out) <- names(x$Original)
    rownames(out) <- paste0(c("L_", "U_"), rep(cl, each = 2))
    out
  })
}
#' @describeIn infer (TODO)
BasicCIResample <- function(.x, .bias_corrected, .probs) {
  # Basic CI 
  # Estimate the distribution of delta_hat = theta_hat - theta by the bootstrap 
  # distribution delta_hat_star = theta_hat_star - theta_hat.
  # Since P(a <= theta_hat - theta <= b) = P(theta_hat - b <= theta <= theta_hat - a)
  # (notice how the limits a and b switch places!!) 
  # Define q_p := p%-quantile of the bootstrap distribution of delta_hat_star
  # and
  #  b = q_(1 - alpha/2) 
  #  a = q_(alpha/2) 
  # then the CI is:
  # CI: [theta_hat - q_(1 - alpha/2); theta_hat - q_(alpha/2)]
  # Note that q_p is just the alpha percentile quantile shifted by theta_hat:
  # q_p = Q_p - theta_hat, where Q_P = F*^-1(p) := the p% quantile of the
  # distribution of theta_hat_star.
  # Therefore:
  # CI: [2*theta_hat - Q_(1 - alpha/2); 2*theta_hat - Q_(alpha/2)]
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  lapply(.x, function(x) {
    out <- t(matrixStats::colQuantiles(x$Resample, probs = .probs, drop = FALSE))
    
    theta_star <- if(.bias_corrected) {
      2*x$Original - colMeans(x$Resampled) 
      # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
    } else {
      x$Original
    }
    
    out <- t(2*theta_star - t(out))
    colnames(out) <- names(x$Original)
    out <- out[1:nrow(out) + rep(c(1, -1), times = nrow(out)/2), , drop = FALSE]
    rownames(out) <- paste0(c("L_", "U_"), rep(cl, each = 2))
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

  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  boot_sd_star <- lapply(.x2, SdResample, .method = .method2, .n = .n) %>% 
    purrr::transpose(.) %>% 
    lapply(function(y) do.call(rbind, y))
  
  boot_sd <- SdResample(.x1, .method = .method, .n = .n)
  
  boot_t_stat <- mapply(function(x, y) t(t(x$Resampled) - x$Original) / y,
                   x = .x1,
                   y = boot_sd_star)

  i <- 1
  y <- boot_sd[[1]]
  z <- .probs[1]
  out <- mapply(function(i, y) {
    lapply(.probs, function(z) {
      qt_boot <- t(matrixStats::colQuantiles(boot_t_stat[[i]], probs = z, drop = FALSE))
      
      theta_star <- if(.bias_corrected) {
        2*.x1[[i]]$Original - colMeans(.x1[[i]]$Resampled) # theta_hat - Bias = 2*theta_hat - mean_theta_hat_star
      } else {
        .x1[[i]]$Original
      }
      ## Confidence interval
      theta_star - qt_boot * y
    })
  }, i = seq_along(.x1), y = boot_sd, SIMPLIFY = FALSE) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) {
      x <- x[1:nrow(x) + rep(c(1, -1), times = nrow(x)/2), , drop = FALSE]
      rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      x
    })

  names(out) <- names(.x1)
  return(out)
}
#' @describeIn infer (TODO)
BcCI <- function(.x, .probs) {
  
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2

  out <- lapply(.x, function(x) {
    p0 <- colMeans(t(t(x$Resampled) <= x$Original))
    z0 <- qnorm(p0)
    
    bc <- lapply(1:length(z0), function(i) {
      bc_quantiles <- c()
      for(j in seq_along(.probs)) {
        q  <- pnorm(2*z0[i] + qnorm(.probs[j]))
        bc_quantiles <- c(bc_quantiles, quantile(x$Resampled[, i], probs = q))
      }
      bc_quantiles
    })
    names(bc) <- names(x$Original)
    bc
  }) %>% 
    lapply(function(x) t(do.call(rbind, x))) %>% 
    lapply(function(x) {
      rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      x
    })
  
  return(out)
}
#' @describeIn infer (TODO)
BcaCI <- function(.object, .x, .probs) {
  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  ## influence values (estimated via jackknife)
  jack <- resamplecSEMResults(
    .object = .object,
    .method = "jackknife"
  )$Resamples$Estimates1
  
  aFun <- function(x) {
    1/6 * (sum(x^3) / sum(x^2)^1.5)
  }
  
  a <- lapply(jack, function(x) t(x$Original - t(x$Resampled))) %>% 
    lapply(function(x) apply(x, 2, aFun)) 
  
  p0 <- lapply(.x, function(x) colMeans(t(t(x$Resampled) <= x$Original)))
  z0 <- lapply(p0, qnorm)

  out <- lapply(seq_along(.x), function(i) {
    bca <- lapply(1:length(z0[[i]]), function(j) {
      q <- pnorm(z0[[i]][j] + (z0[[i]][j] + qnorm(.probs))/ (1 - a[[i]][j]*(z0[[i]][j] + qnorm(.probs)))) 
      quantile(.x[[i]]$Resampled[, j], probs = q)
    })
    names(bca) <- names(.x[[i]]$Original)
    bca
  })%>% 
    lapply(function(x) t(do.call(rbind, x))) %>% 
    lapply(function(x) {
      rownames(x) <- paste0(c("L_", "U_"), rep(cl, each = 2))
      x
    })
  
  names(out) <- names(.x)
  out
}
