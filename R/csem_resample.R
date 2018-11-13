#' Internal: Resample
#'
#' Resample from a given cSEMResults object or data.frame using common
#' resampling methods.
#'
#' @usage resampleData(
#'  .object        = args_default()$.object,
#'  .data          = NULL,
#'  .method        = c("Bootstrap", "Jackknife", "Permutation", "Cross-validation"),
#'  .draws         = 499,
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
resampleData <- function(
  .object = NULL,
  .data = NULL,
  .method = c("Bootstrap", "Jackknife", "Permutation", "Cross-validation"),
  .draws = 499,
  .id = NULL,
  .cv_folds = 10
) {
  ## Match arguments
  method <- match.arg(.method)
  
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
  
  out <- switch (method,
    "Jackknife"   = {
      lapply(1:.draws, function(x) data[-x, ])
    },
    "Bootstrap"   = {
      lapply(1:.draws, function(x) {
        data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]})
    },
    "Permutation" = {
      if(is.null(id)) {
        lapply(1:.draws, function(x) {
          data[sample(1:nrow(data), size = nrow(data), replace = FALSE), ]
        })
      } else {
        lapply(1:.draws, function(x) cbind(data, sample(id)))
      }
    },
    "Cross-validation" = {
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
    }
  )
  ## Return samples
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
  .object = NULL,
  .method = c("Bootstrap", "Jackknife"),
  .runs = 499,
  .handle_inadmissibles = c("drop", "ignore", "replace"),
  .verbose = TRUE,
  .alpha = 0.05,
  .user_fun = NULL,
  ...
) {
  
  method <- match.arg(.method)
  # handle_inadmissibles  <- match.arg(.handle_inadmissibles)
  
  ## Is the object to use the data to resample from admissible at all?
  if(sum(verify(.object)) != 0) {
    warning("Estimation based on the data in `.object` has produced",
            " inadmissible Results.\n", 
            "This may be a sign that something is wrong.",
            call. = FALSE, immediate. = TRUE)
  }
  
  ## Extract relevant quantities
  X    <- .object$Informatio$Data
  args <- .object$Informatio$Arguments
  
  ## Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .runs, style = 3)
  }

  ## Calculate reference distribution
  Est_ls           <- list()
  n_inadmissibles  <- 0
  counter <- 0
  repeat{
    # Counter
    counter <- counter + 1
    
    # Replace the old dataset by a resampled data set (resampleData always returns
    # a list so for just one draw we need to pick the first list element)
    args[[".data"]] <- resampleData(.object, .method = method, .draws = 1)[[1]]
    
    # Estimate model
    Est_temp <- do.call(csem, args)   
    
    # Check status
    status_code <- sum(unlist(verify(Est_temp)))
    
    # Distinguish depending on how inadmissibles should be handled
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
      
      # Compute user defined function if specified
      if(!is.null(.user_fun)) {
        user_funs <- lapply(.user_fun, function(f) c(f(Est_temp)))
      }
      
      x1  <- Est_temp$Estimates
      x2  <- Est_temp$Information
      
      x1$Path_estimates    <- t(x1$Path_estimates)[t(x2$Model$structural) != 0]
      x1$Loading_estimates <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0]
      x1$Weight_estimates  <- t(x1$Weight_estimates)[t(x2$Model$measurement) != 0]
      x1$Inner_weight_estimates <- NULL
      x1$Construct_scores <- NULL
      x1$Indicator_VCV <- NULL
      x1$Proxy_VCV <- c(x1$Proxy_VCV)
      x1$Construct_VCV <- c(x1$Construct_VCV)
      x1$Cross_loadings <- NULL
      
      x1 <- c(x1, user_funs)
      
      Est_ls[[counter]] <- x1
      
    } else if(status_code != 0 & .handle_inadmissibles == "drop") {
      # Set list element to NA if status is not okay and .handle_inadmissibles == "drop"
      Est_ls[[counter]] <- NA
      
    } else {# status is not ok and .handle_inadmissibles == "replace"
      # Reset counter and raise number of inadmissibles by 1
      counter <- counter - 1
      n_inadmissibles <- n_inadmissibles + 1
    }
    
    # Break repeat loop if .runs results have been created.
    if(length(Est_ls) == .runs) {
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
  # Delete NA's
  out <- Filter(Negate(anyNA), Est_ls)

  # Check if at least 3 admissible results were obtained
  if(length(out) < 3) {
    stop("The following error occured in the `infer()` functions:\n",
         "Less than 2 admissible results produced.", 
         " Consider setting .handle_inadmissibles = 'replace' instead.",
         call. = FALSE)
  }
  
  out2 <- purrr::transpose(out) %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) list(
      "Mean"     = colMeans(x),
      "Median"   = matrixStats::colMedians(x),
      "Var"      = matrixStats::colVars(x),
      "Sd"       = matrixStats::colSds(x),
      "Quantile" = t(matrixStats::colQuantiles(x, probs = 1 - .alpha, drop = FALSE))
    )
    )
  
  return(out2)
}
