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

a <- csem(satisfaction, model)
resampleData(.object = a, .method = "Bootstrap", .draws = 1)

.object <- a
.draws <- 20
method <- "Bootstrap"
.handle_inadmissibles <- "replace"
.alpha <- .alpha[order(.alpha)]
.alpha <- c(0.1, 0.2, 0.05)
.verbose = TRUE



infer <- function(
  .object = NULL,
  .method = c("Bootstrap", "Jackknife"),
  .draws = 499,
  .handle_inadmissibles = c("drop", "ignore", "replace"),
  .verbose = TRUE
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
    pb <- txtProgressBar(min = 0, max = .draws, style = 3)
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
    status_code <- sum(verify(Est_temp))
    
    # Distinguish depending on how inadmissibles should be handled
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {

      Est_ls[[counter]] <- Est_temp$Estimates
      
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
  # Delete inadmissibles
  out <- Filter(Negate(anyNA), Est_ls)

  # Check if at least 10 admissible results were obtained
  if(length(out) < 11) {
    stop("The following error occured in the `infer()` functions:\n",
         "Less than 10 admissible results produced.", 
         " Consider setting .handle_inadmissibles = 'replace' instead.",
         call. = FALSE)
  }
  
  out <- purrr::transpose(out)
  
  # Some quantities are matrices and others are vectors. They need to be handled
  # differently so we split them
  out_mat <- out[c("Path_estimates", 
                   "Loading_estimates", 
                   "Weight_estimates", 
                   "Inner_weight_estimates",
                   "Proxy_VCV",
                   "Construct_VCV",
                   "Cross_loadings")] %>% 
    lapply(simplify2array)
  

  funs <- list(
    "Mean"   = function(x) mean(x),
    "Median" = function(x) median(x), 
    "Var"    = function(x) var(x),
    "Sd"     = function(x) sd(x),
    "Quantile" = function(x) lapply(.alpha, function(y) quantile(x, probs = 1 - y))
  )
  
  funs
  ## Apply each function to compute the statistics need 
  tt <- lapply(funs, function(fun) {
    lapply(out_mat, function(x) apply(x, c(1, 2), FUN = fun))
    })
  
 #  tt$Mean
 #  tt$Median
 #  tt$Sd
 #  mm <- tt$Quantile$Path_estimates
 #  str(mm[1,1])
 #  ff <- function(x) lapply(1:length(mm[1,1])mm[i, j]
 # outer(rownames(mm), colnames(mm), FUN = Vectorize())
 #  
 #  lapply(function(x, y) apply(x, c(1, 2), list("mean", "median", "var")), x = out_mat, y = list("mean", "median", "var"))
  
    lapply(function(x) apply(x, c(1, 2), mean))
    
    lapply(function(x) apply(x, c(1, 2), function(y) list(
      "Mean"     = mean(x),
      "Median"   = median(x),
      "Var"      = var(x),
      "Sd"       = sd(x),
      "Quantile" = quantile(x, probs = 1 - .alpha)
      ))
    )
  
  # out_mat$Path_estimates
  # outer(rownames(out_mat$Path_estimates), colnames(out_mat$Path_estimates), 
  #       function(x) )

  out_vec <- out[c("Construct_reliabilities", "R2", "R2adj", "VIF")] %>% 
    lapply(function(x) do.call(rbind, x)) %>% 
    lapply(function(x) list(
      "Mean"     = colMeans(x),
      "Median"   = {tt <- matrixStats::colMedians(x); names(tt) <- colnames(x); tt},
      "Var"      = {tt <- matrixStats::colVars(x); names(tt) <- colnames(x); tt},
      "Sd"       = {tt <- matrixStats::colSds(x); names(tt) <- colnames(x); tt},
      "Quantile" = matrixStats::colQuantiles(x, probs = 1 - .alpha, drop = FALSE)
    ))
  
}

# debugonce(infer)
# infer(a, .method = "Bootstrap", .draws = 20, .handle_inadmissibles = "replace")

# 
# ## Compute csem estimate for each data set
# Est_ls <- lapply(1:.draws, function(x) {
#   
#   ## Replace data
#   args[[".data"]] <- ll[[x]]
#   
#   ## Estimate
#   Est_out <- do.call(csem, args)
#   
#   ## Update progress bar
#   if(.verbose){
#     setTxtProgressBar(pb, x)
#   }
#   
#   ## Return the estimated quantities
#   Est_out
# })
# 
# # Check status (the sum is the number of inadbmissibles)
# ninad <- sum(unlist(lapply(Est_ls, verify)))
# 
# ## Handle inadmissibles
# if(ninad != 0 && handle_inadmissibles == "drop") {
#   
#   # Delete inadmissibles
#   Est_ls <- Filter(function(x) sum(verify(x)) == 0, Est_ls)
# 
# }
# 
# if(ninad != 0 && handle_inadmissibles == "replace") {
#   
#   # Delete inadmissibles
#   Est_ls <- Filter(function(x) sum(verify(x)) == 0, Est_ls)
#   
#   ## Keep on bootstrapping until the number of admissible bootstrapped 
#   ## replicates is equal to .draws
#   while(.draws > .draws - ninad) {
#     
#     ninad <- nrow(boot.out$t[!complete.cases(boot.out$t), , drop = FALSE])
#   }
# }
#   
# 
# }    
#      
# if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore"))
