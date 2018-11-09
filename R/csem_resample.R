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
          split(as.data.frame(data), rep(1:cv_folds, 
                                         each = ceiling(nrow(data)/cv_folds)))
        )
      }
    }
  )
  ## Return samples
  out
}

infer <- function(
  .object = NULL,
  .method = c("Bootstrap", "Jackknife"),
  .draws = 499,
  .handle_inadmissibles = c("drop", "ignore", "replace"),
  .verbose = TRUE
) {
  
  method <- match.arg(.method)
  inad   <- match.arg(.handle_inadmissibles)
  
  ## Is the object to use the data to resample from admissible at all?
  if(sum(verify(.object)) != 0) {
    warning("Estimation based on the data in `.object` has produced inadmissible\n",
            "Results. This may be a sign that something is wrong.")
  }
  
  ## Extract relevant quantities
  args <- .object$Informatio$Arguments
  data <- .object$Informatio$Data
  
  ## Resample
  ll <- resampleData(.object, 
                 .method = method,
                 .draws = .draws)
  
  # Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .draws, style = 3)
  }

  xx <- lapply(1:.draws, function(x) {
    
    args[[".data"]] <- ll[[x]]
    
    Est_out <- do.call(csem, args)
    if(.verbose){
      setTxtProgressBar(pb, x)
    }
    
    Est_out$Estimates
  })
  
  # close progress bar
  if(.verbose){
    close(pb)
  }
  
  ## Process data
  out <- purrr::transpose(xx) %>% 
    lapply(function(x) Reduce("+", x))
  
  
  ### return
  out

  # if(inad == "drop") {
  #   
  #   boot.out$t <- boot.out$t[complete.cases(boot.out$t), , drop = FALSE]
  #   
  # } else {
  #   
  #   ninad <- nrow(boot.out$t[!complete.cases(boot.out$t), , drop = FALSE])
  #   
  #   ## keep on bootstrapping until the number of bootstrapped replicates is equal to .B
  #   while(.draws > .draws - ninad) {
  #     boot.out$t <- boot.out$t[complete.cases(boot.out$t), , drop = FALSE]
  #     tmp <- boot::boot(data = data, 
  #                       statistic = f, 
  #                       R = ninad, 
  #                       names = names,
  #                       inad = inad)
  #     boot.out$t <- rbind(boot.out$t, tmp$t)
  #     # boot.out$R <- boot.out$R + ninad
  #     
  #     ninad <- nrow(boot.out$t[!complete.cases(boot.out$t), , drop = FALSE])
  #   }
  # }
  # 
  # #### Compute quantities ======================================================
  # ## Standard errors
  # if(any(what %in% c("all", "mean", "all"))) {
  #   arith_mean <- matrixStats::colMeans2(boot.out$t)
  # }
  # if(any(what %in% c("all", "se", "t-stat"))) {
  #   se <- matrixStats::colSds(boot.out$t)
  # }
  # if(any(what %in% c("ci", "all"))) {
  #   lapply(1:ncol(boot.out$t), function(x) {
  #     ci <- boot::boot.ci(boot.out, type = c("norm", "basic"), index = x, ...)
  #   })
  # }
  # if(any(what %in% c("all", "t-stat"))) {
  #   t  <- boot.out$t %*% diag(1/se)
  # }
  # 
  # list(
  #   "Mean"   = if(any(what %in% c("mean", "all"))){
  #     arith_mean
  #   } else {
  #     NA
  #   },
  #   "Se"     = if(any(what %in% c("se", "all"))) {
  #     se
  #   } else {
  #     NA
  #   },
  #   # "Ci"   = ifelse(any(what %in% c("ci", "all")), ci, NA),
  #   "T-Stat" = if(any(what %in% c("t-stat", "all"))) {
  #     t
  #   } else {
  #     NA
  #   }
  # )
}
# 
# debugonce(infer)
# yy <- infer(a, .method = "Bootstrap", .draws = 30)
# y1 <- yy[1:4]
# y2 <- lapply(y3, function(x) Reduce("+", x))
# y2
# y3 <- purrr::transpose(yy)
# y3
# debugonce(Reduce)
# Reduce("+", y2)
# 
# listviewer::jsonedit(y3, mode = "view")

