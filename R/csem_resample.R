#' Resample data 
#'
#' Resample from a given dataset using common resampling methods. 
#' For bootstrap or jackknife resampling, package users usually dont need to 
#' call this function but directly use [resamplecSEMResults()] instead.
#'
#' The function `resampleData()` is general purpose. It simply resamples data 
#' from a given dataset according to the resampling method provided 
#' via the `.resample_method` argument. 
#' Currently, `bootstrap`, `jackknife`, `permutation`, and  `cross-validation`
#' (both leave-one-out (LOOCV) and k-fold cross-validation) are implemented. 
#' 
#' The user may provide a dataset to be resampled or a [cSEMResults] 
#' object in which case the original data used in the call to [csem()] is used. 
#' The `.data` argument is `NULL` by default. If both, a [cSEMResults] object and 
#' a dataset via `.data` are provided the latter is ignored. 
#' 
#' As [csem()] accepts a single dataset, a list of datasets as well as datasets
#' that contain a column name used to split the data into groups,
#' the [cSEMResults] object may contain multiple datasets.
#' In this case, resampling is done by dataset or group. Note that depending
#' on the number of datasets provided this computation may be significantly slower
#' as resampling will be repeated for each dataset/group.
#' 
#' If data containing a column to split the data into groups is 
#' provided via the `.data` argument, the column name containing the group levels
#' must be given to `.id`.
#' 
#' The number of bootstrap or permutation runs is given by `.R`. The default is
#' `499` but should be increased in real applications. See e.g.,
#' \insertCite{Hesterberg2015;textual}{cSEM}, p.380 for recommendations.
#' For jackknife and cross-validation `.R` is ignored.
#' 
#' For cross-validation the number of folds defaults to `10` and can be changed
#' via the `.cv_folds` argument. Setting `.cv_folds` to `NULL` produces
#' leave-one-out cross-validation samples.
#' 
#' @usage resampleData(
#'  .object          = NULL,
#'  .resample_method = c("bootstrap", "jackknife", "permutation", "cross-validation"),
#'  .cv_folds        = 10,  
#'  .data            = NULL,
#'  .id              = NULL,
#'  .R               = 499
#' )
#'
#' @param .data A `data.frame`, a `matrix` or a list containing data of either type. 
#'   Possible column types or classes of the data provided are: 
#'   logical, numeric (double or integer), factor (ordered and unordered) 
#'   or a mix of several types. The data may also include
#'   *one* character column whose column name must be given to `.id`. 
#'   This column is assumed to contain group identifiers used to split the data into groups.
#'   If `.object` is provided, `.data` is ignored. Defaults to `NULL`.
#' @param .resample_method Character string. The resampling method to use. One of: 
#'  "*bootstrap*", "*jackknife*", "*permutation*", or "*cross-validation*". 
#'  Defaults to "*bootstrap*".
#' @param .R Integer. The number of bootstrap replications or permutation runs
#'   to use. Defaults to `499`.
#' @inheritParams csem_arguments
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [cSEMResults], [resamplecSEMResults()], [infer()]
#'
#' @examples
#' \dontrun{
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
#' a <- csem(satisfaction, model)
#' 
#' res <- resampleData(a, .resample_method = "bootstrap", .R = 999)
#' str(res)
#' }
#' 
#' @export
#'
resampleData <- function(
  .object          = args_default()$.object,
  .resample_method = c("bootstrap", "jackknife", "permutation", "cross-validation"),
  .cv_folds        = args_default()$.cv_folds,
  .data            = args_default()$.data,
  .id              = args_default()$.id,
  .R               = args_default()$.R
) {
  .resample_method <- match.arg(.resample_method, 
            c("bootstrap", "jackknife", "permutation", "cross-validation"))

  ## Get data set
  if(is.null(.data)) {
    ## Get information according to class of object
    if(any(class(.object) %in% "cSEMResults_default")) {
      data <- as.data.frame(.object$Information$Arguments$.data)
      id   <- NULL
    } else if(any(class(.object) %in% "cSEMResults_multi")) {
      data        <- .object[[1]]$Information$Data_pooled
      data_split  <- lapply(.object, function(x) x$Information$Arguments$.data)
      id          <- .object[[1]]$Information$Arguments$.id
      
    } else if(any(class(.object) %in% "cSEMResults_2ndorder")) {
      data <- .object$First_stage$Information$Arguments$.data
      id   <- NULL
    } else {
      stop2("The following error occured in the `resampleData()` function:\n",
            "`object` must be of class cSEMResults."
            )
    }
  } else {
    
    ## Check data
    if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
      stop2(
        "The following error occured in the `resampleData()` function:\n",
        "Data must be provided as a `matrix`, a `data.frame` or a `list`. ", 
         ".data has class: ", class(.data)
        )
    }
    ## Select cases
    if(!is.null(.id) && !inherits(.data, "list")) {
      
      if(length(.id) != 1) {
        stop2(
          "The following error occured in the `resampleData()` function:\n",
          "`.id` must be a character string or an integer identifying one single column."
          )
      }
      
      if(is.matrix(.data)) {
        .data <- as.data.frame(.data)
      }
      
      data       <- .data 
      data_split <- split(.data, f = .data[, .id])
      id         <- .id

    } else if(any(class(.data) == "list")) {
      
      data <- do.call(rbind, .data)

      if(is.null(names(.data))) {
        data[, "id"] <- rep(paste0("Data_", 1:length(.data)), 
                            times = sapply(.data, nrow))
      } else {
        data[, "id"] <- rep(names(.data), times = sapply(.data, nrow))
      }
      data_split <- .data
      id <- "id"
      
      ##add names here
      
    } else {
      data <- .data
      id   <- NULL
    }
  }
  
  # if(!is.null(id)) {
  #   # extract id column
  #   id_col <- data[, id] # a vector
  #   # data without id column
  #   data   <- data[, -which(colnames(data) == id)] 
  # }
  
  ## Choose resampling method
  out <- switch (.resample_method,
    "jackknife"   = {
      if(exists("data_split")) {
        lapply(data_split, function(y) 
          lapply(1:nrow(y), function(x) y[-x, ]))
      } else {
        lapply(1:nrow(data), function(x) data[-x, ]) 
      }
    },
    "bootstrap"   = {
      if(exists("data_split")) {
        lapply(data_split, function(y) 
          lapply(1:.R, function(x) {
            y[sample(1:nrow(y), size = nrow(y), replace = TRUE), ]}))
      } else {
        lapply(1:.R, function(x) {
          data[sample(1:nrow(data), size = nrow(data), replace = TRUE), ]})
      }
    },
    "permutation" = {
      if(is.null(id)) {
        stop2(
          "The following error occured in the `resampleData()` function:\n",
          "No id column specified to permutate the data with."
        )
      } else {
        lapply(1:.R, function(x) cbind(data[,-which(colnames(data) == id)], "id" = sample(data$id)))
      }
    },
    "cross-validation" = {
      if(is.null(.cv_folds)) {
        # LOOCV (Leave one out CV).
        if(exists("data_split")) {
          lapply(data_split, function(y) 
            lapply(1:nrow(y), function(x) y[x, ]))
        } else {
          lapply(1:nrow(data), function(x) data[x, ])
        }
        
      } else {
        # k-fold cross-validation (=draw k samples of equal size.).
        # Note the last sample may contain less observations if equal sized
        # samples are not possible
        if(exists("data_split")) {
          lapply(data_split, function(y) 
            suppressWarnings(
              split(as.data.frame(y), rep(1:.cv_folds, 
                                          each = ceiling(nrow(y)/.cv_folds)))
            ))
        } else {
          suppressWarnings(
            split(as.data.frame(data), rep(1:.cv_folds, 
                                           each = ceiling(nrow(data)/.cv_folds)))
          )
        }
      }
    } # END cross-validation
  ) # END switch
  ## Return samples
  out
}

#' Resample cSEMResults 
#'
#' The function resamples a [cSEMResults] object using bootstrap or jackknife resampling.
#' 
#' The function essentially calls [csem()] on each of the *M* resamples (created via
#' [resampleData()]) and returns M estimates for each of a subset of practically useful 
#' resampled parameters/statistics computed by [csem()]. Currently, the following 
#' quantities are returned for each resample: 
#' \describe{
#' \item{Parameters}{Path estimates, Weight estimates, Loading estimates}
#' \item{Statistics}{The heterotrait-monotrait ratio (HTMT)}
#' }
#' 
#' If the user needs to resample a statistic that is not returned by default, 
#' this statistic can be provided by a function `f(.object)` via the `.user_fun` argument. 
#' The only accepted argument of this function is `.object` which must be an
#' object of class [cSEMResults]. Internally, the function will be applied on each  
#' cSEMResults resample to produce the desired statistic. As long as this is 
#' the only argument of the function provided, arbitrary complicated statistics
#' may be resampled.
#' 
#' The number of bootstrap runs is given by `.R`. The default is
#' `499` but should be increased in real applications. See e.g.,
#' \insertCite{Hesterberg2015;textual}{cSEM}, p.380 for recommendations.
#' For jackknife `.R` is ignored.
#'
#' @usage resamplecSEMResults(
#'  .object                = NULL,
#'  .resample_method       = c("none", "bootstrap", "jackknife), 
#'  .resample_method2      = c("none", "bootstrap", "jackknife), 
#'  .R                     = 499,
#'  .R2                    = 199,
#'  .handle_inadmissibles  = c("drop", "ignore", "replace"),
#'  .user_funs             = NULL,
#'  .eval_plan             = c("sequential", "multiprocess")
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
  .resample_method       = args_default()$.resample_method,
  .resample_method2      = args_default()$.resample_method2,
  .R                     = args_default()$.R,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .user_funs             = args_default()$.user_funs,
  .eval_plan             = args_default()$.eval_plan
  ) {
  
  ## Set plan on how to resolve futures 
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  future::plan(.eval_plan)
    
  ### Checks, warnings and errors -----------
  match.arg(.resample_method, args_default(.choices = TRUE)$.resample_method)
  match.arg(.resample_method2, args_default(.choices = TRUE)$.resample_method2)
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
  
  ## Has the object to use the data to resample from produced admissible results?
  if(sum(verify(.object)) != 0) {
    warning2(
      "The following issue was encountered in the `resamplecSEMResults()` functions:\n",
      "Estimation based on the original data has produced inadmissible results.\n", 
      "This may be a sign that something is wrong.",
      " Resampling will continue but may not produce valid results.")
  }
  
  ## Check for the minimum number of necessary resamples
  if(.R < 3 | .R2 < 3) {
    stop2("The following error occured in the `resamplecSEMResults()` function:\n",
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
  Est_original$R2                     <- NULL
  Est_original$R2adj                  <- NULL
  Est_original$VIF                    <- NULL
  
  ## Additional statistics to compute by default
  # HTMT
  htmt <- c(HTMT(.object))
  Est_original[[length(Est_original) + 1]] <- htmt
  names(Est_original)[length(Est_original)] <- "HTMT"
  
  ## Add output of the user functions to Est_original
  if(!is.null(.user_funs)) {
    Est_original <- c(Est_original, user_funs)
  }
  
  ## Get argument list
  args <- .object$Information$Arguments
  
  ### Resample and compute -----------------------------------------------------
  
  out <- resamplecSEMResultsCore(
    .object                = .object, 
    .resample_method       = .resample_method,
    .resample_method2      = .resample_method2,
    .R                     = .R,
    .R2                    = .R2,
    .handle_inadmissibles  = .handle_inadmissibles,
    .user_funs             = .user_funs,
    .eval_plan             = .eval_plan
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
  
  if(.resample_method2 != "none") {
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
  #      2.3 Number_of_observations := the number of observations.
  #      2.4 Original_object        := the original cSEMResults object
  #
  # Since resamplecSEMResults is called recursivly within its body it is 
  # difficult to produce the desired output. There is now straightforward way 
  # to determine if a call is a recursive call or not
  # My workaround for now is to simply check if the argument provided for ".object"
  # is "Est_temp" since this is the argument for the recurive call. This
  # is of course not general but a dirty workaround. I guess it is save enough
  # though.
  is_recursive_call <- eval.parent(as.list(match.call()))$.object == "Est_temp"
  
  ## Return
  if(is_recursive_call) {
    out
  } else {
    if(.resample_method2 != "none") {
      out <- list("Estimates1" = out, "Estimates2" = out_2)
      
    } else {
      out <- list("Estimates1" = out, "Estimates2" = NA)
    }
    
    ## Add resamples and additional information to .object
    # Estimates
    estim <- c(.object$Estimates)
    estim[[length(estim) + 1]] <- c(out)
    names(estim)[length(estim)] <- "Estimates_resample"
    
    # Information
    info <- c(.object$Information)
    info[[length(info) + 1]] <- list(
      "Method"                  = .resample_method,
      "Method2"                 = .resample_method2,
      "Number_of_observations"  = nrow(.object$Information$Data)
    )
    names(info)[length(info)] <- "Information_resample"
    
    out <- list(
      "Estimates"   = estim,
      "Information" = info
    )
   
    ## Add/ set class
    class(out) <- c(class(.object), "cSEMResults_resampled")
    return(out)
  }
}

#' Core tasks of the resamplecSEMResults function
#' @noRd
#' 
resamplecSEMResultsCore <- function(
  .object                = args_default()$.object, 
  .resample_method       = args_default()$.resample_method,
  .resample_method2      = args_default()$.resample_method2,
  .R                     = args_default()$.R,
  .R2                    = args_default()$.R2,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .user_funs             = args_default()$.user_funs,
  .eval_plan             = args_default()$.eval_plan
) {
  
  ## Get arguments
  args <- .object$Information$Arguments
  
  ## Resample jackknife
  if(.resample_method == "jackknife") {
    resample_jack <- resampleData(.object, .resample_method = "jackknife") 
    .R <- length(resample_jack)
  }
  
  Est_ls <- future.apply::future_lapply(1:.R, function(i) {
    # Replace the old dataset by a resampled data set (resampleData always returns
    # a list so for just one draw we need to pick the first list element)
    
    data_temp <- if(.resample_method == "jackknife") {
      resample_jack[[i]]
    } else {
      resampleData(.object, .resample_method = "bootstrap", .R = 1)[[1]]
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
      x1$R2                      <- NULL
      x1$R2adj                   <- NULL
      x1$VIF                     <- NULL
      
      ## Additional statistics
      htmt <- c(HTMT(Est_temp))
      x1[[length(x1) + 1]] <- htmt
      names(x1)[length(x1)] <- "HTMT"
      
      ## Add output of the user functions to Est_original
      if(!is.null(.user_funs)) {
        x1 <- c(x1, user_funs)
      }
      
      ## Resampling from a bootstrap sample is required for the
      ## bootstraped t-interval CI, hence the second run
      if(.resample_method2 != "none") {
        
        Est_resamples2 <- resamplecSEMResults(
          .object               = Est_temp,
          .R                    = .R2,
          .handle_inadmissibles = .handle_inadmissibles,
          .resample_method      = .resample_method2,
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
  if(length(out) != .R && .resample_method == "bootstrap" && 
     .handle_inadmissibles == "replace") {
    
    R_new <- .R - length(out)
    
    while (length(out) < .R) {
      Est_replace <- resamplecSEMResultsCore(
        .object               = .object,
        .R                    = R_new,
        .handle_inadmissibles = .handle_inadmissibles,
        .resample_method      = .resample_method,
        .resample_method2     = .resample_method2,
        .R2                   = .R2,
        .user_funs            = .user_funs,
        .eval_plan            = .eval_plan
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
  .alpha           = args_default()$.alpha,
  .bias_corrected  = args_default()$.bias_corrected,
  .statistic       = args_default()$.statistic,
  ...
) {
  
  if(exists(.object$Estimates$Estimates_resamples)) {
    .csem_resample <- resamplecSEMResults(
      .object = .object,
      ...
    )
  }
  
  .alpha <- c(0.1, 0.5, 0.234)
  resamples1 <- a$Estimates$Estimates_resample$Estimates1
  resamples2 <- a$Estimates$Estimates_resample$Estimates2
  resamples1 <- .csem_resample$Resamples$Estimates1 # jackknife or bootstrap
  resamples2 <- .csem_resample$Resamples$Estimates2 # jackknife or bootstrap
  info       <- a$Information$Information_resample
  object     <- a

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
                         .resample_method = info$Method, 
                         .n      = info$Number_of_observations
                         ),
    "Sd"   = SdResample(.x = resamples1, 
                        .resample_method = info$Method, 
                        .n = info$Number_of_observations
                        ),
    "Bias" = BiasResample(.x = resamples1, 
                          .resample_method = info$Method, 
                          .n = info$Number_of_observations
                          ),
    "CI_standard_z" = StandardCIResample(
      .x              = resamples1, 
      .bias_corrected = .bias_corrected,
      .df             = NULL,
      .dist           = "z",
      .resample_method         = info$Method, 
      .n              = info$Number_of_observations,
      .probs          = probs
      ),
    "CI_standard_t" = StandardCIResample(
      .x              = resamples1, 
      .bias_corrected = .bias_corrected,
      .df             = "type1",
      .dist           = "t",
      .resample_method         = info$Method, 
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
    "CI_Bca" = BcaCI(.object = object, resamples1, probs),
    "CI_t_invertall"= if(anyNA(resamples2)) {
      NA
    } else {
      TStatCIResample(
        .x1      = resamples1, 
        .x2      = resamples2, 
        .bias_corrected = .bias_corrected,
        .resample_method  = info$Method, 
        .resample_method2 = info$Method2, 
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
VarResample <- function(.x, .resample_method, .n) {

  lapply(.x, function(x) {
    out        <- matrixStats::colVars(x$Resampled)
    names(out) <- names(x$Original)
    if(.resample_method == "jackknife") {
      out <- out * (.n - 1)^2/.n 
    }
    out
  })
}
#' @describeIn infer (TODO)
SdResample <- function(.x, .resample_method, .n) {
  
  lapply(.x, function(x) {
    out        <- matrixStats::colSds(x$Resampled)
    names(out) <- names(x$Original)
    if(.resample_method == "jackknife") {
      out <- out * (.n - 1)/sqrt(.n)
    }
    out
  })
}
#' @describeIn infer (TODO)
BiasResample <- function(.x, .resample_method, .n) {
  
  lapply(.x, function(x) {
    out        <- colMeans(x$Resampled) - x$Original
    names(out) <- names(x$Original)
    if(.resample_method == "jackknife") {
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
  .resample_method, 
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
  boot_sd <- SdResample(.x, .resample_method = .resample_method, .n = .n)
  
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
  .resample_method, 
  .resample_method2, 
  .n, 
  .probs
  ) {
  # Bootstraped t statistic
  # Bootstrap the t-statisic (since it is roughly pivotal) and compute
  # CI based on bootstraped t-values and bootstraped/jackknife SE
  # CI: [F*^-1(alpha/2) ; F*^-1(1 - alpha/2)]

  ## confidence level (for rownames)
  cl <- 1 - .probs[seq(1, length(.probs), by = 2)]*2
  
  boot_sd_star <- lapply(.x2, SdResample, .resample_method = .resample_method2, .n = .n) %>% 
    purrr::transpose(.) %>% 
    lapply(function(y) do.call(rbind, y))
  
  boot_sd <- SdResample(.x1, .resample_method = .resample_method, .n = .n)
  
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
    .resample_method = "jackknife"
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
