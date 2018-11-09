#' Internal: Resample
#'
#' Resample from a given cSEMResults object or data.frame using common
#' resampling methods.
#'
#' @usage resample(
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
resample <- function(
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
