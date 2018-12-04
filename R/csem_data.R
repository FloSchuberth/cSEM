#' Process data
#'
#' Prepare, standardize, check, and clean data provided via the `.data` argument.
#'
#' @usage processData(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @return A (N x K) data.frame containing the standardized data with columns ordered
#'   according to the order they appear in the measurement model equations provided
#'   via the `.model` argument.
#'
#' @examples
#' dat <- data.frame(x2 = rnorm(100),
#'                   x1 = rnorm(100),
#'                   b2 = rnorm(100),
#'                   b1 = rnorm(100),
#'                   a1 = rnorm(100),
#'                   a2 = rnorm(100)
#'                   )
#'
#' model <- '
#' # Structural model
#' y1 ~ y2 + y3
#'
#' # Measurement model
#' y1 =~ x1 + x2
#' y2 =~ a1 + a2
#' y3 =~ b1 + b2
#' '
#'
#' dat <- processData(dat, model)
#'
#' is.matrix(dat) # TRUE
#' colnames(dat)  # "x1", "x2", "a1", "a2", "b1", "b2"
#'
#' @export

processData <- function(.data, .model) {

  ### Checks, errors and warnings ========
  # Check if any data set is provided
  if(is.null(.data)) {
    stop("No data set provided. Please provide a `data.frame` or a `matrix` of data.",
         call. = FALSE)
  }

  # Check class of the .data object and stop if not of class "data.frame" or "matrix"
  if(!any(class(.data) %in% c("data.frame", "matrix"))) {
    stop("Don't know how to deal with a data object of class: ", class(.data), ".\n",
         "Please provide the data as a `data.frame` or a `matrix`.",
         call. = FALSE)
  }

  # Check if any of the columns is a character. 
  # Allowed types: numeric (double, integer), factor (ordered and unordered), or logical 
  x <- names(which(sapply(.data, is.character)))
  if(length(x) == 1) {
    stop("Column: ",paste0("`", x, "`", collapse = ", "), " of `.data` is of type `character`.\n",
         "Have you forgotten to set `", paste0(".id = '", x, "'`?"), call. = FALSE)
  } else if(length(x) > 1) {
    stop("Columns: ",paste0("`", x, "`", collapse = ", "), "of `.data` are of type `character`.\n",
         "The column type must be one of: `logical`, `numeric`, or `factor`",
         call. = FALSE)
  }

  ### Processing and further checking =========
  # Convert to data.frame if matrix
  # Note we need a data frame to allow for data to have different classes. Namely,
  # factors need to be allowed.
  
  if(is.matrix(.data)) {
    .data <- as.data.frame(.data)
  }

  # Convert .model to cSEMModel format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    .model <- parseModel(.model)
  }
  ## Add indicators to .data if the repeated indicators approach is used
  # Error:
  # Note: the indicators to be added are identified by the string "_2nd_". Hence
  # the string is basically a reserved word. If indicators supplied by the
  # users contain the string this causes and error (unlikely to happen).
  if(any(grepl("_2nd_", colnames(.data)))) {
    stop("Indicator names must not contain the string `_2nd_`.", call. = FALSE)
  }
  
  names_2nd <- colnames(.model$measurement)[grep("_2nd", colnames(.model$measurement))]
  
  if(length(names_2nd) > 0) {
    temp <- do.call(rbind, strsplit(names_2nd, "_2nd_"))
    
    temp <- .data[, temp[, 2]]
    colnames(temp) <- names_2nd
    
    ## extended .data
    .data <- cbind(.data, temp) 
  }
  
  ## Check indicator names
  if(!all(colnames(.model$measurement) %in% colnames(.data))) {
    stop("Unknown indicator(s): ",  
         paste0("`", setdiff(colnames(.model$measurement), colnames(.data)), "`.", collapse = ", "),
         " Please verify your model description.",
         call. = FALSE)
  }
  # Order data according to the ordering of the measurement model
  .data <- .data[, colnames(.model$measurement)]

  ## Return
  return(.data)
}
