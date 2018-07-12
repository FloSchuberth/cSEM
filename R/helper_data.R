#' Process data
#'
#' Prepare, standardize, check, and clean data provided via the `.data` argument.
#'
#' @usage processData(.data, .model)
#'
#' @inheritParams csem_arguments
#'
#' @return A (N x K) matrix containing the standardized data with columns ordered
#'   according to the order they appear in the measurement model equations provided
#'   via the `.model` argument.
#'
#' @examples
#' require(cSEM)
#'
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

processData <- function(.data, .model, .standardize) {

  ### Checks, errors and warnings ========
  # Check if any data set is provided
  if(is.null(.data)) {
    stop("No data set provided. Please provide a data.frame or a matrix of data.")
  }

  # Check class of the .data object and stop if not of class "data.frame" or "matrix"

  if(!(class(.data) %in% c("data.frame", "matrix"))) {
    stop("Don't know how to deal with a data object of class: ", class(.data), ".\n",
         "Please provide the data as a data.frame or a matrix.")
  }

  # Check if all columns are numeric. Stop otherwise
  if(!all(apply(.data, 2, is.numeric))) {
    stop("At least one column of the data is non-numeric.")
  }

  ### Processing =========
  # Convert to matrix if data.frame
  if(is.data.frame(.data)) {
    .data <- as.matrix(.data)
  }

  # Convert to cSEMModel format if not already in this format
  if(!(class(.model) == "cSEMModel")) {
    .model <- parseModel(.model)
  }
  # Order data according to the ordering of the measurement model
  .data <- .data[, colnames(.model$measurement)]

  # Scale data
  if(.standardize==TRUE){
  .data <- scale(.data) # Note: it does not matter if the data is already scaled
                        # Scaling a scaled/standardized data set does not alter
                        # the data.
  }
  ## Set class
  class(.data) <- "cSEMData"

  ## Return
  return(.data)
}

