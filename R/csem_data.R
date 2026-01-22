#' Internal: Process data
#'
#' Prepare, standardize, check, and clean data provided via the `.data` argument.
#'
#' @usage processData(
#'   .data        = NULL, 
#'   .model       = NULL, 
#'   .instruments = NULL
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @return A (N x K) data.frame containing the standardized data with columns ordered
#'   according to the order they appear in the measurement model equations provided
#'   via the `.model` argument.
#'
#' @keywords internal

processData <- function(
  .data        = NULL, 
  .model       = NULL, 
  .instruments = NULL
  ) {

  ### Checks, errors and warnings ========
  # Check if any data set is provided
  if(is.null(.data)) {
    stop2("No data set provided. Please provide a `data.frame` or a `matrix` of data.")
  }

  # Check class of the .data object and stop if not of class "data.frame" or "matrix"
  if(!inherits(.data, c("data.frame", "matrix"))) {
    stop2("Don't know how to deal with a data object of class: ", class(.data), ".\n",
          "Please provide the data as a `data.frame` or a `matrix`.")
  }

  ### Processing and further checking =========
  # Convert to data.frame
  # Note 1: we need a data frame to allow for data to have different classes. Namely,
  #   factors need to be allowed.
  # Note 2: previously as.data.frame() was only called when .data had class
  #    .matrix. However, classes tbl_df, tbl and tibble cause e.g. an error in the
  #    hetcor() function (its basically a programming error on the developers part,
  #    as it checks the class attribute incorrectly). 
  #    Hence, as.data.frame is now always called to make sure .data is
  #    always really a data frame with the single class attribute "data.frame".
    
  .data <- as.data.frame(.data)

  # Convert .model to cSEMModel format if not already in this format
  if(!inherits(.model, "cSEMModel")) {
    .model <- parseModel(.model, .instruments = .instruments)
  }
  
  ## Check if data set is symmetric. This is an indicator that a covariance
  ## matrix has been supplied (which is not supported by cSEM):
  if(
    matrixcalc::is.square.matrix(as.matrix(.data)) &&
    matrixcalc::is.symmetric.matrix(as.matrix(.data))) {
    warning2("Data is symmetric! Did you provide a covariance or correlation matrix to `.data`?\n",
             "Argument `.data` requires a matrix or data.frame of raw data.")
  }
  
  
  ## Check if any of the columns are character and convert them to factors
  # Allowed types: character (converted to factor), numeric (double, integer), 
  # factor (ordered and unordered), or logical 
  x <- names(.data[unlist(lapply(.data, is.character))])
  
  if(length(x) != 0) {
    for(i in x) {
      .data[, i] <- as.factor(.data[, i])
    }
  }

  # if(length(x) == 1) {
  #   stop("Column: ",paste0("`", x, "`", collapse = ", "), " of `.data` is of type `character`.\n",
  #        "Have you forgotten to set `", paste0(".id = '", x, "'`?"), call. = FALSE)
  # } else if(length(x) > 1) {
  #   stop("Columns: ",paste0("`", x, "`", collapse = ", "), "of `.data` are of type `character`.\n",
  #        "The column type must be one of: `logical`, `numeric`, or `factor`",
  #        call. = FALSE)
  # }
  
  ## Add indicators to .data if the repeated indicators approach is used
  # Note: the indicators to be added are identified by the string "_2nd_". Hence
  # the string is basically a reserved word. If indicators supplied by the
  # users contain the string this causes and error (unlikely to happen).
  if(any(grepl("_2nd_", colnames(.data)))) {
    stop2("Indicator names must not contain the string `_2nd_`.")
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
    stop2("The following error occured while processing the data:\n",
      "Unknown indicator(s): ", paste0("`", setdiff(colnames(.model$measurement), 
                                                    colnames(.data)), "`.", collapse = ", "),
         " Please verify your model description.")
  }
  # Order data according to the ordering of the measurement model; delete
  # all columns that are not needed
  .data <- .data[, colnames(.model$measurement)]

  # Check if remaining data set contains NAs
  .data_temp <- .data[!complete.cases(.data), , drop = FALSE]
  
  if(length(rownames(.data_temp)) > 0) {
    stop2("The following error occured while processing the data:\n",
          "Data set contains missing values in rows:", 
          paste0("`", rownames(.data_temp), "`", collapse = ", "),
          "\nRemove NAs or use imputation methods to replace them.")
  }
  ## Return
  return(.data)
}
