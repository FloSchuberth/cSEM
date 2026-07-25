#' Internal: Process data
#'
#' Prepare, standardize, check, and clean data provided via the `.data` argument.
#'
#' @usage processData(
#'   .data            = NULL, 
#'   .model           = NULL, 
#'   .instruments     = NULL,
#'   .handle_missing  = NULL
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
  .data            = NULL, 
  .model           = NULL, 
  .instruments     = NULL,
  .handle_missing  = NULL
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
  
  missing_info <- list(
    "Approach_missing" = .handle_missing,
    "Missing_data" = length(rownames(.data_temp)) > 0,
    "Rows_with_missing" = rownames(.data_temp),
    "Number_of_rows_missing" = length(rownames(.data_temp))
  )
  
  if(nrow(.data_temp) > 0) {
    if (.handle_missing == "none") {
      stop2(
        "The following error occured while processing the data:\n",
        "Data set contains missing values in rows:",
        paste0("`", rownames(.data_temp), "`", collapse = ", "),
        "\nRemove NAs, use imputation methods to replace them, or set `.handle_missing` to \"listwise\", \"mean\", or \"regression\"."
      )
      
    } else if (.handle_missing == "listwise") {
      .data <- .data[complete.cases(.data), , drop = FALSE]
    } else if (.handle_missing == "mean") {
      .data <- imputeMissingMean(.data)
    } else if (.handle_missing == "regression") {
      .data <- imputeMissingRegression(.data)
    }
    missing_info$Number_of_rows_imputed <- nrow(.data_temp)
  }
  attr(.data, "missing_info") <- missing_info
  ## Return
  return(.data)
}

#' Internal: Mean imputation for missing data
#'
#' @keywords internal
#' @noRd
imputeMissingMean <- function(.data) {
  missing_columns <- names(.data)[unlist(lapply(.data, anyNA))]
  nonnumeric_columns <- missing_columns[!unlist(lapply(.data[missing_columns], is.numeric))]
  
  if(length(nonnumeric_columns) != 0) {
    stop2("Mean replacement currently requires numeric indicators. Non-numeric indicators with missing values: ",
          paste0("`", nonnumeric_columns, "`", collapse = ", "), ".")
  }
  
  for(i in missing_columns) {
    replacement <- mean(.data[[i]], na.rm = TRUE)
    if(is.nan(replacement)) {
      stop2("Mean replacement failed because indicator `", i, "` contains only missing values.")
    }
    .data[is.na(.data[[i]]), i] <- replacement
  }
  
  return(.data)
}

#' Internal: Regression imputation for missing data
#'
#' @keywords internal
#' @noRd
imputeMissingRegression <- function(.data) {
  missing_columns <- names(.data)[unlist(lapply(.data, anyNA))]
  nonnumeric_columns <- names(.data)[!unlist(lapply(.data, is.numeric))]
  
  if(length(nonnumeric_columns) != 0) {
    stop2("Regression imputation currently requires all indicators to be numeric. Non-numeric indicators: ",
          paste0("`", nonnumeric_columns, "`", collapse = ", "), ".")
  }
  
  if(ncol(.data) < 2) {
    stop2("Regression imputation requires at least two indicators.")
  }
  
  data_imputed <- imputeMissingMean(.data)
  
  for(i in missing_columns) {
    observed_rows <- !is.na(.data[[i]])
    missing_rows  <- is.na(.data[[i]])
    predictors    <- setdiff(names(.data), i)
    
    if(sum(observed_rows) <= length(predictors)) {
      stop2("Regression imputation for indicator `", i, "` requires more complete observations than predictors.")
    }
    
    model_data <- data_imputed[observed_rows, c(i, predictors), drop = FALSE]
    fit <- stats::lm(stats::reformulate(predictors, response = i), data = model_data)
    predicted <- stats::predict(fit, newdata = data_imputed[missing_rows, predictors, drop = FALSE])
    
    if(anyNA(predicted)) {
      stop2("Regression imputation failed for indicator `", i, "`.")
    }
    
    data_imputed[missing_rows, i] <- predicted
  }
  
  return(data_imputed)
}
