#' Parse lavaan model
#'
#' Turns a model written in [lavaan model syntax][lavaan::model.syntax] into a
#' [cSEMModel] list.
#'
#' @usage parseModel(.model)
#'
#' @inheritParams csem_arguments
#
#' @return A [cSEMModel] list that may be passed to any function requiring
#'   `.csem_model` as a mandatory argument.
#'
#' @examples
#' model <- '
#' # Structural model
#' y1 ~ y2 + y3
#'
#' # Measurement model
#' y1 =~ x1 + x2 + x3
#' y2 =~ x4 + x5
#' y3 =~ x6 + x7
#'
#' # Error correlation
#' x1 ~~ x2
#' '
#'
#' (m <- parseModel(model))
#'
#' # If the model is already a cSEMModel object, the model is returned as is.
#'
#' identical(m, parseModel(m)) # TRUE
#' @export
#'
parseModel <- function(.model) {

  ### Check if already a cSEMModel list  
  if(class(.model) == "cSEMModel") {
    # .model <- model
    return(.model)
  
    ## Check if list contains necessary elements
  } else if(all(c("structural", "measurement") %in% names(.model))) {
    
    x <- setdiff(names(.model), c("structural", "measurement", "error_cor", 
                                  "construct_type", "construct_order", "model_type", "vars_endo", 
                                  "vars_exo", "vars_explana", "explained_by_exo"))
    if(length(x) == 0) {
      
      class(.model) <- "cSEMModel"
      return(.model)
      
    } else {
      
      stop("The model you provided contains element names unknown to cSEM (", 
           paste0("'", x, "'", collapse = ", "), ").\n", 
           "See ?cSEMModel for a list of valid component names.", call. = FALSE)
    }
  } else {
    
    ### Convert to lavaan partable -----------------------------------------------
    m_lav <- lavaan::lavaanify(model = .model, fixed.x = FALSE)
    
    ### Extract relevant information ---------------------------------------------
    tbl_structural    <- m_lav[m_lav$op == "~", ]
    tbl_measurement   <- m_lav[m_lav$op %in% c("=~", "<~"), ]
    tbl_errors        <- m_lav[m_lav$op == "~~" & m_lav$user == 1, ]
    
    names_constructs_all  <- unique(c(tbl_measurement$lhs, tbl_structural$rhs))
    names_constructs      <- unique(unlist(strsplit(names_constructs_all, "\\.")))
    names_constructs_structural <- unique(c(tbl_structural$lhs, unlist(strsplit(tbl_structural$rhs, "\\."))))
    names_indicators      <- unique(tbl_measurement$rhs)
    names_constructs_measurement_rhs <- intersect(names_indicators, names_constructs)
    names_constructs_second_order <- unique(tbl_measurement[tbl_measurement$rhs %in% names_constructs_measurement_rhs, "lhs"])
    
    number_of_constructs_all  <- length(names_constructs_all)
    number_of_constructs      <- length(names_constructs)
    number_of_indicators      <- length(names_indicators)
    
    ## Order
    construct_order <- rep("First order", length(names_constructs))
    names(construct_order) <- names_constructs
    construct_order[names_constructs_second_order] <- "Second order"
    
    ### Checks, errors and warnings ----------------------------------------------
    ## Stop if construct has no obervables/indicators attached
    if(length(setdiff(c(names_constructs_structural, names_constructs_measurement_rhs), 
                      tbl_measurement$lhs)) != 0) {
      if(any(grepl("\\.", tbl_structural$lhs))) {
        stop("Interaction terms cannot appear on the left-hand side of a structural equation.", 
             call. = FALSE)
      } else {
        stop("No measurement equation provided for: ",
             paste0("`", setdiff(names_constructs_structural, 
                                 tbl_measurement$lhs), "`", collapse = ", "),
             call. = FALSE)
      }
    }
    
    ## Stop if construct appears in the measurement but not in the structural model
    if(length(setdiff(tbl_measurement$lhs, 
                      c(names_constructs_structural, 
                        names_constructs_measurement_rhs)) != 0)) {
      ## Stop if user trys to specify a measurement equation for an interaction term
      if(any(grepl("\\.", tbl_measurement$lhs))) {
        stop("Interaction terms cannot appear on the left-hand side of a measurement equation.",
             call. = FALSE)
      } else {
        stop("Construct(s): ",
             paste0("`", setdiff(tbl_measurement$lhs, 
                                 c(names_constructs_structural,
                                   names_constructs_measurement_rhs)), 
                    "`", collapse = ", "),
             " of the measurement model",
             ifelse(length(setdiff(tbl_measurement$lhs, 
                                   c(names_constructs_structural,
                                     names_constructs_measurement_rhs))) == 1, 
                    " does", " do"),
             " not appear in the structural model.",
             call. = FALSE)
      }
    }
    
    ## Stop rhs of second order constructs contains indicators
    if(length(intersect(setdiff(names_indicators, 
                                names_constructs_measurement_rhs), 
                        tbl_measurement[tbl_measurement$lhs 
                                        %in% names_constructs_second_order, "rhs"])) != 0) {
      stop("Second-order constructs must be defined by first-order constructs only.",
           call. = FALSE)
    }
    
    ## Stop if construct has a higher order than 2 (currently not allowed)
    if(length(intersect(names_constructs_second_order, 
                        names_constructs_measurement_rhs)) != 0) {
      stop(paste0("`", unique(tbl_measurement[
        tbl_measurement$rhs == intersect(names_constructs_second_order, 
                                         names_constructs_measurement_rhs), 
        "lhs"]), "`"), " has order > 2. Currently, only first and second-order constructs are supported.",
           call. = FALSE)
    }
    ## Construct type
    
    tbl_measurement$op <- ifelse(tbl_measurement$op == "=~", "Common factor", "Composite")
    construct_type  <- unique(tbl_measurement[, c("lhs", "op")])$op
    names(construct_type) <- names_constructs
    
    ## Type of model (linear or non-linear)
    
    type_of_model <- if(any(grepl("\\.", names_constructs_all))) {
      "Nonlinear"
    } else {
      "Linear"
    }
    ### Construct matrices specifying the relationship between constructs,
    ### indicators and errors ----------------------------------------------------
    model_structural  <- matrix(0,
                                nrow = number_of_constructs,
                                ncol = number_of_constructs_all,
                                dimnames = list(names_constructs, names_constructs_all)
    )
    
    model_measurement <- matrix(0,
                                nrow = number_of_constructs,
                                ncol = number_of_indicators,
                                dimnames = list(names_constructs, names_indicators)
    )
    
    model_error       <- matrix(0,
                                nrow = number_of_indicators,
                                ncol = number_of_indicators,
                                dimnames = list(names_indicators, names_indicators)
    )
    
    ## Structural model
    row_index <- match(tbl_structural$lhs, names_constructs)
    col_index <- match(tbl_structural$rhs, names_constructs_all)
    
    model_structural[cbind(row_index, col_index)] <- 1
    
    ## Measurement model
    row_index <- match(tbl_measurement$lhs, names_constructs)
    col_index <- match(tbl_measurement$rhs, names_indicators)
    
    model_measurement[cbind(row_index, col_index)] <- 1
    
    ## Error model
    row_index <- match(tbl_errors$lhs, names_indicators)
    col_index <- match(tbl_errors$rhs, names_indicators)
    
    model_error[cbind(c(row_index, col_index), c(col_index, row_index))] <- 1
    
    ### Order model ==============================================================
    # Order the structual equations in a way that every equation depends on
    # exogenous variables and variables that have been explained in a previous equation
    # This is necessary for the estimation of models containing non-linear structual
    # relationships.
    
    ### Preparation --------------------------------------------------------------
    temp <- model_structural
    
    ## Extract endogenous and exogenous variables
    var_endo <- rownames(temp)[rowSums(temp) != 0]
    var_exo  <- setdiff(colnames(temp), var_endo)
    
    ## Modify temp to ease computation below
    for(i in var_endo) {
      
      # Extract right-hand side variable names for each structrual equation
      x <- colnames(temp)[temp[i, ] != 0]
      # Split by "." and unlist
      x <- unlist(strsplit(x, "\\."))
      
      # For each structural equation: Get names of single constructs that are already
      # set to 1
      z <- names(temp[i, names_constructs][temp[i, names_constructs] == 1])
      
      if(length(setdiff(x, z)) > 0) {
        stop(paste0("The structural equation for `", i, "` contains interactions between constructs", 
                    " that do not appear individually.\n",
                    "This is almost certainly not a meaningful model."), call. = FALSE)
      }
      
      # If x containes an interaction term, assign 1 to all elements in temp[i, ] whose
      # column names match one or more of the elements/names of the splitted terms
      temp[i, intersect(x, colnames(temp))] <- 1
    }
    
    ## Return error if the structural model contains feedback loops
    if(any(temp[var_endo, var_endo] + t(temp[var_endo, var_endo]) == 2)) {
      stop("The structural model contains feedback loops.",
           " Currently no feedback loops are allowed.",
           call. = FALSE)
    }
    
    # Endo variables that are explained by exo and endo variables
    explained_by_exo_endo <- var_endo[rowSums(temp[var_endo, var_endo, drop = FALSE]) != 0]
    
    # Endo variables explained by exo variables only
    explained_by_exo <- setdiff(var_endo, explained_by_exo_endo)
    
    ### Order =======================
    # First the endo variables that are soley explained by the exo variables
    model_ordered <- temp[explained_by_exo, , drop = FALSE]
    
    # Add variables that have already been ordered/taken care of to a vector
    # (including exogenous variables and interaction terms)
    already_ordered <- c(var_exo, explained_by_exo)
    
    ## Order in a way that the current structural equation does only depend on
    ## exogenous variables and/or variables that have already been ordered
    counter <- 1
    explained_by_exo_endo_temp <- explained_by_exo_endo
    if(length(explained_by_exo_endo) > 0) {
      repeat {
        
        counter <- counter + 1
        
        for(i in explained_by_exo_endo_temp) {
          names_temp <- colnames(temp[i, temp[i, ] == 1, drop = FALSE])
          endo_temp  <- setdiff(names_temp, already_ordered)
          
          if(length(endo_temp) == 0) {
            model_ordered <- rbind(model_ordered, temp[i, , drop = FALSE])
            already_ordered <- c(already_ordered, i)
            explained_by_exo_endo_temp <- setdiff(explained_by_exo_endo_temp, already_ordered)
          }
        } # END for-loop
        if(counter > 50)
          stop("Reordering the structural equations was not succesful. Something is wrong.",
               call. = FALSE)
        if(length(explained_by_exo_endo_temp) == 0) break
      } # END repeat
    } # END if-statement
    
    ## Return a cSEMModel object.
    # A cSEMModel objects contains all the information about the model and its
    # components such as the type of construct used. 
    n <- c(setdiff(names_constructs, rownames(model_ordered)), rownames(model_ordered))
    m <- order(which(model_measurement[n, ] == 1, arr.ind = TRUE)[, "row"])
    structural_ordered <- model_structural[n, c(n, setdiff(colnames(model_ordered), n))]
    
    model_ls <- list(
      "structural"         = structural_ordered,
      # "structural_ordered" = model_ordered, # not needed so far
      "measurement"        = model_measurement[n, m],
      "error_cor"          = model_error[m, m],
      "construct_type"     = construct_type[match(n, names(construct_type))],
      "construct_order"    = construct_order[match(n, names(construct_order))],
      "model_type"         = type_of_model,
      "vars_endo"          = rownames(model_ordered),
      "vars_exo"           = var_exo,
      "vars_explana"       = colnames(structural_ordered)[colSums(structural_ordered) != 0],
      "explained_by_exo"   = explained_by_exo
    )
    class(model_ls) <- "cSEMModel"
    return(model_ls) 
  } # END else
}

#' Convert second order cSEMModel
#'
#' Uses a [cSEMModel] containg second order constructs and turns it into an
#' estimable model using either the "repeated indicators" approach or a 
#' two-step procedure (TODO; link to literature)
#'
#' @usage convertModel(.csem_model)
#'
#' @inheritParams csem_arguments
#
#' @return A [cSEMModel] list that may be passed to any function requiring
#'   `.csem_model` as a mandatory argument.
#'
#' @keywords internal
#'
convertModel <- function(
  .csem_model        = args_default()$.csem_model, 
  .approach_2ndorder = args_default()$.approach_2ndorder
) {
  
  ### Check if a cSEMModel list  
  if(!class(.csem_model) == "cSEMModel") {
    stop("`.model` must be of model of class `cSEMModel`.")
  }
  
  # Get relevant construct/indicator names
  nc_all <- rownames(.csem_model$structural)
  ni_all <- colnames(.csem_model$measurement)
  n1c <- names(.csem_model$construct_order[.csem_model$construct_order == "First order"])
  n2c <- setdiff(nc_all, n1c)
  n1i <- setdiff(ni_all, nc_all)
  n2i <- intersect(ni_all, nc_all)
  
  if(.approach_2ndorder == "repeated_indicators") {
    
    ## Structural .csem_model
    # First order equations
    x1 <- c()
    for(i in nc_all) {
      col_names <- colnames(.csem_model$structural[i, .csem_model$structural[i, , drop = FALSE] == 1, drop = FALSE])
      temp <- if(!is.null(col_names)) {
        paste0(i, "~", paste0(col_names, collapse = "+")) 
      } else {
        "\n"
      }
      x1 <- paste(x1, temp, sep = "\n")
    }
    
    ## Measurement model + second order structural equation 
    # First order constructs
    x2a <- c()
    for(i in n1c) {
      col_names <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
      temp  <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(col_names, collapse = "+"))
      x2a <- paste(x2a, temp, sep = "\n")
    }
    
    # Second order constructs
    x2b <- c()
    for(i in n2c) {
      # i <- n2c[1]
      col_names_1 <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
      col_names_2 <- .csem_model$measurement[col_names_1, colSums(.csem_model$measurement[col_names_1, ,drop = FALSE]) != 0, drop = FALSE]
      temp <- paste0(i, "_2nd_", colnames(col_names_2))
      temp <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(temp, collapse = "+"))
      x2b <- paste(x2b, temp, sep = "\n")
      ## add second order structural equation
      x2b <- paste(x2b, paste0(i, "~", paste0(col_names_1, collapse = "+" )), sep = "\n")
    }
    
    ## Error_cor
    # First order
    x3 <- c()
    for(i in n1i) {
      # - Only upper triagular matrix as lavaan does not allow for double entries such
      #   as x11 ~~ x12 vs x12 ~~ x11
      # - Only 1st order construct indicators are allowed to correlate 
      error_cor <- .csem_model$error_cor
      error_cor[lower.tri(error_cor)] <- 0
      col_names <- colnames(error_cor[i, error_cor[i, , drop = FALSE] == 1, drop = FALSE])
      temp <- if(!is.null(col_names)) {
        paste0(i, "~~", paste0(col_names, collapse = "+"))
      } else {
        "\n"
      }
      x3 <- paste(x3, temp, sep = "\n")
    }
    ## Model to be parsed
    lav_model <- paste(x1, x2a, x2b, x3, sep = "\n")
    
  } else { # First step of the two step approach
    
    ## Structural model
    x1 <- c()
    for(i in 2:length(n1c)) {
      temp <- paste0(n1c[i], "~", paste0(n1c[1:(i-1)], collapse = "+"))
      x1   <- paste(x1, temp, sep = "\n")
    }
    ## Measurement model
    x2 <- c()
    for(i in n1c) {
      col_names <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE] == 1, drop = FALSE])
      temp <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(col_names, collapse = "+"))
      x2   <- paste(x2, temp, sep = "\n")
    }
    ## Error_cor
    x3 <- c()
    for(i in n1i) {
      # only upper triagular matrix as lavaan does not allow for double entries such
      # as x11 ~~ x12 vs x12 ~~ x11
      error_cor <- .csem_model$error_cor
      error_cor[lower.tri(error_cor)] <- 0
      col_names <- colnames(error_cor[i, error_cor[i, , drop = FALSE] == 1, drop = FALSE])
      temp <- if(!is.null(col_names)) {
        paste0(i, "~~", paste0(col_names, collapse = "+"))
      } else {
        "\n"
      }
      x3 <- paste(x3, temp, sep = "\n")
    }
    ## Model to be parsed
    lav_model <- paste(x1, x2, x3, sep = "\n")
  }
  ## Parse model
  model <- parseModel(lav_model)
  return(model)
}