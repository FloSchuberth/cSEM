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
    
    return(.model)
  
    ## Check if list contains necessary elements
  } else if(all(c("structural", "measurement") %in% names(.model))) {
    
    x <- setdiff(names(.model), c("structural", "measurement", "error_cor", 
                                  "construct_type", "model_type", "vars_endo", 
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
    
    number_of_constructs_all  <- length(names_constructs_all)
    number_of_constructs      <- length(names_constructs)
    number_of_indicators      <- length(names_indicators)
    
    ### Checks, errors and warnings ----------------------------------------------
    ## Stop if construct has no obervables/indicators attached
    if(length(setdiff(names_constructs_structural, tbl_measurement$lhs)) != 0) {
      
      stop("No measurement equation provided for: ",
           paste0("`", setdiff(names_constructs_structural, tbl_measurement$lhs), "`", collapse = ", "),
           call. = FALSE)
    }
    
    ## Stop if construct appears in the measurement but not in the structural model
    if(length(setdiff(tbl_measurement$lhs, names_constructs_structural)) != 0) {
      
      stop("Construct(s): ",
           paste0("`", setdiff(tbl_measurement$lhs, names_constructs_structural), "`", collapse = ", "),
           " of the measurement model",
           ifelse(length(setdiff(tbl_measurement$lhs, names_constructs_structural)) == 1, " does", " do"),
           " not appear in the structural model.",
           call. = FALSE)
    }
    
    ## Construct type
    
    tbl_measurement$op <- ifelse(tbl_measurement$op == "=~", "Common factor", "Composite")
    type_of_construct  <- unique(tbl_measurement[, c("lhs", "op")])$op
    names(type_of_construct) <- names_constructs
    
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
        stop(paste0("The structural equation for ", i, " contains interactions between constructs", 
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
    
    model_ls <- list(
      "structural"         = model_structural[n, c(n, setdiff(colnames(model_ordered), n))],
      # "structural_ordered" = model_ordered, # not needed so far
      "measurement"        = model_measurement[n, m],
      "error_cor"          = model_error[m, m],
      "construct_type"     = type_of_construct[match(n, names(type_of_construct))],
      "model_type"         = type_of_model,
      "vars_endo"          = rownames(model_ordered),
      "vars_exo"           = var_exo,
      "vars_explana"       = colnames(model_structural[n, c(n, setdiff(colnames(model_ordered), n))][, colSums(model_structural[n, c(n, setdiff(colnames(model_ordered), n))]) != 0 ]),
      "explained_by_exo"   = explained_by_exo
    )
    class(model_ls) <- "cSEMModel"
    return(model_ls) 
  } # END else
}

