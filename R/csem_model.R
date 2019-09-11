#' Parse lavaan model
#'
#' Turns a model written in [lavaan model syntax][lavaan::model.syntax] into a
#' [cSEMModel] list.
#'
#' Instruments must be supplied seperately as a named list 
#' of vectors of instruments. 
#' The names of the list elements are the names of the dependent constructs of 
#' the structural equation whose explanatory variables are endogenous. 
#' The vectors contain the names of the instruments corresponding to each 
#' equation. Note that exogenous variables of a given equation must be 
#' supplied as instruments for themselves.
#' 
#' By default `parseModel()` attempts to check if the model provided is correct
#' in a sense that all necessary components required to estimate the
#' model are specified (e.g., a construct of the structural model has at least
#' 1 item). To prevent checking for errors use `.check_errors = FALSE`.
#' 
#' @usage parseModel(
#'   .model        = NULL, 
#'   .instruments  = NULL, 
#'   .check_errors = TRUE,
#'   .full_output  = FALSE
#'   )
#'
#' @inheritParams csem_arguments
#' @param .full_output Logical. Should the full output be returned. Defaults to
#'   `FALSE`. Only required for internal purposes.
#' 
#' @inherit csem_model return
#'
#' @examples
#' # ===========================================================================
#' # Providing a model in lavaan syntax 
#' # ===========================================================================
#' model <- "
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
#' "
#'
#' m <- parseModel(model)
#' m
#' 
#' # ===========================================================================
#' # Providing a complete model in cSEM format (class cSEMModel)
#' # ===========================================================================
#' # If the model is already a cSEMModel object, the model is returned as is:
#'
#' identical(m, parseModel(m)) # TRUE
#' 
#' # ===========================================================================
#' # Providing a list 
#' # ===========================================================================
#' # It is possible to provide a list that contains at least the
#' # elements "structural" and "measurement". This is generally discouraged
#' # as this may cause unexpected errors.
#' 
#' m_incomplete <- m[c("structural", "measurement", "construct_type")]
#' parseModel(m_incomplete)
#' 
#' # Providing a list containing list names that are not part of a `cSEMModel`
#' # causes an error:
#' 
#' \dontrun{
#' m_incomplete[c("name_a", "name_b")] <- c("hello world", "hello universe")
#' parseModel(m_incomplete)
#' }
#' 
#' # Failing to provide "structural" or "measurement" also causes an error:
#' 
#' \dontrun{
#' m_incomplete <- m[c("structural", "construct_type")]
#' parseModel(m_incomplete)
#' }
#' 
#' @export
#'
parseModel <- function(
  .model        = NULL, 
  .instruments  = NULL, 
  .check_errors = TRUE,
  .full_output  = FALSE
  ) {

  ### Check if already a cSEMModel list  
  if(class(.model) == "cSEMModel") {

    return(.model)
    
  ### Check if list
  } else if(is.list(.model)) {
    ## Check if list contains minimum necessary elements (structural, measurement)
    if(all(c("structural", "measurement") %in% names(.model))) {
      
      x <- setdiff(names(.model), c("structural", "measurement", "error_cor", 
                                    "construct_type", "construct_order", 
                                    "model_type", "instruments"))
      if(length(x) == 0) {
        
        class(.model) <- "cSEMModel"
        return(.model)
        
      } else {
        
        stop2("The following error occured in the `parseModel()` function:\n", 
              "The list provided contains element names unknown to cSEM: ", 
              paste0("'", x, "'", collapse = ", "), ".\n", 
              "See ?cSEMModel for a list of valid component names.")
      }
    } else {
      stop2("The following error occured in the `parseModel()` function:\n",
            "Structural and measurement matrix required.")
    }
  } else {
    
    ### Convert to lavaan partable ---------------------------------------------
    m_lav <- lavaan::lavaanify(model = .model, fixed.x = FALSE)
    
    ### Extract relevant information -------------------------------------------
    # s := structural
    # m := measurement
    # e := error
    tbl_s  <- m_lav[m_lav$op == "~", ] # structural 
    tbl_m  <- m_lav[m_lav$op %in% c("=~", "<~"), ] # measurement 
    tbl_e  <- m_lav[m_lav$op == "~~" & m_lav$user == 1, ] # error 
    
    ## Check if there are population starting values
    pop_values <- c(tbl_s$ustart, tbl_m$ustart, tbl_e$ustart)
    
    ## Get all relevant subsets of constructs and/or indicators
    # i  := indicators
    # c  := constructs
    # s  := structural model
    # m  := measurement/composite model
    # l  := only linear (terms)
    # [no]_nl := [no] only nonlinear (terms)
    # 
    # Typical name: "name_i/c_s/m_[lhs/rhs]_nl/l"
    
    ### Structural model ---------------------
    # Construct names of the structural model (including nonlinear terms)
    names_c_s_lhs    <- unique(tbl_s$lhs)
    names_c_s_rhs    <- unique(tbl_s$rhs)
    names_c_s        <- union(names_c_s_lhs, names_c_s_rhs)
    
    # Construct names of the structural model including the names of the 
    # individual components of the interaction terms
    names_c_s_lhs_l  <- unique(unlist(strsplit(names_c_s_lhs, "\\.")))
    names_c_s_rhs_l  <- unique(unlist(strsplit(names_c_s_rhs, "\\.")))
    names_c_s_l      <- union(names_c_s_lhs_l, names_c_s_rhs_l)
    
    # Nonlinear construct names of the the structural model 
    names_c_s_lhs_nl <- names_c_s_lhs[grep("\\.", names_c_s_lhs)] # must be empty
    names_c_s_rhs_nl <- names_c_s_rhs[grep("\\.", names_c_s_rhs)]
    
    # Construct names of the structural model without nonlinear terms
    names_c_s_no_nl  <- setdiff(names_c_s, names_c_s_rhs_nl)
    
    ### Measurement/composite model -----------------
    # All indicator names (observables AND constructs that serve as indicators for a 
    # 2nd order construct (linear and nonlinear))
    names_i          <- unique(tbl_m$rhs)
    
    # Indicator names that contain a "." (should only contain 
    # nonlinear 2ndorder terms!)
    names_i_nl       <- names_i[grep("\\.", names_i)] # this catches all terms 
                                                      # with a "."!
    
    # Construct names of the measurement model (including nonlinear terms;
    # Constructs of the rhs of the measurment model are considered first order constructs
    # attached to a second order construct)
    names_c_m_lhs    <- unique(tbl_m$lhs)
    names_c_m_rhs    <- intersect(names_i, c(names_c_m_lhs, names_i_nl))
    names_c_m        <- union(names_c_m_lhs, names_c_m_rhs)
    
    # Construct names of the measurement model including the names of the 
    # individual components of the interaction terms
    names_c_m_lhs_l  <- unique(unlist(strsplit(names_c_m_lhs, "\\.")))
    names_c_m_rhs_l  <- unique(unlist(strsplit(names_c_m_rhs, "\\.")))
    names_c_m_l      <- union(names_c_m_lhs_l, names_c_m_rhs_l)
    
    # Nonlinear construct names of the measurement model 
    names_c_m_lhs_nl <- names_c_m_lhs[grep("\\.", names_c_m_lhs)] # must be empty
    names_c_m_rhs_nl <- names_c_m_rhs[grep("\\.", names_c_m_rhs)]
    
    # Construct names of the measurement model without nonlinear terms
    names_c_m_no_nl  <- setdiff(names_c_s, names_c_s_rhs_nl)
    
    ## Hierachical constructs -------------------------
    # 1st order constructs attached to a second order construct
    names_c_attached_to_2nd <- names_c_m_rhs
    
    # 2nd order construct names
    names_c_2nd      <- unique(tbl_m[tbl_m$rhs %in% names_c_attached_to_2nd, "lhs"])
    
    # 1st order constructs not attached to a second order construct
    names_c_not_attachted_to_2nd <- setdiff(c(names_c_s, names_c_m), c(names_c_attached_to_2nd, names_c_2nd))
    
    # Higher order construct names
    names_c_higher   <- intersect(names_c_2nd, names_c_m_rhs) # must be empty
    
    ## Summary --------------------
    # All construct names (including nonlinear terms)
    names_c_all      <- union(names_c_s, names_c_m) 
    
    # All linear construct names
    names_c          <- names_c_all[!grepl("\\.", names_c_all)] 
    
    # All nonlinear construct names
    names_c_nl       <- names_c_all[grepl("\\.", names_c_all)]

    ## The the number of...
    number_of_constructs_all  <- length(names_c_all)
    number_of_constructs      <- length(names_c)
    number_of_indicators      <- length(names_i)
    
    ## Order
    construct_order <- rep("First order", length(names_c))
    names(construct_order) <- names_c
    construct_order[names_c_2nd] <- "Second order"
    
    ## Instruments
    if(!is.null(.instruments)) {
      names_construct_instruments <- names(.instruments)
      
      ## Note (05/2019): Currently, we only allow linear instruments (i.e. no 
      ##                 interaction terms). We may change that in the future.
      names_instruments_nl <- unlist(.instruments)[grep("\\.", unlist(.instruments))]
      
      if(length(names_instruments_nl) != 0) {
        stop2("The following error occured in the `parseModel()` function:\n",
              "Only linear instruments allowed. Dont know how to handle: ", 
              paste0("`", names_instruments_nl, "`", collapse = ", "))
      }
      names_instruments <- unique(unlist(strsplit(unlist(.instruments),  "\\.")))
    }
    
    ## Sometimes (e.g. for testMGD) we need to parse an incomplete model, hence
    ## errors and warnings should be ignored
    if(.check_errors) {
      ### Checks, errors and warnings --------------------------------------------
      ## Stop if only a subset of starting values is given
      if(!all(is.na(pop_values)) & anyNA(pop_values)) {
        stop2("The following error occured in the `parseModel()` function:\n",
              "Only a subset of population values given. Please specify",
              " all population values or none.")
      }
      if(!is.null(.instruments)) {
        # Note (05/2019): Currently, we only allow instruments from within the model, 
        #                 i.e instruments need to have a structural   
        #                 and a measurement/composite equation. 
        #                 Reason: we need to figure out, how to deal with 
        #                 instruments that have no structural equation. If they are
        #                 not part of the structural model its unclear how to 
        #                 get composites/scores for them.
        
        ## Check if all instruments are part of the structural model (i.e. internal)
        tmp <- setdiff(names_instruments, names_c)
        
        if(!is.null(.instruments) && length(tmp) != 0) {
          stop2("The following error occured in the `parseModel()` function:\n",
                "Currently, only internal instruments allowed. External instruments: ", 
                paste0("`", tmp,  "`", collapse = ", "))
        }
        
        ## Check if construct names for instruments are correct
        tmp <- setdiff(names_construct_instruments, names_c)
        
        if(!is.null(.instruments) && length(tmp) != 0) {
          stop2("The following error occured in the `parseModel()` function:\n",
                "Instruments supplied for unknown constructs: ", 
                paste0("`", tmp,  "`", collapse = ", "))
        } 
      }
      
      ## Stop if one indicator is connected to several constructs
      if(any(duplicated(tbl_m$rhs))) {
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "At least one indicator is connected to several constructs.")
      }
      
      ## Stop if any interaction/nonlinear term is used as an endogenous (lhs) variable in the
      ## structural model 
      if(length(names_c_s_lhs_nl)) {
        
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "Interaction terms cannot appear on the left-hand side of a structural equation.")
      }
      
      ## Stop if any interaction/nonlinear term is used as an endogenous (lhs) variable in the
      ## measurement model 
      if(length(names_c_m_lhs_nl)) {
        
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "Interaction terms cannot appear on the left-hand side of a measurement equation.")
      }
      
      ## Stop if any construct has no observables/indicators attached
      tmp <- setdiff(c(names_c_s_l, names_c_m_rhs_l), names_c_m_lhs)
      
      # Note: code below not required as long as only internal instruments 
      #       are allowed 
      ## Check if any of the individual components of the instruments has no 
      ## observables/indicators attached
      # if(!is.null(.instruments)) {
      #   tmp <- c(tmp, setdiff(names_instruments, names_c_m_lhs))
      # }
      
      if(length(tmp) != 0) {
        
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "No measurement equation provided for: ", 
          paste0("`", tmp,  "`", collapse = ", ")
        )
      } 
      
      ## Stop if a construct appears in the measurement but not in the 
      ## structural model
      tmp <- setdiff(names_c_m_lhs, c(names_c_s_l, names_c_m_rhs_l))
      
      # Note: code below not required as long as only internal instruments 
      #       are allowed 
      # If tmp is non-empty: check if the constructs are instruments
      # (only if not an error should be returned)
      # if(length(tmp) != 0 & !is.null(.instruments)) {
      #   tmp <- setdiff(tmp, names_instruments)
      # }
      
      if(length(tmp) != 0) {
        
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "The following constructs of the measurement model do not appear",
          " in the structural model: ", paste0("`", tmp, "`", collapse = ", ")
        )
      }
      
      ## Stop if any construct has a higher order than 2 (currently not allowed)
      if(length(names_c_higher) != 0) {
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          paste0("`", names_c_higher, "`"), " has order > 2.", 
          " Currently, only first and second-order constructs are supported.")
      }
      
      ## Stop if a nonlinear term is used as a first order construct to build/measure
      ## a second order construct
      if(length(names_i_nl) != 0) {
        stop2("The following error occured in the `parseModel()` function:\n",
              "Only linear first order constructs may be attached to second order constructs.",
              " Dont know how to handle: ", 
              paste0("`", names_i_nl, "`", collapse = ", "))
      }
      ## Stop if at least one of the components of an interaction term does not appear
      ## in any of the structural equations.
      tmp <- setdiff(names_c_s_l, names_c_s_no_nl)
      if(length(tmp) != 0) {
        
        stop2(
          "The following error occured in the `parseModel()` function:\n",
          "The nonlinear terms containing ", paste0("`", tmp, "`", collapse = ", "), 
          " are not embeded in a nomological net.")
      }
    }
    
    ## Construct type
    tbl_m$op <- ifelse(tbl_m$op == "=~", "Common factor", "Composite")
    construct_type  <- unique(tbl_m[, c("lhs", "op")])$op
    names(construct_type) <- unique(tbl_m[, c("lhs", "op")])$lhs
    construct_type <- construct_type[names_c]
    
    ## Type of model (linear or nonlinear)
    
    type_of_model <- if(length(names_c_nl) != 0) {
      "Nonlinear"
    } else {
      "Linear"
    }
    ### Construct matrices specifying the relationship between constructs,
    ### indicators and errors --------------------------------------------------
    # Note: code below not required as long as only internal instruments 
    #       are allowed 
    # model_structural  <- matrix(0,
    #                             nrow = length(names_c_s_l),
    #                             ncol = length(names_c_s),
    #                             dimnames = list(names_c_s_l, names_c_s)
    # )
    model_structural  <- matrix(0,
                                nrow = number_of_constructs,
                                ncol = number_of_constructs_all,
                                dimnames = list(names_c, names_c_all)
    )
    
    model_measurement <- matrix(0,
                                nrow = number_of_constructs,
                                ncol = number_of_indicators,
                                dimnames = list(names_c, names_i)
    )
    
    model_error       <- matrix(0,
                                nrow = number_of_indicators,
                                ncol = number_of_indicators,
                                dimnames = list(names_i, names_i)
    )
    
    ## Structural model
    row_index <- match(tbl_s$lhs, names_c)
    col_index <- match(tbl_s$rhs, names_c_all)
    # Note: code below not required as long as only internal instruments 
    #       are allowed 
    # row_index <- match(tbl_s$lhs, names_c_s_l)
    # col_index <- match(tbl_s$rhs, names_c_s)
    
    model_structural[cbind(row_index, col_index)] <- 1
    
    ## If starting values are given create a supplementary strucutral matrix
    ## that contains the starting values, otherwise assign a 1
    if(!anyNA(pop_values)) {
      model_structural2 <- model_structural
      model_structural2[cbind(row_index, col_index)] <- tbl_s$ustart
    }
    
    ## Measurement model
    row_index <- match(tbl_m$lhs, names_c)
    col_index <- match(tbl_m$rhs, names_i)
    
    model_measurement[cbind(row_index, col_index)] <- 1

    ## If starting values are given create a supplementary strucutral matrix
    ## that contains the starting values, otherwise assign a 1
    if(!anyNA(pop_values)) {
      model_measurement2 <- model_measurement
      model_measurement2[cbind(row_index, col_index)] <- tbl_m$ustart
    }
    
    ## Error covariance matrix
    m_errors   <- tbl_e[tbl_e$lhs %in% names_i, , drop = FALSE]
    con_errors <- tbl_e[tbl_e$lhs %in% names_c, , drop = FALSE] 
    
    row_index <- match(m_errors$lhs, names_i)
    col_index <- match(m_errors$rhs, names_i)
    
    model_error[cbind(c(row_index, col_index), c(col_index, row_index))] <- 1
    
    ## If starting values are given create a supplementary strucutral matrix
    ## that contains the starting values, otherwise assign a 1
    if(!anyNA(pop_values)) {
      
      model_error2 <- model_error
      
      # Extract endogenous and exogenous variables
      vars_endo <- rownames(model_structural)[rowSums(model_structural) != 0]
      vars_exo  <- setdiff(colnames(model_structural), vars_endo)
      
      ## Phi
      Phi <- matrix(0,
                    nrow = length(vars_exo),
                    ncol = length(vars_exo),
                    dimnames = list(vars_exo, vars_exo)
      )
      # Set diagonal elements to 1
      diag(Phi) <- 1
      
      if(length(tbl_e$ustart) != 0) {
        
        model_error2[cbind(c(row_index, col_index), c(col_index, row_index))] <- m_errors$ustart
        
        # Get row and column index for constructs
        row_index <- match(con_errors$lhs, vars_exo)
        col_index <- match(con_errors$rhs, vars_exo)
        
        Phi[cbind(row_index, col_index)] <- con_errors$ustart
        Phi[lower.tri(Phi)] <- t(Phi)[lower.tri(Phi)]
      }
    }
    
    # Currently, composite-based approaches (except GSCA ?) are unable to deal 
    # with measurement errors accros blocks (even if they were, it is not implemented in cSEM).
    # Which errors are across block
    model_error_temp <- model_error
    for(i in names_c) {
      x              <- which(model_measurement[i, ] == 1)
      model_error_temp[x, x] <- NA 
    }
    
    contains_error <- sum(model_error_temp, na.rm = TRUE)

    if(contains_error > 0) {
      warning2("The following warning occured in the `parseModel()` function:\n",
               "Measurement errors across blocks not supported (yet).",
               " Specified error correlation is ignored.")
    }

    
    ### Order model ============================================================
    # Order the structual equations in a way that every equation depends on
    # exogenous variables and variables that have been explained in a previous equation
    # This is necessary for the estimation of models containing nonlinear structual
    # relationships.
    
    ### Preparation ------------------------------------------------------------
    temp <- model_structural
    
    ## Extract endogenous and exogenous variables
    vars_endo <- rownames(temp)[rowSums(temp) != 0]
    vars_exo  <- setdiff(colnames(temp), vars_endo)
    
    # Endo variables that are explained by exo and endo variables
    explained_by_exo_endo <- vars_endo[rowSums(temp[vars_endo, vars_endo, drop = FALSE]) != 0]
    
    # Endo variables explained by exo variables only
    explained_by_exo <- setdiff(vars_endo, explained_by_exo_endo)
    
    ### Order =======================
    # First the endo variables that are soley explained by the exo variables
    model_ordered <- temp[explained_by_exo, , drop = FALSE]
    
    # Add variables that have already been ordered/taken care of to a vector
    # (including exogenous variables and interaction terms)
    already_ordered <- c(vars_exo, explained_by_exo)
    
    ## When there are feedback loops ordering does not work anymore, therefore
    #  ordering is skiped if there are feedback loops. Except for the
    #  the "replace" approach, this should not be a problem.
    
    # Does the structural model contain feedback loops
    if(any(temp[vars_endo, vars_endo] + t(temp[vars_endo, vars_endo]) == 2)) {
      
      model_ordered <- temp[c(already_ordered, explained_by_exo_endo), , drop = FALSE]
      
    } else {
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
            stop2(
              "The following error occured in the `parseModel()` function:\n",
              "Reordering the structural equations was not succesful."
            )
          if(length(explained_by_exo_endo_temp) == 0) break
        } # END repeat
      } # END if-statement
    } # END else
    
    ## Return a cSEMModel object.
    # A cSEMModel objects contains all the information about the model and its
    # components such as the type of construct used. 
    n  <- c(setdiff(names_c, rownames(model_ordered)), rownames(model_ordered))
    n1 <- intersect(n, colnames(model_ordered))
    m <- order(which(model_measurement[n, , drop = FALSE] == 1, arr.ind = TRUE)[, "row"])
    structural_ordered <- model_structural[n, c(n, setdiff(colnames(model_ordered), n)), drop = FALSE]
    # structural_ordered <- model_structural[n1, c(n1, setdiff(colnames(model_ordered), n1))]
    
    model_ls <- list(
      "structural"         = structural_ordered,
      "measurement"        = model_measurement[n, m, drop = FALSE],
      "error_cor"          = model_error[m, m, drop = FALSE],
      "construct_type"     = construct_type[match(n, names(construct_type))],
      "construct_order"    = construct_order[match(n, names(construct_order))],
      "model_type"         = type_of_model
      # "vars_endo"          = rownames(model_ordered),
      # "vars_exo"           = vars_exo,
      # "vars_explana"       = colnames(structural_ordered)[colSums(structural_ordered) != 0],
      # "explained_by_exo"   = explained_by_exo
    )
    
    ## Population values 
    if(!anyNA(pop_values)) {
      model_ls$structural2  <- model_structural2[rownames(structural_ordered), 
                                                 colnames(structural_ordered)]
      model_ls$measurement2 <- model_measurement2[n, m]
      model_ls$error_cor2   <- model_error2[m, m]
      model_ls$Phi          <- Phi
    }
    
    ## Are there instruments? If yes add them
    if(!is.null(.instruments)) {
      # Structural equations are named according to the name of the RHS variable
      # of each equation (i.e. the name of the dependent variables).
      
      # # Name of the equation
      # names_dependent <- names(which((rowSums(structural_ordered) != 0)))
      # 
      # # For which equation are instruments given
      # eq_with_instruments <- intersect(names_dependent, names(.instruments))
      
      for(i in names(.instruments)) {
        
        names_independent <- names(which(structural_ordered[i, ] == 1))
        ## Not sure how to deal with nonlinear terms. For now they are treated
        ## just like any other variable. If there is no instrument for a nonlinear
        ## term, it is treated as exogenous.
        
        names_exogenous   <- intersect(colnames(names_independent), .instruments[[i]])
        names_endogenous  <- setdiff(names_independent, .instruments[[i]])
        
        # First stage relations
        .instruments[[i]] <- matrix(1, nrow = length(names_endogenous), ncol = length(.instruments[[i]]),
                      dimnames = list(names_endogenous, .instruments[[i]]))
      }
      model_ls$instruments <- .instruments
    }
    
    ## Should the full output be returned
    if(.full_output) {
      model_ls$vars_2nd <- names_c_2nd
      model_ls$vars_attached_to_2nd <- names_c_attached_to_2nd
      model_ls$vars_not_attached_to_2nd <- names_c_not_attachted_to_2nd
    }
      
    class(model_ls) <- "cSEMModel"
    return(model_ls) 
  } # END else
}

#' Internal: Convert second order cSEMModel
#'
#' Uses a [cSEMModel] containg second order constructs and turns it into an
#' estimable model using either the "2stage" approach or the "mixed" approach.
#'
#' @usage convertModel(
#'  .csem_model        = NULL, 
#'  .approach_2ndorder = "2stage",
#'  .stage             = "first"
#'  )
#'
#' @inheritParams csem_arguments
#
#' @return A [cSEMModel] list that may be passed to any function requiring
#'   `.csem_model` as a mandatory argument.
#'
#' @keywords internal
#'
convertModel <- function(
  .csem_model        = NULL, 
  .approach_2ndorder = "2stage",
  .stage             = "first"
) {
  
  ## Note: currently we dont include nonlinear relationships in the first
  ##       stage of 2/3stage approach. If we want to change this, we
  ##       need to update the code that subsets the relevant constructs/indicators
  ##       below.
  
  # All linear constructs of the original model
  c_linear_original <- rownames(.csem_model$structural)
  # All constructs used in the first step (= all first order constructs)
  c_linear_1step <- names(.csem_model$construct_order[.csem_model$construct_order == "First order"])
  # All second order constructs
  c_2nd_order <- setdiff(c_linear_original, c_linear_1step)
  # All indicators of the original model (including linear and nonlinear 
  # constructs that form/measure a second order construct)
  i_original <- colnames(.csem_model$measurement)
  i_linear_original <- intersect(c_linear_original, i_original)
  i_nonlinear_original <- grep("\\.", i_original, value = TRUE) 
  # Linear constructs that serve as indicators and need to be replaced
  i_linear_1step <- setdiff(i_original, c(c_linear_original, i_nonlinear_original))
  
  if(.stage %in% c("second")) {
    # Linear constructs that dont form/measure a second order construct
    c_not_attached_to_2nd <- setdiff(c_linear_1step, i_linear_original)
    c_2step <- c(c_not_attached_to_2nd, c_2nd_order)
    
    ## Second step structural model
    x1 <- c()
    for(i in c_2step) {
      col_names <- colnames(.csem_model$structural[i, .csem_model$structural[i, , drop = FALSE] == 1, drop = FALSE])
      # col_names_linear <- intersect(c_linear_original, col_names)
      # col_names_nonlinear <- setdiff(col_names, col_names_linear)
      temp <- if(!is.null(col_names)) {
        ## Modify terms
        temp <- strsplit(x = col_names, split = "\\.")
        temp <- lapply(temp, function(x) {
          x[x %in% c_not_attached_to_2nd] <- paste0(x[x %in% c_not_attached_to_2nd], "_temp")
          paste0(x, collapse = ".")
        })
        # col_names_nonlinear <- unlist(temp)
        # col_names <- c(col_names_linear, col_names_nonlinear)
        col_names <- unlist(temp)
        # col_names[col_names %in% nc_not_to_2nd] <- paste0(col_names[col_names %in% nc_not_to_2nd], "_temp")
        paste0(ifelse(i %in% c_not_attached_to_2nd, paste0(i, "_temp"), i), "~", paste0(col_names, collapse = "+")) 
      } else {
        "\n"
      }
      x1 <- paste(x1, temp, sep = "\n")
    }
    
    ## Measurement model + second order structural equation 
    # Constructs that dont form/measure a second order construct
    x2a <- c()
    for(i in c_not_attached_to_2nd) {
      temp <- paste0(paste0(i, "_temp"), ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), i)
      x2a <- paste(x2a, temp, sep = "\n")
    }
    # Second order constructs
    x2b <- c()
    for(i in c_2nd_order) {
      col_names <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
      temp  <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(col_names, collapse = "+"))
      x2b <- paste(x2b, temp, sep = "\n")
    }
    
    ## Model to be parsed
    lav_model <- paste(x1, x2a, x2b, sep = "\n")
  } else { # BEGIN: first step
    
    if(.approach_2ndorder %in% c("mixed")) {
      
      ## Structural model
      # First order equations
      x1 <- c()
      for(i in c_linear_original) {
        col_names <- colnames(.csem_model$structural[i, .csem_model$structural[i, , drop = FALSE] == 1, drop = FALSE])
        temp <- if(!is.null(col_names)) {
          paste0(i, "~", paste0(col_names, collapse = "+")) 
        } else {
          "\n"
        }
        x1 <- paste(x1, temp, sep = "\n")
      }
      
      # Add indirect effect for the extended repeated indicators approach
      # if(.approach_2ndorder == "RI_extended") {
      #   # 1. Which constructs are attached to which 2nd order construct
      #   # 2. Which antecedent constructs does the second order have
      #   # 3. Add structural equations for all antecedent constructs of second order
      #   #    construct j on all first order constructs building/measuring the
      #   #    second order construct
      #   # 4. Add path between all constructs attached to second construct j
      #   #    (done at the end of the functions)
      #   
      #   for(j in c_2nd_order) {
      #     attached <-  names(.csem_model$measurement[j, .csem_model$measurement[j, ] == 1])
      #     antecedents <- names(which(.csem_model$structural[j, ] == 1))
      #     
      #     for(i in attached) {
      #       x1 <- paste(x1, paste0(i, "~", antecedents, collapse = "+"), sep = "\n")
      #     } 
      #   }
      # }
      
      ## Measurement model + second order structural equation 
      # First order constructs
      x2a <- c()
      for(i in c_linear_1step) {
        col_names <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
        temp  <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(col_names, collapse = "+"))
        x2a <- paste(x2a, temp, sep = "\n")
      }
      
      # Second order constructs
      x2b <- c()
      for(i in c_2nd_order) {
        # i <- c_2nd_order[1]
        col_names_1 <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE ] == 1, drop = FALSE])
        col_names_1_nonlinear <- grep("\\.", col_names_1, value = TRUE) 
        col_names_1_linear <- setdiff(col_names_1, col_names_1_nonlinear)
        col_names_2 <- .csem_model$measurement[col_names_1_linear, colSums(.csem_model$measurement[col_names_1_linear, ,drop = FALSE]) != 0, drop = FALSE]
        temp <- paste0(i, "_2nd_", colnames(col_names_2))
        temp <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(temp, collapse = "+"))
        x2b <- paste(x2b, temp, sep = "\n")
        ## add second order structural equation
        if(.csem_model$construct_type[i] == "Composite") {
          x2b <- paste(x2b, paste0(i, "~", paste0(col_names_1, collapse = "+" )), sep = "\n")
        } else {
          x2ba <- c()
          for(j in i_linear_original) {
            x2ba <- paste(x2ba, paste0(j, "~", i), sep = "\n")
          }
          x2b <- paste(x2b, x2ba, sep = "\n")
        }
      }
      
      ## Error_cor
      # First order
      x3 <- c()
      for(i in i_linear_1step) {
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
      for(i in 2:length(c_linear_1step)) {
        temp <- paste0(c_linear_1step[i], "~", paste0(c_linear_1step[1:(i-1)], collapse = "+"))
        x1   <- paste(x1, temp, sep = "\n")
      }
      ## Measurement model
      x2 <- c()
      for(i in c_linear_1step) {
        col_names <- colnames(.csem_model$measurement[i, .csem_model$measurement[i, , drop = FALSE] == 1, drop = FALSE])
        temp <- paste0(i, ifelse(.csem_model$construct_type[i] == "Composite", "<~", "=~"), paste0(col_names, collapse = "+"))
        x2   <- paste(x2, temp, sep = "\n")
      }
      ## Error_cor
      x3 <- c()
      for(i in i_linear_1step) {
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
    } # END first step of the 2/3 stage approach
  } # END first step

  model <- parseModel(lav_model)
  
  ## Add path between all constructs attached to second construct j
  ## if the extended repeated indicators approach was used
  # if(.approach_2ndorder == "RI_extended") {
  #   # Note there must not be feedback loops, so we make sure model$structural is
  #   # lower triangular
  #   
  #   m <- model$structural[.csem_model$vars_attached_to_2nd, .csem_model$vars_attached_to_2nd]
  #   m[lower.tri(m)] <- 1
  #   
  #   model$structural[.csem_model$vars_attached_to_2nd, .csem_model$vars_attached_to_2nd] <- m
  #   }
  
  ## add
  return(model)
}