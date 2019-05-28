#' Estimate the structural coefficients
#'
#' Estimates the coefficients of the structural model (nonlinear and linear) using
#' OLS or 2SLS.
#'
#' @usage estimatePath(
#'  .approach_nl    = args_default()$.approach_nl,
#'  .approach_paths = args_default()$.approach_paths,
#'  .csem_model     = args_default()$.csem_model,
#'  .H              = args_default()$.H,
#'  .instruments    = args_default()$.instruments,
#'  .normality      = args_default()$.normality,
#'  .P              = args_default()$.P,
#'  .Q              = args_default()$.Q
#'  )
#'   
#' @inheritParams csem_arguments
#'
#' @return A named list containing the estimated structural coefficients, the
#'   R2, the adjusted R2, and the VIF's for each regression.
#'

estimatePath <- function(
  .approach_nl    = args_default()$.approach_nl,
  .approach_paths = args_default()$.approach_paths,
  .csem_model     = args_default()$.csem_model,
  .H              = args_default()$.H,
  .instruments    = args_default()$.instruments,
  .normality      = args_default()$.normality,
  .P              = args_default()$.P,
  .Q              = args_default()$.Q
  ) {
  
  ## Check approach_path argument:
  if(!any(.approach_paths %in% c("OLS", "2SLS", "3SLS"))) {
    stop2("The following error occured in the `estimatePath()` function:\n",
          paste0("'", .approach_paths, "'"), 
          " is an unknown approach to estimate the path model.")
  }
  
  ## Warning if instruments are given but .approach_path = "OLS"
  if(!is.null(.instruments) & .approach_paths == "OLS") {
    warning2("Instruments supplied but path approach is 'OLS'.\n",
             "Instruments are ignored.", 
             " Consider setting `.approach_path = '2SLS'.")
  }
  
  ## Error if no instruments are given but .approach_path = "2SLS"
  if(is.null(.instruments) & (.approach_paths == "2SLS" | .approach_paths == "3SLS")) {
    stop2("`.approach_path = '2SLS' '3SLS' requires instruments.")
  }
  
  m         <- .csem_model$structural
  vars_endo <- rownames(m)[rowSums(m) != 0]
  vars_exo  <- setdiff(colnames(m), vars_endo)
  explained_by_exo_endo <- vars_endo[rowSums(m[vars_endo, vars_endo, drop = FALSE]) != 0]
  vars_ex_by_exo <- setdiff(vars_endo, explained_by_exo_endo)
  vars_explana   <- colnames(m)[colSums(m) != 0]

  # Number of observations (required for the adjusted R^2)
  n <- dim(.H)[1]
  
  if(.csem_model$model_type == "Linear") {

    res <- lapply(vars_endo, function(x) {
      # Which of the variables in vars_endo have instruments specified, i.e.
      # have endogenous variables on the RHS. By default: FALSE.
      endo_in_RHS <- FALSE
      
      if(!is.null(.instruments)) {
        endo_in_RHS <- x %in% names(.instruments)
      }
      
      ## Independent variables of the structural equation of construct x
      indep_var <-  colnames(m[x, m[x, ] != 0, drop = FALSE])
      
      # Compute "OLS" if endo_in_RHS is FALSE, i.e no instruments are 
      # given for this particular equation or .approach_path is "OLS"
      if(!endo_in_RHS | .approach_paths == "OLS") {
        
        # Coef = (X'X)^-1X'y = V(eta_indep)^-1 Cov(eta_indep, eta_dep)
        coef <- solve(.P[indep_var, indep_var, drop = FALSE]) %*% .P[indep_var, x, drop = FALSE]
        
        # Since Var(dep_Var) = 1 we have R2 = Var(X coef) = t(coef) %*% X'X %*% coef
        r2   <- c(t(coef) %*% .P[indep_var, indep_var, drop = FALSE] %*% coef)
        # names(r2) <- x
        
        # Calculation of the adjusted R^2
        r2adj <- c(1 - (1 - r2)*(n - 1)/(n - length(indep_var)-1))
        # names(r2adj) <- x
        
        # Calculation of the VIF values (VIF_k = 1 / (1 - R^2_k)) where R_k is
        # the R^2 from a regression of the k'th explanatory variable on all other
        # explanatory variables of the same structural equation.
        # VIF's require at least two explanatory variables to be meaningful
        vif <- if(length(indep_var) > 1) {
          diag(solve(cov2cor(.P[indep_var, indep_var, drop = FALSE])))
        } else {
          NA
        } 
      } # END OLS
      

      # Compute "2SLS" if endo_in_RHS is TRUE, i.e instruments are 
      # given for this particular equation and .approach_path is "2SLS".
      
      ## Two stage least squares (2SLS) and three stage least squares (3SLS)
      if(endo_in_RHS & (.approach_paths == "2SLS" | .approach_paths == "3SLS")) {
        
        ## First stage
        # Note: Regress the P endogenous variables (X) on the L instruments 
        #       and the K exogenous independent variables (which must be part of Z).
        #       Therefore: X (N x P) and Z (N x (L + K)) and
        #       beta_1st = (Z'Z)^-1*(Z'X)
        names_X <- rownames(.csem_model$instruments[[x]])
        names_Z <- colnames(.csem_model$instruments[[x]])
        
        ## Error if the number of instruments (including the K exogenous variables)
        ## is less than the number of independent variables in the original 
        ## structural equation for construct "x"
        if(length(names_Z) < length(indep_var)) {
          stop2("The following error occured in the `estimatePath()` function:\n",
                "The number of instruments for the structural equation of construct ",
                paste0("'", x, "'"), " is less than the number of independent ",
                "variables.\n", "Make sure all exogenous variables correctly ",
                " supplied as instruments to `.instruments`.")
        }
        
        # Assuming that .P (the construct correlation matrix) also contains 
        # the instruments (ensured if only internal instruments are allowed)
        # we can use .P.
        
        beta_1st <- solve(.P[names_Z, names_Z, drop = FALSE], 
                          .P[names_Z, indep_var, drop = FALSE])
        
        ## Second stage
        # Note: X_hat = beta_1st*Z --> X_hat'X_hat = beta_1st' (Z'Z) beta_1st
        
        coef <- solve(t(beta_1st) %*% .P[names_Z, names_Z, drop = FALSE] %*% beta_1st, 
                      t(beta_1st) %*% .P[names_Z, x, drop = FALSE])
        
        
        # Although the r^2 can be calculated in case of 2SLS,
        # the r^2 and all corresponding statistics are not correct. 
        # Hence, I suggest to overwrite it with NA. This might help to detect potential problems.
        # 
        r2    = NA
        r2adj = NA
        
        # The VIF should be based on the second-stage equation 
        vif <- if(length(names_Z) > 1) {
          diag(solve(cov2cor(.P[names_Z, names_Z, drop = FALSE])))
        } else {
          NA
        }
      } # END 2SLS
      
      ## Collect results
      list("coef" = coef, "r2" = r2, "r2adj" = r2adj, "vif" = vif)
    }) # END lapply
    
    names(res) <- vars_endo
    res <- purrr::transpose(res)

    if(.approach_paths == "3SLS"){
      # Variance covariance matrix of the error term
      VCVresid=matrix(0,nrow=length(vars_endo),ncol=length(vars_endo),
                      dimnames=list(vars_endo,vars_endo))
      
      # Fill the VCV of the error terms
      for(i in   vars_endo){
        for(j in  vars_endo){
          coefsi <- res$coef[[i]]
          coefsj <- res$coef[[j]] 
          VCVresid[i,j] <- .P[i,j] - t(coefsi) %*%  .P[m[i,]!=0,j,drop = FALSE] -
            .P[i, m[j,]!=0, drop = FALSE] %*% coefsj +
            t(coefsi) %*% .P[m[i,]!=0, m[j,]!=0, drop=FALSE] %*% coefsj

      part = lapply(vars_endo,function(x){
        independents = colnames(m)[m[x,]!=0]
        indendo = intersect(independents,vars_endo)
        indexog = intersect(independents,vars_exo)
        InvOfVCVresid=solve(VCVresid)

      LHSpart=sapply(vars_endo,function(mue){ 
          as.numeric(InvOfVCVresid[x,mue])* 
            .P[c(indendo,indexog),vars_exo , drop = FALSE]%*%
            solve(.P[vars_exo,vars_exo,drop=FALSE])%*%
            .P[vars_exo,mue,drop=FALSE]
        })
      
      # sum up all elements
      if(is.matrix(LHSpart)){
        LHSpart=matrix(rowSums(LHSpart),ncol=1)
      }else{ 
        LHSpart=sum(LHSpart)
      }
      
      RHSpart=lapply(vars_endo, function(mue){
        InvOfVCVresid[x,mue, drop=TRUE]* 
          .P[c(indendo,indexog),vars_exo,drop=FALSE]%*%
          solve(.P[vars_exo,vars_exo])%*%
          .P[vars_exo,c(intersect(colnames(m)[m[mue,]!=0],vars_endo),
                        intersect(colnames(m)[m[mue,]!=0],vars_exo))]
        })
      
      RHSpart = do.call(cbind,RHSpart)
        
      list(LHSpart = LHSpart,RHSpart = RHSpart)
      })
      
      
      part = purrr::transpose(part)
      LHS = do.call(rbind,part[["LHSpart"]])
      RHS = do.call(rbind,part[["RHSpart"]])
      
      # solve equation
      allparas=solve(RHS,LHS)
      
      # parameters need to be sorted back, i.e., overwrite the res object of 2SLS
      
      stop2("3SLS is not implemented yet.")
      
    }
    
    
  } else {
    ## Error if approach_paths is not "OLS"
    # Note (05/2019): Currently, only "OLS" is allowed for nonlinear models
    if(.approach_paths != "OLS") {
      stop2("The following error occured in the `estimatePath()` function:\n",
           "Currently, ", .approach_paths, " is only applicable to linear models.")
    }
    
    ### Preparation ============================================================
    # Implementation and notation is based on:
    # Dijkstra & Schermelleh-Engel (2014) - PLSc for nonlinear structural
    #                                       equation models
    
    ### Calculation ============================================================
    ## Calculate elements of the VCV matrix of the explanatory variables -------
    if(.normality == TRUE) {
      # For the sequential approach normality = TRUE requires all 
      # explanatory variables to be exogenous!
      if(length(setdiff(vars_explana, vars_exo)) != 0 & .approach_nl == "sequential") {
        
        stop("The following error was encountered while calculating the path coefficients:\n",
             "The sequential approach can only be used in conjunction with `normality = TRUE`", 
             " if all explanatory variables are exogenous.", call. = FALSE)
      } else {
        vcv_explana <- outer(vars_explana,
                             vars_explana,
                             FUN = Vectorize(f3, vectorize.args = c(".i", ".j")),
                             .Q  = .Q,
                             .H  = .H)
      }

      # It can happen that this matrix is not symmetric
      vcv_explana[lower.tri(vcv_explana)] = t(vcv_explana)[lower.tri(vcv_explana)]
      
    } else {
      
      # Define the type/class of the moments in the VCV matrix of the explanatory
      # variables 
      class_explana <- outer(vars_explana, vars_explana, FUN = Vectorize(f1))
      rownames(class_explana) <- colnames(class_explana) <- vars_explana
      
      # Calculate
      vcv_explana <- outer(vars_explana,
                           vars_explana,
                           FUN = Vectorize(f2, vectorize.args = c(".i", ".j")),
                           .select_from = class_explana,
                           .Q = .Q,
                           .H = .H)
    
      # It can happen that this matrix is not symmetric
      vcv_explana[lower.tri(vcv_explana)] = t(vcv_explana)[lower.tri(vcv_explana)]
      
      }
    
    # Set row- and colnames for matrix
    rownames(vcv_explana) <- colnames(vcv_explana) <- vars_explana
    
    # Create list with each list element holding the VCV matrix of the
    # explanatory variables of one endogenous variable
    vcv_explana_ls <- lapply(vars_endo, function(x) {
      res <- colnames(m[x, m[x, , drop = FALSE] == 1, drop = FALSE])
      vcv_explana[res, res, drop = FALSE]
    })
    names(vcv_explana_ls) <- vars_endo
    
    ## Check if all vcv matrices are semi positive-definite and warn if not
    semidef <- lapply(vcv_explana_ls, function(x) {
      matrixcalc::is.positive.semi.definite(x)
    })
    
    if(any(!unlist(semidef))) {
      warning("The following issue was encountered while calculating the path coefficients:\n",
              "The variance-covariance matrix of the explanatory variables for ",
              "at least one of the structural equations is not positive semi-definite.",
              call. = FALSE, immediate. = TRUE)
    }
    ## Calculate covariances between explanatory and endogenous variables ------
    
    # Define the class of the moments in the VCV matrix between explanatory
    # and endogenous variables
    class_endo_explana <- outer(vars_endo, vars_explana, FUN = Vectorize(f1))
    rownames(class_endo_explana) <- vars_endo
    colnames(class_endo_explana) <- vars_explana
    
    # Calculate
    cv_endo_explana <- outer(vars_endo, vars_explana,
                             FUN = Vectorize(f2, vectorize.args = c(".i", ".j")),
                             .select_from = class_endo_explana,
                             .Q = .Q,
                             .H = .H)
    rownames(cv_endo_explana) <- vars_endo
    colnames(cv_endo_explana) <- vars_explana
    
    # Create list with each list element holding the covariances between one
    # endogenous variable and its corresponding explanatory variables
    cv_endo_explana_ls <- lapply(vars_endo, function(x) {
      res <- colnames(m[x, m[x, , drop = FALSE] == 1, drop = FALSE])
      cv_endo_explana[x, res, drop = FALSE]
    })
    names(cv_endo_explana_ls) <- vars_endo
    
    ## Calculate path coef, R2 and VIF ----------------------------------------------
    # Path coefficients
    coef <- mapply(function(x, y) solve(x) %*% t(y),
                   x = vcv_explana_ls,
                   y = cv_endo_explana_ls,
                   SIMPLIFY = FALSE)
    
    # Coefficient of determinaten (R^2)
    r2 <- mapply(function(x, y) t(y) %*% x %*% y,
                 x = vcv_explana_ls,
                 y = coef,
                 SIMPLIFY = FALSE)
    
    # Adjusted R^2 
    r2adj = mapply(function(x,y) 1-(1-x)*(n-1)/(n-nrow(y)),
                   x = r2,
                   y = coef)
    
    # Variance inflation factor
    vif = lapply(vcv_explana_ls, function(x) diag(solve(cov2cor(x))))
    
    ##==========================================================================
    # Replacement approach
    ### ========================================================================
    if(.approach_nl == "replace") {
      # warning("Something is wrong here!")
      ### Preparation ==========================================================
      if(.normality == FALSE) {
        
        stop("The following error was encountered while calculating the path coefficients:\n",
             "The replacement approach is only implemented for `normality = TRUE`.",
             call. = FALSE)
      }
      # Create list with each list element holding one structural equation
      struc_coef_ls <- lapply(coef, function(x) {
        a <- c(x)
        names(a) <- rownames(x)
        a
      })
      
      # Add a "structural equation" for all exogenous constructs
      temp <- intersect(rownames(.csem_model$structural), vars_exo)
      
      if(length(temp) > 0 ) {
        
        struc_coef_ls <- lapply(temp, function(x) {
          struc_coef_ls[[x]] <- 1
          names(struc_coef_ls[[x]]) <- x
          struc_coef_ls
        })[[1]] # there is a problem here 
      }
      
      ### Calculation ==========================================================
      ## Calculate variance of the structural errors
      var_struc_error <- 1 - unlist(r2)
      
      ## Preallocate
      vcv  <- list()

      ## Loop over each endogenous variable
      for(k in vars_endo) {
        
        if(k %in% vars_ex_by_exo) {
          # If the endogenous variable is only explained by exogenous variables:
          # add an error term (zeta)
          
          struc_coef_ls[[k]][paste0("zeta_", k)] <- 1
          
        } else {
          # If the endogenous variable is explained by at least one other
          # endogenous variable the covariances between all explanatory variables
          # needs to be computed in order to compute path coefficients later on
          
          ## Preallocate
          temp <- list()
          explana_k <- names(struc_coef_ls[[k]])
          
          ## Loop over each explanatory variable of structural equation k
          for(m in explana_k) {
            
            # Split term
            a <- strsplit(m, "\\.")[[1]]
            
            # Insert corresponding equation for the first componenent of a
            temp[[m]] <- struc_coef_ls[[a[1]]]
            
            if(length(a) > 1) {
              
              ## Insert the (previously build) corresponding equation for each
              ## component of the splitted term
              for(l in 1:(length(a) - 1)) {
                
                rr             <- temp[[m]] %o% struc_coef_ls[[a[l + 1]]]
                rr_vec         <- c(rr)
                names(rr_vec)  <- c(outer(rownames(rr),
                                          colnames(rr),
                                          FUN = paste, sep = "."))
                
                temp[[m]] <- rr_vec
              } # END for l in 1:(length(a) - 1)
            } # END if
          } # END for m in explana_k
          
          ## Calculate vcv matrix of the explana variables ---------------------
          vcv[[k]] <- outer(explana_k, explana_k,
                            FUN = Vectorize(f4, vectorize.args = c(".i", ".j")),
                            .Q  = .Q,
                            .H  = .H,
                            .var_struc_error = var_struc_error,
                            .temp = temp)
          
          # Set row- and colnames for vcv matrix
          rownames(vcv[[k]]) <- colnames(vcv[[k]]) <- explana_k
          
          ## Calculate path coefs, R^2, adjusted R^2, VIF and update "struc_coef_ls" (= matrix of
          ## structural equations) and "var_struc_error" (= vector of
          ## structural error variances) ---------------------------------------
          
          coef[[k]] <- solve(vcv[[k]]) %*% t(cv_endo_explana_ls[[k]])
          r2[[k]]   <- t(coef[[k]]) %*% vcv[[k]] %*% coef[[k]]
          r2adj[[k]] = 1-(1-r2[[k]])*(n-1)/(n-nrow(coef[[k]]))
          vif[[k]] = diag(solve(cov2cor(vcv[[k]])))
          var_struc_error[k]    <- 1 - r2[[k]]
          
          temp <- mapply(function(x, y) x * y,
                         x = temp,
                         y = coef[[k]],
                         SIMPLIFY = FALSE)
          
          struc_coef_ls[[k]]        <- unlist(temp)
          names(struc_coef_ls[[k]]) <- unlist(lapply(temp, names), use.names = FALSE)
          struc_coef_ls[[k]][paste0("zeta_", k)] <- 1
          
        } # END else
      } # END for k in vars_endo
    } # END if(.approach_nlhod = replace)
    res <- list("coef" = coef, "r2" = r2, "r2adj" = r2adj, "vif" = vif)
  } # END if nonlinear
  ### Structure results --------------------------------------------------------
  tm <- t(.csem_model$structural)
  tm[which(tm == 1)] <- do.call(rbind, res$coef)
  
  ## Delete VIF's that are set to NA
  res$vif <- Filter(Negate(anyNA), res$vif)
  
  ## Return result -------------------------------------------------------------
  list("Path_estimates" = t(tm), "R2" = unlist(res$r2),"R2adj" = unlist(res$r2adj), "VIF" = res$vif)
}