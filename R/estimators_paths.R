#' Estimate the structural coefficients
#'
#' Estimates the coefficients of the structural model (non-linear and linear)
#'
#' More details here
#'
#' @inheritParams csem_arguments
#'
#' @return A named list containing the estimated structural coefficients and the
#'   r^2 for each regression.
#'
#' @export
#'

estimatePathOLS <- function(
  .H                   = NULL,
  .W                   = NULL,
  .Q                   = NULL,
  .csem_model          = NULL,
  .normal              = NULL,
  .approach_nl         = NULL
) {
  ### Preparation ==============================================================
  # Implementation and notation is based on:
  # Dijkstra & Schermelleh-Engel (2014) - PLSc for nonlinear structural
  #                                       equation models

  m              <- .csem_model$structural
  # m              <- .csem_model$structural_ordered
  vars_explana   <- .csem_model$vars_explana
  vars_endo      <- .csem_model$vars_endo
  vars_exo       <- .csem_model$vars_exo
  vars_ex_by_exo <- .csem_model$explained_by_exo

  # Evaluate Proxies H and Q's
  .H <- .H
  .Q <- .Q

  ### Calculation ==============================================================
  ## Calculate elements of the VCV matrix of the explanatory variables ---------
  if(.normal == TRUE) {

    vcv_explana <- outer(vars_explana,
                         vars_explana,
                         FUN = Vectorize(f3, vectorize.args = c(".i", ".j")),
                         .Q  = .Q,
                         .H  = .H)
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

  ## Calculate covariances between explanatory and endogenous variables --------

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

  ## Calculate path coef and R2 ------------------------------------------------
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

  ##============================================================================
  # Replacement approach
  ### ==========================================================================
  if(.approach_nl == "replace") {
    ### Preparation ============================================================
    if(.normal == FALSE) {
      stop("Only implemented for normal = TRUE")
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
        })[[1]]
    }

    ### Calculation ============================================================
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

        ## Calculate vcv matrix of the explana variables -----------------------
        vcv[[k]] <- outer(explana_k, explana_k,
                          FUN = Vectorize(f4, vectorize.args = c(".i", ".j")),
                          .Q  = .Q,
                          .H  = .H,
                          .var_struc_error = var_struc_error,
                          .temp = temp)

        # Set row- and colnames for vcv matrix
        rownames(vcv[[k]]) <- colnames(vcv[[k]]) <- explana_k

        ## Calculate path coefs, R^2 and update "struc_coef_ls" (= matrix of
        ## structural equations) and "var_struc_error" (= vector of
        ## structural error variances) -----------------------------------------

        coef[[k]] <- solve(vcv[[k]]) %*% t(cv_endo_explana_ls[[k]])
        r2[[k]]   <- t(coef[[k]]) %*% vcv[[k]] %*% coef[[k]]
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

  ## Structure results
  names_path <- paste0(rep(names(coef), times = sapply(coef, length)), " ~ ",
                       unlist(sapply(coef, rownames)))

  path_estimates  <- data.frame("Path" = names_path,
                                "Estimate" = unlist(coef, use.names = FALSE))

  ## Return result -------------------------------------------------------------
  list("Path_estimates" = path_estimates, R2 = r2)
}

# estimatePath2SLS <- function(
#   .H                   = NULL,
#   .W                   = NULL,
#   .Q                   = NULL,
#   .csem_model          = NULL,
# ) {
# 
# }