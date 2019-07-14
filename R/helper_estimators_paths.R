#' Internal: Utility functions for the estimation of non-linear models
#' 
#' 
#' @param .i Row index
#' @param .j Column index
#' @param .select_from matrix to select from
#' @inheritParams csem_arguments
#' 
#' @name nonlinear_estimation_utilities
#' @rdname nonlinear_estimation_utilities
#' @keywords internal
f1 <- function(.i, .j) {

  tab_i <- classifyConstructs(.i)[[1]]
  tab_j <- classifyConstructs(.j)[[1]]

  class_ij <- paste(unique(tab_i$Term_class),
                    unique(tab_j$Term_class), sep = "_")

  ## Rename/reorder joint class so that "class1_class2" = "class2_class1"
  switch (class_ij,
          # Class "Single_x"
          "Quadratic_Single"      = {class_ij <- "Single_Quadratic"},
          "Cubic_Single"          = {class_ij <- "Single_Cubic"},
          "TwInter_Single"        = {class_ij <- "Single_TwInter"},
          "ThrwInter_Single"      = {class_ij <- "Single_ThrwInter"},
          "QuadTwInter_Single"    = {class_ij <- "Single_QuadTwInter"},
          # Class "Quadratic_x"
          "Cubic_Quadratic"       = {class_ij <- "Quadratic_Cubic"},
          "TwInter_Quadratic"     = {class_ij <- "Quadratic_TwInter"},
          "ThrwInter_Quadrtic"    = {class_ij <- "Quadratic_ThrwInter"},
          "QuadTwInter_Quadratic" = {class_ij <- "Quadratic_QuadTwInter"},
          # Class "Cubic_x"
          "TwInter_Cubic"         = {class_ij <- "Cubic_TwInter"},
          "ThrwInter_Cubic"       = {class_ij <- "Cubic_ThrwInter"},
          "QuadTwInter_Cubic"     = {class_ij <- "Cubic_QuadTwInter"},
          # Class "TwInter_x"
          "ThrwInter_TwInter"     = {class_ij <- "TwInter_ThrwInter"},
          "QuadTwInter_TwInter"   = {class_ij <- "TwInter_QuadTwInter"},
          # Class "ThrwInter_x"
          "QuadTwInter_ThrwInter" = {class_ij <- "ThrwInter_QuadTwInter"},
          # Dont change the object if none of the above apply
          # Note: If none of the alternatives apply,
          # switch evaluates the last alternative...not sure anymore if thats
          # true. Needs to be checked!
          "None_of_the_above"     = {class_ij}
  )
  ## Return
  class_ij
}
#' @rdname nonlinear_estimation_utilities
f2 <- function(.i, .j, .select_from, .Q, .H) {

  x <- .select_from[.i, .j, drop = FALSE]

  ## Select estimator according to joint class of the ij'th element of
  ## the VCV matrix of the explanatory variables

  args <- list(.i = .i, .j = .j, .Q = .Q, .H = .H)

  switch (x,
          # Class "Single_x"
          "Single_Single"           = {res <- do.call(SingleSingle, args)},
          "Single_Quadratic"        = {res <- do.call(SingleQuadratic, args)},
          "Single_Cubic"            = {res <- do.call(SingleCubic, args)},
          "Single_TwInter"          = {res <- do.call(SingleTwInter, args)},
          "Single_ThrwInter"        = {res <- do.call(SingleThrwInter, args)},
          "Single_QuadTwInter"      = {res <- do.call(SingleQuadTwInter, args)},
          # Class "Quadratic_x"
          "Quadratic_Quadratic"     = {res <- do.call(QuadraticQuadratic, args)},
          "Quadratic_Cubic"         = {res <- do.call(QuadraticCubic, args)},
          "Quadratic_TwInter"       = {res <- do.call(QuadraticTwInter, args)},
          "Quadratic_ThrwInter"     = {res <- do.call(QuadraticThrwInter, args)},
          "Quadratic_QuadTwInter"   = {res <- do.call(QuadraticQuadTwInter, args)},
          # Class "Cubic_x"
          "Cubic_Cubic"             = {res <- do.call(CubicCubic, args)},
          "Cubic_TwInter"           = {res <- do.call(CubicTwInter, args)},
          "Cubic_ThrwInter"         = {res <- do.call(CubicThrwInter, args)},
          "Cubic_QuadTwInter"       = {res <- do.call(CubicQuadTwInter, args)},
          # Class "TwInter_x"
          "TwInter_TwInter"         = {res <- do.call(TwInterTwInter, args)},
          "TwInter_ThrwInter"       = {res <- do.call(TwInterThrwInter, args)},
          "TwInter_QuadTwInter"     = {res <- do.call(TwInterQuadTwInter, args)},
          # Class "ThrwInter_x"
          "ThrwInter_ThrwInter"     = {res <- do.call(ThrwInterThrwInter, args)},
          "ThrwInter_QuadTwInter"   = {res <- do.call(ThrwInterQuadTwInter, args)},
          # Class "QuadTwInter_x"
          "QuadTwInter_QuadTwInter" = {res <- do.call(QuadTwInercQuadTwInter, args)}
  )
  ## Retun result
  res
}
#' @rdname nonlinear_estimation_utilities
f3 <- function(.i, .j, .Q, .H) {

  Q = .Q
  H = .H

  tab_i <- classifyConstructs(.i)[[1]]
  tab_j <- classifyConstructs(.j)[[1]]

  x1 <- merge(tab_i, tab_j, "Component", all = TRUE)
  x2 <- rowSums(x1[, c("Component_freq.x", "Component_freq.y")], na.rm = TRUE)

  component <- list(x1$Component, tab_i$Component, tab_j$Component)
  freq      <- list(x2, tab_i$Component_freq, tab_j$Component_freq)

  res <- c()
  for (i in 1:3) {

    invisible(utils::capture.output(out <-symmoments::callmultmoments(freq[[i]])))

    if (identical(out, -1)) {

      res[i] <- 0

    } else {

      x  <- outer(component[[i]], component[[i]],
                  FUN = Vectorize(SingleSingle, vectorize.args = c(".i", ".j")),
                  .H = H,
                  .Q = Q)

      xVect   <- x[lower.tri(x, diag = TRUE)]

      out_rep <- out$representation
      xVect_m <- matrix(rep(xVect, times = nrow(out_rep)),
                        nrow = nrow(out_rep), byrow = TRUE)^out_rep
      xVect_m <- cbind(out$coefficients, xVect_m)

      res[i] <- sum(matrixStats::rowProds(xVect_m))
    }
  }
  ## Calculate covariance: E(XY) - E(X)*E(Y)
  res <- res[1] - res[2] * res[3]

  ## Return
  res
}
#' @rdname nonlinear_estimation_utilities
f4 <- function(.i, .j, .Q, .H, .var_struc_error, .temp = NULL) {

  Q               <- .Q
  H               <- .H
  var_struc_error <- .var_struc_error
  temp            <- .temp

  eq_i <- temp[[.i]]
  eq_j <- temp[[.j]]

  eq_ij            <- eq_i %o% eq_j
  eq_ij_vec        <- c(eq_ij)
  names(eq_ij_vec) <- c(outer(rownames(eq_ij),
                              colnames(eq_ij),
                              FUN = paste, sep = "."))

  E <- lapply(list(eq_ij_vec, eq_i, eq_j), function(.x) {

    name_i <- names(.x)

    moment <- c()
    for (k in 1:length(name_i)) {

      ## Split term
      term_split <- strsplit(name_i[k], "\\.")

      ## Count instances of each construct name (used for classifying)
      tab_i <- term_split %>%
        table(.) %>%
        as.data.frame(., stringsAsFactors = FALSE)
      colnames(tab_i) <- c("Component", "Component_freq")

      freq      <- tab_i$Component_freq
      component <- tab_i$Component
      select    <- tab_i[tab_i$Component_freq == 1 && grepl("zeta", tab_i$Component)]

      if(length(select) != 0) {

      moment[k] <- 0

      } else {

        invisible(utils::capture.output(out <-symmoments::callmultmoments(freq)))

        if (identical(out, -1)) {

          moment[k] <- 0

        } else {

          x  <- outer(component, component,
                      FUN = Vectorize(f5, vectorize.args = c(".i", ".j")),
                      .H = H,
                      .Q = Q,
                      .var_struc_error = var_struc_error)

          xVect   <- x[lower.tri(x, diag = TRUE)]

          out_rep <- out$representation
          xVect_m <- matrix(rep(xVect, times = nrow(out_rep)),
                            nrow = nrow(out_rep), byrow = TRUE)^out_rep
          xVect_m <- cbind(out$coefficients, xVect_m)

          moment[k]  <- sum(matrixStats::rowProds(xVect_m))
        } # END else
      } # END else
    } # END for k in name_i

    sum(moment * .x)
  }) # END lapply

  ## Calculate covariance: E(XY) - E(X)*E(Y)
  Cov_ij <- E[[1]] - E[[2]] * E[[3]]

  ## Return
  Cov_ij
}
#' @rdname nonlinear_estimation_utilities
f5 <- function(.i, .j, .H, .Q, .var_struc_error) {

  err_i <- grepl("zeta", .i)
  err_j <- grepl("zeta", .j)

  if(err_i == TRUE && err_j == TRUE && .i == .j) {

    x <- .var_struc_error[unlist(strsplit(.i, "_"))[2]]

  } else if (any(err_i == TRUE, err_j == TRUE)) {

    x <- 0

  } else {

    ## Calculate M
    M <- mean(.H[, .i] * .H[, .j])

    ## Calculate denominator D
    D <- prod(.Q[c(.i, .j)])

    ## Calculate VCV element
    if(.i != .j) { # Covariance (rowindex != columnindex)
      x <- M / D
    } else { # Variance: (rowindex == columnindex)
      x <- 1
    }
  }
  ## Return
  x
}

#' Internal: Calculate consistent moments of a non-linear model 
#' 
#' Collection of various moment estimators. See [classifyConstructs] for a list of
#' possible moments.
#' 
#' M is the matrix of the sample counterparts (estimates) of the
#' left-hand side terms in Equation (21) - (24) \insertCite{Dijkstra2014}{cSEM}.
#' The label "M" did not appear in the paper and is only used in the package.
#' Similar is suggested by \insertCite{Wall2000;textual}{cSEM} using classical factor scores.
#' 
#' @param .i Row index
#' @param .j Column index
#' @inheritParams csem_arguments
#'  
#' @name moments
#' @rdname moments
#' @keywords internal

SingleSingle <- function(.i, .j, .Q, .H) {
  
  ## Calculate M
  M <- mean(.H[, .i] * .H[, .j])
  
  ## Calculate denominator D
  D <- prod(.Q[c(.i, .j)])
  
  ## Calculate VCV element
  if(.i != .j) { # Covariance (rowindex != columnindex)
    x <- M / D
  } else { # Variance: (rowindex == columnindex)
    x <- 1
  }
  ## Return
  x
}
#' @rdname moments
SingleQuadratic <- function(.i, .j, .Q, .H) {
  
  ## Label
  i_split <- unlist(strsplit(.i, "\\."))
  j_split <- unlist(strsplit(.j, "\\."))
  ij      <- c(i_split, j_split)
  
  ## Calculate M
  M <- mean(matrixStats::rowProds(.H[, ij]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  x <- M / D
  
  ## Return
  x
}
#' @rdname moments
SingleCubic <- function(.i, .j, .Q, .H) {
  
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Single" term and which one contains
  ## the "Cubic" (=other) term?
  single_tab <- if(term_class_i == "Single") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Single") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Single" term
  # j := name of the "Cubic" term
  # j_split   := the "Cubic" term split into its components
  # j_single  := the name of the component(s) of the "Cubic" term
  # ij        := all components
  
  i           <- names(single_tab)
  j           <- names(other_tab)
  j_split     <- unlist(strsplit(j, "\\."))
  j_single    <- unique(j_split)
  ij          <- c(i, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, c(i, j_single)]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i, j_split)) == 0) {
    # 1.) The "Single" term does not equal the components of the "Cubic" term
    
    x <- (M1 - 3 * (1 - .Q[j_single]^2) * M2) / D
    
  } else {
    # 2.) The "Single" term is equal to the components of the "Cubic" term
    
    x <- (M1 - 3 * (1 - .Q[j_single]^4)) / D
    
  }
  ## Return
  x
}
#' @rdname moments
SingleTwInter <- function(.i, .j, .Q, .H) {
  
  ## Label
  i_split <- unlist(strsplit(.i, "\\."))
  j_split <- unlist(strsplit(.j, "\\."))
  ij      <- c(i_split, j_split)
  
  ## Calculate M
  M <- mean(matrixStats::rowProds(.H[, ij]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  x <- M / D
  
  ## Return
  x
}
#' @rdname moments
SingleThrwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Single" term and which one contains
  ## the "ThrwInter" (=other) term?
  single_tab <- if(term_class_i == "Single") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Single") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Single" term
  # j := name of the "ThrwInter term
  # j_split     := the "ThrwInter" term split into its components
  # ij          := all components
  # ij_match    := the name of the component of the "ThrwInter" term that matches
  #                the "Single" term
  
  i        <- names(single_tab)
  j        <- names(other_tab)
  j_split  <- unlist(strsplit(j, "\\."))
  ij       <- c(i, j_split)
  ij_match <- intersect(i, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, c(setdiff(j_split, ij_match))]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- M1 / D
    
  } else {
    # 2.) The (component of) "Single" matches one of the components of the
    #     "ThrwInter" term
    
    x <- (M1 - (1 - .Q[ij_match]^2) * M2) / D
    
  }
  ## Return
  x
}
#' @rdname moments
SingleQuadTwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Single" term and which one contains
  ## the "QuadTwInter" (=other) term?
  single_tab <- if(term_class_i == "Single") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Single") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Single" term
  # j := name of the "QuadTwInter" term
  # ij          := all components
  # j_split     := the "QuadTwInter" term split into its components
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  
  i           <- names(single_tab)
  j           <- names(other_tab)
  j_split     <- unlist(strsplit(j, "\\."))
  ij          <- c(i, j_split)
  j_single    <- other_tab[[1]][other_tab[[1]]$Component_type == "Single", ]$Component
  j_quadratic <- other_tab[[1]][other_tab[[1]]$Component_type == "Quadratic", ]$Component
  
  ## Calculation ---------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, c(i, j_single)]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i, j_split)) == 0) {
    # 1.) None of the components match
    
    x <- (M1 - (1 - .Q[j_quadratic]^2) * M2) / D
    
  } else if (i == j_single) {
    # 2.) The "Single" term matches the "Single" component of the
    #     "QuadTwInter" term
    
    x <- (M1 - 1) / D + 1
    
  } else {
    # 3.) The "Single" term matches the "Quadratic" component of the
    #     "QuadTwInter" term
    
    x <- (M1 - 3 * (1 - .Q[j_quadratic]^2) * M2) / D
    
  }
  ## Return
  x
}
#' @rdname moments
QuadraticQuadratic <- function(.i, .j, .Q, .H) {
  
  ## Label
  i_split  <- unlist(strsplit(.i, "\\."))
  j_split  <- unlist(strsplit(.j, "\\."))
  i_single <- unique(i_split)
  ij       <- c(i_split, j_split)
  
  ## Calculate M
  M <- mean(matrixStats::rowProds(.H[, ij]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(.i != .j) {
    # 1.) None of the components match
    
    x <- (M - 1) / D
    
  } else {
    # 2.) All components match (terms are identical)
    
    x <- (M - 3 * (1 - .Q[i_single]^4)) / D - 1
    
  }
  ## Return
  x
}
#' @rdname moments
QuadraticCubic <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Quadratic" term and which
  ## one contains the "Cubic" (=other) term?
  quad_tab   <- if(term_class_i == "Quadratic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Quadratic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Quadratic" term
  # j := name of the "Cubic" term
  # i_split   := the "Quadratic" term split into its components
  # j_split   := the "Cubic" term split into its components
  # i_single  := the name of the component(s) of the "Quadratic" term
  # j_single  := the name of the component(s) of the "Cubic" term
  # ij        := all components
  
  i         <- names(quad_tab)
  j         <- names(other_tab)
  i_split   <- unlist(strsplit(i, "\\."))
  j_split   <- unlist(strsplit(j, "\\."))
  i_single  <- unique(i_split)
  j_single  <- unique(j_split)
  ij        <- c(i_split, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M3 <- mean(matrixStats::rowProds(.H[, j_split]))
  
  ## Calculate denominator (D)
  D0 <- prod(.Q[j_split])
  D1 <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(i_single != j_single) {
    # 1.) Non of the components match
    
    x <- (M1 - 3*(1 - .Q[j_single]^2) * M2 - (1 - .Q[i_single]^2) * M3) / D1 -
      M3 / D0
    
  } else {
    # 2.) The component(s) of the "Quadratic" term matche those of the
    #    "Cubic" term
    
    x <- (M1 - 10 * (1 - .Q[j_single]^2) * M3) / D1 - M3/D0
    
  }
  ## Return
  x
}
#' @rdname moments
QuadraticTwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Quadratic" term and which
  ## one contains the "TwInter" (=other) term?
  quad_tab   <- if(term_class_i == "Quadratic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Quadratic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Quadratic" term
  # j := name of the "TwInter" term
  # i_split   := the "Quadratic" term split into its components
  # j_split   := the "TwInter" term split into its components
  # i_single  := the name of the component(s) of the "Quadratic" term
  # ij        := all components
  
  i         <- names(quad_tab)
  j         <- names(other_tab)
  i_split   <- unlist(strsplit(i, "\\."))
  j_split   <- unlist(strsplit(j, "\\."))
  i_single  <- unique(i_split)
  ij        <- c(i_split, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, j_split]))
  
  ## Calculate denominator (D)
  D0 <- prod(.Q[j_split])
  D  <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i_split, j_split)) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - (1 - .Q[i_single]^2) * M2) / D - M2 / D0
    
  } else {
    # 2.) The component(s) of the "Quadratic" term match one of the
    #     component(s) of the "TwInter" term
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2) / D - M2 / D0
    
  }
  ## Return
  x
}
#' @rdname moments
QuadraticThrwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Quadratic" term and which one contains
  ## the "ThrwInter" (=other) term?
  quad_tab   <- if(term_class_i == "Quadratic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Quadratic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Quadratic" term
  # j := name of the "ThrwInter" term
  # i_split   := the "Quadratic" term split into its components
  # j_split   := the "ThrwInter" term split into its components
  # i_single  := the name of the component(s) of the "Quadratic" term
  # ij        := all components
  
  i         <- names(quad_tab)
  j         <- names(other_tab)
  i_split   <- unlist(strsplit(i, "\\."))
  j_split   <- unlist(strsplit(j, "\\."))
  i_single  <- unique(i_split)
  ij        <- c(i_split, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, j_split]))
  
  ## Calculate denominator (D)
  D0 <- prod(.Q[j_split])
  D  <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i_split, j_split)) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - (1 - .Q[i_single]^2) * M2) / D - M2 / D0
    
  } else {
    # 2.) The component(s) of the "Quadratic" match one of the component(s)
    #     of the "ThrwInter" term
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2) / D - M2 / D0
    
  }
  ## Return
  x
}
#' @rdname moments
QuadraticQuadTwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Quadratic" term and which
  ## one contains the "TwInter" (=other) term?
  quad_tab   <- if(term_class_i == "Quadratic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Quadratic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Quadratic" term
  # j := name of the "QuadTwInter" term
  # i_split     := the "Quadratic" term split into its components
  # j_split     := the "QuadTwInter" term split into its components
  # i_single    := the name of the component(s) of the "Quadratic" term
  
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  
  i           <- names(quad_tab)
  j           <- names(other_tab)
  i_split     <- unlist(strsplit(i, "\\."))
  j_split     <- unlist(strsplit(j, "\\."))
  i_single    <- unique(i_split)
  ij          <- c(i_split, j_split)
  j_single    <- other_tab[[1]][other_tab[[1]]$Component_type == "Single", ]$Component
  j_quadratic <- other_tab[[1]][other_tab[[1]]$Component_type == "Quadratic", ]$Component
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, j_split]))
  M3 <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  
  ## Calculate denominator (D)
  D0 <- prod(.Q[j_split])
  D  <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i_split, j_split)) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - (1 - .Q[i_single]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3) / D -
      M2 / D0
    
  } else if (i_single == j_single) {
    # 2.) The component(s) of the "Quadratic" term match the "Single" component
    #     of the "QuadTwInter" term
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3) / D -
      M2 / D0
    
  } else {
    # 3.) The component(s) of the "Quadratic" term match the "Quadratic"
    #     component of the "QuadTwInter" term
    
    x <- (M1 - 6 * (1 - .Q[i_single]^2) * M2) / D -
      M2 / D0
    
  }
  ## Return
  x
}
#' @rdname moments
CubicCubic <- function(.i, .j, .Q, .H) {
  
  ## Label
  # i_split  := the term .i split into its components
  # j_split  := the term .j split into its components
  # i_single := the name of the component of the "Cubic" term .i
  # j_single := the name of the component of the "Cubic" term .j
  # ij       := all components
  
  i_split  <- unlist(strsplit(.i, "\\."))
  j_split  <- unlist(strsplit(.j, "\\."))
  i_single <- unique(i_split)
  j_single <- unique(j_split)
  ij       <- c(i_split, j_split)
  
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M3  <- mean(matrixStats::rowProds(.H[, c(j_split, i_single)]))
  M4  <- mean(matrixStats::rowProds(.H[, c(i_single, j_single)]))
  
  
  ## Calculate denominator (D)
  D0a <- prod(.Q[i_single])
  D0b <- prod(.Q[j_single])
  D  <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(i_single != j_single) {
    # 1.) Non of the components match
    
    x <- (M1 - 3 * (1 - .Q[j_single]^2) * M2 - 3 * (1 - .Q[i_single]^2) * M3 +
            9 * (1 - .Q[i_single]^2) * (1 - .Q[j_single]) * M4) / D -
      (M0a / D0a) * (M0b / D0b)
    
  } else {
    # 2.) All components match (terms are identical)
    
    x <- (M1 - 15 * (1 - .Q[i_single]^2) * M2 + 15 * .Q[i_single]^6 -
            45 * .Q[i_single]^2 + 30) / D -
      (M0a / D0a)^2
    
  }
  ## Return
  x
}
#' @rdname moments
CubicTwInter <- function(.i, .j, .Q, .H) {
  ### Prepartion ---------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Cubic" term and which one contains
  ## the "TwInter" (=other) term?
  cubic_tab  <- if(term_class_i == "Cubic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Cubic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Cubic" term
  # j := name of the "TwInter" term
  # i_split   := the "Cubic" term split into its components
  # j_split   := the "TwInter" term split into its components
  # i_single  := the name of the component(s) of the "Cubic" term
  # ij        := all components
  
  i         <- names(cubic_tab)
  j         <- names(other_tab)
  i_split   <- unlist(strsplit(i, "\\."))
  j_split   <- unlist(strsplit(j, "\\."))
  i_single  <- unique(i_split)
  ij        <- c(i_split, j_split)
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, c(i_single, j_split)]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i_split, j_split)) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2) / D - (M0a * M0b) / D
    
  } else {
    # 2.) The component(s) of the "Cubic" term match one of the component(s)
    #     of the "TwInter" term
    
    x <- (M1 - 6 * (1 - .Q[i_single]^2) * M2) / D - (M0a * M0b) / D
    
  }
  ## Return
  x
}
#' @rdname moments
CubicThrwInter <- function(.i, .j, .Q, .H) {
  ### Prepartion ---------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Cubic" term and which one contains
  ## the "ThrwInter" (=other) term?
  cubic_tab  <- if(term_class_i == "Cubic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Cubic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Cubic" term
  # j := name of the "ThrwInter" term
  # i_split     := the "Cubic" term split into its components
  # j_split     := the "ThrwInter" term split into its components
  # i_single    := the name of the component(s) of the "Cubic" term
  # ij          := all components
  # ij_match    := the name of the component of the "ThrwInter" term that matches
  #                the component of the "Cubic" term
  
  i         <- names(cubic_tab)
  j         <- names(other_tab)
  i_split   <- unlist(strsplit(i, "\\."))
  j_split   <- unlist(strsplit(j, "\\."))
  i_single  <- unique(i_split)
  ij        <- c(i_split, j_split)
  ij_match  <- intersect(i_single, j_split)
  
  ### Prepartion ---------------------------------------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, c(i_single, j_split)]))
  M3  <- mean(matrixStats::rowProds(.H[, setdiff(j_split, ij_match)]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2) / D - (M0a * M0b) / D
    
  } else {
    # 2.) The component(s) of the "Cubic" term match one of the component(s)
    #     of the "ThrwInter" term
    
    x <- (M1 - 6 * (1 - .Q[i_single]^2) * M2 + 3 *
            (1 - .Q[i_single]^2)^2 * M3) / D - (M0a * M0b) / D
    
  }
  ## Return
  x
}
#' @rdname moments
CubicQuadTwInter <- function(.i, .j, .Q, .H) {
  ### Prepartion ---------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "Cubic" term and which one contains
  ## the "QuadTwInter" (=other) term?
  cubic_tab  <- if(term_class_i == "Cubic") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "Cubic") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "Cubic" term
  # j := name of the "QuadTwInter" term
  # i_split     := the "Cubic" term split into its components
  # j_split     := the "QuadTwInter" term split into its components
  # i_single    := the name of the component(s) of the "Cubic" term
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  # ij          := all components
  
  i           <- names(cubic_tab)
  j           <- names(other_tab)
  i_split     <- unlist(strsplit(i, "\\."))
  j_split     <- unlist(strsplit(j, "\\."))
  i_single    <- unique(i_split)
  ij          <- c(i_split, j_split)
  j_single    <- other_tab[[1]][other_tab[[1]]$Component_type == "Single", ]$Component
  j_quadratic <- other_tab[[1]][other_tab[[1]]$Component_type == "Quadratic", ]$Component
  
  ### Prepartion ---------------------------------------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, c(i_single, j_split)]))
  M3  <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M4  <- mean(matrixStats::rowProds(.H[, c(i_single, j_single)]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(intersect(i_split, j_split)) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3 +
            3 * (1 -  .Q[i_single]^2) * (1 - .Q[j_single]^2) * M4) / D -
      (M0a * M0b) / D
    
  } else if (i_single == j_single) {
    # 2.) The component(s) of the "Cubic" term match the "Single" component
    #     of the "QuadTwInter" term
    
    x <- (M1 - 6 * (1 - .Q[i_single]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3 -
            3 * (1 - .Q[i_single]^2) * (.Q[i_single]^2 * .Q[j_quadratic]^2 +
                                          .Q[j_quadratic]^2 - 2)) / D -
      (M0a * M0b) / D
    
  } else {
    # 3.) The component(s) of the "Cubic" term match the "Quadratic" component
    #     of the "QuadTwInter" term
    
    x <- (M1 - 10 * (1 - .Q[i_single]^2) * M2 + 15 *
            (1 - .Q[i_single]^2)^2 * M4) / D - (M0a * M0b) / D
    
  }
  ## Return
  x
}
#' @rdname moments
TwInterTwInter <- function(.i, .j, .Q, .H) {
  ## Label
  # i_split  := the term .i split into its components
  # j_split  := the term .j split into its components
  # i_single := the names of the components of the "TwInter" term .i
  # j_single := the names of the components of the "TwInter" term .j
  # ij       := all components
  # ij_match := the name of the component of .i that match those of .j
  
  i_split  <- unlist(strsplit(.i, "\\."))
  j_split  <- unlist(strsplit(.j, "\\."))
  i_single <- unique(i_split)
  j_single <- unique(j_split)
  ij       <- c(i_split, j_split)
  ij_match <- intersect(i_split, j_split)
  
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_single]))
  M0b <- mean(matrixStats::rowProds(.H[, j_single]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)]]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- M1 / D - M0a * M0b / D
    
  } else if (length(ij_match) == 1){
    # 2.) One of the components of .i matches one of .j
    
    x <- (M1 - (1 - .Q[ij_match]^2) * M2) / D - (M0a * M0b) / D
    
  } else {
    # 3.) All components of .i match those of .j
    
    x <- (M1 - 1) / D + 1 - (M0a * M0b) / D
    
  }
  ## Return
  x
}
#' @rdname moments
TwInterThrwInter <- function(.i, .j, .Q, .H) {
  
  ## Label
  # i_split  := the term .i split into its components
  # j_split  := the term .j split into its components
  # ij       := all components
  # ij_match := the name of the component of .i that match those of .j
  
  i_split  <- unlist(strsplit(.i, "\\."))
  j_split  <- unlist(strsplit(.j, "\\."))
  ij       <- c(i_split, j_split)
  ij_match <- intersect(i_split, j_split)
  
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)], drop = FALSE]))
  M3  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[1])], drop = FALSE]))
  M4  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[2])], drop = FALSE]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- M1 / D - (M0a * M0b) / D
    
  } else if (length(ij_match) == 1){
    # 2.) One of the components of .i matches one of .j
    
    x <- (M1 - (1 - .Q[ij_match]^2) * M2) / D - (M0a * M0b) / D
    
  } else {
    # 3.) All components of .i match with any two of the components of .j
    
    x <- (M1 - (1 - .Q[ij_match[1]]^2) * M3 - (1 - .Q[ij_match[2]]^2) * M4) / D - 
      (M0a * M0b) / D
    
  }
  ## Return x
  x
}
#' @rdname moments
TwInterQuadTwInter <- function(.i, .j, .Q, .H) {
  
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "TwInter" term and which
  ## one contains the "QuadTwInter" (=other) term?
  tw_tab    <- if(term_class_i == "TwInter") {tab_i} else {tab_j}
  other_tab <- if(term_class_i != "TwInter") {tab_i} else {tab_j}
  
  ## Label
  # i := name of the "TwInter" term
  # j := name of the "QuadTwInter" term
  # i_split  := the "TwInter" term split into its components
  # j_split  := the "QuadTwInter" term split into its components
  # ij       := all components
  # ij_match := the names of the components of the "TwInter" term that match
  #             those of the "QuadTwInter" term
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  
  i           <- names(tw_tab)
  j           <- names(other_tab)
  i_split     <- unlist(strsplit(j, "\\."))
  j_split     <- unlist(strsplit(j, "\\."))
  ij          <- c(i_split, j_split)
  ij_match    <- intersect(i_split, j_split)
  j_single    <- other_tab[[1]][other_tab[[1]]$Component_type == "Single", ]$Component
  j_quadratic <- other_tab[[1]][other_tab[[1]]$Component_type == "Quadratic", ]$Component
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% j_quadratic)]]))
  M3 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)]]))
  M4 <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M5 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% j_single)]]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) None of the components of the "TwInter" term match those of the
    #     "QuadTwInter" term
    
    x <- (M1 - (1 - .Q[j_quadratic]^2) * M2) / D - (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && ij_match == j_single) {
    # 2.) One of the components of the "TwInter" term matches the "Single"
    #     component of the "QuadTwInter" term
    
    x <- (M1 - (1 - .Q[j_quadratic]^2) * M2 - (1 - .Q[ij_match]^2) * M3) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && ij_match == j_quadratic) {
    # 3.) One of the components of the "TwInter" term matches the "Quadratic"
    #     component of the "QuadTwInter" term
    
    x <- (M1 - 3 * (1 - .Q[ij_match]^2) * M4) / D -
      (M0a * M0b) / D
    
  } else {
    # 4.) One of the components of the "TwInter" term matches the "Single"
    #     and the other matches the "Quadratic" component  of the
    #     "QuadTwInter" term
    
    x <- (M1 - 3* (1 - .Q[j_quadratic]^2) * M4 - (1 - .Q[j_single]^2) * M5) / D -
      (M0a * M0b) / D
    
  }
  ## Return
  x
}
#' @rdname moments
ThrwInterThrwInter <- function(.i, .j, .Q, .H){
  ## Label
  # i_split     := the "ThrwInter" term .i split into its components
  # j_split     := the "ThrwInter" term .j split into its components
  # ij          := all terms
  # ij_match    := the name of the component of the "ThrwInter" terms that match
  
  i_split    <- unlist(strsplit(.i, "\\."))
  j_split    <- unlist(strsplit(.j, "\\."))
  ij         <- c(i_split, j_split)
  ij_match   <- intersect(i_split, j_split)
  
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)]]))
  M3  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[1])]]))
  M4  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[2])]]))
  M5  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[3])]]))
  
  ## Calculate denominator
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the terms matches
    
    x <- M1 / D - (M0a * M0b) / D
    
  } else if (length(ij_match) == 1) {
    # 2.) Exactly one of the components in the first "ThrwInter" term
    #     matches one in the second term
    
    x <- (M1 - (1 - .Q[ij_match]^2) * M2) / D - (M0a * M0b) / D
    
  } else if (length(ij_match) == 2) {
    # 3.) Two components of the first "ThrwInter" term match with two components
    #     of the other "ThrwInter" term.
    
    x <- (M1 - (1 - .Q[ij_match[1]]^2) * M3 - (1 - .Q[ij_match[2]]^2) * M4 -
            (1 - .Q[ij_match[1]]^2) * (1 - .Q[ij_match[2]]^2) * M2) / D -
      (M0a * M0b) / D
  } else {
    # 4.) All components of the first "ThrwInter" term match those of the second
    
    x <- (M1 - (1 - .Q[ij_match[1]]^2) * M3 - (1 - .Q[ij_match[2]]^2) * M4 -
            (1 - .Q[ij_match[3]]^2) * M5) / D +
      (.Q[ij_match[1]]^2 * .Q[ij_match[2]]^2 * .Q[ij_match[3]]^2 -
         .Q[ij_match[1]]^2 - .Q[ij_match[2]]^2 - .Q[ij_match[3]]^2 + 2) / D -
      (M0a * M0b) / D
    
  }
  ## Return x
  x
}
#' @rdname moments
ThrwInterQuadTwInter <- function(.i, .j, .Q, .H){
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)
  tab_j <- classifyConstructs(.j)
  
  ## Get the classification for term .i
  term_class_i <- unique(tab_i[[.i]]$Term_class)
  
  ## Which classification table contains the "ThrwInter" term and which
  ## one contains the "QuadTwInter" (=other) term?
  thrw_tab   <- if(term_class_i == "ThrwInter") {tab_i} else {tab_j}
  other_tab  <- if(term_class_i != "ThrwInter") {tab_i} else {tab_j}
  
  ## Label
  # i_split     := the "ThrwInter" term .i split into its components
  # j_split     := the "QuadTwInter" term .j split into its components
  # ij          := all component names
  # ij_match    := the name of the component of the two terms that match
  # ij_nomatch  := the name of the component of the two terms that dont match
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  
  i           <- names(thrw_tab)
  j           <- names(other_tab)
  i_split     <- unlist(strsplit(i, "\\."))
  j_split     <- unlist(strsplit(j, "\\."))
  j_single    <- other_tab[[1]][other_tab[[1]]$Component_type == "Single", ]$Component
  j_quadratic <- other_tab[[1]][other_tab[[1]]$Component_type == "Quadratic", ]$Component
  
  ij          <- c(i_split, j_split)
  ij_match    <- intersect(i_split, j_split)
  ij_nomatch  <- setdiff(ij, ij_match)
  
  ### Calculation ---------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1 <- mean(matrixStats::rowProds(.H[, ij]))
  M2 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% j_quadratic)]]))
  M3 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% j_single)]]))
  M4 <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% c(j_quadratic, ij_match))]]))
  M5 <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M6 <- mean(matrixStats::rowProds(.H[, c(ij_nomatch, j_quadratic)]))
  
  ## Calculate denominator
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - (1 - .Q[j_quadratic]^2) * M2) / D - (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && ij_match == j_single) {
    # 2.) One of the components of the "ThrwInter" term matches the "Single"
    #     component of the "QuadTwInter" term.
    
    x <- (M1 - (1 - .Q[j_quadratic]^2) * M2 - (1 - .Q[j_single]^2) * M3 -
            (1 - .Q[j_single]^2) * (1 - .Q[j_quadratic]^2) * M4) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && ij_match == j_quadratic) {
    # 3.) One of the components of the "ThrwInter" term matches the "Quadratic"
    #     component of the "QuadTwInter" term.
    
    x <- (M1 - 3 * (1 - .Q[j_quadratic]^2) * M5) / D - (M0a * M0b) / D
    
  } else {
    # 4.) One component of the "ThrwInter" term matches the "Single" and another
    #     the "Quadratic" component of the "QuadTwInter" term.
    
    x <- (M1 - 3 * (1 - .Q[j_quadratic]^2) * M5 - (1 - .Q[j_single]^2) * M3 -
            3 * (1 - .Q[j_quadratic]^2) * ( 1 - .Q[j_single]^2) *  M6) / D -
      (M0a * M0b) / D
    
  }
  ## Return x
  x
}
#' @rdname moments
QuadTwInercQuadTwInter <- function(.i, .j, .Q, .H) {
  ### Preparation --------------------------------------------------------------
  ## Get the classification table
  tab_i <- classifyConstructs(.i)[[1]]
  tab_j <- classifyConstructs(.j)[[1]]
  
  ## Label
  # i := name of the "Quadratic" term
  # j := name of the "QuadTwInter" term
  # i_split     := the "Quadratic" term split into its components
  # j_split     := the "QuadTwInter" term split into its components
  # i_single    := the name of the component(s) of the "Quadratic" term
  # j_single    := name of the "Single" component of the "QuadTwInter" term
  # j_quadratic := name of the "Quadratic" component of the "QuadTwInter" term
  
  i_split     <- unlist(strsplit(.i, "\\."))
  j_split     <- unlist(strsplit(.j, "\\."))
  ij          <- c(i_split, j_split)
  ij_match    <- intersect(i_split, j_split)
  
  i_single    <- tab_i[tab_i$Component_type == "Single", ]$Component
  i_quadratic <- tab_i[tab_i$Component_type == "Quadratic", ]$Component
  j_single    <- tab_j[tab_j$Component_type == "Single", ]$Component
  j_quadratic <- tab_j[tab_j$Component_type == "Quadratic", ]$Component
  
  ### Calculation --------------------------------------------------------------
  ## Calculate M
  M0a <- mean(matrixStats::rowProds(.H[, i_split]))
  M0b <- mean(matrixStats::rowProds(.H[, j_split]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% i_quadratic)]]))
  M3  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% j_quadratic)]]))
  M4  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% c(i_quadratic, j_quadratic))]]))
  M5  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)]]))
  M6  <- mean(matrixStats::rowProds(.H[, c(i_split, j_single)]))
  M7  <- mean(matrixStats::rowProds(.H[, c(j_split, i_single)]))
  M8  <- mean(matrixStats::rowProds(.H[, c(i_single, j_single)]))
  M9  <- mean(matrixStats::rowProds(.H[, ij[ij == i_quadratic]]))
  
  ## Calculate denominator (D)
  D <- prod(.Q[ij])
  
  ## Calculate VCV element
  if(length(ij_match) == 0) {
    # 1.) Non of the components match
    
    x <- (M1 - (1 - .Q[i_quadratic]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3 -
            (1 - .Q[i_quadratic]^2) * (1 - .Q[j_quadratic]^2) * M4) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && i_single == j_single) {
    # 2.) The "Single" component of the term .i matches the "Single"
    #     component of the term .j
    
    x <- (M1 - (1 - .Q[i_quadratic]^2) * M2 - (1 - .Q[j_quadratic]^2) * M3 -
            (1 - .Q[ij_match]^2) * M5) / D +
      (.Q[i_quadratic]^2 * .Q[j_quadratic]^2 * .Q[ij_match]^2 -
         .Q[i_quadratic]^2 - .Q[j_quadratic]^2 - .Q[ij_match]^2 + 2) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && i_single == j_quadratic) {
    # 3.) The "Single" component of the term .i matches the "Quadratic"
    #     component of the term .j
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M6 - (1 - .Q[i_quadratic]) * M7 -
            3 * (1 - .Q[i_single]^2) * (1 - .Q[i_quadratic]^2) * M8) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && i_quadratic == j_single) {
    # 4.) The "Single" component of the term .j matches the "Quadratic"
    #     component of the term .i
    
    x <- (M1 - 3 * (1 - .Q[j_single]^2) * M7 - (1 - .Q[j_quadratic]) * M6 -
            3 * (1 - .Q[j_single]^2) * (1 - .Q[j_quadratic]^2) * M8) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 1 && i_quadratic == j_quadratic) {
    # 5.) The "Quadratic" component of the term .i matches the "Quadratic"
    #     component of the term .j
    
    x <- (M1 - 6 * (1 - .Q[i_quadratic]^2) * M6 +
            3 * (1 - .Q[i_quadratic]^2)^2 * M8) / D -
      (M0a * M0b) / D
    
  } else if (length(ij_match) == 2 &&
             i_single         == j_single &&
             i_quadratic      == j_quadratic) {
    # 6.) Both components match. The "Single" component of .i matches the "Single"
    #     component of .j and the "Quadratic" component of .i matches the
    #     "Quadratic" component of .j
    
    x <- (M1 - (1 - .Q[i_single]^2) * M9 - 6 * (1 - .Q[i_quadratic]^2) * M6 -
            3 * (1 - .Q[i_quadratic]^2) *
            (.Q[i_single]^2 * .Q[i_quadratic]^2 + .Q[i_single]^2 - 2)) / D -
      (M0a * M0b) / D
    
  } else {
    # 7.) Both components match. The "Single" component of .i matches the "Quadratic"
    #     component of .j and the "Quadratic" component of .i matches the
    #     "Single" component of .j
    
    x <- (M1 - 3 * (1 - .Q[i_single]^2) * M6 - 3 * (1 - .Q[i_quadratic]^2) * M7 +
            9 * (1 - .Q[i_quadratic]^2) * (1 - .Q[i_quadratic]^2) * M8) -
      (M0a * M0b) / D
    
  }
  ## Return
  x
}

