# Internal functions used by the estimator function

f1 <- function(.i, .j) {

  tab_i <- classifyConstructs(.i)[[1]]
  tab_j <- classifyConstructs(.j)[[1]]

  class_ij <- paste(unique(tab_i$Term_class),
                    unique(tab_j$Term_class), sep = "_")

  ## Rename joint class so "class1_class2" = "class2_class1"
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
          # Dont change the object if non of the above apply
          # Note: If non of the alternatives apply,
          # switch evaluates the last alternative...not sure anymore if thats
          # true. Needs to be checked!
          "None_of_the_above"     = {class_ij}
  )
  ## Return
  class_ij
}

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

    invisible(capture.output(out <-symmoments::callmultmoments(freq[[i]])))

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

        invisible(capture.output(out <-symmoments::callmultmoments(freq)))

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
