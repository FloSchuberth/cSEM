# Calculate VCV elements
#
# Does the estimation
#
# More on that
#
# Note: M is the matrix of the sample counterparts (estimates) of the
# left-hand side terms in equation (21) - (24) in the
# Dijkstra & Schermelleh-Engel (2014) paper. The label "M" did not appear in
# the paper and is only used in the package.
#
# @param .i xx
# @param .j xx
# @param .Q xx
# @param .H xx
#
# @return to be written
#

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
  D <- prod(.Q[c(i, j_split)])

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

    x <- (M1 - (1 - .Q[i_single]^2) * M2 - (1 - .Q[j_single]^2) * M3) / D -
      M3 / D0

  } else {
    # 2.) The component(s) of the "Quadratic" term matche those of the
    #    "Cubic" term

    x <- (M1 - 10 * (1 - .Q[i_single]^2) * M3) / D - M3/D0

  }
  ## Return
  x
}

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
            (1 - .Q[i_single]^2) * M4) / D - (M0a * M0b) / D

  }
  ## Return
  x
}

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
  M0a <- mean(matrixStats::rowProds(.H[, i_single]))
  M0b <- mean(matrixStats::rowProds(.H[, j_single]))
  M1  <- mean(matrixStats::rowProds(.H[, ij]))
  M2  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match)]]))
  M3  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[1])]]))
  M4  <- mean(matrixStats::rowProds(.H[, ij[!(ij %in% ij_match[2])]]))

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

    x <- (M1 - (1 - .Q[ij_match[1]]^2) * M3 - .Q[ij_match[2]]^2 * M4) / D + 1 -
      (M0a * M0b) / D

  }
  ## Return x
  x
}

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
    # 2.) Exactly one of the components in the first "ThrwInter" term matches
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
  # j_split     := the "ThrwInter" term .j split into its components
  # ij          := all terms
  # ij_match    := the name of the component of the "ThrwInter" terms that match

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

