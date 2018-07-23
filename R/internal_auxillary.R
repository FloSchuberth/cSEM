#' Internal: Calculate proxies
#'
#' Calculates proxy values for the J constructs of the structural model.
#'
#' More on that
#'
#' @inheritParams csem_arguments
#'
#' @return The (N x J) matrix of proxy values with column names equal to
#'   the names of the J constructs.
#'
calculateProxies <- function(.X = NULL, .W = NULL) {

  ## Proxies for the linear terms/latent variables
  H <- .X %*% t(.W)
  
  ## Alternative 
  # Note: functions like var, sd, scale, and cov use n-1. 
  # H <- apply(.X, 2, function(x) {(x - mean(x)) / (sd(x) * sqrt((length(x) - 1) / length(x)))}) %*% t(.W)
  return(H)
}

#' Internal: Calculate proxy variance-covariance matrix
#'
#' Calculates the variance-covariance (VCV) matrix of the proxies.
#'
#' More on that
#'
#' @inheritParams csem_arguments
#'
#' @return A (J x J) proxy VCV matrix.
#'
calculateProxyVCV <- function(.S = NULL, .W = NULL) {

  x <- .W %*% .S %*% t(.W)

  # Due to floting point errors may not be symmetric anymore. In order
  # prevent that replace the lower triangular elements by the upper
  # triangular elements

  x[lower.tri(x)] <- t(x)[lower.tri(x)]

  ## Return
  x
}

#' Internal: Calculate proxy-construct covariance
#'
#' Calculates the proxy-covariance covariance (Q)
#' (also called reliability coefficient).
#'
#' More on that
#'
#' @inheritParams csem_arguments
#'
#' @return A vector of proxy-construct covariances.
#'
calculateProxyConstructCV <- function(
  .W             = NULL,
  .csem_model    = NULL,
  .disattenuate  = NULL,
  .modes         = NULL,
  .correction_factors = NULL,
  .reliabilities = NULL
  ) {

  x <- rep(1, times = nrow(.W))
  names(x) <- rownames(.W)

  if(is.null(.reliabilities) & .disattenuate == TRUE) {
    ## Get names of constructs modeled as composites
    names_c  <- names(.csem_model$construct_type[.csem_model$construct_type == "Composite"])
    ## Get names of constructs modeled as common factors
    names_cf <- setdiff(rownames(.csem_model$structural), names_c)
    ## Get names of the common factors whose weights where estimated with "ModeA"
    names_modeA <- intersect(names(.modes[.modes == "ModeA"]), names_cf)
    ## Get names of the common factors whose weights where estimated with "ModeB"
    names_modeB <- intersect(names(.modes[.modes == "ModeB"]), names_cf)

    if(length(names_modeA) > 0) {
      
      x_modeA <- c(diag(.W[names_modeA, ] %*% t(.W[names_modeA, ])) %*%
                     diag(.correction_factors[names_modeA]))
    }
    
    if(length(names_modeB) > 0) {
      x_modeB <- .correction_factors[names_modeB]
    } else {
      x_modeB <- 1
    }
    
    x[names_modeA] <- x_modeA
    x[names_modeB] <- x_modeB
    
  
  } else if(is.null(.reliabilities) & .disattenuate == FALSE) {
    
    return (x)
    
  } else if(!is.null(.reliabilities)) {
    
    ## Check construct names:
    # Do all construct names in .reliabilities match the construct
    # names used in the model?
    tmp <- setdiff(names(.reliabilities), rownames(.W))
    
    if(length(tmp) != 0) {
      stop("Construct name(s): ", paste0("`", tmp, "`", collapse = ", "), 
           " provided to `.reliabilities`", 
           ifelse(length(tmp) == 1, " is", " are"), " unknown.", call. = FALSE)
    }
    
    # Check whether defined external reliabilities are correctly defined
    if(any(.reliabilities > 1 | .reliabilities < 0)) {
      stop('Reliabilities must be between 0 and 1.', call. = FALSE)
    }
    
    x[names(.reliabilities)] <- sqrt(.reliabilities)
  } # END if
  
  return(x)
}

#' Internal: Calculate construct VCV matrix
#'
#' Calculates the covariance matrix of the constructs.
#'
#' More description here
#'
#' @inheritParams
#'
#' @return The (J x J) construct VCV matrix. Disattenuated if requested.
#'
calculateConstructVCV <- function(.C, .Q, .csem_model) {

  f <- function(.i, .j) { .C[.i, .j] / (.Q[.i] * .Q[.j]) }
  m <- rownames(.csem_model$measurement)

  x <- outer(m, m, FUN = Vectorize(f))
  diag(x) <- 1
  rownames(x) <- colnames(x) <- m

  ## Return
  return(x)
}

#' Internal: Classify interaction terms by type
#'
#' Classifies interaction between constructs according to their type.
#'
#' Classification is required to estimate nonlinear structural relationships.
#' Note that exponential terms are modeled as "interactions with itself"
#' (i.e. x^3 = x.x.x).  The following terms and classifications are currently
#' supported:
#'
#' \itemize{
#' \item Single term: e.g., `eta1`
#' \item Quadratic term: e.g., `eta1.eta1`
#' \item Cubic term: e.g., `eta1.eta1.eta1`
#' \item Two-way interactions term: e.g., `eta1.eta2`
#' \item Three-way interactions term: e.g., `eta1.eta2.eta3`
#' \item Quadratic and two-way interaction term: e.g., `eta1.eta1.eta3`"
#' }
#'
#' @usage classifyConstructs(.terms)
#'
#' @inheritParams csem_arguments
#'
#' @return A named list of length equal to the number of terms provided containing
#'   a data frame with columns "*Term_class*", "*Component*",
#'   "*Component_type*", and "*Component_freq*".
#'
#' @examples
#' x <- c("eta1.eta2", "eta1", "eta3.eta1.eta1")
#'
#' classifyConstructs(x)
#'
classifyConstructs <- function(.terms) {
  ## Split term
  terms_split <- strsplit(.terms, "\\.")

  ## Count instances of each construct name (used for classifying)
  terms_classified <- lapply(terms_split, function(.x) {
    x <- .x %>%
      table(.) %>%
      as.data.frame(., stringsAsFactors = FALSE)

    ## To save typing
    a <- sum(x$Freq)
    b <- length(unique(x$.))

    ## Do the actual classification --------------------------------------------

    if(a > 3) {
      stop("The nonlinear term(s): ", paste0("`", .terms, "`", collapse = ", "), 
           ifelse(length(.terms == 1, " is", " are")), " currently not supported.\n",
           "Please see ?classifyConstructs for a list of supported terms.",
           call. = FALSE)
    } else {
      switch(a,
             "1" = {
               x <- data.frame("Term_class"     = "Single",
                               "Component"      = x$.,
                               "Component_type" = "Single",
                               "Component_freq" = x$Freq,
                               stringsAsFactors = FALSE)
             },
             "2" = {
               if(b == 1) {
                 x <- data.frame("Term_class"     = "Quadratic",
                                 "Component"      = x$.,
                                 "Component_type" = "Quadratic",
                                 "Component_freq" = x$Freq,
                                 stringsAsFactors = FALSE)
               } else if (b == 2){
                 x <- data.frame("Term_class"     = "TwInter",
                                 "Component"      = x$.,
                                 "Component_type" = c("Single", "Single"),
                                 "Component_freq" = x$Freq,
                                 stringsAsFactors = FALSE)
               }
             },
             "3" = {
               if(b == 1) {
                 x <- data.frame("Term_class"     = "Cubic",
                                 "Component"      = x$.,
                                 "Component_type" = "Cubic",
                                 "Component_freq" = x$Freq,
                                 stringsAsFactors = FALSE)
               } else if (b == 2) {
                 x <- data.frame("Term_class"     = "QuadTwInter",
                                 "Component"      = x$.,
                                 "Component_type" =
                                   ifelse(x$Freq == 1, "Single", "Quadratic"),
                                 "Component_freq" = x$Freq,
                                 stringsAsFactors = FALSE)
               } else if (b == 3) {

                 x <- data.frame("Term_class"     = "ThrwInter",
                                 "Component"      = x$.,
                                 "Component_type" =
                                   c("Single", "Single", "Single"),
                                 "Component_freq" = x$Freq,
                                 stringsAsFactors = FALSE)
               }
             }
      ) # END switch
    } # END else
  }) # END lapply

  names(terms_classified) <- unlist(.terms)
  terms_classified
}

#' Internal: Scale weights
#'
#' Scales weights to ensure that the proxies have variance of one.
#'
#' More description here
#'
#' @usage scaleWeights(.S = NULL, .W = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of scaled weights.

scaleWeights <- function(.S = NULL, .W = NULL) {

  ## Calculate the variance of the proxies:
  var_proxies <- diag(.W %*% .S %*% t(.W))

  ## Scale the weights to ensure that the proxies have a variance of one
  W_scaled <- solve(diag(sqrt(var_proxies))) %*% .W

  ## Assign rownames and colnames to the scaled weights and return
  rownames(W_scaled) <- rownames(.W)
  colnames(W_scaled) <- colnames(.W)

  return(W_scaled)
}

#' Internal: A symmetric version of setdiff()
#'
#' A symmetric version of setdiff()
#'
#' More description here
#'
#' @usage setdiff2(x, y)
#'
#' @param x,y elements to compare.
#'
#' @return A vector containing the non-overlapping elements.

setdiff2 <- function(x,  y) { 
  
  setdiff(union(x, y), intersect(x, y))
}
