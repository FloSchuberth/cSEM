#' Internal: Calculate PLSc correction factors
#'
#' Calculates the correction factor used by PLSc.
#'
#' Currently, seven approaches are available:
#'
#' \itemize{
#' \item "dist_squared_euclid" (default)
#' \item "dist_euclid_weighted"
#' \item "fisher_transformed"
#' \item "mean_geometric"
#' \item "mean_harmonic"
#' \item "mean_arithmetic"
#' \item "geo_of_harmonic" (not yet implemented)
#' }
#'
#' See \insertCite{Dijkstra2013}{cSEM} for details.
#' @usage calculateCorrectionFactors(
#'  .S               = args_default()$.S,
#'  .W               = args_default()$.W,
#'  .modes           = args_default()$.modes,
#'  .csem_model      = args_default()$.csem_model,
#'  .PLS_approach_cf = args_default()$.PLS_approach_cf
#'  )
#' @inheritParams csem_arguments
#'
#' @return A numeric vector of correction factors with element names equal
#'   to the names of the J constructs used in the measurement model.
#' @references
#'   \insertAllCited{}
#' @keywords internal

calculateCorrectionFactors <- function(
  .S               = args_default()$.S,
  .W               = args_default()$.W,
  .modes           = args_default()$.modes,
  .csem_model      = args_default()$.csem_model,
  .PLS_approach_cf = args_default()$.PLS_approach_cf
) {
  
  ### Compute correction factors  ----------------------------------------------
  correction_factors <- vector(mode = "double", length = nrow(.W))
  names(correction_factors) <- rownames(.W)
  
  L <- .W %*% .S * .csem_model$measurement
  
  for(j in rownames(.W)) {
    
    ## Depending on the mode: extract vector of weights or indicator-proxy
    # correlations (composite loadings) of block j 
    if(.modes[j] == "modeA") {
      w_j <- .W[j, ] %>%
        .[. != 0] %>%
        as.matrix(.)
    } else {
      w_j <- L[j, ] %>%
        .[. != 0] %>%
        as.matrix(.)
    }
    
    ## Check if single indicator block or composite; If yes, set cf to 1
    if(nrow(w_j) == 1 | .csem_model$construct_type[j] == "Composite") {
      correction_factors[j] <- 1
    } else {
      ## Extract relevant objects
      E_jj <- .csem_model$error_cor[rownames(w_j), rownames(w_j)]
      S_jj <- .S[rownames(w_j), rownames(w_j)]
      W_jj <- w_j %*% t(w_j)
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of S_jj (result is vectorized)
      S_vect <- replace(S_jj, which(E_jj == 1), 0) %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[. != 0]
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of W_jj (vectorized)
      W_vect <- replace(W_jj, which(E_jj == 1), 0)  %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[. != 0]
      
      if(length(S_vect) == 0) {
        stop("At least one pair of indicators with uncorrelated measurement
             errors required.\n Please revise your model.",
             call. = FALSE)
      }
      
      ## Do the actual computation ---------------------------------------------
      switch (.PLS_approach_cf,
              "dist_squared_euclid"          = {
                cf <- sum(W_vect * S_vect) / sum(W_vect^2)
              },
              "dist_euclid_weighted" = {
                weights <- 1 / (1 - S_vect^2)
                cf      <- sum(W_vect * S_vect * weights) / sum(W_vect^2 * weights)
              },
              "fisher_transformed"   = {
                # Function to be minimized
                temp_fun <- function(.c, .W_vect, .S_vect){
                  sum((0.5*log((1 + .S_vect) / (1 - .S_vect)) -
                         0.5*log((1 + .c*.W_vect) / (1 - .c*W_vect)))^2)
                }
                
                # Optimaziation
                temp_optim <- optim(fn = temp_fun, par = 0.5, method = "BFGS",
                                    .W_vect = W_vect, .S_vect = S_vect)
                cf <- temp_optim$par
              },
              "mean_geometric"       = {
                cf <- prod(S_vect/W_vect)^(1/length(S_vect))
              },
              "mean_arithmetic"      = {
                cf <- mean(S_vect/W_vect)
              },
              "mean_harmonic"        = {
                cf <- 1/mean(1/(S_vect/W_vect))
              },
              "geo_of_harmonic"      = {stop("not implemented yet")}
      )
      
      ## Compute absolute value and take the sqrt since cf = c^2
      correction_factors[j] <- sqrt(abs(cf))
    }
    }
  return(correction_factors)
  ### For maintenance: ### -------------------------------------------------------------------------
  # w_j (K_j x 1)      := Column vector of indicator weights of block j
  # W_vect (1 x 2*K_j) := Vector of off-diagonal elements of w_jw'_j
  # S_vect (1 x 2*K_j) := Vector of off-diagonal elements of S_jj
}

#' Internal: Calculate composite variance-covariance matrix
#'
#' Calculate the variance-covariance (VCV) matrix of the composites/proxies.
#'
#' @usage calculateCompositeVCV(
#'  .S  = args_default()$.S,
#'  .W  = args_default()$.W
#'  )
#'  
#' @inheritParams csem_arguments
#'
#' @return A (J x J) composite VCV matrix.
#' @keywords internal

calculateCompositeVCV <- function(
  .S = args_default()$.S, 
  .W = args_default()$.W
  ){

  x <- .W %*% .S %*% t(.W)

  # Due to floting point errors may not be symmetric anymore. In order
  # prevent that replace the lower triangular elements by the upper
  # triangular elements

  x[lower.tri(x)] <- t(x)[lower.tri(x)]

  ## Return
  x
}

#' Internal: Calculate construct variance-covariance matrix
#'
#' Calculate the variance-covariance matrix (VCV) of the constructs.
#'
#' @usage calculateConstructVCV(
#'  .C          = args_default()$.C, 
#'  .Q          = args_default()$.Q, 
#'  .csem_model = args_default()$.csem_model
#'  )
#'  
#' @inheritParams csem_arguments
#'
#' @return The (J x J) construct VCV matrix. Disattenuated if requested.
#' @keywords internal

calculateConstructVCV <- function(
  .C          = args_default()$.C, 
  .Q          = args_default()$.Q, 
  .csem_model = args_default()$.csem_model
  ) {

  f <- function(.i, .j) { .C[.i, .j] / (.Q[.i] * .Q[.j]) }
  m <- rownames(.csem_model$measurement)

  x <- outer(m, m, FUN = Vectorize(f))
  diag(x) <- 1
  rownames(x) <- colnames(x) <- m

  ## Return
  return(x)
}

#' Internal: Calculate indicator correlation matrix
#' 
#' Calculate the indicator correlation matrix using conventional or robust methods.
#' 
#' If `"none"` the method depends on the type of column of `.X_cleaned`:
#' \describe{
#'   \item{`Numeric-numeric`}{Bravais-Pearson product-moment correlation}
#'   \item{`Numeric-factor`}{Polyserial correlation}
#'   \item{`Factor-factor`}{Polychoric correlation}
#' }
#' Note: logical input is treated as a 0-1 factor variable.
#' 
#' If  `"mcd"` (= minimum covariance determinant) a robust covariance is estimated
#' via `MASS::cov.rob()`. See `?MASS::cov.rob()` for details.
#'
#' @usage calculateIndicatorCor(
#'   .X_cleaned           = args_default()$.X_cleaned, 
#'   .approach_cor_robust = args_default()$.approach_cor
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @return A list containing the (K x K) indicator correlation matrix S. 
#' @keywords internal

calculateIndicatorCor <- function(
  .X_cleaned           = args_default()$.X_cleaned,
  .approach_cor_robust = args_default()$.approach_cor
){
  
  only_numeric_cols <- all(unlist(lapply(.X_cleaned, is.numeric)))
  
  if(.approach_cor_robust != "none" && !only_numeric_cols) {
    stop2("Setting `.approach_cor_robust = ", .approach_cor_robust, "` requires all",
          " columns of .data to be numeric.")
  }
  
  ## polycor::hetcor() is relatively slow. If all columns are numeric use cor
  ## directly
  switch (.approach_cor_robust,
          "none" = {
            if(only_numeric_cols) {
              S <- stats::cor(.X_cleaned)
              cor_type <- "Bravais-Pearson" 
            } else {
              # Pd is TRUE by default. See ?hetcor for details
              temp <- polycor::hetcor(.X_cleaned, std.err = FALSE, pd = TRUE)
              S    <- temp$correlations
              cor_type <- unique(c(temp$type))
              cor_type <- cor_type[which(nchar(cor_type) != 0)] # delete '""' 
            }
          },

          "mcd" = {
            S <- MASS::cov.rob(.X_cleaned, cor = TRUE, method = "mcd")$cor
            S[upper.tri(S) == TRUE] = t(S)[upper.tri(S) == TRUE]

            cor_type <-  "Robust (MCD)"
          },
          "spearman" = {
            S <- cor(.X_cleaned, method = "spearman")

            cor_type <-  "Robust (Spearman)"
          }
  )
  # (TODO) not sure how to name the "type" yet and what to do with it. Theoretically,
  # a polycoric correlation could also be used with GSCA or some other non-PLS-PM method.
  list(S = S, cor_type = cor_type)
}

#' Internal: Calculate Reliabilities
#'
#' (TODO)
#'  
#' @inheritParams csem_arguments
#'
#' @keywords internal

calculateReliabilities <- function(
  .X                = args_default()$.X,
  .S                = args_default()$.S,
  .W                = args_default()$.W,
  .approach_weights = args_default()$.approach_weights,
  .csem_model       = args_default()$.csem_model,
  .disattenuate     = args_default()$.disattenuate,
  .PLS_approach_cf  = args_default()$.PLS_approach_cf,
  .reliabilities    = args_default()$.reliabilities
){
  modes   <- .W$Modes
  W        <- .W$W
  names_cf <- names(.csem_model$construct_type[.csem_model$construct_type == "Common factor"])
  names_c  <- setdiff(names(.csem_model$construct_type), names_cf)
  Q        <- rep(1, times = nrow(W))
  names(Q) <- rownames(W)
  
  Lambda   <- W %*% .S * .csem_model$measurement # composite loadings
  # These are the defaults. If disattenuation is requested all loadings/Q's
  # that belong to a construct modeled as a common factor will be replaced now.
  
  if(is.null(.reliabilities)) {
    ## Congeneric reliability: Q^2 = (w' * lambda)^2 where lambda is a consistent estimator
    ## of the true factor loading.
    ## Approaches differ in the way the loadings are calculated but in the end
    ## it is always: Q = w' lambda.
    
    if(.approach_weights == "PLS-PM") {
      
      if(.disattenuate) {
        ## Consistent loadings are obtained using PLSc, which uses the fact that
        ## lambda = c * w, where c := correction factor
        
        ## 1. Compute correction factor (c)
        correction_factors <- calculateCorrectionFactors(
          .S               = .S,
          .W               = W,
          .modes           = modes,
          .csem_model      = .csem_model,
          .PLS_approach_cf = .PLS_approach_cf
        )
        
        ## 2. Compute consistent loadings and Q (Composite/proxy-construct correlation)
        ## for constructs modeled as common factors.
        for(j in names_cf) {
          Lambda[j, ] <- correction_factors[j] * W[j, ]
          Q[j]        <- c(W[j, ] %*% Lambda[j, ])
        }
      } 
    } else if(.approach_weights == "GSCA") {
      
      if(.disattenuate) {
        # Currently, GSCAm only supports pure common factor models. This may change
        # in the future.
        
        # Compute consistent loadings and Q (Composite/proxy-construct correlation)
        # Consistent factor loadings are obtained from GSCAm.
        
        for(j in rownames(Lambda)) {
          Lambda[j, ] <- .W$C[j, ]
          Q[j]        <- c(W[j, ] %*% Lambda[j, ]) 
        }
      }
    } else if(.approach_weights %in% c("unit", "bartlett", "regression", "PCA", 
                "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR")) {
    
      #   Note: "bartlett" and "regression" weights are obtained AFTER running a CFA. 
      #   Therefore loadings are consistent and as such "disattenuation"
      #   has already implicitly happend. Hence, it is impossible to 
      #   specify "bartlett" or "regression" when .disattenuate = FALSE.
      if(!.disattenuate & .approach_weights %in% c("bartlett", "regression")) {
        stop2("Unable to use `bartlett` or `regression` weights when `.disattenuate = FALSE`.")
      }
      
      #   Note: Only necessary if at least one common factor is in the model
      if(.disattenuate & any(.csem_model$construct_type == "Common factor")) {
        ## "Croon" approach
        # The Croon reliabilities assume that the scores are built by sumscores 
        # (i.e. unit weights).
        # Moreover, all constructs must be modeled as common factors.
        
        if(any(.csem_model$construct_type == "Composite")) {
          stop2("The following error occured in the `calculateReliabilities()` function:\n",
                "Disattenuation only applicable if all constructs are modeled as common factors.")
        }
        
        for(j in names_cf) {
          indicator_names <- colnames(Lambda[j, Lambda[j, ] != 0, drop = FALSE])
          model <- paste0(j, '=~', paste0(indicator_names, collapse = '+'))
          
          # Run CFA for each construct separately to obtain the consistent loadings
          res_lavaan <- lavaan::cfa(model, .X, std.lv = TRUE, std.ov = TRUE)
          
          if(.approach_weights %in% c("bartlett", "regression")) {
            
            method <- ifelse(.approach_weights == "bartlett", "Bartlett", "regression")
            
            # Compute factor scores using "methode"
            scores <- lavaan::lavPredict(res_lavaan, fsm = TRUE, method = method)
            
            # Get weights and scale
            w <- c(attr(scores, 'fsm')[[1]])
            # w <- w /  c(sqrt(w %*% .S[indicator_names, indicator_names, drop = FALSE] %*% w))
            W[j, indicator_names] <- w
          }
          
          # Extract loading estimates
          lambda           <- res_lavaan@Model@GLIST$lambda
          dimnames(lambda) <- res_lavaan@Model@dimNames[[1]]
          
          # Extract weights for construct j
          w <- W[j, indicator_names, drop = FALSE]
          
          # Reorder indicator names to ensure they are identical to the order in cSEM
          lambda <- lambda[colnames(w), ]
          Lambda[j, indicator_names] <- lambda
          
          # Compute composite/proxy-construct correlation: Q 
          Q[j] <- w %*% lambda
        }
      }
    } else {
      stop2("The following error occured in the `calculateReliabilities()` function:\n",
            .approach_weights, " is an unknown weight scheme.")
    }
  } else {
    
    ## Check if weighting scheme is PLS:
    #  Note: In PLS we have lambda = c * w and Q = (w'w) *c --> c = Q / (w'w) 
    #        and therefore lambda = w'Q / (w'w). For other approaches we need to
    #        figure out how weights, loadings and reliabilities are related s.t.
    #        we can solve for lambda depending on w and Q.
    #        Currently we are only aware of such an approach for PLS/PLSc
    
    if(.approach_weights != "PLS-PM") {
      stop2("The following error occured in the `calculateReliabilities()` function:\n",
            "User-provided reliabilities are not supported for weighting approach: ", 
            .approach_weights)
    }
    
    if(!.disattenuate) {
      .disattenuate <- TRUE
      warning2("`.disattenuate = FALSE` is set to `TRUE`.")
    }
    
    # Do all construct names in .reliabilities match the construct
    # names used in the model?
    tmp1 <- setdiff(names(.reliabilities), rownames(W))
    
    if(length(tmp1) != 0) {
      stop2("Construct name(s): ", paste0("`", tmp1, "`", collapse = ", "), 
           " provided to `.reliabilities`", 
           ifelse(length(tmp1) == 1, " is", " are"), " unknown.")
    }
    
    ## Compute reliabilities not supplied by the user 
    
    #  Only for those common factors, since for composites they are either 
    #  supplied or 1.
    Q[names(.reliabilities)] <- sqrt(.reliabilities) # replace the one already supplied
    
    names_cf_not_supplied <- setdiff(names_cf, names(.reliabilities))
    # Get names of the common factors whose weights where estimated with "modeA"
    names_modeA <- intersect(names(modes[modes == "modeA"]), names_cf_not_supplied)
    # Get names of the common factors whose weights where estimated with "modeB"
    names_modeB <- intersect(names(modes[modes == "modeB"]), names_cf_not_supplied)
    
    correction_factors <- calculateCorrectionFactors(
      .S               = .S,
      .W               = W,
      .modes           = modes,
      .csem_model      = .csem_model,
      .PLS_approach_cf = .PLS_approach_cf
    )
    
    if(length(names_cf_not_supplied) != 0) {
      if(length(names_modeA) > 0) {
        
        Q_modeA <- c(diag(W[names_modeA, , drop = FALSE] %*% t(W[names_modeA, ,drop = FALSE])) %*%
                       diag(correction_factors[names_modeA], nrow = length(names_modeA)))
        Q[names_modeA] <- Q_modeA
      }
      
      if(length(names_modeB) > 0) {
        Q_modeB <- correction_factors[names_modeB]
        Q[names_modeB] <- Q_modeB
      }
    }

    ## Compute Loadings based on reliabilities
    # Loadings are calculated for common factors and for those composites
    # that a reliability is give.
    names_cf_rel <- union(names(.reliabilities), names_cf)
    names_modeA  <- intersect(names(modes[modes == "modeA"]), names_cf_rel)
    names_modeB  <- intersect(names(modes[modes == "modeB"]), names_cf_rel)
    # if(length(names_cf) > 0) {
    if(length(names_modeA) > 0) {
      ## Disattenuate loadings and cross-loadings 
      for(i in names_modeA) {
        w  <- W[i, ] # becomes a vector!
        w1 <- w[which(w != 0)]
        w2 <- w[which(w == 0)]
        
        Lambda[i, names(w1)] <- Q[i] * w1 / c(t(w1) %*% w1)
        Lambda[i, names(w2)] <- Lambda[i, names(w2)] / Q[i]
        
      }
    }
    if(length(names_modeB) > 0) {
      # Composite loading estimates 
      L <- Lambda*(W != 0)
      
      for(i in names_modeB){
        temp  <- L[i, ] # becomes a vector!
        temp1 <- temp[which(temp != 0)]
        temp2 <- temp[which(temp == 0)]
        
        Lambda[i, names(temp1)] <- temp1*Q[i] 
        Lambda[i, names(temp2)] <- Lambda[i, names(temp2)] / Q[i]
      } # END for i in names_modeB
    } # END if(length(names_modeB))
    # } # END if(length(names_cf))
  } # END if(is.null(.reliabilities))
  
  ## Return Loadings, Reliabilities, and Cross loadings
  out <- list("Lambda" = Lambda, "Q2" = Q^2, "W" = W)
  return(out)
}

#' Internal: Set the dominant indicator
#' 
#' Set the dominant indicator for each construct.
#'
#' @usage setDominantIndicator(
#'  .W                   = W,
#'  .dominant_indicators = .dominant_indicators
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @return The (J x K) matrix of weights with the dominant indicator set.
#' @keywords internal
#'
setDominantIndicator <- function(
  .W                   = W,
  .dominant_indicators = .dominant_indicators
) {
  ## Check construct names:
  # Do all construct names in .dominant_indicators match the construct
  # names used in the model?
  tmp <- setdiff(names(.dominant_indicators), rownames(.W))
  
  if(length(tmp) != 0) {
    stop("Construct name(s): ", paste0("`", tmp, "`", collapse = ", "), 
         " provided to `.dominant_indicators`", 
         ifelse(length(tmp) == 1, " is", " are"), " unknown.", call. = FALSE)
  }
  
  ## Check indicators
  # Do all indicators names in .dominant_indicators match the indicator
  # names used in the model?
  tmp <- setdiff(.dominant_indicators, colnames(.W))
  
  if(length(tmp) != 0) {
    stop("Indicator name(s): ", paste0("`", tmp, "`", collapse = ", "), 
         " provided to `.dominant_indicators`", 
         ifelse(length(tmp) == 1, " is", " are"), " unknown.", call. = FALSE)
  }
  
  for(i in names(.dominant_indicators)) {
    .W[i, ] = .W[i, ] * sign(.W[i, .dominant_indicators[i]])
  }
  return(.W)
}