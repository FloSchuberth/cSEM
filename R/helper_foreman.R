<<<<<<< HEAD
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
    } else if(.modes[j] == "modeB") {
      w_j <- L[j, ] %>%
        .[. != 0] %>%
        as.matrix(.)
    } else {
      w_j <- .W[j, ] %>% 
        as.matrix(.)
    }
    
    ## Check if single indicator block or composite; If yes, set cf to 1
    if(!(.modes[j] %in% c("modeA", "modeB")) |
       nrow(w_j) == 1 | .csem_model$construct_type[j] == "Composite") {
      correction_factors[j] <- 1
    } else {
      ## Extract relevant objects
      E_jj <- .csem_model$error_cor[rownames(w_j), rownames(w_j)]
      S_jj <- .S[rownames(w_j), rownames(w_j)]
      W_jj <- w_j %*% t(w_j)
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of S_jj (result is vectorized)
      S_vect <- replace(S_jj, which(E_jj == 1), NA) %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[!is.na(.)]
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of W_jj (vectorized)
      W_vect <- replace(W_jj, which(E_jj == 1), NA)  %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[!is.na(.)]
      
      if(length(S_vect) == 0) {
        stop2(
          "The following error occured while calculating the correction factor:\n",
          "At least one pair of indicators with uncorrelated measurement errors",
          " required in each measurement equation.\n", 
          "Measurement equation: `", j, "` has none.")
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
#' Calculate the sample variance-covariance (VCV) matrix of the composites/proxies.
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
#' Calculate the variance-covariance matrix (VCV) of the constructs, i.e., correlations 
#' that involve common factors/latent variables are diattenuated.
#'
#' @usage calculateConstructVCV(
#'  .C          = args_default()$.C, 
#'  .Q          = args_default()$.Q
#'  )
#'  
#' @inheritParams csem_arguments
#'
#' @return The (J x J) construct VCV matrix. Disattenuated if requested.
#' @keywords internal

calculateConstructVCV <- function(
  .C          = args_default()$.C, 
  .Q          = args_default()$.Q
  ) {

  f <- function(.i, .j) { .C[.i, .j] / (.Q[.i] * .Q[.j]) }
  m <- rownames(.C)
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
#' If `.approach_cor_robust = "none"` (the default) the type of correlation computed
#' depends on the types of the columns of `.X_cleaned` (i.e., the indicators) 
#' involved in the computation. 
#' \describe{
#'   \item{`Numeric-numeric`}{If both columns (indicators) involved are numeric, the
#'      Bravais-Pearson product-moment correlation is computed (via [stats::cor()][stats::cor()]).}
#'   \item{`Numeric-factor`}{If any of the columns is a factor variable, the 
#'     polyserial correlation \insertCite{Drasgow1988}{cSEM} is computed (via 
#'     [polycor::polyserial()][polycor::polyserial()]).}
#'   \item{`Factor-factor`}{If both columns are factor variables, the 
#'     polychoric correlation \insertCite{Drasgow1988}{cSEM} is computed (via 
#'     [polycor::polychor()][polycor::polychor()]).}
#' }
#' Note: logical input is treated as a 0-1 factor variable.
#' 
#' If  `"mcd"` (= minimum covariance determinant), the MCD estimator 
#' \insertCite{Rousseeuw1999}{cSEM}, a robust covariance estimator, is applied
#' (via [MASS::cov.rob()][MASS::cov.rob()]).
#' 
#' If `"spearman"`, the Spearman rank correlation is used (via [stats::cor()][stats::cor()]).
#'
#' @usage calculateIndicatorCor(
#'   .X_cleaned           = NULL, 
#'   .approach_cor_robust = "none"
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @return A list with elements:
#' \describe{
#'   \item{`$S`}{The (K x K) indicator correlation matrix}
#'   \item{`$cor_type`}{The type(s) of indicator correlation computed ( 
#'    "Pearson", "Polyserial", "Polychoric")}
#'    \item{`$thre_est`}{Currently ignored (NULL)}
#' }
#' @keywords internal

calculateIndicatorCor <- function(
  .X_cleaned           = NULL,
  .approach_cor_robust = "none"
){
  
  is_numeric_indicator <- lapply(.X_cleaned, is.numeric)
  
  only_numeric_cols <- all(unlist(is_numeric_indicator))
  
  if(.approach_cor_robust != "none" && !only_numeric_cols) {
    stop2("Setting `.approach_cor_robust = ", .approach_cor_robust, "` requires all",
          " columns of .data to be numeric.")
  }
  
  ## polycor::hetcor() is relatively slow. If all columns are numeric use cor
  ## directly
  switch (.approach_cor_robust,
          "none" = {
            if(only_numeric_cols) {
              S <- cor(.X_cleaned)
              cor_type <- "Pearson" 
              thres_est = NULL
            } else {
              
              # Indicator's correlation matrix
              S <- matrix(0, ncol = ncol(.X_cleaned), nrow = ncol(.X_cleaned),
                             dimnames = list(colnames(.X_cleaned), colnames(.X_cleaned)))
              # matrix containing the type of correlation 
              cor_type <- S
              
              # list for the thresholds
              thres_est <- NULL
              
              # temp is used to only calculate the correlations between two 
              # indicators once (upper triangular matrix)
              temp <- colnames(.X_cleaned)
              for(i in colnames(.X_cleaned)){
                temp <- temp[temp!=i]
                for(j in temp){
                  # If both indicators are not continous, the polychoric 
                  # correlation is calculated
                  if (is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == FALSE){
                    # The polycor package gives a list with the polychoric correlation and
                    # the thresholds estimates
                    cor_temp <- polycor::polychor(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$row.cuts
                    thres_est[[j]] <- cor_temp$col.cuts
                    
                    # If one indicator is continous, the polyserial correlation 
                    # is calculated.Note: polyserial needs the continous 
                    # indicator as the first argument.
                  }else if(is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == TRUE){
                    # The polycor package gives the polyserial correlation and the thresholds
                    cor_temp <- polycor::polyserial(.X_cleaned[,j], .X_cleaned[,i], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$cuts
                    thres_est[[j]] <- NA
                  }else if(is_numeric_indicator[[i]] == TRUE & is_numeric_indicator[[j]] == FALSE){
                    cor_temp <- polycor::polyserial(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[j]] <- cor_temp$cuts
                    thres_est[[i]] <- NA
                    
                    # If both indicators are continous, the Pearson correlation
                    # is calculated.
                  }else{
                    S[i,j] <- cor(.X_cleaned[,i], .X_cleaned[,j])
                    cor_type[i,j] <- "Pearson"
                    thres_est[[i]] <- NA
                    thres_est[[j]] <- NA
                  }
                }
              }
              S <- S + t(S)
              diag(S) <- 1
              cor_type <- unique(c(cor_type))
              cor_type <- cor_type[cor_type != "0"]
              
            
              
              # The lavCor function does no smoothing in case of empty cells, which creates problems during bootstrap
              # # Use lavCor function from the lavaan package for the calculation of the polychoric and polyserial correlation 
              # # No smoothing is conducted to ensure positive definiteness of the correlation matrix
              # S <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "cor")
              # 
              # # Estimate thresholds
              # thres_est <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "th")
              # 
              # # Define type of correlation, that can either be polyserial or polychoric
              # type_var=unlist(sapply(.X_cleaned,class))
              # 
              # # if at least one numeric variable is included the polyserial correlation is applied
              # if('numeric' %in% type_var){
              #   cor_type = "Polyserial"
              # } else { #only if all variables are categorical the type of correlation is set to polychoric
              #   cor_type = "Polychoric"
              # }

            }
          },

          "mcd" = {
            S <- MASS::cov.rob(.X_cleaned, cor = TRUE, method = "mcd")$cor
            S[upper.tri(S) == TRUE] = t(S)[upper.tri(S) == TRUE]

            cor_type <-  "Robust (MCD)"
            
            thres_est = NULL
          },
          "spearman" = {
            S <- cor(.X_cleaned, method = "spearman")

            cor_type <-  "Robust (Spearman)"
            
            thres_est = NULL
          }
  )
  # (TODO) not sure how to name the "type" yet and what to do with it. Theoretically,
  # a polycoric correlation could also be used with GSCA or some other non-PLS-PM method.
  list(S = S, cor_type = cor_type, thres_est = thres_est)
}

#' Internal: Calculate Reliabilities
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
        ## Currently only done for Mode A and Mode B. For unit, modeBNNLS, PCA, it is not clear how the Qs are calculated. 
        for(j in names_cf) {
          
          if(modes[j]=='modeA'){
          Lambda[j, ] <- correction_factors[j] * W[j, ]
          }
          
          if(modes[j]=='modeB'){
            Lambda[j, ] <- correction_factors[j] * Lambda[j, ]
          }
          
          Q[j]        <- c(W[j, ] %*% Lambda[j, ])
        }
      } 
    } else if(.approach_weights == "GSCA") {
      
      if(.disattenuate & all(.csem_model$construct_type == "Common factor")) {
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
    
      #  Note: "bartlett" and "regression" weights are obtained AFTER running a CFA. 
      #  Therefore loadings are consistent and as such "disattenuation"
      #  has already implicitly happend. Hence, it is impossible to 
      #  specify "bartlett" or "regression" when .disattenuate = FALSE.
      if(!.disattenuate & .approach_weights %in% c("bartlett", "regression")) {
        stop2(
          "The following error occured in the `calculateReliabilities()` function:\n",
          "Unable to obtain `bartlett` or `regression` weights when `.disattenuate = FALSE`.")
      }
      
      if(any(.csem_model$construct_type == "Composite") & .approach_weights %in% c("bartlett", "regression")) {
        stop2(
          "The following error occured in the `calculateReliabilities()` function:\n",
          "Unable to obtain `bartlett` or `regression` weights ",
          " for models containing constructs modeled as composites.")
      }
      
      # Note: Only necessary if at least one common factor is in the model
      if(.disattenuate & any(.csem_model$construct_type == "Common factor")) {
        ## "Croon" approach
        # The Croon reliabilities assume that the scores are built by sum scores 
        # (i.e. unit weights).
        # Moreover, all constructs must be modeled as common factors.
        
        if(any(.csem_model$construct_type == "Composite")) {
          stop2(
            "The following error occured in the `calculateReliabilities()` function:\n",
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
            
            
            # Standardization does not affect the path coefficients, only the weights and the corresponding 
            # proxy scores are affected. In the manual we should refer to standardized regression weights
            # We decided to use standardized weights as it might be problematic for other function if the scores are not standardized.
            w <- w /  c(sqrt(w %*% .S[indicator_names, indicator_names, drop = FALSE] %*% w))
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
    
    
    # If reliabilities are given by the user, a correction for attenuation is always conducted. 
    # Consequently, in case of the two/three-stage approach for models containing second-order
    # constructs always a correction for attenuation is conducted in the second stage 
    # because in the first stage we calculate the reliabilities and therefore in the 
    # second stage reliabilities are given which leads to correction for attenuation. 
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
  
  # Return Loadings, Reliabilities, and Cross loadings
  # Additionally, the disattenuate argument is returned because it can change in
  # the calculate reliabilities function.
  out <- list("Lambda" = Lambda, "Q2" = Q^2, "W" = W,".disattenuate" = .disattenuate)
  return(out)
}

#' Internal: Set the dominant indicator
#' 
#' Set the dominant indicator for each construct. Since the sign of the weights, 
#' and thus the loadings is often not determined, a dominant indicator can be chosen
#' per block. The sign of the weights are chosen that the correlation between the 
#' dominant indicator and the composite is positive. 
#'
#' @usage setDominantIndicator(
#'  .W                   = args_default()$.W,
#'  .dominant_indicators = args_default()$.dominant_indicators, 
#'  .S                   = args_default()$.S
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @return The (J x K) matrix of weights with the dominant indicator set.
#' @keywords internal
#'
setDominantIndicator <- function(
  .W                   = args_default()$.W,
  .dominant_indicators = args_default()$.dominant_indicators,
  .S                   = args_default()$.S  
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
  
  # Calculate the loadings
  L = .W%*%.S[colnames(.W),colnames(.W)] * abs(sign(.W))
  
  for(i in names(.dominant_indicators)) {
    # ensure that the dominant indicator of a block has positive correlation/loading with the composite
    .W[i, ] = .W[i, ] * sign(L[i, .dominant_indicators[i]])
    
    # Old version which ensures that the weight for this indicator has a positive sign.
    # .W[i, ] = .W[i, ] * sign(.W[i, .dominant_indicators[i]])
  }
  return(.W)
=======
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
    } else if(.modes[j] == "modeB") {
      w_j <- L[j, ] %>%
        .[. != 0] %>%
        as.matrix(.)
    } else {
      w_j <- .W[j, ] %>% 
        as.matrix(.)
    }
    
    ## Check if single indicator block or composite; If yes, set cf to 1
    if(!(.modes[j] %in% c("modeA", "modeB")) |
       nrow(w_j) == 1 | .csem_model$construct_type[j] == "Composite") {
      correction_factors[j] <- 1
    } else {
      ## Extract relevant objects
      E_jj <- .csem_model$error_cor[rownames(w_j), rownames(w_j)]
      S_jj <- .S[rownames(w_j), rownames(w_j)]
      W_jj <- w_j %*% t(w_j)
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of S_jj (result is vectorized)
      S_vect <- replace(S_jj, which(E_jj == 1), NA) %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[!is.na(.)]
      
      ## Set indicator pairs whose measurement errors are correlated to zero and
      ## extract non-zero off-diagonal elements of W_jj (vectorized)
      W_vect <- replace(W_jj, which(E_jj == 1), NA)  %>%
        .[lower.tri(.) | upper.tri(.)] %>%
        .[!is.na(.)]
      
      if(length(S_vect) == 0) {
        stop2(
          "The following error occured while calculating the correction factor:\n",
          "At least one pair of indicators with uncorrelated measurement errors",
          " required in each measurment equation.\n", 
          "Measurement equation: `", j, "` has none.")
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
#' Calculate the sample variance-covariance (VCV) matrix of the composites/proxies.
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
#' Calculate the variance-covariance matrix (VCV) of the constructs, i.e., correlations 
#' that involve common factors/latent variables are diattenuated.
#'
#' @usage calculateConstructVCV(
#'  .C          = args_default()$.C, 
#'  .Q          = args_default()$.Q
#'  )
#'  
#' @inheritParams csem_arguments
#'
#' @return The (J x J) construct VCV matrix. Disattenuated if requested.
#' @keywords internal

calculateConstructVCV <- function(
  .C          = args_default()$.C, 
  .Q          = args_default()$.Q
  ) {

  f <- function(.i, .j) { .C[.i, .j] / (.Q[.i] * .Q[.j]) }
  m <- rownames(.C)
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
#' If `.approach_cor_robust = "none"` (the default) the type of correlation computed
#' depends on the types of the columns of `.X_cleaned` (i.e., the indicators) 
#' involved in the computation. 
#' \describe{
#'   \item{`Numeric-numeric`}{If both columns (indicators) involved are numeric, the
#'      Bravais-Pearson product-moment correlation is computed (via [stats::cor()][stats::cor()]).}
#'   \item{`Numeric-factor`}{If any of the columns is a factor variable, the 
#'     polyserial correlation \insertCite{Drasgow1988}{cSEM} is computed (via 
#'     [polycor::polyserial()][polycor::polyserial()]).}
#'   \item{`Factor-factor`}{If both columns are factor variables, the 
#'     polychoric correlation \insertCite{Drasgow1988}{cSEM} is computed (via 
#'     [polycor::polychor()][polycor::polychor()]).}
#' }
#' Note: logical input is treated as a 0-1 factor variable.
#' 
#' If  `"mcd"` (= minimum covariance determinant), the MCD estimator 
#' \insertCite{Rousseeuw1999}{cSEM}, a robust covariance estimator, is applied
#' (via [MASS::cov.rob()][MASS::cov.rob()]).
#' 
#' If `"spearman"`, the Spearman rank correlation is used (via [stats::cor()][stats::cor()]).
#'
#' @usage calculateIndicatorCor(
#'   .X_cleaned           = NULL, 
#'   .approach_cor_robust = "none"
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @return A list with elements:
#' \describe{
#'   \item{`$S`}{The (K x K) indicator correlation matrix}
#'   \item{`$cor_type`}{The type(s) of indicator correlation computed ( 
#'    "Pearson", "Polyserial", "Polychoric")}
#'    \item{`$thre_est`}{Currently ignored (NULL)}
#' }
#' @keywords internal

calculateIndicatorCor <- function(
  .X_cleaned           = NULL,
  .approach_cor_robust = "none"
){
  
  is_numeric_indicator <- lapply(.X_cleaned, is.numeric)
  
  only_numeric_cols <- all(unlist(is_numeric_indicator))
  
  if(.approach_cor_robust != "none" && !only_numeric_cols) {
    stop2("Setting `.approach_cor_robust = ", .approach_cor_robust, "` requires all",
          " columns of .data to be numeric.")
  }
  
  ## polycor::hetcor() is relatively slow. If all columns are numeric use cor
  ## directly
  switch (.approach_cor_robust,
          "none" = {
            if(only_numeric_cols) {
              S <- cor(.X_cleaned)
              cor_type <- "Pearson" 
              thres_est = NULL
            } else {
              
              # Indicator's correlation matrix
              S <- matrix(0, ncol = ncol(.X_cleaned), nrow = ncol(.X_cleaned),
                             dimnames = list(colnames(.X_cleaned), colnames(.X_cleaned)))
              # matrix containing the type of correlation 
              cor_type <- S
              
              # list for the thresholds
              thres_est <- NULL
              
              # temp is used to only calculate the correlations between two 
              # indicators once (upper triangular matrix)
              temp <- colnames(.X_cleaned)
              for(i in colnames(.X_cleaned)){
                temp <- temp[temp!=i]
                for(j in temp){
                  # If both indicators are not continous, the polychoric 
                  # correlation is calculated
                  if (is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == FALSE){
                    # The polycor package gives a list with the polychoric correlation and
                    # the thresholds estimates
                    cor_temp <- polycor::polychor(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$row.cuts
                    thres_est[[j]] <- cor_temp$col.cuts
                    
                    # If one indicator is continous, the polyserial correlation 
                    # is calculated.Note: polyserial needs the continous 
                    # indicator as the first argument.
                  }else if(is_numeric_indicator[[i]] == FALSE & is_numeric_indicator[[j]] == TRUE){
                    # The polycor package gives the polyserial correlation and the thresholds
                    cor_temp <- polycor::polyserial(.X_cleaned[,j], .X_cleaned[,i], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[i]] <- cor_temp$cuts
                    thres_est[[j]] <- NA
                  }else if(is_numeric_indicator[[i]] == TRUE & is_numeric_indicator[[j]] == FALSE){
                    cor_temp <- polycor::polyserial(.X_cleaned[,i], .X_cleaned[,j], thresholds = TRUE)
                    S[i,j] <- cor_temp$rho
                    cor_type[i,j] <- cor_temp$type
                    thres_est[[j]] <- cor_temp$cuts
                    thres_est[[i]] <- NA
                    
                    # If both indicators are continous, the Pearson correlation
                    # is calculated.
                  }else{
                    S[i,j] <- cor(.X_cleaned[,i], .X_cleaned[,j])
                    cor_type[i,j] <- "Pearson"
                    thres_est[[i]] <- NA
                    thres_est[[j]] <- NA
                  }
                }
              }
              S <- S + t(S)
              diag(S) <- 1
              cor_type <- unique(c(cor_type))
              cor_type <- cor_type[cor_type != "0"]
              
            
              
              # The lavCor function does no smoothing in case of empty cells, which creates problems during bootstrap
              # # Use lavCor function from the lavaan package for the calculation of the polychoric and polyserial correlation 
              # # No smoothing is conducted to ensure positive definiteness of the correlation matrix
              # S <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "cor")
              # 
              # # Estimate thresholds
              # thres_est <- lavaan::lavCor(.X_cleaned, se = 'none', estimator = "two.step", output = "th")
              # 
              # # Define type of correlation, that can either be polyserial or polychoric
              # type_var=unlist(sapply(.X_cleaned,class))
              # 
              # # if at least one numeric variable is included the polyserial correlation is applied
              # if('numeric' %in% type_var){
              #   cor_type = "Polyserial"
              # } else { #only if all variables are categorical the type of correlation is set to polychoric
              #   cor_type = "Polychoric"
              # }

            }
          },

          "mcd" = {
            S <- MASS::cov.rob(.X_cleaned, cor = TRUE, method = "mcd")$cor
            S[upper.tri(S) == TRUE] = t(S)[upper.tri(S) == TRUE]

            cor_type <-  "Robust (MCD)"
            
            thres_est = NULL
          },
          "spearman" = {
            S <- cor(.X_cleaned, method = "spearman")

            cor_type <-  "Robust (Spearman)"
            
            thres_est = NULL
          }
  )
  # (TODO) not sure how to name the "type" yet and what to do with it. Theoretically,
  # a polycoric correlation could also be used with GSCA or some other non-PLS-PM method.
  list(S = S, cor_type = cor_type, thres_est = thres_est)
}

#' Internal: Calculate Reliabilities
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
        ## Currently only done for Mode A and Mode B. For unit, modeBNNLS, PCA, it is not clear how the Qs are calculated. 
        for(j in names_cf) {
          
          if(modes[j]=='modeA'){
          Lambda[j, ] <- correction_factors[j] * W[j, ]
          }
          
          if(modes[j]=='modeB'){
            Lambda[j, ] <- correction_factors[j] * Lambda[j, ]
          }
          
          Q[j]        <- c(W[j, ] %*% Lambda[j, ])
        }
      } 
    } else if (.approach_weights == "GSCA") {
      
        if ((.disattenuate &
            all(.csem_model$construct_type == "Common factor")) | any(.csem_model$construct_type == "Common factor")) {
          
          # Currently, GSCAm only supports pure common factor models. This may change
          # in the future.
          
          # Compute consistent loadings and Q (Composite/proxy-construct correlation)
          # Consistent factor loadings are obtained from GSCAm.
          
          for (j in rownames(Lambda)) {
            Lambda[j, ] <- .W$C[j, ]
            Q[j]        <- c(W[j, ] %*% Lambda[j, ])
          }
        }
      
    } else if(.approach_weights %in% c("unit", "bartlett", "regression", "PCA", 
                "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR")) {
    
      #  Note: "bartlett" and "regression" weights are obtained AFTER running a CFA. 
      #  Therefore loadings are consistent and as such "disattenuation"
      #  has already implicitly happend. Hence, it is impossible to 
      #  specify "bartlett" or "regression" when .disattenuate = FALSE.
      if(!.disattenuate & .approach_weights %in% c("bartlett", "regression")) {
        stop2(
          "The following error occured in the `calculateReliabilities()` function:\n",
          "Unable to obtain `bartlett` or `regression` weights when `.disattenuate = FALSE`.")
      }
      
      if(any(.csem_model$construct_type == "Composite") & .approach_weights %in% c("bartlett", "regression")) {
        stop2(
          "The following error occured in the `calculateReliabilities()` function:\n",
          "Unable to obtain `bartlett` or `regression` weights ",
          " for models containing constructs modeled as composites.")
      }
      
      # Note: Only necessary if at least one common factor is in the model
      if(.disattenuate & any(.csem_model$construct_type == "Common factor")) {
        ## "Croon" approach
        # The Croon reliabilities assume that the scores are built by sum scores 
        # (i.e. unit weights).
        # Moreover, all constructs must be modeled as common factors.
        
        if(any(.csem_model$construct_type == "Composite")) {
          stop2(
            "The following error occured in the `calculateReliabilities()` function:\n",
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
            
            
            # Standardization does not affect the path coefficients, only the weights and the corresponding 
            # proxy scores are affected. In the manual we should refer to standardized regression weights
            # We decided to use standardized weights as it might be problematic for other function if the scores are not standardized.
            w <- w /  c(sqrt(w %*% .S[indicator_names, indicator_names, drop = FALSE] %*% w))
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
    
    
    # If reliabilities are given by the user, a correction for attenuation is always conducted. 
    # Consequently, in case of the two/three-stage approach for models containing second-order
    # constructs always a correction for attenuation is conducted in the second stage 
    # because in the first stage we calculate the reliabilities and therefore in the 
    # second stage reliabilities are given which leads to correction for attenuation. 
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
  
  # Return Loadings, Reliabilities, and Cross loadings
  # Additionally, the disattenuate argument is returned because it can change in
  # the calculate reliabilities function.
  out <- list("Lambda" = Lambda, "Q2" = Q^2, "W" = W,".disattenuate" = .disattenuate)
  return(out)
}

#' Internal: Set the dominant indicator
#' 
#' Set the dominant indicator for each construct. Since the sign of the weights, 
#' and thus the loadings is often not determined, a dominant indicator can be chosen
#' per block. The sign of the weights are chosen that the correlation between the 
#' dominant indicator and the composite is positive. 
#'
#' @usage setDominantIndicator(
#'  .W                   = args_default()$.W,
#'  .dominant_indicators = args_default()$.dominant_indicators, 
#'  .S                   = args_default()$.S
#'  )
#'
#' @inheritParams csem_arguments
#' 
#' @return The (J x K) matrix of weights with the dominant indicator set.
#' @keywords internal
#'
setDominantIndicator <- function(
  .W                   = args_default()$.W,
  .dominant_indicators = args_default()$.dominant_indicators,
  .S                   = args_default()$.S  
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
  
  # Calculate the loadings
  L = .W%*%.S[colnames(.W),colnames(.W)] * abs(sign(.W))
  
  for(i in names(.dominant_indicators)) {
    # ensure that the dominant indicator of a block has positive correlation/loading with the composite
    .W[i, ] = .W[i, ] * sign(L[i, .dominant_indicators[i]])
    
    # Old version which ensures that the weight for this indicator has a positive sign.
    # .W[i, ] = .W[i, ] * sign(.W[i, .dominant_indicators[i]])
  }
  return(.W)
>>>>>>> 503f2dec (Pre Kronecker work)
}