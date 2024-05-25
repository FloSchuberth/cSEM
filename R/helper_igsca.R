#' R Implementation of igsca_sim.m
#' 
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#' 
#' In the example section, the specified model is based on the tutorial I-GSCA model associated with GSCA Pro \insertCite{hwangetal2023StructuralEquationModelingAMultidisciplinaryJournal}{cSEM}.
#' 
#' @param Z0 Data matrix of N cases (measurements) x J indicators with named
#'   columns, unstandardized.
#' @param W0 Indicator matrix of weights: J indicators (rows) and their
#'   corresponding Gamma construct variables (columns).
#' @param C0 Indicator matrix of loadings: J indicators (rows) and their
#'   corresponding Gamma construct variable (columns).
#' @param B0 Square indicator matrix of path coefficients:
#'   from-construct-variable (rows) and to-construct-variable (columns). The
#'   order of Gamma construct variables should match the order in C0 and W0.
#' @param con_type A vector that denotes whether each construct variable
#'   (columns in W0 and C0) is a common factor or composite. Its length should
#'   be equal to the number of columns of W0 and C0.
#' @param indicator_type An indicator vector that indices whether a j indicator
#'   (rows of W0 and C0) corresponds to a common factor variable (1) or a
#'   composite variable (0). This vector is important for computing the
#'   uniqueness terms (D) because it zeros the entries for composite indicators.
#' @param .dominant_indicators A named vector that indices the dominant
#'   indicator for each construct variable. *It is to be clarified whether this
#'   should only apply to factor latent variables or also composite latent
#'   variables.* This is important for ensuring that the signs of the
#'   path-coefficients and loadings are consistent. It is sometimes the case in
#'   composite-based structural equation modelling methods that
#'   loadings/path-coefficients may have the opposite sign. The length of this
#'   vector should be equal to the number of construct variables and each value
#'   should represent the row number of the dominant indicator for that
#'   construct variable.
#' @param .iter_max Maximum number of iterations of the Alternating Least Squares
#'   (ALS) algorithm.
#' @param .tolerance Minimum amount of absolute change in the estimates of the
#'   path-coefficients (if B0 is non-zero) or the loadings (if B0 is all zero,
#'   meaning there are no path-coefficients) between ALS iterations before ending
#'   the optimization.
#' @inheritParams csem
#' 
#' @author Michael S. Truong
#' @return List of 4 matrices that make up a fitted I-GSCA Model: 
#' * (1) Weights 
#' * (2) Loadings 
#' * (3) Uniqueness Terms D^2 
#' * (4) Path Coefficients.
#'  
#' @importFrom MASS ginv 
#' 
#' @references
#'   \insertAllCited{}
#' @examples
#' 
#' # Specify the model according to GSCA Pro's example
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' data(LeDang2022)
#' 
#' csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "IGSCA",
#' .dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion = "sum_diff_absolute")
igsca <-
  function(Z0, W0, C0, B0, con_type, indicator_type, .dominant_indicators, .iter_max = 100, .tolerance = 0.0001, .conv_criterion) {
  
## Initialize Computational Variables -----------------------------------------------------
  n_case <- nrow(Z0)
  n_indicators <- ncol(Z0)
  n_constructs <- ncol(W0)
  n_total_var <- n_indicators + n_constructs
  
  w_index <- which(c(W0) == 1) 
  c_index <- which(c(C0) == 1)
  b_index <- which(c(B0) == 1)
  

### Initial Estimates and Preparation -------------------------------------
  prepared_for_ALS <- prepare_for_ALS(
    Z0 = Z0,
    W0 = W0,
    B0 = B0,
    n_indicators = n_indicators,
    n_case = n_case,
    n_constructs = n_constructs,
    indicator_type = indicator_type
  )
  
  list2env(prepared_for_ALS[c("W", "C", "B", "V", "Z", "D", "U", "Gamma")],
           envir = environment())
  
  
## Alternating Least Squares Algorithm -------------------------------------

### Optimization Preparation -----------------------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (length(b_index) > 0) {
      est <- B[b_index] 
    } else {
      est <- C[c_index]
    }    
    est0 <- est + 1 
    it <- 0
    
### Optimization Loop: Alternating Between Weights and Loadings ---------------

    while ((!checkConvergence(.W_new = est, .W_old = est0, .tolerance = .tolerance, .conv_criterion = .conv_criterion)) && (it <= .iter_max)) {

      # Update Counter Variables
      it <- it + 1
      est0 <- est
  
#### Compute Pseudo-Weights --------------------------------------------------
      updated_X_WW_pseudo_weights <- update_X_WW_pseudo_weights(
        Z = Z,
        U = U,
        D = D,
        C = C,
        n_constructs = n_constructs,
        B = B
      )
      # Creates X (for updating Composites) and WW (for updating Factors)
      list2env(updated_X_WW_pseudo_weights[c("X", "WW")], envir = environment())

#### Sequential Updating of Constructs ---------------------------------------

      A <- cbind(C, B) 
      
      # After each cycle, the Gamma, W and V matrices are updated
      for (gamma_idx in seq_len(n_constructs)) {
      
        tot <- n_indicators + gamma_idx
        windex_gamma_idx <- (W0[, gamma_idx] == 1)
        X_gamma_idx <- X[, windex_gamma_idx]
        
        if (con_type[gamma_idx] == "Composite") {
          
          theta <-
            update_composite(
              n_total_var = n_total_var,
              tot = tot, # Changes per gamma_idx iteration
              n_constructs = n_constructs,
              gamma_idx = gamma_idx, # Changes per gamma_idx iteration
              W = W, # Changes per gamma_idx iteration
              A = A,
              V = V, # Changes per gamma_idx iteration
              X = X,
              windex_gamma_idx = windex_gamma_idx # Changes per gamma_idx iteration
            )
          
        } else if (con_type[gamma_idx] == "Common factor") {
          
          theta <-
             update_common_factor(
               WW = WW, 
               windex_gamma_idx = windex_gamma_idx, # Changes per gamma_idx iteration
               gamma_idx = gamma_idx # Changes per gamma_idx iteration
               )
          
        } else {
          
          stop("con_type should only either be `Composite` or `Common factor`")
          
        }
        
        # This is where Gamma, Weights and V are updated based on which gamma_idx
        theta <- theta / norm(X_gamma_idx %*% theta, "2")
        Gamma[, gamma_idx] <- X_gamma_idx %*% theta
        W[windex_gamma_idx, gamma_idx] <- theta
        V[windex_gamma_idx, tot] <- theta
    }
      
#### Update Loadings, Path Coefficients and Uniqueness Terms ----------
      updated_C_B_D_U <- update_C_B_D_U(
        X = X,
        n_indicators = n_indicators,
        Gamma = Gamma,
        c_index = c_index,
        C = C,
        n_constructs = n_constructs,
        b_index = b_index,
        B = B,
        n_case = n_case,
        D = D,
        Z = Z,
        indicator_type = indicator_type
      )
      
      list2env(updated_C_B_D_U[c("D", "U", "C", "B")], envir = environment())

#### Update estimates for while-loop condition testing -----------------------
      if (length(b_index) > 0) {
        est <- B[b_index]
      } else {
        est <- C[c_index]
      }
    }

## Flip Signs for Factors and Composites Based on Dominant Indicators --------
    if (!is.null(.dominant_indicators)) {
      flipped_signs <-
        flip_signs_ind_domi(
          n_constructs = n_constructs,
          Z = Z,
          .dominant_indicators = .dominant_indicators,
          Gamma = Gamma,
          C = C,
          B = B
        )
      
      list2env(flipped_signs[c("Gamma", "C", "B")], envir = environment())
    }

## Output Formatting -------------------------------------------------------
  Weights <- W
  rownames(Weights) <- rownames(W0) 
  
  Loadings <- t(C)
  rownames(Loadings) <- rownames(C0)
  
  uniqueD <- diag(D) ^ 2
  names(uniqueD) <- rownames(Loadings) 
  
  PathCoefficients <- B
  
  
    return(
      list(
        "Weights" = Weights,
        "Loadings" =  Loadings,
        "Path Coefficients" =  PathCoefficients,
        "Uniqueness Terms" = uniqueD,
        "Construct Scores" = Gamma,
        "Iterations" = it
      )
    )

}

#' R Implementation of gsca_inione.m from Heungsun Hwang
#' 
#' Internal I-GSCA function that is a slightly modified implementation of ordinary Generalised Structured Component Analysis (GSCA).
#' 
#' Initializes the values for I-GSCA
#' @inheritParams igsca
#' 
#' @author Michael S. Truong
#' @returns Returns a list of starting values for:
#' * Weights (W)
#' * Loadings (C)
#' * Path Coefficients (B)  
#'
gsca_inione <- function(Z0, W0, B0) {
  
  N <- nrow(Z0)
  J <- nrow(W0)
  P <- ncol(W0)
  TRep <- J + P
  C0 <- t(W0)
  A0 <- cbind(C0, B0)
  
  w_index <- which(W0 != 0)
  aindex <- which(A0 != 0)
  
  Z <- scale(Z0, center = TRUE, scale = TRUE) * sqrt(N) / sqrt(N - 1)
  # Random Values to W and A
  W <- W0
  A <- A0
  
  W[w_index] <- rep(1, length(w_index))
  A[aindex] <- runif(length(aindex), min = 0, max = 1)
  
  Gamma <- Z %*% W
  W <-
    W / (t(sqrt(diag(t(
      Gamma
    ) %*% Gamma))) |> rep(each = nrow(W)))
  V <- cbind(diag(J), W)
  Gamma <- Z %*% W
  
  Psi <- Z %*% V
  vecPsi <- c(Psi)
  it <- 0
  f0 <- 100000000
  imp <- 100000000
  while ((it <= 300) && (imp > 0.0001)) {
    it <- it + 1
    Phi = kronecker(diag(TRep), Gamma)
    Phi <- Phi[, aindex]
    A[aindex] <- solve(t(Phi) %*% Phi, t(Phi)) %*% vecPsi
    
    for (p in seq_len(P)) {
      t_lil <- J + p
      windex_p <- which(W0[, p] != 0)
      m <- matrix(0, nrow = 1, ncol = TRep)
      m[t_lil] <- 1
      a <- A[p,]
      beta <- m - a
      H1 <- diag(P)
      H2 <- diag(TRep)
      H1[p, p] <- 0
      H2[t_lil, t_lil] <- 0
      Delta <- (W %*% H1 %*% A) - (V %*% H2)
      vecZDelta <- c(Z %*% Delta)
      
      XI <- kronecker(t(beta), Z)
      XI <- XI[, windex_p]
      theta <- MASS::ginv(t(XI) %*% XI) %*% t(XI) %*% vecZDelta
      zw <- Z[, windex_p] %*% theta
      
      
      theta <- sqrt(N) * theta / norm(zw, "2")
      W[windex_p, p] <- theta
      V[windex_p, t_lil] <- theta
    }
    Gamma <- Z %*% W
    Psi <- Z %*% V
    dif <- Psi - Gamma %*% A
    f <- sum(diag(t(dif) %*% dif))
    imp <- f0 - f
    f0 <- f
    vecPsi <- c(Psi)
  }
  C <- A[, 1:J]
  B <- A[, (J + 1):ncol(A)]
  return(list(
    "W" = W,
    "C" = C,
    "B" = B,
    "it" = it
  ))
}


#' Prepare for ALS Algorithm
#'
#' Internal I-GSCA function
#' 
#' @inheritParams igsca
#' @param n_indicators Number of indicators
#' @param n_case Number of measurements
#' @param n_constructs Number of constructs
#'
#' @returns List of matrices to put through the Alternating Least Squared (ALS) algorithm: 
#' * Estimated Weights matrix (W)
#' * Estimated Loadings matrix (C)
#' * Estimated Path Coefficients matrix (B)
#' * Computational Helper Matrix (V)
#' * Estimated Uniqueness Errors vector (D)
#' * Estimated Related to Uniqueness Errors vector (U)
#' * Estimated Construct Scores matrix (Gamma)
prepare_for_ALS <- function(Z0, W0, B0, n_indicators, n_case, n_constructs, indicator_type) {
  
  # Initial estimates using GSCA
  initial_est <-
    gsca_inione(
      Z0 = Z0,
      W0 = apply(W0 != 0, 2, as.numeric),
      B0 = apply(B0 != 0, 2, as.numeric)
    )
  
  list2env(initial_est[c("W", "C", "B")], envir = environment())
  
  # Initialize matrix to hold V
  V <- cbind(diag(n_indicators), W)
  # Standardize Data Matrix
  Z <- scale(Z0, center = TRUE, scale = TRUE) / sqrt(n_case - 1)
  # Create Initial Construct Scores Matrix
  Gamma <- Z %*% W
  # Initialize Unique Error Matrix
  D <- diag(n_indicators)
  
  # 
  # Solution for Q is copied from estimators_weights.R
  Q <- qr.Q(qr(Gamma), complete =  TRUE)
  F_o <- Q[, (n_constructs + 1):n_case, drop = FALSE]
  
  # svd between R and Matlab by Ahmed Fasih on February 1/2017 https://stackoverflow.com/a/41972818
  svd_out <- (D %*% t(Z) %*% F_o) |>
    {
      \(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))
    }()
  u <- svd_out$u
  v <- svd_out$v
  
  # Utilde deviates from Matlab because of the SVD
  Utilde <- v[, 1:n_indicators] %*% t(u)
  U <- F_o %*% Utilde
  D <- diag(diag(t(U) %*% Z))
  D[indicator_type == "Composite", indicator_type == "Composite"] <- 0
  
  
  return(list(
    "W" = W,
    "C" = C,
    "B" = B,
    "V" = V,
    "Z" = Z,
    "D" = D,
    "U" = U,
    "Gamma" = Gamma
  ))
}

#' Update Pseudo Weights for Composites and Common-Factors
#' 
#' Computation of WW matrix was converted between Matlab to R based on page 44 of \insertCite{Hiebeler2015;textual}{cSEM}.  
#' @param Z Standardized data matrix
#' @param D Matrix of estimated unique error
#' @param U Matrix of estimates related to unique error
#' @param C Matrix of estimated loadings
#' @param n_constructs  Number of constructs
#' @param B Matrix of path coefficients
#'
#' @returns Two matrices:
#' * X: Remaining part of data (Z) after accounting for uniqueness terms (U) and (D), used for estimating composite loadings. Also used for standardizing theta when updating Gamma, W and V
#' * WW: Weights after accounting for current Loading and Path-Coefficients values, used for estimating common-factor loadings
#' 
#' @references
#'   \insertAllCited{}
#' 
#'
update_X_WW_pseudo_weights <- function(Z, U, D, C, n_constructs, B) {
  # X deviates from Matlab because it is an offspring of svd_out
  X <- Z - U %*% D

  WW <-
    t(C) %*% solve((C %*% t(C) + diag(n_constructs) - 2 * B + (B %*% t(B))))
  return(list("X" = X, "WW" = WW))
}


#' Update Common Factor Variable
#'
#' @param WW Pseudo-weights for Common Factors
#' @param windex_gamma_idx Index of weights related to the indicators for the construct of interest
#' @param gamma_idx Index of which construct we are examining
#'
#' @return theta: Used to update factor latent variables -- after accounting for loadings and path-coefficients.
#' 
#' 
#'
update_common_factor <- function(WW, windex_gamma_idx, gamma_idx) {
  theta <- WW[windex_gamma_idx, gamma_idx]
  return(theta)
}

#' Update Composite Variables
#'
#' @param n_total_var Number of indicators and constructs
#' @param tot Index dependent on which construct variable we are examining
#' @param n_constructs Number of constructs
#' @param gamma_idx Index of which construct we are examining
#' @param W Weights matrix
#' @param A Stacked matrix of loadings and path coefficients
#' @param V Unclear meaning
#' @param X Pseudo-weights for composite variables
#' @param windex_gamma_idx Index of weights related to the indicators for the construct of interest
#'
#' @return theta: A matrix that will later be used to update the weights for the composite variable.
#' 
#'
update_composite <-
  function(n_total_var, tot, n_constructs, gamma_idx, W, A, V, X, windex_gamma_idx) {
    
    e <- matrix(0, nrow = 1, ncol = n_total_var)
    e[tot] <- 1
    H1 <- diag(n_total_var)
    H2 <- diag(n_constructs)
    H1[tot, tot] <- 0
    H2[gamma_idx, gamma_idx] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)
    
    
    vecZDelta <- c(X %*% Delta)
    beta <- e - A[gamma_idx, ]
    XI <- kronecker(t(beta), X)
    XI <- XI[, windex_gamma_idx]
    
    theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta
    return(theta)
  }

#' Update Loadings, Path-Coefficients and Uniqueness Terms After Updating Latent Variables
#' 
#' @inheritParams prepare_for_ALS
#' @inheritParams update_X_WW_pseudo_weights
#' @param X Pseudo-weights for composites
#' @param Gamma Construct Scores
#' @param c_index Index of loadings
#' @param b_index Index of Path Coefficients
#' @param n_case Number of Cases
#' @param indicator_type Vector of whether each indicator corresponds to a common factor or composite
#'
#' @return List of matrices:
#'
#' * (1) Uniqueness terms (D) and (U)
#' * (2) Estimated Path Coefficients matrix (B)
#' * (3) Estimated Loadings matrix (C)
#'
update_C_B_D_U <-
  function(X,
           n_indicators,
           Gamma,
           c_index,
           C,
           n_constructs,
           b_index,
           B,
           n_case,
           D,
           Z,
           indicator_type) {
    t1 <- c(X)
    M1 <- kronecker(diag(n_indicators), Gamma)
    M1 <- M1[, c_index]
    C[c_index] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)
    
    t2 <- c(Gamma)
    M2 <- kronecker(diag(n_constructs), Gamma)
    M2 <- M2[, b_index]
    B[b_index] <- MASS::ginv(t(M2) %*% M2) %*% (t(M2) %*% t2)
    
    # Solution for Q is copied from estimators_weights.R
    Q <- qr.Q(qr(Gamma), complete =  TRUE)
    F_o <- Q[, (n_constructs + 1):n_case]
    
    # svd between R and Matlab by Ahmed Fasih on February 1/2017 https://stackoverflow.com/a/41972818
    svd_out2 <- (D %*% t(Z) %*% F_o) |>
      {
        \(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))
      }()
    u <- svd_out2$u
    v <- svd_out2$v
    Utilde <- v[, 1:n_indicators] %*% t(u)
    U <- F_o %*% Utilde
    D <- diag(diag(t(U) %*% Z))
    D[indicator_type == "Composite", indicator_type == "Composite"] <- 0
    
    return(
      list(
        "D" = D,
        "U" = U,
        "B" = B,
        "C" = C
      )
    )
  }

#' Flip signs of Gamma, Loadings and Path-Coefficients Cells Based on Dominant Indicator
#' 
#' @inheritParams igsca
#' @inheritParams update_X_WW_pseudo_weights
#' @inheritParams update_C_B_D_U
#'
#' @return List of matrices:
#' * Estimated Construct Scores (Gamma)
#' * Estimated Loadings matrix (C)
#' * Estimated Path-Coefficients matrix (B)
#'
flip_signs_ind_domi <- function(n_constructs, Z, .dominant_indicators, Gamma, C, B) {
  for (gamma_idx in seq_len(n_constructs)) {
    if (exists(.dominant_indicators[gamma_idx])) {
      if ((t(Z[, .dominant_indicators[gamma_idx]]) %*% Gamma[, gamma_idx]) < 0) {
        Gamma[, gamma_idx] <- (-1 * Gamma[, gamma_idx])
        C[gamma_idx,] <- (-1 * C[gamma_idx,])
        B[gamma_idx,] <- (-1 * B[gamma_idx,])
        B[, gamma_idx] <- (-1 * B[, gamma_idx])
      }
    }
  }
  return(list(
    "Gamma" = Gamma,
    "C" = C,
    "B" = B
  ))
}

#' A parseModel extractor function for the purposes of running I-GSCA code example
#' 
#' In the context of igsca, this function prepares: (1) the initial indicators
#' (Z0), weights (W0), structural (B0), loadings(C0) matrices; (2) whether a
#' construct is a latent or composite variable (con_type); (3) whether an
#' indicator corresponds to a latent or composite variable (indicator_type); and
#' (4) the dominant indicator of each construct (.dominant_indicators).
#' 
#' @inheritParams csem
#' 
#'
#' @return Returns a list of matrices/vectors required for igsca_sim() to run:
#' * Z0
#' * W0
#' * B0
#' * C0
#' * con_type
#' * indicator_type
extract_parseModel <-
  function(.data, .model) {
    # Note: parseModel is from cSEM internal
    
    csemify <- parseModel(.model = .model)
    
    Z0 <- .data[, csemify$indicators]
    
    # B0 <- csemify$structural
    B0 <- t(csemify$structural)
    W0 <- t(csemify$measurement)
    
    # con_type <- csemify$construct_type == "Common factor"
    con_type <- csemify$construct_type
    
    # Constructing indicator_type
    indicator_type <- vector(mode = "character", length = ncol(csemify$measurement))
    indicator_type <- rep("Composite", length(indicator_type)) # Default value must be "Composite" because of how the measurement matrix returned by parseModel uses 0 to denote both composite variables and the absence of any corresponding construct variable.
    names(indicator_type) <- colnames(csemify$measurement)
    
    
    for (gamma_idx in rownames(csemify$measurement)) {
      for (indicator in colnames(csemify$measurement)) {
        if (csemify$construct_type[gamma_idx] == "Common factor") {
          
          if (csemify$measurement[gamma_idx, indicator] == 1) {
            indicator_type[indicator] <- "Common factor"
          }
          
        } else if (csemify$construct_type[gamma_idx] == "Composite") {
          # The reason why the following code breaks behavior is because
          # it makes indicator_type turn into all "Composite because in the parseModel's return
          # the value 0 in the measurement matrix denotes both a composite and
          #  no correspondence to any latent variable
          # 
          # if(csemify$measurement[gamma_idx, indicator] == 0) {
          #   indicator_type[indicator] <- "Composite"
          # }
        } else {
          warning("Indicator does not correspond to either a Composite or Common factor. Unsupported behavior may ensue.")
        }
      }
    }
    
    C0 <- W0
    
    
    stopifnot(
      "Data matrix (Z0) does not correctly correspond with Weights matrix (W0)" = identical(colnames(Z0), rownames(W0)),
      "Weights matrix (W0) does not correctly correspond with Structural matrix (B0)" = identical(colnames(W0), colnames(B0)),
      "Construct indicator does not correctly correspond with Weights Matrix (W0)" = identical(names(csemify$construct_type), colnames(W0)),
      "All indicators for composite and common factor should have loadings in I-GSCA" = identical(W0, C0),
      "Every indicator should only have loadings from one construct" = identical(unique(rowSums(C0)), 1),
      "Every indicator should only have weights to one construct" = identical(unique(rowSums(W0)), 1)
    )
    
    return(
      list(
        "Z0" = Z0,
        "B0" = B0,
        "W0" = W0,
        "con_type" = con_type,
        "indicator_type" = indicator_type,
        "C0" = C0
      )
    )
  }

#' Converts Output of igsca functions into a table to facilitate comparisons
#' 
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#' @inheritParams extract_parseModel
#' @param weights Weights matrix
#' @param loadings Loadings matrix
#' @param uniqueD Vector of Uniqueness for each indicator of a common factor
#' @param paths Path coefficients matrix
#' @importFrom lavaan lavaanify
#' @return Table of Weights, Loadings, Path-Coefficients and Uniqueness terms from i-gsca algorithms in Matlab or R.
#' @export
#'
get_lavaan_table_igsca_matrix <- function(model, weights, loadings, uniqueD, paths) {
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"),]
  # Pre-allocate Columns
  table <-
    cbind(table, list(
      "weights" = 0,
      "loadings" = 0,
      "uniqueD" = 0,
      "paths" = 0
    ))
  
  # Slide in weights
  for (indicator in rownames(weights)) {
    for (lv in colnames(weights)) {
      table[((table$lhs == lv &
              table$rhs == indicator) &
              table$op %in% c("<~", "=~")), "weights"] <- weights[indicator, lv]
    }
  }
  
  # Slide in loadings
  for (indicator in rownames(loadings)) {
    for (lv in colnames(loadings)) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "loadings"] <- loadings[indicator, lv]
    }
  }
  
  # Slide in uniqueD
  for (indicator in names(uniqueD)) {
    # This assumes that every indicator only loads onto one factor
    # Cross-factor loadings will not work with this
      table[((table$rhs == indicator) &
               (table$op == "=~")), "uniqueD"] <- uniqueD[indicator]
  }
  
  # Slide in Paths
  for (lv_from in rownames(paths)) {
    for (lv_to in colnames(paths)) {
      table[((table$rhs == lv_from &
                table$lhs == lv_to) &
               table$op == "~"), "paths"] <- paths[lv_from, lv_to]
    }
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  
  return(table)
  
}

igsca_warning <- function() {
  # TODO: Might want to do something about this during the simulations...
  # TODO: Maybe only call with the print Method?
  warning("IGSCA in cSEM currently only supports...")
}