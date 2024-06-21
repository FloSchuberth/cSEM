#' R Implementation of igsca_sim.m
#' 
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#' 
#' @param Z0 Data matrix of N cases (measurements) x J indicators with named columns, unstandardized.
#' @param W0 Indicator matrix of weights: J indicators (rows) and their corresponding Gamma construct variables (columns).
#' @param C0 Indicator matrix of loadings: J indicators (rows) and their corresponding Gamma construct variable (columns).
#' @param B0 Square indicator matrix of path coefficients: from-construct-variable (rows) and to-construct-variable (columns). The order of Gamma construct variables should match the order in C0 and W0.
#' @param con_type A vector that denotes whether each construct variable (columns in W0 and C0) is a common factor or composite. Its length should be equal to the number of columns of W0 and C0. 
#' @param indicator_type An indicator vector that indices whether a j indicator (rows of W0 and C0) corresponds to a common factor variable (1) or a composite variable (0). This vector is important for computing the uniqueness terms (D) because it zeros the entries for composite indicators. 
#' @param ind_domi A numeric vector that indices the dominant indicator for each construct variable. *It is to be clarified whether this should only apply to factor latent variables or also composite latent variables.* This is important for ensuring that the signs of the path-coefficients and loadings are consistent. It is sometimes the case in composite-based structural equation modelling methods that loadings/path-coefficients may have the opposite sign. The length of this vector should be equal to the number of construct variables and each value should represent the row number of the dominant indicator for that construct variable. 
#' @param itmax Maximum number of iterations of the Alternating Least Squares (ALS) algorithm.
#' @param ceps Minimum amount of absolute change in the estimates of the path-coefficients (if B0 is non-zero) or the loadings (if B0 is all zero, meaning ther are no path-coefficients) between ALS iterations before ending the optimization.
#' 
#' @author Michael S. Truong
#' 
#' @return List of 4 matrices that make up a fitted I-GSCA Model: (1) Weights, (2) Loadings, (3) Uniqueness Terms D^2, and (4) Path Coefficients.
#' TODO: Cite the GSCA Pro SEM reference instead of GSCA Pro's Example
#' @importFrom MASS ginv 
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
#' # Load the data
#' data(LeDang2022)
#' 
#' # Pre-generate the required matrices for the algorithm and assume that the first indicator for each latent variable is the dominant indicator.
#' igsca_in <- extract_parseModel(model = tutorial_igsca_model,
#'                                    data = LeDang2022,
#'                                    ind_domi_as_first = TRUE)
#'
#' # Fit the I-GSCA model
#' (igsca_out <- with(igsca_in, igsca(Z0 = Z0,
#'                             W0 = W0,
#'                             C0 = C0,
#'                             B0 = B0,
#'                             con_type = con_type,
#'                             indicator_type = indicator_type,
#'                             ind_domi = ind_domi)
#'                             ))
igsca <-
  function(Z0, W0, C0, B0, con_type, indicator_type, ind_domi, itmax = 100, ceps = 0.001) {
  
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
    while ((sum(abs(est0 - est)) > ceps) && (it <= itmax)) {

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
    flipped_signs <-
      flip_signs_ind_domi(
        n_constructs = n_constructs,
        Z = Z,
        ind_domi = ind_domi,
        Gamma = Gamma,
        C = C,
        B = B
      )
    
    list2env(flipped_signs[c("Gamma", "C", "B")], envir = environment())
    

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
        "Uniqueness Terms" = uniqueD
      )
    )

}

#' R Implementation of gsca_inione.m from Heungsun Hwang
#' 
#' Internal I-GSCA Function
#' 
#' Initializes the values for I-GSCA
#' @param Z0 Data Matrix
#' @param W0 Indicator matrix of weights
#' @param B0 Indicator matrix of Path Coefficients.
#' 
#' @author Michael S. Truong
#' @return When used in the context of igsca_sim(), it returns a list of the starting values for Weights, Loadings and Path Coefficients. In principle, otherwise, it is a slightly modified implementation of ordinary Generalised Structured Component Analysis (GSCA). 
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
#' @param Z0 
#' @param W0 
#' @param B0 
#' @param n_indicators 
#' @param n_case 
#' @param n_constructs 
#' @param indicator_type 
#'
#' @return List of matrices to put through the Alternating Least Squared (ALS) algorithm: 
#'
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
#' @param Z Standardized Data
#' @param U 
#' @param D Unique Error
#' @param C Loadings
#' @param n_constructs  Number of constructs
#' @param B Path Coefficients
#'
#' @return Two matrices:
#' - *X*: Remaining part of data (Z) after accounting for uniqueness terms (U) and (D), used for estimating composite loadings. Also used for standardizing theta when updating Gamma, W and V
#' - *WW*: Weights after accounting for current Loading and Path-Coefficients values, used for estimating common-factor loadings
#' 
#' @export
#'
update_X_WW_pseudo_weights <- function(Z, U, D, C, n_constructs, B) {
  # X deviates from Matlab because it is an offspring of svd_out
  X <- Z - U %*% D
  # FIXME: warning("I don't quite understand why this is the equivalent of the Matlab expression, revisit page 44 of Hiebeler 2015 R and Matlab")
  # TODO: Should find alternative to solve() to avoid matrix inversions
  WW <-
    t(C) %*% solve((C %*% t(C) + diag(n_constructs) - 2 * B + (B %*% t(B))))
  return(list("X" = X, "WW" = WW))
}


#' Update Common Factor Variable
#'
#' @param WW 
#' @param windex_gamma_idx 
#' @param gamma_idx 
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
#' @param n_total_var 
#' @param tot 
#' @param n_constructs 
#' @param gamma_idx 
#' @param W 
#' @param A 
#' @param V 
#' @param X 
#' @param windex_gamma_idx 
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
#' @param X 
#' @param n_indicators 
#' @param Gamma 
#' @param c_index 
#' @param C 
#' @param n_constructs 
#' @param b_index 
#' @param B 
#' @param n_case 
#' @param D 
#' @param Z 
#' @param indicator_type 
#'
#' @return List of matrices:
#'
#' 1) Uniqueness terms (D) and (U)
#' 2) Path-Coefficients (B) and Loadings (C)
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
    # FIXME: warning("It's unclear to me whether this SVD is safe. The Matlab code seems to either do a 'normal' svd or a economy svd depending on the dimensionality of the input matrix. Should look into this more.")
    # https://www.mathworks.com/help/matlab/ref/double.svd.html?searchHighlight=svd&s_tid=srchtitle_support_results_1_svd#d126e1597915
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
#' @param n_constructs 
#' @param Z 
#' @param ind_domi 
#' @param gamma_idx 
#' @param Gamma 
#' @param C 
#' @param B 
#'
#' @return List of matrices: Gamma, Loadings (C) and Path-Coefficients (B)
#'
flip_signs_ind_domi <- function(n_constructs, Z, ind_domi, Gamma, C, B) {
  for (gamma_idx in seq_len(n_constructs)) {
    if ((t(Z[, ind_domi[gamma_idx]]) %*% Gamma[, gamma_idx]) < 0) {
      Gamma[, gamma_idx] <- (-1 * Gamma[, gamma_idx])
      C[gamma_idx, ] <- (-1 * C[gamma_idx, ])
      B[gamma_idx, ] <- (-1 * B[gamma_idx, ])
      B[, gamma_idx] <- (-1 * B[, gamma_idx])
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
#' In the context of igsca, this function prepares: (1) the initial indicators (Z0), weights (W0), structural (B0), loadings(C0) matrices; (2) whether a construct is a latent or composite variable (con_type); (3) whether an indicator corresponds to a latent or composite variable (indicator_type); and (4) the dominant indicator of each construct (ind_domi). 
#' 
#' @param model Specified Model in lavaan style
#' @param data Dataframe 
#' @param ind_domi_as_first Boolean for whether the first indicator for each latent factor should be chosen as the dominant indicator
#'
#' @return Returns a list of matrices required for igsca_sim() to run: Z0, W0, B0, C0, con_type, indicator_type, ind_domi. 
#' @export
#' @examples
#'
#' 
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
#' data("LeDang2022")
#' 
#' extract_parseModel(model = tutorial_igsca_model, data = LeDang2022, ind_domi_as_first = TRUE)
extract_parseModel <-
  function(model, data, ind_domi_as_first = TRUE) {
    # Note: parseModel is from cSEM internal
    csemify <- parseModel(.model = model)
    
    Z0 <- data[, csemify$indicators]
    
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
      "Construct indicator does not correctly correspond with Weights Matrix (W0)" = identical(names(csemify$construct_type), colnames(W0))
    )
    
    if (isTRUE(ind_domi_as_first)) {
      # Row Indices of W0 that correspond to the dominant indicator for each factor/composite
      # TODO: Change to W0 to W0[, con_type] if it turns out that the correction should only apply to latent factors
      ind_domi <- apply(W0, 2, as.logical) |>
        apply(2, which, TRUE) |>
        lapply(FUN = \(x) x[[1]]) |>
        unlist()
    } else {
      ind_domi <- NA # FIXME: The length/data-structure might need to be adjusted to match ind_domi when isTRUE(ind_domi_as_first)
    }
    
    return(
      list(
        "Z0" = Z0,
        "B0" = B0,
        "W0" = W0,
        "con_type" = con_type,
        "indicator_type" = indicator_type,
        "C0" = C0,
        "ind_domi" = ind_domi
      )
    )
  }

#' Converts Output of igsca functions into a table to facilitate comparisons
#' 
#' Assumes that indicators only load onto one factor and that there are no cross-factor loadings
#'
#' @param model 
#' @param weights 
#' @param loadings 
#' @param uniqueD 
#' @param paths 
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


#' Takes GSCAPro input and Creates a Lavaan-style Table
#'
#' TODO: Rename the functions here
#' 
#' Assumes that every indicator loads onto only one latent variable (composite/factor)
#' 
#' Expects output from parse_GSCAPro_FullResults.
#' 
#' @param gscapro_in 
#' @param model 
#'
#' @return Table of GSCA Pro results of Weights, Loadings, Path Coefficients and Uniqueness Terms in lavaan style. 
#' @export
#'
get_lavaan_table_igsca_gscapro <- function(gscapro_in, model) {
  
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
  
  
  # Correct names of latent variables
  gscapro_in$Weights$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  gscapro_in$Weights$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  
  
  lv_idx <- list()
  span <- list() 
  
  lv_idx$weights <- which(is.na(gscapro_in$Weights$Estimate))
  # Number of rows after the row of the LV that correspond to that LV's indicators
  span$weights <- diff(lv_idx$weights)
  gscapro_in$Weights$lv <- ""
  gscapro_in$Weights$lv <-
    c(with(gscapro_in$Weights, rep(V1[is.na(Estimate) &
                                      (nchar(V1) > 0)], times = span$weights)), "")
  
  
  # Retrieving indicators from gscapro
  gscapro_indicators <-
    with(gscapro_in$Weights, V1[!(V1 %in% lv) & (nchar(V1) > 0)])
  gscapro_lv <-
    with(gscapro_in$Weights, unique(lv)[nchar(unique(lv)) > 0])
  
  rownames(gscapro_in$Weights)<-gscapro_in$Weights$V1
  
  for (indicator in gscapro_indicators) {
    for (lv in gscapro_lv) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "weights"] <- gscapro_in$Weights[indicator, "Estimate"]
    }
  }
  
  
  # Loadings Parser Shortcut
  gscapro_in$Loadings$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  gscapro_in$Loadings$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  if (identical(gscapro_in$Weights$V1, gscapro_in$Loadings$V1)) {
    rownames(gscapro_in$Loadings) <- gscapro_in$Loadings$V1
    
    for (indicator in gscapro_indicators) {
      for (lv in gscapro_lv) {
        table[((table$lhs == lv &
                  table$rhs == indicator) &
                 table$op %in% c("<~", "=~")), "loadings"] <-
          gscapro_in$Loadings[indicator, "Estimate"]
      }
    }
  } else {
    stop("Cannot take the shortcut of weight indicators for loadings")
  }
  
  # Paths need their own parser 
  
  gscapro_in$`Path coefficients` <-
    gscapro_in$`Path coefficients`[!is.na(gscapro_in$`Path coefficients`$Estimate), ]
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "->",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = " ")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  # ephpostfacto on November 6, 2009; editted by Jilber Urbina on Jan 1/2014 https://stackoverflow.com/a/1690753
  paths <-strsplit(gscapro_in$`Path coefficients`$V1, " ")
  gscapro_in$`Path coefficients`$lvfrom <-sapply(paths, FUN = \(.) .[1])
  gscapro_in$`Path coefficients`$lvto <-sapply(paths, FUN = \(.) .[2])
  
  for (lv_to_iter in gscapro_in$`Path coefficients`$lvto) {
    for (lv_from_iter in gscapro_in$`Path coefficients`$lvfrom) {
      table[((table$rhs == lv_from_iter &
                table$lhs == lv_to_iter) &
               table$op == "~"), "paths"] <- with(gscapro_in$`Path coefficients`, Estimate[(lvfrom == lv_from_iter) & (lvto == lv_to_iter)])
                 
    }
  }
  
  # Parsing for Uniqueness Terms
  for (indicator in names(gscapro_in$`Unique D^2`)) {
    table[((table$rhs == indicator) &
             (table$op == "=~")), "uniqueD"] <- gscapro_in$`Unique D^2`[indicator]
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  return(table)
}

#' Title
#'
#' @param igsca_results 
#'
#' @return FIT index statistic
#' @export
#'
#' @examples summarize_FIT_idx(igsca_results)
summarize_FIT_idx <- function(igsca_results) {
  # TODO: Compute FIT statistic
  return(FIT)
}

#' Title
#'
#' @param model 
#' @param group 
#'
#' @return Multi-group igsca model
#' @export 
#'
#' @examples model_multigroup_igsca(model, group)
model_multigroup_igsca <- function(model, group) {
  
  return(mutligroup_igsca_result)
}