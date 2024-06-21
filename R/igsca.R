#' R Implementation of igsca_sim.m
#' 
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#' 
#' @param Z0 Data matrix of N cases (measurements) x J indicators with named columns
#' @param W0 Indicator matrix of weights: J indicators (rows) and their corresponding Gamma construct variables (columns).
#' @param C0 Indicator matrix of loadings: J indicators (rows) and their corresponding Gamma construct variable (columns).
#' @param B0 Square indicator matrix of path coefficients: from-construct-variable (rows) and to-construct-variable (columns). The order of Gamma construct variables should match the order in C0 and W0.
#' @param con_type A vector that denotes whether each construct variable (columns in W0 and C0) is a common factor or composite. Its length should be equal to the number of columns of W0 and C0. 
#' @param ov_type An indicator vector that indices whether a j indicator (rows of W0 and C0) corresponds to a common factor variable (1) or a composite variable (0). This vector is important for computing the uniqueness terms (D) because it zeros the entries for composite indicators. 
#' @param ind_domi A numeric vector that indices the dominant indicator for each construct variable. *It is to be clarified whether this should only apply to factor latent variables or also composite latent variables.* This is important for ensuring that the signs of the path-coefficients and loadings are consistent. It is sometimes the case in composite-based structural equation modelling methods that loadings/path-coefficients may have the opposite sign. The length of this vector should be equal to the number of construct variables and each value should represent the row number of the dominant indicator for that construct variable. 
#' @param itmax Maximum number of iterations of the Alternating Least Squares (ALS) algorithm.
#' @param ceps Minimum amount of absolute change in the estimates of the path-coefficients (if B0 is non-zero) or the loadings (if B0 is all zero, meaning ther are no path-coefficients) between ALS iterations before ending the optimization.
#' 
#' @author Michael S. Truong
#' 
#' @return List of 4 matrices that make up a fitted I-GSCA Model: (1) Weights, (2) Loadings, (3) Uniqueness Terms D^2, and (4) Path Coefficients.
#' TODO: Cite the GSCA Pro SEM reference instead of GSCA Pro's Example
#' @export
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
#'                             ov_type = ov_type,
#'                             ind_domi = ind_domi)
#'                             ))
igsca <-
  function(Z0,
           W0,
           C0,
           B0,
           con_type,
           ov_type,
           ind_domi,
           itmax = 100,
           ceps = 0.001) {
    
  

  
# Initialize Auxiliary Computational Variables -----------------------------------------------------
  n_case <- nrow(Z0)
  n_indicators <- ncol(Z0)
  n_constructs <- ncol(W0)
  n_total_var <- n_indicators + n_constructs
  
  windex <- which(c(W0) == 1) 
  cindex <- which(c(C0) == 1)
  bindex <- which(c(B0) == 1)
  
  # column and row names are kept aside because swapping with matlab's environment removes the column and row names. Arrays in matlabs don't have row or column names.
  special_names <- list()
  special_names$weights$row <- rownames(W0)
  special_names$weights$col <- colnames(W0)
  special_names$loadings$row <- rownames(C0)
  special_names$loadings$col <- colnames(C0)
  special_names$path$row <- rownames(B0)
  special_names$path$col <- colnames(B0)
  




# Initial Estimates and Preparation -------------------------------------
  prepared_for_ALS <- prepare_for_ALS(
    Z0 = Z0,
    W0 = W0,
    B0 = B0,
    n_indicators = n_indicators,
    n_case = n_case,
    n_constructs = n_constructs,
    ov_type = ov_type
  )
  ### Updates W, C, B, V, Z, D, U, Gamma
  list2env(prepared_for_ALS, envir = environment())
  
  
# Alternating Least Squares Algorithm -------------------------------------

## While Counters and Initial Estimates ----------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (length(bindex) > 0) {
      est <- B[bindex] 
    } else {
      est <- C[cindex]
    }    
    est0 <- est + 1 
    it <- 0
    
## Alternate Between Updating Weights and Loadings -----------------------
    while ((sum(abs(est0 - est)) > ceps) && (it <= itmax)) {

      # Counter Things
      it <- it + 1
      est0 <- est
      
### Compute pseudo-weights ----------------------------------------------------------
      updated_X_weights <- update_X_weights(
        Z = Z,
        U = U,
        D = D,
        C = C,
        n_constructs = n_constructs,
        B = B
      )
      # Creates X (for updating Composites) and WW (for updating Factors)
      list2env(updated_X_weights, envir = environment())
      
### Iterative Update of LVs -------------------------------------------------
      A <- cbind(C, B) 
      
      for (j in seq_len(n_constructs)) {
      # After each cycle, the Gamma, W and V matrices are updated
        tot <- n_indicators + j
        windex_j <- (W0[, j] == 1)
        Xj <- X[, windex_j]
        
        
        if (con_type[j] == "Composite") {
          # Update Composite LV
          theta <-
            update_composite_LV(
              n_total_var = n_total_var,
              tot = tot, # Changes per lv_update iteration
              n_constructs = n_constructs,
              j = j, # Changes per lv_update iteration
              W = W, # Changes per lv_update iteration
              A = A,
              V = V, # Changes per lv_update iteration
              X = X,
              windex_j = windex_j # Changes per lv_update iteration
            )
          
          
          
        } else if (con_type[j] == "Common factor") {
          # Update Factor LV
          theta <-
             update_factor_LV(
               WW = WW, 
               windex_j = windex_j, # Changes per iteration
               j = j # Changes per iteration
               )
          
        } else {
          stop("con_type should either be `Composite` or `Common factor`")
        }
        # This is where the 'actual' updating occurs, in terms of the Gamma matrix, Weights and V(?)
        theta <- theta / norm(Xj %*% theta, "2")
        
        Gamma[, j] <- Xj %*% theta
        W[windex_j, j] <- theta
        V[windex_j, tot] <- theta
        
        
    }
      
### Update Loadings, Path Coefficients and Uniqueness Terms ----------

      updated_C_B_D <- update_C_B_D(
        X = X,
        n_indicators = n_indicators,
        Gamma = Gamma,
        cindex = cindex,
        C = C,
        n_constructs = n_constructs,
        bindex = bindex,
        B = B,
        n_case = n_case,
        D = D,
        Z = Z,
        ov_type = ov_type
      )
      # Updates C, B, D, uniqueD, est, and U
      list2env(updated_C_B_D, envir = environment())
    }
    

# Flip Signs for Factors and Composites Based on Dominant Indicators --------

    flipped_signs <-
      flip_signs_ind_domi(
        n_constructs = n_constructs,
        Z = Z,
        ind_domi = ind_domi,
        j = j,
        Gamma = Gamma,
        C = C,
        B = B
      )
    # Updates Gamma, C and B
    list2env(flipped_signs, envir = environment())
    
    
    
    mW <- W[windex] 
    # This was commented out in the Matlab version because the loadings dimensions
    # are transposed relative to the weights matrix
    # mC <- c[cindex] 
    mC <- t(C)[t(C0) == 1] 
    mB <- B[bindex] 
    mD <- uniqueD
  
  
  # Format Output
  Weights <- W
  rownames(Weights) <- special_names$weights$row
  colnames(Weights) <- special_names$weights$col
  
  Loadings <- t(C)
  rownames(Loadings) <- special_names$loadings$row
  colnames(Loadings) <- special_names$loadings$col
  
  PathCoefficients <- B
  rownames(PathCoefficients) <- special_names$path$row
  colnames(PathCoefficients) <- special_names$path$col
  
  # browser()
  
  names(uniqueD) <- rownames(Loadings)
  
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
  
  windex <- which(W0 != 0)
  aindex <- which(A0 != 0)
  
  Z <- scale(Z0, center = TRUE, scale = TRUE) * sqrt(N) / sqrt(N - 1)
  # Random Values to W and A
  W <- W0
  A <- A0
  
  W[windex] <- rep(1, length(windex))
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
#' @param ov_type 
#'
#' @return List of matrices to put through the Alternating Least Squared (ALS) algorithm: 
#'
prepare_for_ALS <- function(Z0, W0, B0, n_indicators, n_case, n_constructs, ov_type) {
  
  
  
  initial_est <-
    gsca_inione(
      Z0 = Z0,
      W0 = apply(W0 != 0, 2, as.numeric),
      B0 = apply(B0 != 0, 2, as.numeric)
    )
  # TODO: Update using list2env()
  W <- initial_est$W
  C <- initial_est$C
  B <- initial_est$B
  
  
  V <- cbind(diag(n_indicators), W)
  Z <- scale(Z0, center = TRUE, scale = TRUE) / sqrt(n_case - 1)
  Gamma <- Z %*% W
  D <- diag(n_indicators)
  
  
  # TODO: Does the LAPACK argument for qr() still matter? Why does `complete` matter?
  # Solution for Q is copied from estimators_weights.R
  Q <- qr.Q(qr(Gamma), complete =  TRUE)
  F_o <- Q[, (n_constructs + 1):n_case, drop = FALSE]
  # FIXME: Should compare the svd of matlab and R, unsure if it matters that they don't match
  # FIXME: Come back to this stack overflow thing for proper citation
  # Ahmed Fasih https://stackoverflow.com/a/41972818
  svd_out <- (D %*% t(Z) %*% F_o) |>
    {
      \(mx) svd(mx, nu = nrow(mx),  nv = ncol(mx))
    }()
  u <- svd_out$u
  v <- svd_out$v
  
  # TODO: Utilde deviates from Matlab because of the SVD
  Utilde <- v[, 1:n_indicators] %*% t(u)
  U <- F_o %*% Utilde
  D <- diag(diag(t(U) %*% Z))
  D[ov_type == 0, ov_type == 0] <- 0
  
  
  return(list(
    "W" = W,
    "C" = C,
    "B" = B,
    "V" = V,
    "Z" = Z,
    "D" = D,
    "U" = U,
    "Gamma" = Gamma
    # "Utilde" = Utilde,
    # "v" = v,
    # "u" = u,
    # "F_o" = F_o,
    # "Q" = Q,
  ))
}

#' Update Weights and X (?)
#'
#' @param Z 
#' @param U 
#' @param D 
#' @param C 
#' @param n_constructs 
#' @param B 
#'
#' @return Two matrices:
#' - *X*: Remaining part of data (Z) after accounting for uniqueness terms (U) and (D)(?)
#' - *WW*: *presumably* Weights after accounting for current Loading and Path-Coefficients values.
#' TODO: Double check that I've got this right.
#' @export
#'
update_X_weights <- function(Z, U, D, C, n_constructs, B) {
  # X deviates from Matlab because it is an offspring of svd_out
  X <- Z - U %*% D
  # FIXME: warning("I don't quite understand why this is the equivalent of the Matlab expression, revisit page 44 of Hiebeler 2015 R and Matlab")
  # TODO: Should find alternative to solve() to avoid matrix inversions
  WW <-
    t(C) %*% solve((C %*% t(C) + diag(n_constructs) - 2 * B + (B %*% t(B))))
  return(list("X" = X, "WW" = WW))
}


#' Update Factor Latent Variable
#'
#' @param WW 
#' @param windex_j 
#' @param j 
#'
#' @return theta: Used to update factor latent variables -- after accounting for loadings and path-coefficients.
#' 
#' TODO: Double check that this is correct.
#' @export
#'
update_factor_LV <- function(WW, windex_j, j) {
  theta <- WW[windex_j, j]
  return(theta)
}





#' Update Composite Latent Variables
#'
#' @param n_total_var 
#' @param tot 
#' @param n_constructs 
#' @param j 
#' @param W 
#' @param A 
#' @param V 
#' @param X 
#' @param windex_j 
#'
#' @return theta: A matrix that will later be used to update the weights for the composite variable.
#' TODO: Double check that this is for weights only?
#' @export
#'
update_composite_LV <-
  function(n_total_var, tot, n_constructs, j, W, A, V, X, windex_j) {
    # This updates the composite
    e <- matrix(0, nrow = 1, ncol = n_total_var)
    e[tot] <- 1
    H1 <- diag(n_total_var)
    H2 <- diag(n_constructs)
    H1[tot, tot] <- 0
    H2[j, j] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)
    
    
    vecZDelta <- c(X %*% Delta)
    beta <- e - A[j, ]
    XI <- kronecker(t(beta), X)
    XI <- XI[, windex_j]
    
    theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta
    return(theta)
  }

#' Update Loadings, Path-Coefficients and Uniqueness Terms After Updating Latent Variables
#'
#' @param X 
#' @param n_indicators 
#' @param Gamma 
#' @param cindex 
#' @param C 
#' @param n_constructs 
#' @param bindex 
#' @param B 
#' @param n_case 
#' @param D 
#' @param Z 
#' @param ov_type 
#'
#' @return List of matrices:
#'
#' 1) Uniqueness terms (D) and (\eqn{D^2})
#' 2) Non-zero Path-Coefficients (if Paths are specified in B0) or non-zero Loadings 
#' 3) Path-Coefficients (B) and Loadings (C)
#' 4) U, Utilde, v, F_o, Q, U *something*
#' TODO: Name these accessory matrices
#'
update_C_B_D <-
  function(X,
           n_indicators,
           Gamma,
           cindex,
           C,
           n_constructs,
           bindex,
           B,
           n_case,
           D,
           Z,
           ov_type) {
    t1 <- c(X)
    M1 <- kronecker(diag(n_indicators), Gamma)
    M1 <- M1[, cindex]
    C[cindex] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)
    
    t2 <- c(Gamma)
    M2 <- kronecker(diag(n_constructs), Gamma)
    M2 <- M2[, bindex]
    B[bindex] <- MASS::ginv(t(M2) %*% M2) %*% (t(M2) %*% t2)
    
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
    D[ov_type != 1, ov_type != 1] <- 0
    
    
    if (length(bindex) > 0) {
      est <- B[bindex]
    } else {
      est <- C[cindex]
    }
    # browser()
    uniqueD <- diag(D) ^ 2
    
    return(
      list(
        "D" = D,
        "uniqueD" = uniqueD,
        "est" = est,
        "U" = U,
        "Utilde" = Utilde,# TODO: This is used somewhere?
        "v" = v,
        "u" = u,
        "F_o" = F_o,
        "Q" = Q, # TODO: This is used somewhere else?
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
#' @param j 
#' @param Gamma 
#' @param C 
#' @param B 
#'
#' @return List of matrices: Gamma, Loadings (C) and Path-Coefficients (B)
#'
flip_signs_ind_domi <- function(n_constructs, Z, ind_domi, j, Gamma, C, B) {
  for (j in seq_len(n_constructs)) {
    if ((t(Z[, ind_domi[j]]) %*% Gamma[, j]) < 0) {
      Gamma[, j] <- (-1 * Gamma[, j])
      C[j, ] <- (-1 * C[j, ])
      B[j, ] <- (-1 * B[j, ])
      B[, j] <- (-1 * B[, j])
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
#' In the context of igsca, this function prepares: (1) the initial indicators (Z0), weights (W0), structural (B0), loadings(C0) matrices; (2) whether a construct is a latent or composite variable (con_type); (3) whether an indicator corresponds to a latent or composite variable (ov_type); and (4) the dominant indicator of each construct (ind_domi). 
#' 
#' @param model Specified Model in lavaan style
#' @param data Dataframe 
#' @param ind_domi_as_first Boolean for whether the first indicator for each latent factor should be chosen as the dominant indicator
#'
#' @return Returns a list of matrices required for igsca_sim() to run: Z0, W0, B0, C0, con_type, ov_type, ind_domi. 
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
    # browser()
    # B0 <- csemify$structural
    B0 <- t(csemify$structural)
    W0 <- t(csemify$measurement)
    
    # con_type <- csemify$construct_type == "Common factor"
    con_type <- csemify$construct_type
    
    # Constructing ov_type
    ov_type <- vector(length = ncol(csemify$measurement))
    names(ov_type) <- colnames(csemify$measurement)
    
    
    for (lv in rownames(csemify$measurement)) {
      for (indicator in colnames(csemify$measurement)) {
        if (csemify$construct_type[lv] == "Common factor") {
          
          if (csemify$measurement[lv, indicator] == 1) {
            ov_type[indicator] <- TRUE
          }
          
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
        "ov_type" = ov_type,
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