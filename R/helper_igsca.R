#' R Implementation of igsca_sim.m
#'
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#'
#' In the example section, the specified model is based on the tutorial I-GSCA model associated with GSCA Pro \insertCite{hwangetal2023StructuralEquationModelingAMultidisciplinaryJournal}{cSEM}.
#' 
#' **Note**: Here, we assume that there is only one unique loading per indicator. 
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
#' @inheritParams csem_arguments
#'
#' @author Michael S. Truong
#' @return List of 4 matrices that make up a fitted I-GSCA Model:
#' * (1) Weights
#' * (2) Loadings
#' * (3) Squared Unique Loadings D^2
#' * (4) Path Coefficients.
#' * (5) Unique Component of Indicators (for common factors) DU
#'
#' @importFrom MASS ginv
#' @import RcppArmadillo
#'
#' @references
#'   \insertAllCited{}
#' @examples
#' \dontrun{
#' # Specify the model according to GSCA Pro's example
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 +
#'                       Behavior7 + Behavior8 +  Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2
#'
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 +
#'                     Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 +
#'                 Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 +
#'                 Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 +
#'                  Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 +
#'                      Conscientious6 + Conscientious7 + Conscientious8 +
#'                      Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 +
#'                         Openness7 + Openness8 + Openness9 + Openness10
#'
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion +
#'                      Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#'
#' data(LeDang2022)
#'
#' csem(.data = LeDang2022, tutorial_igsca_model, .approach_weights = "GSCA",
#' .dominant_indicators = NULL, .tolerance = 0.0001, .conv_criterion =
#' "sum_diff_absolute")
#' }
igsca <-
  function(
    Z0,
    W0,
    C0,
    B0,
    con_type,
    indicator_type,
    .dominant_indicators,
    .iter_max = 100,
    .tolerance = 0.0001,
    .conv_criterion,
    .S = args_default()$.S
  ) {
    ## Initialize Computational Variables -----------------------------------------------------
    n_case <- nrow(Z0)
    n_indicators <- ncol(Z0)
    n_constructs <- ncol(W0)
    n_total_var <- n_indicators + n_constructs

    w_index <- which(c(W0) == 1)
    c_index <- which(c(C0) == 1)
    b_index <- which(c(B0) == 1)

    # Initialize Bindpoints for list2env
    Z <- matrix()
    C <- matrix()
    W <- matrix()
    B <- matrix()
    V <- matrix()
    D <- matrix()
    U <- matrix()
    Gamma <- matrix()
    X <- matrix()
    WW <- matrix()

    ### Initial Estimates and Preparation -------------------------------------
    prepared_for_ALS <- initializeAlsEstimates(
      Z0 = Z0,
      W0 = W0,
      B0 = B0,
      n_indicators = n_indicators,
      n_case = n_case,
      n_constructs = n_constructs,
      indicator_type = indicator_type,
      .S = .S
    )

    list2env(
      prepared_for_ALS[c("W", "C", "B", "V", "Z", "D", "U", "Gamma")],
      envir = environment()
    )

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

    while (
      (!checkConvergence(
        .W_new = est,
        .W_old = est0,
        .tolerance = .tolerance,
        .conv_criterion = .conv_criterion
      )) &&
        (it <= .iter_max)
    ) {
      # Update Counter Variables
      it <- it + 1
      est0 <- est

      #### Compute Pseudo-Weights --------------------------------------------------
      updated_X_WW_pseudo_weights <- updateXAndWW(
        Z = Z,
        U = U,
        D = D,
        C = C,
        B = B,
        n_constructs = n_constructs
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
          Theta <-
            updateCompositeTheta(
              W = W, # Changes per gamma_idx iteration
              A = A,
              X = X, # X = Z- UD, as per updateXAndWW()
              V = V, # Changes per gamma_idx iteration
              gamma_idx = gamma_idx, # Changes per gamma_idx iteration
              tot = tot, # Changes per gamma_idx iteration
              n_constructs = n_constructs,
              n_total_var = n_total_var,
              windex_gamma_idx = windex_gamma_idx, # Changes per gamma_idx iteration
              .S = .S
            )
        } else if (con_type[gamma_idx] == "Common factor") {
          Theta <-
            updateCommonFactorTheta(
              WW = WW,
              windex_gamma_idx = windex_gamma_idx, # Changes per gamma_idx iteration
              gamma_idx = gamma_idx # Changes per gamma_idx iteration
            )
        } else {
          stop("con_type should only either be `Composite` or `Common factor`")
        }

        # This is where Gamma, Weights and V are updated based on which gamma_idx
        Theta <- Theta / norm(X_gamma_idx %*% Theta, "2")
        Gamma[, gamma_idx] <- X_gamma_idx %*% Theta
        W[windex_gamma_idx, gamma_idx] <- Theta
        V[windex_gamma_idx, tot] <- Theta
      }

      #### Update Loadings, Path Coefficients and Uniqueness Terms ----------
      updated_C_B_D_U <- updateCBDU(
        X = X,
        Gamma = Gamma,
        C = C,
        B = B,
        D = D,
        Z = Z,
        n_indicators = n_indicators,
        indicator_type = indicator_type,
        n_constructs = n_constructs,
        n_case = n_case,
        c_index = c_index,
        b_index = b_index,
        con_type = con_type
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
        flipGammaCBSigns(
          Z = Z,
          C = C,
          B = B,
          Gamma = Gamma,
          n_constructs = n_constructs,
          .dominant_indicators = .dominant_indicators
        )

      list2env(flipped_signs[c("Gamma", "C", "B")], envir = environment())
    }

    ## Output Formatting -------------------------------------------------------
    # identical(diag(D^2), diag(D)^2)
    D_diag <- diag(D)
    names(D_diag) <- colnames(Z)
    # colnames(D_squared) <- colnames(Z)
    # rownames(D_squared) <- colnames(Z)

    Unique_scores <- U
    colnames(Unique_scores) <- colnames(Z)

    return(
      list(
        "W" = t(W), # W is J X P, so W^T is P \times J. As shown in the examples ?csem, res$Estimates$Weight_estimates should be P \times J
        "C" = C, # C is P \times J. As shown in the examples ?csem, res$Estimates$Loading_estimates should be P \times J
        "B" = t(B), # B is From \times To; so t(B) is To \times From. As shown in the examples ?csem, res$Estimates$Path_estimates should be To \times From
        "Construct_scores" = (Z - (U %*% D)) %*% W,
        "Unique_loading_estimates" = D_diag,
        "Unique_scores" = Unique_scores,
        "Modes" = "gsca", 
        "Conv_status" = ifelse(it > .iter_max, FALSE, TRUE),
        "Iterations" = it,
        "Data" = Z # Z is N \times P
      )
    )
  }

#' R Implementation of gsca_inione.m from Heungsun Hwang
#'
#' Initializes the values for I-GSCA by using modified Generalized Structured
#' Component Analysis (GSCA) to estimate the model.
#'
#' @inheritParams igsca
#' @inheritParams csem_arguments
#' @importFrom MASS ginv
#' @author Michael S. Truong
#' @returns Returns a list of starting values for:
#' * Weights (W)
#' * Loadings (C)
#' * Path Coefficients (B)
#'
initializeIgscaEstimates <- function(Z0, W0, B0, .S = args_default()$.S) {
  N <- nrow(Z0)
  J <- nrow(W0)
  P <- ncol(W0)
  TRep <- J + P
  C0 <- t(W0)
  A0 <- cbind(C0, B0)

  w_index <- which(W0 != 0)
  aindex <- which(A0 != 0)

  # The original algorithm may be doing this to use a biased consistent estimator of the standard deviation instead of unbiased, not sure why.
  Z <- scale(Z0, center = TRUE, scale = TRUE) * sqrt(N) / sqrt(N - 1)
  # Random Values to W and A
  W <- W0
  A <- A0

  W[w_index] <- rep(1, length(w_index))
  A[aindex] <- runif(length(aindex), min = 0, max = 1)

  Gamma <- Z %*% W
  W <-
    W /
    (t(sqrt(diag(
      t(
        Gamma
      ) %*%
        Gamma
    ))) |>
      rep(each = nrow(W)))
  V <- cbind(diag(J), W)
  Gamma <- Z %*% W

  Psi <- Z %*% V
  vecPsi <- c(Psi)
  it <- 0
  f0 <- 100000000
  imp <- 100000000
  while ((it <= 300) && (imp > 0.0001)) {
    it <- it + 1
    # Phi = kronecker(diag(TRep), Gamma)
    # Phi <- Phi[, aindex]
    Phi <- kroneckerC(diag(TRep), Gamma, aindex)
    A[aindex] <- solve(t(Phi) %*% Phi, t(Phi)) %*% vecPsi

    for (p in seq_len(P)) {
      t_lil <- J + p
      windex_p <- which(W0[, p] != 0)
      m <- matrix(0, nrow = 1, ncol = TRep)
      m[t_lil] <- 1
      a <- A[p, ]
      beta <- m - a
      H1 <- diag(P)
      H2 <- diag(TRep)
      H1[p, p] <- 0
      H2[t_lil, t_lil] <- 0
      Delta <- (W %*% H1 %*% A) - (V %*% H2)
      # vecZDelta <- c(Z %*% Delta)

      # XI <- kronecker(t(beta), Z)
      # XI <- XI[, windex_p]

      # XI <- kroneckerC(t(beta), Z, windex_p)
      # Theta <- MASS::ginv(t(XI) %*% XI) %*% t(XI) %*% vecZDelta

      # Kronecker bypass -- as shown in calculateWeightsGSCA.R
      Theta <- MASS::ginv(
        as.numeric(beta %*% t(beta)) * .S[windex_p, windex_p]
      ) %*%
        t(beta %*% t(Delta) %*% .S[, windex_p])

      zw <- Z[, windex_p] %*% Theta

      Theta <- sqrt(N) * Theta / norm(zw, "2")
      W[windex_p, p] <- Theta
      V[windex_p, t_lil] <- Theta
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
#' @inheritParams csem_arguments
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
initializeAlsEstimates <- function(
  Z0,
  W0,
  B0,
  n_indicators,
  n_case,
  n_constructs,
  indicator_type,
  .S = args_default()$.S
) {
  # Initialize Bindpoints for list2env
  W <- matrix()
  B <- matrix()

  # Initial estimates using GSCA
  initial_est <-
    initializeIgscaEstimates(
      Z0 = Z0,
      W0 = W0,
      B0 = B0,
      .S = .S
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
  Q <- qr.Q(qr(Gamma), complete = TRUE)
  F_o <- Q[, (n_constructs + 1):n_case, drop = FALSE]

  if (n_case <= n_indicators) {
    # I don't know if the following warning is necessary given that it should
    # still work -- this is said to come from the old gsca_m implementation in
    # estimators_weights.R
    #
    # This also comes from the Supplementary Material to GSCA_m, so I think it
    # should work
    #
    # warning("Model may not be correctly estimated because the number of cases is less than the number of indicators and the handling is experimental.")

    # From estimator_weights.R, which is from one of the older GSCA_m implementations
    # Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (n_constructs+1):n_case, drop = FALSE]
    s <- svd(D %*% t(Z) %*% F_o)
    Utilde <- s$v %*% t(s$u)
    U <- F_o %*% Utilde # Gamma_orth = F_o
  } else if (n_case > n_indicators) {
    # svd between R and Matlab by Ahmed Fasih on February 1/2017
    # https://stackoverflow.com/a/41972818

    U <- tryCatch(
      {
        svd_out <- (D %*% t(Z) %*% F_o) |>
          {
            \(mx) svd(mx, nu = nrow(mx), nv = ncol(mx))
          }()
        u <- svd_out$u
        v <- svd_out$v
        # Utilde deviates from Matlab because of the SVD
        Utilde <- v[, 1:n_indicators] %*% t(u)
        U <- F_o %*% Utilde
      },
      error = function(e) {
        s <- svd(D %*% t(Z) %*% F_o)
        Utilde <- s$v %*% t(s$u)
        U <- F_o %*% Utilde # Gamma_orth = F_o
        return(U)
      }
    )
  }
  # Zeroing the initialized Unique_scores of composites is unique to our R implementation and not found in the original Matlab version
  # AFAICT, this shouldn't affect the parameter estimates because U only affects the other parameter estimates via U %*% D, and D will
  #  already zero the effects of the initialization
  U[, indicator_type == "Composite"] <- 0
  D <- diag(diag(t(U) %*% Z)) # TODO: This forces the D matrix to be diagonal-entries only.
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

#' Update X and WW (Pseudo Weights) for Composites and Common-Factors
#'
#' Computation of WW matrix was converted between Matlab to R based on page 44 of \insertCite{Hiebeler2015;textual}{cSEM}.
#' @param Z Standardized data matrix
#' @param D Matrix of estimated unique error
#' @param U Matrix of estimates related to unique error
#' @param C Matrix of estimated loadings
#' @param B Matrix of path coefficients
#' @param n_constructs  Number of constructs
#'
#' @returns Two matrices:
#' * X: Remaining part of data (Z) after accounting for uniqueness terms (U) and (D), used for estimating composite loadings. Also used for standardizing Theta when updating Gamma, W and V
#' * WW: Weights after accounting for current Loading and Path-Coefficients values, used for estimating common-factor loadings
#'
#' @references
#'   \insertAllCited{}
#'
#'
updateXAndWW <- function(Z, U, D, C, B, n_constructs) {
  # X deviates from Matlab because it is an offspring of svd_out
  X <- Z - U %*% D

  WW <-
    t(C) %*% solve((C %*% t(C) + diag(n_constructs) - 2 * B + (B %*% t(B))))
  return(list("X" = X, "WW" = WW))
}


#' Update Theta for Common Factor Variable
#'
#' @param WW Pseudo-weights for Common Factors
#' @param windex_gamma_idx Index of weights related to the indicators for the construct of interest
#' @param gamma_idx Index of which construct we are examining
#'
#' @return Theta: Used to update factor latent variables -- after accounting for loadings and path-coefficients.
#'
#'
#'
updateCommonFactorTheta <- function(WW, windex_gamma_idx, gamma_idx) {
  Theta <- WW[windex_gamma_idx, gamma_idx]
  return(Theta)
}

#' Update Theta for Composite Variables
#' 
#' It is unintuitive that X is used here, seeing as how X = Z-UD; and we use X to update composite variables. 
#' However, from a non-computational point of view, it shouldn't matter because for the composite indicators, X_comp = Z_comp
#'
#' @param W Weights matrix
#' @param A Stacked matrix of loadings and path coefficients
#' @param V Unclear meaning
#' @param X The matrix X is equal to Z - UD
#' @param windex_gamma_idx Index of weights related to the indicators for the construct of interest
#' @param n_total_var Number of indicators and constructs
#' @param tot Index dependent on which construct variable we are examining
#' @param n_constructs Number of constructs
#' @param gamma_idx Index of which construct we are examining
#' @inheritParams csem_arguments
#' @importFrom MASS ginv
#' @return Theta: A matrix that will later be used to update the weights for the composite variable.
#'
#'
updateCompositeTheta <-
  function(
    W,
    A,
    V,
    X,
    windex_gamma_idx,
    n_total_var,
    tot,
    n_constructs,
    gamma_idx,
    .S = args_default()$.S
  ) {
    e <- matrix(0, nrow = 1, ncol = n_total_var)
    e[tot] <- 1
    H1 <- diag(n_total_var)
    H2 <- diag(n_constructs)
    H1[tot, tot] <- 0
    H2[gamma_idx, gamma_idx] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)

    # vecZDelta <- c(X %*% Delta) # Commented out because no longer computing kronecker product
    beta <- e - A[gamma_idx, ]
    # XI <- kronecker(t(beta), X)
    # XI <- XI[, windex_gamma_idx]

    # XI <- kroneckerC(t(beta), X, which(windex_gamma_idx))
    # Theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta

    # Kronecker bypass -- as shown in calculateWeightsGSCA.R
    # Bypassing the kronecker product this way leads to minor differences in results,
    # but the original matlab solution may already have been inaccurate due to matrix inversion.
    #
    # The best way to identify whether this difference does or does not matter is by a combination of mathematics
    #  (theoretical equivalence, ignoring computational details) and simulation study
    Theta <- MASS::ginv(
      as.numeric(beta %*% t(beta)) * .S[windex_gamma_idx, windex_gamma_idx]
    ) %*%
      t(beta %*% t(Delta) %*% .S[, windex_gamma_idx])

    return(Theta)
  }

#' Update Loadings, Path-Coefficients and Uniqueness Terms After Updating Latent Variables
#'
#' @inheritParams initializeAlsEstimates
#' @inheritParams updateXAndWW
#' @inheritParams igsca
#' @param X Pseudo-weights for composites
#' @param Gamma Construct Scores
#' @param c_index Index of loadings
#' @param b_index Index of Path Coefficients
#' @param n_case Number of Cases
#' @param indicator_type Vector of whether each indicator corresponds to a common factor or composite
#' @importFrom MASS ginv
#' @return List of matrices:
#'
#' * (1) Uniqueness terms (D) and (U)
#' * (2) Estimated Path Coefficients matrix (B)
#' * (3) Estimated Loadings matrix (C)
#'
updateCBDU <-
  function(
    Z,
    X,
    Gamma,
    C,
    B,
    D,
    indicator_type,
    n_indicators,
    c_index,
    n_constructs,
    b_index,
    n_case,
    con_type
  ) {
    ## Loading Update ----------------------------------------------------------

    t1 <- c(X)
    # M1 <- kronecker(diag(n_indicators), Gamma)
    # M1 <- M1[, c_index]
    M1 <- kroneckerC(diag(n_indicators), Gamma, c_index)
    C[c_index] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)
    #
    # Kronecker bypass as shown in calculateWeightsGSCAm
    # vcv_gamma <- t(Gamma) %*% Gamma # Came from path coefficients update in calculateWeightsGSCAm
    # TODO: This will have to be adjusted to differentiate between canonical and
    # nomological composites for IGSCA. Right now, it assumes that all the
    # constructs need loadings. If there are only canonical composites, then
    # there needs to be a bypass
    # vars_cf_nomocomp <- which(con_type %in% c("Common factor", "Composite"))
    # cov_gamma_indicators <- t(Gamma) %*% Z
    # Y <- which(colSums(C[vars_cf_nomocomp, ]) != 0)
    # loadings <- lapply(Y, function(y) {
    #   x    <-  which(C[vars_cf_nomocomp, y] != 0)
    # FIXME: The difference between inverse via solve and pseudo-inverse does
    # not account for differences between this approach and kronecker
    # coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*% cov_gamma_indicators[x, y, drop = FALSE]
    #   coef <- solve(vcv_gamma[x, x, drop = FALSE]) %*% cov_gamma_indicators[x, y, drop = FALSE]
    # })
    #
    # Ct<-t(C)
    # Ct[c_index] <- unlist(loadings, use.names = FALSE)
    # C<-t(Ct)

    # Path Coefficients Update ------------------------------------------------

    # t2 <- c(Gamma)
    # M2 <- kronecker(diag(n_constructs), Gamma)
    # M2 <- M2[, b_index]
    #
    # M2 <- kroneckerC(diag(n_constructs), Gamma, b_index)
    # B[b_index] <- MASS::ginv(t(M2) %*% M2) %*% (t(M2) %*% t2)
    #
    # Kronecker bypass as shown in calculateWeightsGSCAm
    # TODO: Remember to comment vcv_gamma assignment out if I restore kronecker by-pass to loadings
    vcv_gamma <- t(Gamma) %*% Gamma
    vars_endo <- which(colSums(B) != 0)

    beta <- lapply(vars_endo, function(y) {
      x <- which(B[, y, drop = FALSE] != 0)
      coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*%
        vcv_gamma[x, y, drop = FALSE]
    })
    B[b_index] <- unlist(beta, use.names = FALSE)

    ## Uniqueness Component Update ---------------------------------------------
    
    # TODO: The following procedure creates non-zero columns of U, even for composite indicators. It would be better if this could be avoided.
    # Solution for Q is copied from estimators_weights.R
    Q <- qr.Q(qr(Gamma), complete = TRUE)
    F_o <- Q[, (n_constructs + 1):n_case]

    if (n_case <= n_indicators) {
      # From estimator_weights.R, which is from one of the older GSCA_m implementations
      # Gamma_orth <- qr.Q(qr(Gamma), complete = TRUE)[, (n_constructs+1):n_case, drop = FALSE]
      s <- svd(D %*% t(Z) %*% F_o)
      Utilde <- s$v %*% t(s$u)
      U <- F_o %*% Utilde # Gamma_orth = F_o
    } else if (n_case > n_indicators) {
      # svd between R and Matlab by Ahmed Fasih on February 1/2017
      # https://stackoverflow.com/a/41972818

      U <- tryCatch(
        {
          svd_out <- (D %*% t(Z) %*% F_o) |>
            {
              \(mx) svd(mx, nu = nrow(mx), nv = ncol(mx))
            }()
          u <- svd_out$u
          v <- svd_out$v
          # Utilde deviates from Matlab because of the SVD
          Utilde <- v[, 1:n_indicators] %*% t(u)
          U <- F_o %*% Utilde
        },
        error = function(e) {
          s <- svd(D %*% t(Z) %*% F_o)
          Utilde <- s$v %*% t(s$u)
          U <- F_o %*% Utilde # Gamma_orth = F_o
          return(U)
        }
      )
    }

    U[, indicator_type == "Composite"] <- 0
    D <- diag(diag(t(U) %*% Z))
    D[indicator_type == "Composite", indicator_type == "Composite"] <- 0

    return(
      list(
        "C" = C,
        "B" = B,
        "D" = D,
        "U" = U
      )
    )
  }

#' Flip signs of Gamma, Loadings and Path-Coefficients Cells Based on Dominant Indicator
#'
#' @inheritParams igsca
#' @inheritParams updateXAndWW
#' @inheritParams updateCBDU
#'
#' @return List of matrices:
#' * Estimated Construct Scores (Gamma)
#' * Estimated Loadings matrix (C)
#' * Estimated Path-Coefficients matrix (B)
#'
flipGammaCBSigns <- function(
  Z,
  Gamma,
  C,
  B,
  n_constructs,
  .dominant_indicators
) {
  for (gamma_idx in seq_len(n_constructs)) {
    if (.dominant_indicators[gamma_idx] %in% colnames(Z)) {
      if (
        (t(Z[, .dominant_indicators[gamma_idx]]) %*% Gamma[, gamma_idx]) < 0
      ) {
        Gamma[, gamma_idx] <- (-1 * Gamma[, gamma_idx])
        C[gamma_idx, ] <- (-1 * C[gamma_idx, ])
        B[gamma_idx, ] <- (-1 * B[gamma_idx, ])
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

#' Gets the Relevant Inputs for IGSCA
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
getIgscaInputs <-
  function(.data, .model) {
    # Note: parseModel is from cSEM internal

    csemify <- parseModel(.model = .model)

    Z0 <- .data[, csemify$indicators]

    # Igsca assumes \eta \times B, so the rows should be from
    # and the columns should be to. This is in contrast to 
    # the rest of cSEM
    B0 <- t(csemify$structural) 

    C0 <- csemify$measurement
    W0 <- t(csemify$measurement)

    # con_type <- csemify$construct_type == "Common factor"
    con_type <- csemify$construct_type

    # Constructing indicator_type
    indicator_type <- vector(
      mode = "character",
      length = ncol(csemify$measurement)
    )
    # Default value must be "Composite" because of how the measurement matrix
    #  returned by parseModel uses 0 to denote both composite variables and the 
    # absence of any corresponding construct variable.
    indicator_type <- rep("Composite", length(indicator_type))
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
          warning(
            "Indicator does not correspond to either a Composite or Common factor. Unsupported behavior may ensue."
          )
        }
      }
    }

    return(
      list(
        "Z0" = Z0,
        "B0" = B0,
        "W0" = W0,
        "C0" = C0,
        "con_type" = con_type,
        "indicator_type" = indicator_type
      )
    )
  }

#' Block Diagonalize Estimated Parameter Matrices to Facilitate Computation of FIT Statistics
#'
#' Block diagonalizes the estimated paramater matrices as shown on Equations
#' 3.28-3.29 on page 111 of \insertCite{Hwang2014;textual}{cSEM}. Should only be used on multi-group models
#'
#' @inheritParams tidy.cSEMResults
#'
#' @return cSEMResults in single-group data structure with block diagonalized parameter estimates
#' @export
#' @importFrom Matrix bdiag
bdiagonalizeMultiGroupIgscaEstimates <- function(x) {
  if (!identical(names(x), c("Estimates", "Information"))) {
    # Multi-Group Code

    ## Extract Matrices ------------------------------------------------------
    # Extract Estimated Matrices
    # TODO: Test whether the following code works

    estimates_to_be_extracted <- list(
      "Path_estimates",
      "Loading_estimates",
      "Weight_estimates",
      "Construct_scores",
      "Unique_scores"
    )

    names(estimates_to_be_extracted) <- unlist(estimates_to_be_extracted)

    extraction <- lapply(
      X = estimates_to_be_extracted,
      function(matrix_name, multigroup_output) {
        extraction <- lapply(
          multigroup_output,
          function(onegroup_output, matrix_name) {
            return(onegroup_output[["Estimates"]][[matrix_name]])
          },
          matrix_name = matrix_name
        )
        return(extraction)
      },
      multigroup_output = x
    )

    # Extract Data
    extraction$Data <- lapply(x, function(onegroup_output) {
      return(onegroup_output[["Information"]][["Data"]])
    })

    # Remove Null Matrices
    extracts_to_remove <- lapply(extraction, \(x) {
      lapply(x, is.null) |> unlist() |> all()
    })

    extraction <- extraction[which(
      !unlist(lapply(extraction, \(x) lapply(x, is.null) |> unlist() |> all()))
    )]

    extraction <- mapply(
      function(extract, extract_name) {
        # We don't keep it as a sparse matrix because that might break
        # functionality with other functions unless much more of Matrix is
        # imported

        bdiaged <- Matrix::bdiag(extract)
        colnames(bdiaged) <- rep(
          colnames(extract[[1]]),
          times = length(extract)
        )
        if (
          !(extract_name %in% c("Data", "Construct_scores", "Unique_scores"))
        ) {
          rownames(bdiaged) <- rep(
            rownames(extract[[1]]),
            times = length(extract)
          )
        }

        return(as.matrix(bdiaged))
      },
      extract = extraction,
      extract_name = names(extraction),
      SIMPLIFY = FALSE
    )

    ## Create Surrogate Output -----------------------------------------------
    surrogate_out <- list()
    # Insert the diagonalized matrices into the surrogate
    ## It's not simple to take x[[1]] as the surrogate structure because it has
    ## many estimatates (such as reliabilities) that are specific to the group
    ## model

    for (extract_name in names(extraction)[which(
      names(extraction) != "Data"
    )]) {
      surrogate_out[["Estimates"]][[extract_name]] <- extraction[[extract_name]]
    }

    surrogate_out[["Information"]][["Data"]] <- extraction[["Data"]]

    return(surrogate_out)
  } else if (identical(names(x), c("Estimates", "Information"))) {
    # Single-Group Code
    stop("This function is only meant for multi-group models.")
  }
}