#' Integrated Generalized Structured Component Analysis
#'
#' This R implementation of I-GSCA is based on the Matlab implementation in igsca_sim.m by Dr. Heungsun Hwang.
#'
#' In the example section, the specified model is based on the tutorial I-GSCA model associated with GSCA Pro \insertCite{hwangetal2023StructuralEquationModelingAMultidisciplinaryJournal}{cSEM}.
#'
#' **Note**: Here, we assume that there is only one unique loading per indicator.
#'
#'
#' @param .Z0 Data matrix of N cases (measurements) x J indicators with named
#'   columns, unstandardized.
#' @param .W0 Indicator matrix of weights: J indicators (rows) and their
#'   corresponding Gamma construct variables (columns).
#' @param .C0 Indicator matrix of loadings: J indicators (rows) and their
#'   corresponding Gamma construct variable (columns).
#' @param .B0 Square indicator matrix of path coefficients:
#'   from-construct-variable (rows) and to-construct-variable (columns). The
#'   order of Gamma construct variables should match the order in C0 and W0.
#' @param .con_type A vector that denotes whether each construct variable
#'   (columns in W0 and C0) is a common factor or composite. Its length should
#'   be equal to the number of columns of W0 and C0.
#' @param .indicator_type An indicator vector that indices whether a j indicator
#'   (rows of W0 and C0) corresponds to a common factor variable (1) or a
#'   composite variable (0). This vector is important for computing the
#'   uniqueness terms (D) because it zeros the entries for composite indicators.
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
    .X = args_default()$.X,
    .S = args_default()$.S,
    .csem_model = args_default()$.csem_model,
    .conv_criterion = args_default()$.conv_criterion,
    .GSCA_modes = args_default()$.GSCA_modes,
    .iter_max = args_default()$.iter_max,
    .starting_values = args_default()$.starting_values,
    .tolerance = args_default()$.tolerance,
    .Z0,
    .W0,
    .C0,
    .B0,
    .con_type,
    .indicator_type
  ) {
    ## Initialize Computational Variables -----------------------------------------------------
    n_case <- nrow(.Z0)
    n_indicators <- ncol(.Z0)
    n_constructs <- ncol(.W0)
    n_total_var <- n_indicators + n_constructs
    
    # Normalize data matrix to make normalized Z
    normalization_factor <- sqrt(n_case - 1)
    Z <- .Z0 / normalization_factor
    w_index <- which(c(.W0) == 1)
    b_index <- which(c(.B0) == 1)

    ### Initial Values ------------------------------------
    GSCA_starting_values <- calculateWeightsGSCA(
      .X = .X,
      .S = .S,
      .csem_model = .csem_model,
      .conv_criterion = .conv_criterion,
      .GSCA_modes = .GSCA_modes,
      .iter_max = 10,
      .tolerance = .tolerance,
      .starting_values = .starting_values
    )

    W <- t(GSCA_starting_values$W)
    B <- t(GSCA_starting_values$B)
    C <- GSCA_starting_values$C
    V <- cbind(diag(n_indicators), W)

    # Create initial values for U and D, using normalized data (Z) and construct scores (Eta) --------------


    # Create Initial Normalized Construct Scores Matrix
    Eta <- Z %*% W

    # Initalize unique loadings

    list_UD <- updateUD(
      D = diag(n_indicators),
      Eta_normed = Eta,
      .indicator_type = .indicator_type,
      n_constructs = n_constructs,
      n_indicators = n_indicators,
      n_case = n_case,
      Z_normed = Z
    )

    "U" <- list_UD[["U"]]
    "D" <- list_UD[["D"]]

    # if starting values are provided
    if (!is.null(.starting_values)) {
      W <- setStartingValues(.W = t(W), .starting_values = .starting_values) |>
        t()
    }

    # Define and Check Modes: 'Common factor', 'NCMP', and 'CCMP'
    modes <- ifelse(
      .csem_model$construct_type == "Composite",
      yes = "CCMP",
      no = .csem_model$construct_type
    )
    if (!is.null(.GSCA_modes)) {
      stopifnot(
        "Invalid .GSCA_modes selected. Only NCMP or CCMP are supported" = all(sapply(
          .GSCA_modes,
          function(x) x %in% c('NCMP', 'CCMP')
        ))
      )

      if (length(names(.GSCA_modes)) != 0) {
        stopifnot(
          "Invalid construct name listed in .GSCA_modes" = length(setdiff(
            names(.GSCA_modes),
            names(.csem_model$construct_type)
          )) ==
            0
        )
        modes[names(.GSCA_modes)] <- .GSCA_modes
      } else if (length(names(.GSCA_modes)) == 0) {
        stopifnot(
          "When passing a global setting to .GSCA_modes, only one choice may be passed" = length(
            .GSCA_modes
          ) ==
            1
        )

        # If only an unnamed element is given (NCMP or CCMP), then set all composites to .GSCA_modes
        if (length(.GSCA_modes) == 1) {
          modes <- ifelse(modes == "CCMP", yes = .GSCA_modes, no = modes)
        }
      }
    }
    # Set Loadings of canonical composites to 0 just in case it wasn't handled properly in the initial values
    C[which(modes == "CCMP"), ] <- 0
    # Create c_index which should only exist for nomological composites and common factors
    .C0[which(modes == "CCMP"), ] <- 0
    c_index <- which(c(.C0) == 1)

    ## Alternating Least Squares Algorithm -------------------------------------

    ### Optimization Preparation -----------------------------------------------
    # Set the initial estimates based on either the structural model or the loadings
    # if there's no structural model
    if (length(b_index) > 0) {
      est <- B[b_index]
    } else {
      # Because canonical composites do not estimate loadings, here the weights are generally used to determine convergenece, instead of loadings
      est <- W[w_index]
      # est <- C[c_index]
    }
    est0 <- est + 1
    it <- 0

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

      ### Update Weights --------------------------------------------------

      # X is the indicators Z with measurement error removed UD.
      X <- Z - (U %*% D)

      # WW is important for updating the Theta for common factors
      WW <-
        t(C) %*%
        solve((C %*% t(C) + diag(n_constructs) - (2 * B) + (B %*% t(B))))

      #### for-loop update per construct ---------------------------------------
      A <- cbind(C, B)

      # After each cycle, the Gamma, W and V matrices are updated
      for (eta_idx in seq_len(n_constructs)) {
        tot <- n_indicators + eta_idx
        windex_eta_idx <- (.W0[, eta_idx, drop = FALSE] == 1)
        X_eta_idx <- X[, windex_eta_idx, drop = FALSE]

        if (.con_type[eta_idx] == "Composite") {
          Theta <-
            updateCompositeTheta(
              W = W,
              A = A,
              X = X,
              V = V,
              eta_idx = eta_idx,
              tot = tot,
              n_constructs = n_constructs,
              n_total_var = n_total_var,
              windex_eta_idx = windex_eta_idx,
              .S = .S
            )
        } else if (.con_type[eta_idx] == "Common factor") {
          Theta <- WW[windex_eta_idx, eta_idx, drop = FALSE]
        } else {
          stop(".con_type should only either be `Composite` or `Common factor`")
        }

        normed_weights <- Theta / norm(X_eta_idx %*% Theta, "2")
        W[windex_eta_idx, eta_idx] <- normed_weights
        V[windex_eta_idx, tot] <- normed_weights
      }
      Eta <- X %*% W # Trying to save on compute time

      ### Update Loadings, Path Coefficients and Uniqueness Terms ----------
      list_CB <- updateCB(
        X = X,
        Eta = Eta,
        Lambda = C,
        B = B,
        n_indicators = n_indicators,
        .indicator_type = .indicator_type,
        n_constructs = n_constructs,
        n_case = n_case,
        lambda_index = c_index,
        b_index = b_index,
        .con_type = .con_type,
        modes = modes
      )

      list2env(list_CB[c("C", "B")], envir = environment())

      ## Uniqueness Scores and Loadings Update ------------------------------------
      list_UD <- updateUD(
        D = D,
        Eta_normed = Eta,
        .indicator_type = .indicator_type,
        n_constructs = n_constructs,
        n_indicators = n_indicators,
        n_case = n_case,
        Z_normed = Z
      )

      list2env(list_UD[c("U", "D")], envir = environment())

      if (length(b_index) > 0) {
        est <- B[b_index]
      } else {
        est <- W[w_index]
        # Because canonical composites do not estimate loadings, here the weights are generally used to determine convergenece, instead of loadings
        # est <- C[c_index]
      }
    }

    ## Output Formatting -------------------------------------------------------

    # Compute loadings for Canonical Composites
    if (any(modes %in% "CCMP")) {
      CCMP_C <- t(W) %*% .S * .csem_model$measurement
      C[names(modes)[modes == "CCMP"], ] <- CCMP_C[
        names(modes)[modes == "CCMP"],
        ,
        drop = FALSE
      ]
    }

    D_diag <- diag(D)
    names(D_diag) <- colnames(Z)

    # Get the standardized unique scores back from normalized U and zero out the unique scores for composite indicators
    U[, .indicator_type == "Composite"] <- 0
    Unique_scores <- U * normalization_factor
    colnames(Unique_scores) <- colnames(Z)

    return(
      list(
        "W" = t(W), # W is J X P, so W^T is P \times J. As shown in the examples ?csem, res$Estimates$Weight_estimates should be P \times J
        "C" = C, # C is P \times J. As shown in the examples ?csem, res$Estimates$Loading_estimates should be P \times J
        "B" = t(B), # B is From \times To; so t(B) is To \times From. As shown in the examples ?csem, res$Estimates$Path_estimates should be To \times From
        # Recall that Z0 is the original standardized data
        "Construct_scores" = (.Z0 - (Unique_scores %*% D)) %*% W,
        "Unique_loading_estimates" = D_diag,
        "Unique_scores" = Unique_scores,
        "Modes" = "gsca (igsca)",
        "Conv_status" = ifelse(it > .iter_max, FALSE, TRUE),
        "Iterations" = it
      )
    )
  }

#' Update Theta for Composite Variables
#' 
#' It is unintuitive that X is used here, seeing as how X = Z-UD; and we use X to update composite variables. 
#' However, from a non-computational point of view, it shouldn't matter because for the composite indicators, X_comp = Z_comp
#'
#' @param W Weights matrix
#' @param A Stacked matrix of loadings and path coefficients \eqn{\left[\Lambda \mid B \right]}
#' @param V Stacked matrix of identity matrix and weights \eqn{\left[I \mid W \right]}
#' @param X The matrix X is equal to \eqn{Z - UD}
#' @param windex_eta_idx Index of weights related to the indicators for the construct of interest
#' @param n_total_var Number of indicators and constructs
#' @param tot Index dependent on which construct variable we are examining
#' @param n_constructs Number of constructs
#' @param eta_idx Index of which construct we are examining
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
    windex_eta_idx,
    n_total_var,
    tot,
    n_constructs,
    eta_idx,
    .S = args_default()$.S
  ) {
    e <- matrix(0, nrow = 1, ncol = n_total_var)
    e[tot] <- 1
    H1 <- diag(n_total_var)
    H2 <- diag(n_constructs)
    H1[tot, tot] <- 0
    H2[eta_idx, eta_idx] <- 0
    Delta <- (W %*% H2 %*% A) - (V %*% H1)

    beta <- e - A[eta_idx, , drop = FALSE]

    Theta <- MASS::ginv(
      as.numeric(beta %*% t(beta)) * .S[windex_eta_idx, windex_eta_idx, drop = FALSE]
    ) %*%
      t(beta %*% t(Delta) %*% .S[, windex_eta_idx, drop = FALSE])

    return(Theta)

    # Kronecker Method
    # vecZDelta <- c(X %*% Delta) 
    # XI <- kronecker(t(beta), X)
    # XI <- XI[, windex_eta_idx]
    # XI <- kroneckerC(t(beta), X, which(windex_eta_idx))
    # Theta <- solve((t(XI) %*% XI), t(XI)) %*% vecZDelta
  }

#' Update Loadings and Path-Coefficients
#'
#' @inheritParams igsca
#' @param X Indicators with measurement error removed
#' @param Eta Construct Scores
#' @param Lambda Loadings matrix
#' @param B Path coefficients matrix
#' @param n_indicators Number of indicators
#' @param n_constructs Number of oncstructs
#' @param lambda_index Index of loadings
#' @param b_index Index of Path Coefficients
#' @param n_case Number of Cases
#' @param .indicator_type Vector of whether each indicator corresponds to a common factor or composite
#' @param modes Named vector of whether the construct is a Common factor, nomological composite or canonical composite.
#' @importFrom MASS ginv
#' @return List of matrices:
#'
#' * (1) Estimated Loadings matrix (C)
#' * (2) Estimated Path Coefficients matrix (B)
#'
updateCB <-
  function(
    X,
    Eta,
    Lambda,
    B,
    .indicator_type,
    n_indicators,
    lambda_index,
    n_constructs,
    b_index,
    n_case,
    .con_type,
    modes
  ) {
    # Loading Update ----------------------------------------------------------

    # Kronecker bypass
    # browser()
    vars_cf_ncmp <- names(modes)[modes %in% c("Common factor", "NCMP")]
    # cov_eta_indicators <- t(Eta) %*% X # Interestingly, this is not the same as cor(Eta, X)?
    cov_eta_indicators <- crossprod(Eta, X) # TODO: Check that this is equivalent to t(Eta) %*% X
    # cor_eta <- t(Eta) %*% Eta 
    cor_eta <- crossprod(Eta) # TODO: Check that this is equivalent to t(Eta) %*% Eta
    

    dep_vars <- (colSums(Lambda[vars_cf_ncmp, , drop = FALSE]) != 0) |> 
        which() |> 
        names()
    # This approach assumes that every factor/NCMP loads onto one indicator: no cross-loadings
    loadings <- lapply(dep_vars, function(y) {
      x <- (rowSums(Lambda[vars_cf_ncmp, y, drop = FALSE]) != 0) |> 
          which() |> 
          names()
      coef <- MASS::ginv(cor_eta[x, x, drop = FALSE]) %*% cov_eta_indicators[x, y, drop = FALSE]
    })
    # A future approach should consider avoiding c_index and using explicit names, for safety.
    Lambda[lambda_index] <- unlist(loadings, use.names =  FALSE)

    # Kronecker Approach and Assumes All Composites are Nomological
    # t1 <- c(X)
    # M1 <- kroneckerC(diag(n_indicators), Eta, c_index)
    # C[c_index] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)

    # Path Coefficients Update ------------------------------------------------
    vars_endo <- colnames(B)[colSums(B) != 0]
    beta <- lapply(vars_endo, function(y) {
      x <- (rowSums(B[, y, drop = FALSE]) != 0) |> 
        which() |> 
        names()
      coef <- MASS::ginv(cor_eta[x, x, drop = FALSE]) %*%
        cor_eta[x, y, drop = FALSE]
    })
    B[b_index] <- unlist(beta, use.names = FALSE)

    return(
      list(
        "C" = Lambda,
        "B" = B
      )
    )
  }


#' Update unique scores and unique loadings
#' 
#' Intended to be used within the alternating least squares algorithm for either GSCA_M or IGSCA. Assumes that the construct scores and data are normalized.
#'
#' @param D Unique loadings
#' @param Eta_normed Normalized data
#' @param Z_normed Normalized data
#' @param n_constructs Number of constructs
#' @param n_case Number of cases
#' @param n_indicators Number of indicators
#'
#' @returns List of 2 elements, normalized unique scores (`U`) and normalized unique loadings (`D`)
#' 
#' @inheritParams igsca
#' @inheritParams csem_arguments
#'
#' @export
updateUD <- function(D, Eta_normed, .indicator_type, n_constructs, n_case, n_indicators, Z_normed) {

  # Update Unique Scores---------------------
  # U <- tryCatch(
  #   {
  #     # TODO: Need to understand why these two approaches differ from one-another. The old approach may have assumed something about the data or the product that I'm not aware of...
  #     # Implementation based on the updated MATLAB code provided by Heungsun in
  #     # private communication to deal with big dataset (20.10.2021)
  #     svd_etaprod <- svd(t(Eta_normed) %*% Eta_normed)
  #     gd2 <- diag(svd_etaprod$d)
  #     gv <- svd_etaprod$v
  #     GU <- Eta_normed %*% gv %*% solve(sqrt(gd2))
  #     M3 <- Z_normed %*% D - GU %*% (t(GU) %*% Z_normed) %*% t(D)

  #     if (n_case > n_indicators) {
  #       svd_M3prod <- svd(t(M3) %*% M3)
  #       d2 <- diag(svd_M3prod$d)
  #       v <- svd_M3prod$v

  #       u <- M3 %*% v %*% solve(sqrt(d2))
  #       U <- u[, 1:n_indicators, drop = FALSE] %*% t(v)
  #     } else {
  #       svd_M3 = svd(M3)
  #       u <- svd_M3$u
  #       v <- svd_M3$v
  #       U <- u[, 1:n_indicators, drop = FALSE] %*% t(v)
  #     }
  #     U
  #   },
  #   error = function(e) {
  #     # Old method based on Hwang et al. (2017)
  #     Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
  #       (n_constructs + 1):n_case,
  #       drop = FALSE
  #     ]
  #     svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
  #     Utilde <- svd_mx$v %*% t(svd_mx$u)
  #     U <- Eta_Q2 %*% Utilde
  #     return(U)
  #   }
  # )

  # Claude's Approach --- Efficient projection method: O(NJ^2 + NPJ) instead of O(N^2 J)
  # Uses thin QR (N x P, not N x N) to avoid forming the null space basis.
  # Mathematically equivalent to the Hwang et al. (2017) complete QR approach
  # because P_perp Z D = Gamma_perp (Gamma_perp' Z D), preserving SVD structure.
  # See dev/igsca/benchmarking_R/unique_scores_optimization.md for full proof.
  # TODO: Investigate Claude's approach more deeply.
  # Q_thin <- qr.Q(qr(Eta_normed))
  # ZD <- Z_normed %*% D
  # M_proj <- ZD - Q_thin %*% crossprod(Q_thin, ZD)
  # svd_mx <- svd(M_proj)
  # U <- svd_mx$u %*% t(svd_mx$v)

  # Old method based on Hwang et al. (2017) — O(N^2) memory and computation
  Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
    (n_constructs + 1):n_case,
    drop = FALSE
  ]
  svd_mx <- svd(D %*% t(Z_normed) %*% Eta_Q2)
  Utilde <- svd_mx$v %*% t(svd_mx$u)
  U <- Eta_Q2 %*% Utilde
  
  # U[, .indicator_type == "Composite"] <- 0

  # Update Unique Loadings

  # D <- diag(diag(t(U) %*% Z_normed))
  D <- diag(diag(crossprod(U, Z_normed)))
  D[.indicator_type == "Composite", .indicator_type == "Composite"] <- 0

  # Return output
  return(list("U" = U, "D" = D))
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