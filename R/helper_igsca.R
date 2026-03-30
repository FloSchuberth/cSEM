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
    normalization_factor <- sqrt(n_case - 1)

    w_index <- which(c(.W0) == 1)
    b_index <- which(c(.B0) == 1)

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

    ### Starting Values ------------------------------------
    prepared_for_ALS <- initializeAlsEstimates(
      .X = .X,
      .S = .S,
      .csem_model = .csem_model,
      .conv_criterion = .conv_criterion,
      .GSCA_modes = .GSCA_modes,
      .iter_max = .iter_max,
      .starting_values = .starting_values,
      .tolerance = .tolerance,
      Z0 = .Z0,
      W0 = .W0,
      B0 = .B0,
      n_indicators = n_indicators,
      n_case = n_case,
      n_constructs = n_constructs,
      .indicator_type = .indicator_type,
      normalization_factor = normalization_factor
    )

    list2env(
      prepared_for_ALS[c("W", "C", "B", "V", "Z", "D", "U", "Gamma")],
      envir = environment()
    )

    # if starting values are provided
    if (!is.null(.starting_values)) {
      W <- setStartingValues(.W = t(W), .starting_values = .starting_values) |>
        t()
    }

    # Define and Check Modes: 'Common factor', 'NCMP', and 'CCMP'
    modes <- ifelse(.csem_model$construct_type == "Composite", yes = "CCMP", no = .csem_model$construct_type)
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
        stopifnot("When passing a global setting to .GSCA_modes, only one choice may be passed" = length(.GSCA_modes) == 1)

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
        t(C) %*% solve((C %*% t(C) + diag(n_constructs) - (2 * B) + (B %*% t(B))))

      #### for-loop Updating of Constructs ---------------------------------------
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
      Gamma <- X %*% W # Trying to save on compute time

      ### Update Loadings, Path Coefficients and Uniqueness Terms ----------
      updated_C_B_D_U <- updateCBDU(
        X = X,
        Gamma = Gamma,
        C = C,
        B = B,
        D = D,
        Z = Z,
        n_indicators = n_indicators,
        .indicator_type = .indicator_type,
        n_constructs = n_constructs,
        n_case = n_case,
        c_index = c_index,
        b_index = b_index,
        .con_type = .con_type,
        modes = modes
      )

      list2env(updated_C_B_D_U[c("D", "U", "C", "B")], envir = environment())

      #### Update estimates for while-loop condition testing -----------------------
      if (length(b_index) > 0) {
        est <- B[b_index]
      } else {
        # Because canonical composites do not estimate loadings, here the weights are generally used to determine convergenece, instead of loadings
        est <- W[w_index]
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
  
    # Get the standardized unique scores back from normalized U
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


#' Prepare for ALS Algorithm
#'
#' Internal I-GSCA function
#'
#' @inheritParams igsca
#' @inheritParams csem_arguments
#' @param n_indicators Number of indicators
#' @param n_case Number of measurements
#' @param n_constructs Number of constructs
#' @param normalization_factor Factor to normalize the data by
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
  .X                           = args_default()$.X,
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .GSCA_modes                  = args_default()$.GSCA_modes,
  .iter_max                    = args_default()$.iter_max,
  .starting_values             = args_default()$.starting_values,
  .tolerance                   = args_default()$.tolerance,
  Z0,
  W0,
  B0,
  n_indicators,
  n_case,
  n_constructs,
  .indicator_type,
  normalization_factor
) {
  
  # Initial weights, path-coefficients and loadings using GSCA
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
  
  # Create initial values for U and D, using normalized data (Z) and construct scores (Gamma) --------------

  # Normalize data matrix to make normalized Z
  Z <- Z0 / normalization_factor

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

  return(list(
    "W" = W,
    "C" = C,
    "B" = B,
    "V" = cbind(diag(n_indicators), W),
    "Z" = Z,
    "D" = list_UD[["D"]],
    "U" = list_UD[["U"]],
    "Gamma" = Eta
  ))
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

#' Update Loadings, Path-Coefficients and Uniqueness Terms After Updating Latent Variables
#'
#' @inheritParams initializeAlsEstimates
#' @inheritParams igsca
#' @param X Pseudo-weights for composites
#' @param Gamma Construct Scores
#' @param c_index Index of loadings
#' @param b_index Index of Path Coefficients
#' @param n_case Number of Cases
#' @param .indicator_type Vector of whether each indicator corresponds to a common factor or composite
#' @param modes Named vector of whether the construct is a Common factor, nomological composite or canonical composite.
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
    .indicator_type,
    n_indicators,
    c_index,
    n_constructs,
    b_index,
    n_case,
    .con_type,
    modes
  ) {
    ## Loading Update ----------------------------------------------------------

    # Kronecker Approach and Assumes All Composites are Nomological
    # t1 <- c(X)
    # M1 <- kroneckerC(diag(n_indicators), Gamma, c_index)
    # C[c_index] <- MASS::ginv(t(M1) %*% M1) %*% (t(M1) %*% t1)
    
    # Kronecker bypass
    # browser()
    vars_cf_ncmp <- names(modes)[modes %in% c("Common factor", "NCMP")]
    cov_gamma_indicators <- t(Gamma) %*% Z
    vcv_gamma <- t(Gamma) %*% Gamma
    dep_vars <- (colSums(C[vars_cf_ncmp, , drop = FALSE]) != 0) |> 
        which() |> 
        names()
    # This approach assumes that every factor/NCMP loads onto one indicator: no cross-loadings
    loadings <- lapply(dep_vars, function(y) {
      x <- (rowSums(C[vars_cf_ncmp, y, drop = FALSE]) != 0) |> 
          which() |> 
          names()
      coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*% cov_gamma_indicators[x, y, drop = FALSE]
    })
    # A future approach should consider avoiding c_index and using explicit names, for safety.
    C[c_index] <- unlist(loadings, use.names =  FALSE)

    # Path Coefficients Update ------------------------------------------------
    vars_endo <- colnames(B)[colSums(B) != 0]

    beta <- lapply(vars_endo, function(y) {
      x <- (rowSums(B[, y, drop = FALSE]) != 0) |> 
        which() |> 
        names()
      coef <- MASS::ginv(vcv_gamma[x, x, drop = FALSE]) %*%
        vcv_gamma[x, y, drop = FALSE]
    })
    B[b_index] <- unlist(beta, use.names = FALSE)

    ## Uniqueness Component Update ---------------------------------------------

    list_UD <- updateUD(
      D = D,
      Eta_normed = Gamma,
      .indicator_type = .indicator_type,
      n_constructs = n_constructs,
      n_indicators = n_indicators,
      n_case = n_case,
      Z_normed = Z
    )

    return(
      list(
        "C" = C,
        "B" = B,
        "D" = list_UD[["D"]],
        "U" = list_UD[["U"]]
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
#'
#' @returns List of 2 elements, normalized unique scores (`U`) and normalized unique loadings (`D`)
#' 
#' @inheritParams igsca
#' @inheritParams csem_arguments
#'
#' @export
updateUD <- function(D, Eta_normed, .indicator_type, n_constructs, n_case, n_indicators, Z_normed) {


  # Update Unique Scores
  Eta_Q2 <- qr.Q(qr(Eta_normed), complete = TRUE)[,
    (n_constructs + 1):n_case,
    drop = FALSE
  ]

  # svd between R and Matlab by Ahmed Fasih on February 1/2017
  # https://stackoverflow.com/a/41972818

  mx <- D %*% t(Z_normed) %*% Eta_Q2

  U <- tryCatch(
    {
      svd_mx <- svd(mx, nu = nrow(mx), nv = ncol(mx))
      u <- svd_mx$u
      v <- svd_mx$v
      Utilde <- v[, 1:n_indicators, drop = FALSE] %*% t(u)
      U <- Eta_Q2 %*% Utilde
    },
    error = function(e) {
      svd_mx <- svd(mx)
      Utilde <- svd_mx$v %*% t(svd_mx$u)
      U <- Eta_Q2 %*% Utilde
      return(U)
    }
  )
  
  U[, .indicator_type == "Composite"] <- 0

  # Update Unique Loadings
  D <- diag(diag(t(U) %*% Z_normed))
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