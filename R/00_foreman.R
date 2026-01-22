#' Internal: Composite-based SEM
#'
#' The central hub of the \pkg{cSEM} package. It acts like a
#' foreman by collecting all (estimation) tasks, distributing them to lower
#' level package functions, and eventually recollecting all of their results.
#' It is called by [csem()] to manage the actual calculations.
#' It may be called directly by the user, however, in most cases it will likely
#' be more convenient to use [csem()] instead.
#'
#' @usage foreman(
#'   .data                        = args_default()$.data,
#'   .model                       = args_default()$.model,
#'   .approach_cor_robust         = args_default()$.approach_cor_robust,
#'   .approach_nl                 = args_default()$.approach_nl,
#'   .approach_paths              = args_default()$.approach_paths,
#'   .approach_weights            = args_default()$.approach_weights,
#'   .conv_criterion              = args_default()$.conv_criterion,
#'   .disattenuate                = args_default()$.disattenuate,
#'   .dominant_indicators         = args_default()$.dominant_indicators,
#'   .estimate_structural         = args_default()$.estimate_structural,
#'   .GSCA_modes                  = args_default()$.GSCA_modes,
#'   .id                          = args_default()$.id,
#'   .instruments                 = args_default()$.instruments,
#'   .iter_max                    = args_default()$.iter_max,
#'   .normality                   = args_default()$.normality,
#'   .PLS_approach_cf             = args_default()$.PLS_approach_cf,
#'   .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
#'   .PLS_modes                   = args_default()$.PLS_modes,
#'   .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
#'   .reliabilities               = args_default()$.reliabilities,
#'   .starting_values             = args_default()$.starting_values,
#'   .tolerance                   = args_default()$.tolerance
#'   )
#'
#' @inheritParams csem_arguments
#'
#' @inherit csem_results return
#'
#' @seealso [csem], [cSEMResults]
#'
#' @keywords internal

foreman <- function(
  .data = args_default()$.data,
  .model = args_default()$.model,
  .approach_cor_robust = args_default()$.approach_cor_robust,
  .approach_nl = args_default()$.approach_nl,
  .approach_paths = args_default()$.approach_paths,
  .approach_weights = args_default()$.approach_weights,
  .conv_criterion = args_default()$.conv_criterion,
  .disattenuate = args_default()$.disattenuate,
  .dominant_indicators = args_default()$.dominant_indicators,
  .estimate_structural = args_default()$.estimate_structural,
  .GSCA_modes = args_default()$.GSCA_modes,
  .id = args_default()$.id,
  .instruments = args_default()$.instruments,
  .iter_max = args_default()$.iter_max,
  .normality = args_default()$.normality,
  .PLS_approach_cf = args_default()$.PLS_approach_cf,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner = args_default()$.PLS_weight_scheme_inner,
  .reliabilities = args_default()$.reliabilities,
  .starting_values = args_default()$.starting_values,
  .tolerance = args_default()$.tolerance
) {
  args_used <- c(as.list(environment(), all.names = TRUE))

  ### Preprocessing ============================================================
  ## Parse and order model to "cSEMModel" list
  csem_model <- parseModel(.model, .instruments = .instruments)

  ## Prepare, check, and clean data (a data.frame)
  X_cleaned <- processData(
    .data = .data,
    .model = csem_model,
    .instruments = NULL
  )

  ### Computation ==============================================================
  ## Calculate empirical indicator covariance/correlation matrix
  Cor <- calculateIndicatorCor(
    .X_cleaned = X_cleaned,
    .approach_cor_robust = .approach_cor_robust
  )

  # Extract the correlation matrix
  S <- Cor$S

  ## Standardize
  X <- scale(data.matrix(X_cleaned), center = TRUE, scale = TRUE)

  ## Calculate weights
  if (.approach_weights != "GSCA") {
    if (.approach_weights == "PLS-PM") {
      W <- calculateWeightsPLS(
        .S = S,
        .csem_model = csem_model,
        .iter_max = .iter_max,
        .PLS_modes = .PLS_modes,
        # Arguments passed on to calculateInnerWeightsPLS
        .PLS_ignore_structural_model = .PLS_ignore_structural_model,
        .PLS_weight_scheme_inner = .PLS_weight_scheme_inner,
        # Arguments passed on to calcuateOuterWeightsPLS
        .data = X,
        # Arguments passed to checkConvergence
        .conv_criterion = .conv_criterion,
        .tolerance = .tolerance,
        # starting values
        .starting_values = .starting_values
      )
    } else if (
      .approach_weights %in%
        c("SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR")
    ) {
      W <- calculateWeightsKettenring(
        .S = S,
        .csem_model = csem_model,
        .approach_gcca = .approach_weights
      )
    } else if (.approach_weights == "unit") {
      W <- calculateWeightsUnit(
        .S = S,
        .csem_model = csem_model,
        .starting_values = .starting_values
      )
    } else if (.approach_weights %in% c("bartlett", "regression")) {
      # Note:  1. "bartlett" and "regression" weights are calculated later in the
      #        calculateReliabilities() function. Here only placeholders for
      #        the weights are set in order to keep the list structure of W.

      W <- list(
        "W" = csem_model$measurement,
        "E" = NULL,
        "Modes" = NULL,
        "Conv_status" = NULL,
        "Iterations" = 0
      )
    } else if (.approach_weights == "PCA") {
      W <- calculateWeightsPCA(
        .S = S,
        .csem_model = csem_model
      )
    }

    ## Dominant indicators:
    if (!is.null(.dominant_indicators)) {
      W$W <- setDominantIndicator(
        .W = W$W,
        .dominant_indicators = .dominant_indicators,
        .S = S,
        .approach_weights = .approach_weights
      )
    }

    LambdaQ2W <- calculateReliabilities(
      .X = X,
      .S = S,
      .W = W,
      .approach_weights = .approach_weights,
      .csem_model = csem_model,
      .disattenuate = .disattenuate,
      .PLS_approach_cf = .PLS_approach_cf,
      .reliabilities = .reliabilities
    )

    Weights <- LambdaQ2W$W
    Lambda <- LambdaQ2W$Lambda
    Q <- sqrt(LambdaQ2W$Q2)

    # Overwrite the .disattenuate argument
    # The .disattenuate argument can changed in the calculate Reliabilities function
    # in the case that reliabilities are
    # given by the user or that the two stage approach is used.
    args_used$.disattenuate <- LambdaQ2W$.disattenuate

    ## Calculate measurement error correlation
    # Compute theta
    Theta <- S - t(Lambda) %*% Lambda

    ## Calculate proxies/scores
    H <- X %*% t(Weights)

    ## Calculate proxy covariance matrix and assumes that weights are for constructs with measurement error
    C <- calculateCompositeVCV(.S = S, .W = Weights)

    ## Calculate construct correlation matrix and corrects for attenuation
    P <- calculateConstructVCV(.C = C, .Q = Q)

    ## Estimate structural coef
    if (sum(rowSums(csem_model$structural)) == 0) {
      .estimate_structural <- FALSE
    }

    if (.estimate_structural) {
      estim_results <- estimatePath(
        .approach_nl = .approach_nl,
        .approach_paths = .approach_paths,
        .csem_model = csem_model,
        .H = H,
        .normality = .normality,
        .P = P,
        .Q = Q
      )
    } else {
      estim_results <- NULL
    }
    # To ensure compatibility with GSCA
    Unique_loading_estimates <- NULL
    Unique_scores <- NULL
  } else if (.approach_weights == "GSCA") {
    if (csem_model$model_type == "Nonlinear") {
      stop2(
        "cSEM currently does not support GSCA and GSCAm for models containing nonlinear terms."
      )
    }
    if (all(csem_model$construct_type == "Common factor")) {
      if (isTRUE(.disattenuate)) {
        W <- calculateWeightsGSCAm(
          .X = X,
          .csem_model = csem_model,
          .conv_criterion = .conv_criterion,
          .iter_max = .iter_max,
          .tolerance = .tolerance,
          .starting_values = .starting_values
        )
      } else {
        stop2(
          "The csem argument .disattenuate should be TRUE for GSCA_M to fit with common factors"
        )
      }
    } else if (all(csem_model$construct_type == "Composite")) {
      W <- calculateWeightsGSCA(
        .X = X,
        .S = S,
        .csem_model = csem_model,
        .conv_criterion = .conv_criterion,
        .iter_max = .iter_max,
        .tolerance = .tolerance,
        .starting_values = .starting_values
      )
    } else if (
      any(csem_model$construct_type == "Composite") &
        any(csem_model$construct_type == "Common factor")
    ) {
      if (isTRUE(.disattenuate)) {
        # IGSCA Algorithm for a mixture of Common factor and Composite constructs
        W <- calculateWeightsIGSCA(
          .data = X,
          .csem_model = csem_model,
          .tolerance = .tolerance,
          .iter_max = .iter_max,
          .dominant_indicators = .dominant_indicators,
          .conv_criterion = .conv_criterion,
          .S = S
        )
      } else {
        stop2(
          "The csem argument .disattenuate should be TRUE for IGSCA to fit with both common factors and composites."
        )
      }
    }

    if (!is.null(.dominant_indicators)) {
      flipped_W_Eta_L_B <- setDominantIndicator(
        .W = W$W,
        .dominant_indicators = .dominant_indicators,
        .approach_weights = "GSCA",
        .X = X,
        .L = W$C,
        .B = W$B,
        Construct_scores = W$Construct_scores,
        .csem_model = csem_model
      )
      W$W <- flipped_W_Eta_L_B$W
      W$Construct_scores <- flipped_W_Eta_L_B$Eta
      W$C <- flipped_W_Eta_L_B$L
      W$B <- flipped_W_Eta_L_B$B
    }

    LambdaQ2W <- calculateReliabilities(
      .X = X,
      .S = S,
      .W = W,
      .approach_weights = .approach_weights,
      .csem_model = csem_model,
      .disattenuate = .disattenuate,
      .PLS_approach_cf = .PLS_approach_cf,
      .reliabilities = .reliabilities
    )

    Q <- sqrt(LambdaQ2W$Q2)

    Theta <- NULL
    # Calculate Construct Scores
    H <- W$Construct_scores
    ## Calculate proxy covariance matrix
    C <- cov(H)
    ## Calculate construct correlation matrix
    P <- cor(H)

    ## Estimate structural coef
    if (sum(rowSums(csem_model$structural)) == 0) {
      .estimate_structural <- FALSE
    }

    if (.estimate_structural) {
      estim_results <- estimatePath(
        .approach_nl = .approach_nl,
        .approach_paths = .approach_paths,
        .csem_model = csem_model,
        .H = H,
        .normality = .normality,
        .P = P,
        .Q = Q
      )
      # TODO: Discuss this section with Flo
      estim_results$Path_estimates <- W$B
      
      # estim_results$R2 <- NULL # FIXME: Required by calculateGoF() and calculatef2
      estim_results$R2adj <- NULL
      estim_results$VIF <- NULL
      estim_results$SE <- NA
    } else {
      estim_results <- NULL
    }

    # Format Output
    Weights <- W$W
    Lambda <- W$C
    Unique_loading_estimates <- W$Unique_loading_estimates
    Unique_scores <- W$Unique_scores
  }

  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates" = list(
      "Path_estimates" = if (.estimate_structural) {
        estim_results$Path_estimates
      } else {
        estim_results
      },
      "Loading_estimates" = Lambda,
      "Weight_estimates" = Weights,
      "Inner_weight_estimates" = W$E,
      "Residual_correlation" = Theta,
      "Construct_scores" = H,
      "Indicator_VCV" = S,
      "Proxy_VCV" = C,
      "Construct_VCV" = P,
      "Unique_loading_estimates" = Unique_loading_estimates,
      "Unique_scores" = Unique_scores,
      "Reliabilities" = Q^2,
      "R2" = if (.estimate_structural) {
        estim_results$R2
      } else {
        estim_results
      },
      "R2adj" = if (.estimate_structural) {
        estim_results$R2adj
      } else {
        estim_results
      },
      "VIF" = if (.estimate_structural) {
        estim_results$VIF
      } else {
        estim_results
      },
      "SE" = if (.estimate_structural) {
        estim_results$SE
      } else {
        estim_results
      }
    ),
    "Information" = list(
      "Data" = X,
      "Model" = csem_model,
      "Arguments" = args_used,
      "Type_of_indicator_correlation" = Cor$cor_type,
      "Threshold_parameter_estimates" = Cor$thres_est,
      "Weight_info" = list(
        "Modes" = W$Modes,
        "Number_iterations" = W$Iterations,
        "Convergence_status" = W$Conv_status
      ),
      "Approach_2ndorder" = NA
    )
  )

  class(out) <- c("cSEMResults", "cSEMResults_default")
  invisible(out)

  ### For maintenance: ---------------------------------------------------------
  # X (N x K)   := Matrix of indicator values (=data)
  # S (K x K)   := Empirical indicator covariance/correlation matrix
  # W$W (J x K) := (Block-)diagonal matrix of weights with the same dimension as
  #                .modellist$measurement
  # H (N x J)   := Matrix of proxy values for the J constructs
  # Lambda (K x J) := Blockdiagonal matrix of factor loadings.
  # Q (1 x J)   := Vector of estimated construct-proxy correlations: Q := R(eta; eta_bar)
  # C (J x J)   := Proxy covariance matrix: R(eta_bar, eta_bar)
  # P (J x J)   := Construct correlation matrix (possibly disattenuated)
  }
