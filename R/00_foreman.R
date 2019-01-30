#' Composite based SEM
#'
#' The function is the central hub of the `cSEM`package. It acts like a 
#' foreman by collecting all (estimation) tasks, distributing them to lower 
#' level package functions, and eventually recollecting all of their results. 
#' It is called by [csem()] to manage the actual calculations.
#' It may be called directly by the user, however, in most cases it will likely
#' be more convenient to use [csem()] or [cca()] instead.
#'
#' More details here (TODO).
#' 
#' @usage foreman(
#'     .data                        = args_default()$.data,
#'     .model                       = args_default()$.model,
#'     .approach_cor_robust         = args_default()$.approach_cor_robust,
#'     .approach_nl                 = args_default()$.approach_nl,
#'     .approach_paths              = args_default()$.approach_paths,
#'     .approach_weights            = args_default()$.approach_weights,
#'     .disattenuate                = args_default()$.disattenuate,
#'     .dominant_indicators         = args_default()$.dominant_indicators,
#'     .estimate_structural         = args_default()$.estimate_structural,
#'     .iter_max                    = args_default()$.iter_max,
#'     .normality                   = args_default()$.normality,
#'     .PLS_approach_cf             = args_default()$.PLS_approach_cf,
#'     .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
#'     .PLS_modes                   = args_default()$.PLS_modes,
#'     .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
#'     .reliabilities               = args_default()$.reliabilities,
#'     .tolerance                   = args_default()$.tolerance
#'     )
#'
#' @inheritParams csem_arguments
#'
#' @inherit csem_results return
#'
#' @seealso [csem], [cca], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export
#'

foreman <- function(
  .data                        = args_default()$.data,
  .model                       = args_default()$.model,
  .approach_cor_robust         = args_default()$.approach_cor_robust,
  .approach_nl                 = args_default()$.approach_nl,
  .approach_paths              = args_default()$.approach_paths,
  .approach_weights            = args_default()$.approach_weights,
  .conv_criterion              = args_default()$.conv_criterion,
  .disattenuate                = args_default()$.disattenuate,
  .dominant_indicators         = args_default()$.dominant_indicators,
  .estimate_structural         = args_default()$.estimate_structural,
  .id                          = args_default()$.id,
  .iter_max                    = args_default()$.iter_max,
  .normality                   = args_default()$.normality,
  .PLS_approach_cf             = args_default()$.PLS_approach_cf,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes                   = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
  .reliabilities               = args_default()$.reliabilities,
  .tolerance                   = args_default()$.tolerance
  ) {

  ### Preprocessing ============================================================
  ## Parse and order model to "cSEMModel" list
  csem_model <- parseModel(.model)

  ## Prepare, check, and clean data (a data.frame)
  X_cleaned <- processData(.data = .data, .model = csem_model) 
  
  ### Computation ==============================================================
  ## Calculate empirical indicator covariance/correlation matrix
  Cor <- calculateIndicatorCor(.X_cleaned = X_cleaned, 
                               .approach_cor_robust = .approach_cor_robust)
  
  # Extract the correlation matrix
  S <- Cor$S
  
  ## Standardize
  X <- scale(data.matrix(X_cleaned))
  
  ## Calculate weights
  if(.approach_weights == "PLS-PM") {
    W <- calculateWeightsPLS(
      .S                        = S,
      .csem_model               = csem_model,
      .iter_max                 = .iter_max,
      .PLS_modes                = .PLS_modes,
      .tolerance                = .tolerance,
      # Arguments passed on to calculateInnerWeightsPLS
      .PLS_ignore_structural_model  = .PLS_ignore_structural_model,
      .PLS_weight_scheme_inner      = .PLS_weight_scheme_inner,
      # Arguments passed to checkConvergence
      .conv_criterion           = .conv_criterion
    )
  } else if(.approach_weights %in% c("SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR")) {
    W <- calculateWeightsKettenring(
      .S                        = S,
      .csem_model               = csem_model,
      .approach                 = .approach_weights
    )
  } else if(.approach_weights == "GSCA") {
    W <- calculateWeightsGSCA(
      .S                        = S,
      .csem_model               = csem_model,
      .conv_criterion           = .conv_criterion,
      .iter_max                 = .iter_max,
      .tolerance                = .tolerance
    )
  } else if(.approach_weights == "fixed") {
    W <- calculateWeightsFixed(
      .data                     = NULL,
      .model                    = csem_model
    )
  } else if(.approach_weights == "unit") {
    W <- calculateWeightsUnit(
      .S                        = S,
      .csem_model               = csem_model
    )
  }

  ## Dominant indicators:
  # Use the dominant indicators approach (Henseler et al. (2016)) for PLS-PM. 
  # Perhaps this applicable to weights obtained from algorithms other than PLS-PM.
  # Currently only PLS-PM is supported.
  
  if(!is.null(.dominant_indicators)) {
    if(.approach_weights %in% c("PLS-PM",'MAXVAR','MINVAR','SUMCORR','SSQCORR','GENVAR')) {
      
      ## Check construct names:
      # Do all construct names in .dominant_indicators match the construct
      # names used in the model?
      tmp <- setdiff(names(.dominant_indicators), rownames(W$W))
      
      if(length(tmp) != 0) {
        stop("Construct name(s): ", paste0("`", tmp, "`", collapse = ", "), 
             " provided to `.dominant_indicators`", 
             ifelse(length(tmp) == 1, " is", " are"), " unknown.", call. = FALSE)
      }
      
      ## Check indicators
      # Do all indicators names in .dominant_indicators match the indicator
      # names used in the model?
      tmp <- setdiff(.dominant_indicators, colnames(W$W))
      
      if(length(tmp) != 0) {
        stop("Indicator name(s): ", paste0("`", tmp, "`", collapse = ", "), 
             " provided to `.dominant_indicators`", 
             ifelse(length(tmp) == 1, " is", " are"), " unknown.", call. = FALSE)
      }

      for(i in names(.dominant_indicators)) {
        W$W[i, ] = W$W[i, ] * sign(W$W[i, .dominant_indicators[i]])
      }
    } # END if PLS-PM
  } # END if 
  
  ## Calculate proxies/scores
  H <- calculateComposites(
    .X          = X,
    .W          = W$W
  )

  ## Calculate PLSc-type correction factors if no or only a subset of reliabilities
  # are given and disattenuation is requested. Otherwise use only the reliabilities
  # to disattenuate both weights and proxy correlations to obtain consistent
  # estimators for the loadings, cross-loadings and path coefficients.

  if(.approach_weights == "PLS-PM" & (is.null(.reliabilities) | length(.reliabilities) != ncol(H))  & .disattenuate == TRUE) {
    correction_factors <- calculateCorrectionFactors(
      .S               = S,
      .W               = W$W,
      .modes           = W$Modes,
      .csem_model      = csem_model,
      .PLS_approach_cf = .PLS_approach_cf)
  } else {
    correction_factors <- NULL
  }
  
  ## Calculate Q's (correlation between construct and proxy)
  # Note: Q_i^2 := R^2(eta_i; eta_bar_i) is also called the reliability coefficient
  # rho_A in Dijkstra & Henseler (2015) - Consistent partial least squares path modeling
  
  if(.approach_weights == "GSCA") {
    Q <- rep(1, nrow(csem_model$structural))
  } else {
    Q <- calculateCompositeConstructCV(
      .W                  = W$W,
      .csem_model         = csem_model,
      .modes              = W$Modes,
      .disattenuate       = .disattenuate,
      .correction_factors = correction_factors,
      .reliabilities      = .reliabilities
    ) 
  }

  ## Calculate loadings and cross-loadings (covariance between construct and indicators)
  Lambda <- calculateLoadings(
    .S             = S,
    .W             = W$W,
    .Q             = Q,
    .csem_model    = csem_model,
    .modes         = W$Modes,
    .disattenuate  = .disattenuate
  )

  ## Calculate proxy covariance matrix
  C <- calculateCompositeVCV(.S = S, .W = W$W)
  
  ## Calculate construct correlation matrix
  P <- calculateConstructVCV(.C = C, .Q = Q, .csem_model = csem_model)

  ## Estimate structural coef
  if(.estimate_structural) {
    estim_results <- estimatePathOLS(
      .H            = H,
      .Q            = Q,
      .P            = P,
      .csem_model   = csem_model,
      .normality    = .normality,
      .approach_nl  = .approach_nl
    )
  } else {
    estim_results <- NULL
  }

  ### Output -------------------------------------------------------------------
  out <- list(
    "Estimates"   = list(
      "Path_estimates"         = if(.estimate_structural) {
        estim_results$Path_estimates
      } else {
        estim_results
      },
      "Loading_estimates"      = Lambda * csem_model$measurement,
      "Weight_estimates"       = W$W,
      "Inner_weight_estimates" = W$E,
      "Construct_scores"       = H,
      "Indicator_VCV"          = S,
      "Proxy_VCV"              = C,
      "Construct_VCV"          = P,
      "Cross_loadings"         = Lambda,
      "Construct_reliabilities"= Q^2,
      "Correction_factors"     = correction_factors,
      "R2"                     = if(.estimate_structural) {
        estim_results$R2
      } else {
        estim_results
      }, 
      "R2adj"                     = if(.estimate_structural) {
        estim_results$R2adj
      } else {
        estim_results
      }, 
      "VIF"                     = if(.estimate_structural) {
        estim_results$VIF
      } else {
        estim_results
      }
    ),
    "Information" = list(
      "Data"          = X,
      "Model"         = csem_model,
      "Arguments"     = as.list(match.call())[-1],
      "Type_of_indicator_correlation" = Cor$cor_type,
      "Weight_info"   = list(
        "Modes"              = W$Modes,
        "Number_iterations"  = W$Iterations,
        "Convergence_status" = W$Conv_status
      )
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
