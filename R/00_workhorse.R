#' Composite based estimation
#'
#' This function is the backbone of the `cSEM` package. It is called by [csem],
#' [cca], and test functions such as [testMICOM] to do the actual calculations.
#' It may be called directly by the user, however, in most cases it will likely
#' be more convenient to use [csem] or [cca] instead.
#'
#' @usage workhorse(
#'   .data                    = NULL,
#'   .model                   = NULL,
#'   .approach_cf             = NULL,
#'   .approach_nl             = NULL,
#'   .approach_paths          = NULL,
#'   .approach_weights        = NULL,
#'   .disattenuate            = NULL,
#'   .estimate_structural     = NULL,
#'   .ignore_structural_model = NULL,
#'   .iter_max                = NULL,
#'   .normal                  = NULL,
#'   .PLS_mode                = NULL,
#'   .PLS_weight_scheme_inner = NULL,
#'   .tolerance               = NULL,
#'   .reliabilities           = NULL,
#'   .dominant_indicators = NULL, 
#'   .standardize             = NULL
#'    )
#'
#' @inheritParams csem_arguments
#'
#' @inherit csem_results return
#'
#' @seealso [csem], [cca], [testMICOM]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export
#'

workhorse <- function(
  .data                    = NULL,
  .model                   = NULL,
  .approach_cf             = c("dist_euclid", "dist_euclid_weighted", "fisher_transformed",
                               "mean_geometric", "mean_harmonic", "mean_arithmetic",
                               "geo_of_harmonic"),
  .approach_nl             = c("none", "replace"),
  .approach_paths          = c("OLS", "2SLS", "3SLS"),
  .approach_weights        = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR",
                               "GSCA", "fixed", "unit"),
  .disattenuate            = TRUE,
  .estimate_structural     = TRUE,
  .ignore_structural_model = FALSE,
  .iter_max                = 100,
  .normal                  = TRUE,
  .PLS_mode                = NULL,
  .PLS_weight_scheme_inner = c("centroid", "factorial", "path"),
  .tolerance               = 1e-05,
  .reliabilities           = NULL, 
  .dominant_indicators = NULL
  ) {

  ### Preprocessing ============================================================
  ## Parse and order model to "cSEMModel" list
  csem_model <- parseModel(.model)

  ## Prepare, standardize, check, and clean data
  X <- processData(.data = .data, .model = csem_model) 
  X = scale(X)
  ### Computation ==============================================================
  ## Calculate empirical indicator covariance/correlation matrix
  S <- stats::cov(X)

  ## Calculate weights
  if(.approach_weights == "PLS") {
    W <- calculateWeightsPLS(
      .data                     = S,
      .model                    = csem_model,
      .PLS_mode                 = .PLS_mode,
      .tolerance                = .tolerance,
      .iter_max                 = .iter_max,
      # Arguments passed on to calculateInnerWeightsPLS
      .PLS_weight_scheme_inner  = .PLS_weight_scheme_inner,
      .ignore_structural_model  = .ignore_structural_model
    )
    
    # use dominant indicators approach Henseler et al. (2016), perhaps this
    # applicaple to other weight function.
    # If yes then put at the end
    # Input is a Vektor with indicator names, where the names contain the constructs
    
    if(!is.null(.dominant_indicators)){
      # Check the provided dominant indicators
      if(FALSE %in% (.dominant_indicators %in% colnames(W$W))){
        stop('At least one defined dominant indicator does is not used in the model')}
      if(FALSE %in% names(.dominant_indicators)%in%rownames(W$W) | is.null(names(.dominant_indicators))){
        stop('At least one dominant indicator is defined for a non-defined construct')}
      
     for(i in names(.dominant_indicators)){
       W$W[i,]=W$W[i,]*sign(W$W[i,.dominant_indicators[i]])
     } 
    }
    
    
    
    
    
  } else if(.approach_weights %in% c("SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR")) {
    W <- calculateWeightsKettenring(
      .data                     = data,
      .model                    = csem_model,
      .criteria                 = .approach_weights
    )
  } else if(.approach_weights == "GSCA") {
    W <- calculateWeightsGSCA(
      .data                     = NULL,
      .model                    = csem_model
    )
  } else if(.approach_weights == "fixed") {
    W <- calculateWeightsFixed(
      .data                     = NULL,
      .model                    = csem_model
    )
  } else if(.approach_weights == "unit") {
    W <- calculateWeightsUnit(
      .data                     = NULL,
      .model                    = csemmodel
    )
  } else {
    stop("Unknown weighting approach.")
  }

  ## Calculate proxies/scores
  H <- calculateProxies(
    .X          = X,
    .W          = W$W
  )

  ## Calculate correction factors
  correction_factors <- calculateCorrectionFactors(
    .S            = S,
    .W            = W$W,
    .csem_model   = csem_model,
    .approach_cf  = .approach_cf
  )

  ## Calculate loadings
  Lambda <- calculateLoadings(
    .S             = S,
    .W             = W$W,
    .csem_model    = csem_model,
    .modes         = W$Modes,
    .disattenuate  = .disattenuate,
    .correction_factors = correction_factors
  )

  ## Calculate Q's (correlation between construct and proxy)
  # Note: Q_i^2 := R^2(eta_i; eta_bar_i) is also called the reliability coefficient
  # rho_A in Dijkstra (2015) - Consistent partial least squares path modeling
  # The the Q presented in matrixpls differ from our Q since matrixpls prints the squared Q  

    Q <- calculateProxyConstructCV(
      .S             = S,
      .W             = W$W,
      .csem_model    = csem_model,
      .modes         = W$Modes,
      .disattenuate  = .disattenuate,
      .correction_factors = correction_factors
    )
    
  if(!is.null(.reliabilities)){
    # Check whether defined external reliabilities are correctly defined
    if(TRUE %in% (.reliabilities > 1)){stop('External reliabilities must be smaller equal 1.')}
    if(FALSE %in% (.reliabilities %in% names(Q)) | is.null(.reliabilities)){
      stop('The provided reliability refers to a non-defined construct.')}
    
    Q[names(.reliabilities)] = sqrt(.reliabilities)
  }  
    
 

    

  ## Calculate proxy correlation matrix
  # C <- calculateProxyVCV(.S = S, .W = W$W) 
  C = cov(H)
  ## Calculate construct correlation matrix
  P <- calculateConstructVCV(.C = C, .Q = Q, .csem_model = csem_model)

  ## Estimate structural coef
  if(.estimate_structural) {
    estim_results <- estimatePathOLS(
      .H            = H,
      .W            = W$W,
      .Q            = Q,
      .csem_model   = csem_model,
      .normal       = .normal,
      .approach_nl  = .approach_nl
    )
  } else {
    estim_results <- NA
  }

  ## Set class for Output
  out <- list(
    "Estimates"           = list("Path_estimates"         = if(.estimate_structural) {
      estim_results$Path_estimates
    } else {
      estim_results
    },
                                 "Loading_estimates"      = Lambda,
                                 "Weight_estimates"       = W$W,
                                 "Inner_weight_estimates" = W$E,
                                 "Construct_scores"                = H,
                                 "Indicator_VCV"          = S,
                                 "Proxy_VCV"              = C,
                                 "Construct_VCV"          = P,
                                 "Construct_Reliability"     = Q^2,
                                 "Correction_factors"     = correction_factors
                                 ),
    "Tests"               = list("Overall_model_fit"      = c("still to implement"),
                                 "Some_other_test"        = c("still to implement")
                                 ),
    "Meta_information"    = list("Number_of_observations" = nrow(X),
                                 "Weight_approach"        = .approach_weights,
                                 "Path_approach"          = .approach_paths,
                                 "Construct_types"        = csem_model$construct_type,
                                 "PLS_Modes"              = W$Modes,
                                 "PLS_Inner_Weightning_scheme" = .PLS_weight_scheme_inner
                                 )
    )

  class(out) <- "cSEMResults"
  return(out)


  ### For maintenance: ---------------------------------------------------------
  # X (N x K) := Matrix of indicator values (=data)
  # S (K x K) := Empirical indicator covariance/correlation matrix
  # W (J x K) := (Block-)diagonal matrix of weights with the same dimension as
  #              .modellist$measurement
  # H (N x J) := Matrix of proxy values for the J constructs
  # Lambda (K x J) := Blockdiagonal matrix of factor loadings.
  # Q (1 x J) := Vector of reliability coefficients: Q := R(eta; eta_bar)
  # C (J x J) := Proxy correlation matrix: R(eta_bar, eta_bar)
  # P (J x J) := Construct correlation matrix
}
