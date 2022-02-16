#' Assess model
#'
#' \lifecycle{maturing}
#' 
#' Assess a model using common quality criteria.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article on the
#' \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} for details.
#' 
#' The function is essentially a wrapper around a number of internal functions
#' that perform an "assessment task" (called a **quality criterion** in \pkg{cSEM}
#' parlance) like computing reliability estimates,
#' the effect size (Cohen's f^2), the heterotrait-monotrait ratio of correlations (HTMT) etc.
#' 
#' By default every possible quality criterion is calculated (`.quality_criterion = "all"`). 
#' If only a subset of quality criteria are needed a single character string
#' or a vector of character strings naming the criteria to be computed may be 
#' supplied to [assess()] via the `.quality_criterion` argument. Currently, the
#' following quality criteria are implemented (in alphabetical order):
#' 
#' \describe{
#' \item{Average variance extracted (AVE); "ave"}{An estimate of the 
#'   amount of variation in the indicators that is due to the underlying latent variable. 
#'   Practically, it is calculated as the ratio of the (indicator) true score variances 
#'   (i.e., the sum of the squared loadings)
#'   relative to the sum of the total indicator variances. The AVE is inherently
#'   tied to the common factor model. It is therefore unclear how to meaningfully 
#'   interpret AVE results for constructs modeled as composites. 
#'   It is possible to report the AVE for constructs modeled as composites by setting 
#'   `.only_common_factors = FALSE`, however, result should be interpreted with caution 
#'   as they may not have a conceptual meaning. Calculation is done
#'   by [calculateAVE()].}
#' \item{Congeneric reliability; "rho_C", "rho_C_mm", "rho_C_weighted", "rho_C_weighted_mm"}{
#'   An estimate of the reliability assuming a congeneric measurement model (i.e., loadings are
#'   allowed to differ) and a test score (proxy) based on unit weights.
#'   There are four different versions implemented. See the 
#'   \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#methods}{Methods and Formulae} section
#'   of the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#'   article on the
#'   \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} for details.
#'   Alternative but synonemmous names for `"rho_C"` are: 
#'   composite reliability, construct reliablity, reliability coefficient, 
#'   Joereskog's rho, coefficient omega, or Dillon-Goldstein's rho. 
#'   For `"rho_C_weighted"`: (Dijkstra-Henselers) rhoA. `rho_C_mm` and `rho_C_weighted_mm`
#'   have no corresponding names. The former uses unit weights scaled by (w'Sw)^(-1/2) and
#'   the latter weights scaled by (w'Sigma_hat w)^(-1/2) where Sigma_hat is 
#'   the model-implied indicator correlation matrix.
#'   The Congeneric reliability is inherently
#'   tied to the common factor model. It is therefore unclear how to meaningfully 
#'   interpret congeneric reliability estimates for constructs modeled as composites. 
#'   It is possible to report the congeneric reliability for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateRhoC()].}
#' \item{Distance measures; "dg", "dl", "dml"}{Measures of the distance
#'   between the model-implied and the empirical indicator correlation matrix.
#'   Currently, the geodesic distance (`"dg"`), the squared Euclidian distance
#'   (`"dl"`) and the the maximum likelihood-based distance function are implemented 
#'   (`"dml"`). Calculation is done by [calculateDL()], [calculateDG()], 
#'   and [calculateDML()].}
#' \item{Degrees of freedom, "df"}{
#'   Returns the degrees of freedom. Calculation is done by [calculateDf()].
#'   }
#' \item{Effects; "effects"}{Total and indirect effect estimates. Additionally, 
#'   the variance accounted for (VAF) is computed. The VAF is defined as the ratio of a variables
#'   indirect effect to its total effect. Calculation is done
#'   by [calculateEffects()].}
#' \item{Effect size; "f2"}{An index of the effect size of an independent
#'   variable in a structural regression equation. This measure is commonly 
#'   known as Cohen's f^2. The effect size of the k'th
#'   independent variable in this case
#'   is definied as the ratio (R2_included - R2_excluded)/(1 - R2_included), where 
#'   R2_included and R2_excluded are the R squares of the 
#'   original structural model regression equation (R2_included) and the
#'   alternative specification with the k'th variable dropped (R2_excluded).
#'   Calculation is done by [calculatef2()].}
#' \item{Fit indices; "chi_square", "chi_square_df", "cfi", "cn", "gfi", "ifi", "nfi", 
#'       "nnfi",  "rmsea", "rms_theta", "srmr"}{
#'   Several absolute and incremental fit indices. Note that their suitability
#'   for models containing constructs modeled as composites is still an
#'   open research question. Also note that fit indices are not tests in a 
#'   hypothesis testing sense and
#'   decisions based on common one-size-fits-all cut-offs proposed in the literature 
#'   suffer from serious statistical drawbacks. Calculation is done by [calculateChiSquare()],
#'   [calculateChiSquareDf()], [calculateCFI()], 
#'   [calculateGFI()], [calculateIFI()], [calculateNFI()], [calculateNNFI()], 
#'   [calculateRMSEA()], [calculateRMSTheta()] and [calculateSRMR()].}
#' \item{Fornell-Larcker criterion; "fl_criterion"}{A rule suggested by \insertCite{Fornell1981;textual}{cSEM}
#'   to assess discriminant validity. The Fornell-Larcker
#'   criterion is a decision rule based on a comparison between the squared
#'   construct correlations and the average variance extracted. FL returns
#'   a matrix with the squared construct correlations on the off-diagonal and 
#'   the AVE's on the main diagonal. Calculation is done by `calculateFLCriterion()`.}
#' \item{Goodness of Fit (GoF); "gof"}{The GoF is defined as the square root 
#'   of the mean of the R squares of the structural model times the mean 
#'   of the variances in the indicators that are explained by their 
#'   related constructs (i.e., the average over all lambda^2_k).
#'   For the latter, only constructs modeled as common factors are considered
#'   as they explain their indicator variance in contrast to a composite where 
#'   indicators actually build the construct.
#'   Note that, contrary to what the name suggests, the GoF is **not** a 
#'   measure of model fit in a Chi-square fit test sense. Calculation is done 
#'   by [calculateGoF()].}
#' \item{Heterotrait-monotrait ratio of correlations (HTMT); "htmt"}{
#'   An estimate of the correlation between latent variables assuming tau equivalent
#'   measurement models. The HTMT is used 
#'   to assess convergent and/or discriminant validity of a construct. 
#'   The HTMT is inherently tied to the common factor model. If the model contains
#'   less than two constructs modeled as common factors and 
#'   `.only_common_factors = TRUE`, `NA` is returned.
#'   It is possible to report the HTMT for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateHTMT()].}
#' \item{HTMT2; "htmt2"}{
#'   An estimate of the correlation between latent variables assuming congeneric
#'   measurement models. The HTMT2 is used 
#'   to assess convergent and/or discriminant validity of a construct. 
#'   The HTMT is inherently tied to the common factor model. If the model contains
#'   less than two constructs modeled as common factors and 
#'   `.only_common_factors = TRUE`, `NA` is returned.
#'   It is possible to report the HTMT for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateHTMT()].}
#' \item{Model selection criteria: "aic", "aicc", "aicu", "bic", "fpe", "gm", 
#'       "hq", "hqc", "mallows_cp"}{
#'   Several model selection criteria as suggested by \insertCite{Sharma2019;textual}{cSEM}
#'   in the context of PLS. See: [calculateModelSelectionCriteria()] for details.}
#' \item{Reliability: "reliability"}{
#'   As described in the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#methods}{Methods and Formulae} 
#'   section of the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#'   article on the \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} 
#'   there are many different estimators for the (internal consistency) reliability.
#'   Choosing `.quality_criterion = "reliability"` computes the three most common
#'   measures, namely: "Cronbachs alpha" (identical to "rho_T"), "Jöreskogs rho" (identical to "rho_C_mm"),
#'   and "Dijkstra-Henselers rho A" (identical to "rho_C_weighted_mm").
#'   Reliability is inherently
#'   tied to the common factor model. It is therefore unclear how to meaningfully 
#'   interpret reliability estimates for constructs modeled as composites. 
#'   It is possible to report the three common reliability estimates for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   }
#' \item{R square and R square adjusted; "r2", "r2_adj"}{The R square and the adjusted
#'   R square for each structural regression equation.
#'   Calculated when running [csem()].}
#' \item{Tau-equivalent reliability; "rho_T"}{An estimate of the
#'   reliability assuming a tau-equivalent measurement model (i.e. a measurement
#'   model with equal loadings) and a test score (proxy) based on unit weights.
#'   Tau-equivalent reliability is the preferred name for reliability estimates
#'   that assume a tau-equivalent measurement model such as Cronbach's alpha.
#'   The tau-equivalent
#'   reliability (Cronbach's alpha) is inherently
#'   tied to the common factor model. It is therefore unclear how to meaningfully 
#'   interpret tau-equivalent
#'   reliability estimates for constructs modeled as composites. 
#'   It is possible to report tau-equivalent
#'   reliability estimates for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateRhoT()].}
#' \item{Variance inflation factors (VIF); "vif"}{An index for the amount of (multi-) 
#'   collinearity between independent variables of a regression equation. Computed
#'   for each structural equation. Practically, VIF_k is defined
#'   as the ratio of 1 over (1 - R2_k) where R2_k is the R squared from a regression
#'   of the k'th independent variable on all remaining independent variables.
#'   Calculated when running [csem()].}
#' \item{Variance inflation factors for PLS-PM mode B (VIF-ModeB); "vifmodeB"}{An index for 
#'   the amount of (multi-) collinearity between independent variables (indicators) in
#'   mode B regression equations. Computed only if `.object` was obtained using
#'   `.weight_approach = "PLS-PM"` and at least one mode was mode B. 
#'   Practically, VIF-ModeB_k is defined as the ratio of 1 over (1 - R2_k) where 
#'   R2_k is the R squared from a regression of the k'th indicator of block j on
#'   all remaining indicators of the same block.
#'   Calculation is done by [calculateVIFModeB()].}
#' }
#' 
#' For details on the most important quality criteria see the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#methods}{Methods and Formulae} section
#' of the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article on the on the
#' \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website}.
#' 
#' Some of the quality criteria are inherently tied to the classical common
#' factor model and therefore only meaningfully interpreted within a common
#' factor model (see the 
#' \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article for details). 
#' It is possible to force computation of all quality criteria for constructs 
#' modeled as composites by setting `.only_common_factors = FALSE`, however, 
#' we explicitly warn to interpret quality criteria in analogy to the common factor 
#' model in this case, as the interpretation often does not carry over to composite models.
#'
#' \subsection{Resampling}{
#' To resample a given quality criterion supply the name of the function
#' that calculates the desired quality criterion to [csem()]'s `.user_funs` argument.
#' See [resamplecSEMResults()] for details.
#' }
#' 
#' @usage assess(
#'   .object              = NULL, 
#'   .quality_criterion   = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
#'                            "hqc", "mallows_cp", "ave",
#'                            "rho_C", "rho_C_mm", "rho_C_weighted", 
#'                            "rho_C_weighted_mm", "dg", "dl", "dml", "df",
#'                            "effects", "f2", "fl_criterion", "chi_square", "chi_square_df",
#'                            "cfi", "cn", "gfi", "ifi", "nfi", "nnfi", 
#'                            "reliability",
#'                            "rmsea", "rms_theta", "srmr",
#'                            "gof", "htmt", "htmt2", "r2", "r2_adj",
#'                            "rho_T", "rho_T_weighted", "vif", 
#'                            "vifmodeB"),
#'   .only_common_factors = TRUE, 
#'   ...
#' )
#' 
#' @inheritParams csem_arguments
#' @param ... Further arguments passed to functions called by [assess()].
#'   See [args_assess_dotdotdot] for a complete list of available arguments.
#'
#' @seealso [csem()], [resamplecSEMResults()], [exportToExcel()]
#'
#' @return A named list of quality criteria. Note that if only a single quality
#'   criteria is computed the return value is still a list!
#'   
#' @example inst/examples/example_assess.R
#' 
#' @export

assess <- function(
  .object              = NULL, 
  .quality_criterion   = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                           "hqc", "mallows_cp", "ave",
                           "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "dg", "dl", "dml", "df",
                           "effects", "f2", "fl_criterion", "chi_square", "chi_square_df",
                           "cfi", "cn", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "htmt2", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB"),
  .only_common_factors = TRUE, 
  ...
){
  
  # Note to authors:
  #   If you add a new argument to .quality_criterion, make sure to add it
  #   to the exportToExcel() function as well!
  
  ## Check arguments
  match.arg(.quality_criterion, 
            args_default(.choices = TRUE)$.quality_criterion, several.ok = TRUE)
  
  ## 
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, assess, 
           .only_common_factors = .only_common_factors,
           .quality_criterion = .quality_criterion,
           ...
    )
    return(out)
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    x11 <- .object$First_stage$Estimates
    x12 <- .object$First_stage$Information
    
    x21 <- .object$Second_stage$Estimates
    x22 <- .object$Second_stage$Information
    
  } else if(inherits(.object, "cSEMResults_default")) {

    x21 <- .object$Estimates
    x22 <- .object$Information
    
  } else {
    stop2(
      "The following error occured in the assess() function:\n",
      "`.object` must be a `cSEMResults` object."
    )
  }
  
  if(x22$Model$model_type != "Linear") {
    stop2("Currently, `assess()` does not support models containing nonlinear terms.",
          "Use the individual `calculateXXX()` functions instead.")
  }
  ## Set up empty list
  out <- list()
  out[["Information"]] <- list()
  
  ## Select quality criteria
  if(any(.quality_criterion %in% c("all", "ave"))) {
    # AVE
    out[["AVE"]] <- calculateAVE(
      .object, 
      .only_common_factors = .only_common_factors
    )
  }
  if(any(.quality_criterion %in% c("all", "aic"))) {
    # Akaike information criterion (AIC)
    out[["AIC"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "aic",
      .by_equation = TRUE
      )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "aicc"))) {
    # Corrected AIC (AICc)
    out[["AICc"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "aicc",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "aicu"))) {
    # Unbiased AIC (AICu)
    out[["AICu"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "aicu",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "bic"))) {
    # Bayesian information criterion (BIC)
    out[["BIC"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "bic",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "fpe"))) {
    # Bayesian information criterion (BIC)
    out[["FPE"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "fpe",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "gm"))) {
    # Bayesian information criterion (BIC)
    out[["GM"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "gm",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "hq"))) {
    # Bayesian information criterion (BIC)
    out[["HQ"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "hq",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "hqc"))) {
    # Bayesian information criterion (BIC)
    out[["HQc"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "hqc",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "mallows_cp"))) {
    # Bayesian information criterion (BIC)
    out[["Mallows_Cp"]]  <- calculateModelSelectionCriteria(
      .object,
      .ms_criterion = "mallows_cp",
      .by_equation = TRUE
    )[[1]]
  }
  if(any(.quality_criterion %in% c("all", "rho_C"))) {
    # RhoC
    out[["RhoC"]]  <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors,
      .weighted = FALSE,
      .model_implied = TRUE
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_C_mm"))) {
    # RhoC
    out[["RhoC_mm"]]  <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors,
      .weighted = FALSE,
      .model_implied = FALSE
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_C_weighted"))) {
    # RhoC weighted
    out[["RhoC_weighted"]]  <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors, 
      .weighted = TRUE,
      .model_implied = FALSE
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_C_weighted_mm"))) {
    # RhoC weighted
    out[["RhoC_weighted_mm"]]  <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors, 
      .weighted = TRUE,
      .model_implied = TRUE
    )
  }
  if(any(.quality_criterion %in% c("all", "dg"))) {
    # dG
    out[["DG"]]    <- calculateDG(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "dl"))) {
    # dL
    out[["DL"]]    <- calculateDL(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "dml"))) {
    # dML
    out[["DML"]]   <- calculateDML(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "df"))) {
    # dML
    out[["Df"]]   <- calculateDf(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "effects")) && 
     !all(x22$Model$structural == 0)) {
    # Direct, total and indirect effects
    out[["Effects"]] <- if(inherits(.object, "cSEMResults_2ndorder")) {
      summarize(.object)$Second_stage$Estimates$Effect_estimates
    } else {
      summarize(.object)$Estimates$Effect_estimates 
    }
  }
  if(any(.quality_criterion %in% c("all", "f2")) && !all(x22$Model$structural == 0) 
     && x22$Arguments$.approach_paths == "OLS") {
    
    # Effect size (f2)
    out[["F2"]] <- calculatef2(.object)
  }
  if(any(.quality_criterion %in% c("all", "chi_square"))) {
    out[["Chi_square"]] <- calculateChiSquare(.object)
  }
  if(any(.quality_criterion %in% c("all", "chi_square_df"))) {
    out[["Chi_square_df"]] <- calculateChiSquareDf(.object)
  }
  if(any(.quality_criterion %in% c("all", "cfi"))) {
    out[["CFI"]] <- calculateCFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "gfi"))) {
    out[["GFI"]] <- calculateGFI(
      .object,
      ...
      )
  }
  if(any(.quality_criterion %in% c("all", "cn"))) {
    # Hoelter's (critical) N (CN)
    out[["CN"]] <- calculateCN(
      .object,
      ...
    )
  }
  if(any(.quality_criterion %in% c("all", "ifi"))) {
    out[["IFI"]] <- calculateIFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "nfi"))) {
    out[["NFI"]] <- calculateNFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "nnfi"))) {
    out[["NNFI"]] <- calculateNNFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "rmsea"))) {
    out[["RMSEA"]] <- calculateRMSEA(.object)
  }
  if(any(.quality_criterion %in% c("all", "rms_theta"))) {
    if(inherits(.object, "cSEMResults_default")) {
      out[["RMS_theta"]] <- calculateRMSTheta(.object)
    } else {
      warning("Computation of the RMS_theta",
              " not supported for models containing second-order constructs:\n",
              "Argument 'rms_theta' is ignored.", call. = FALSE) 
    }
  }
  if(any(.quality_criterion %in% c("all", "srmr"))) {
    out[["SRMR"]] <- calculateSRMR(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "fl_criterion"))) {
    if(inherits(.object, "cSEMResults_default")) {
      out[["Fornell-Larcker"]] <- calculateFLCriterion(
        .object, 
        .only_common_factors = .only_common_factors)
    } else {
      warning("Computation of the Fornell-Larcker criterion",
              " not supported for models containing second-order constructs:\n",
              "Argument 'fl_criterion' is ignored.", call. = FALSE) 
    }
  }
  if(any(.quality_criterion %in% c("all", "gof")) && !all(x22$Model$structural == 0)) {
    # GoF
    out[["GoF"]]   <- calculateGoF(.object)
  }
  if(any(.quality_criterion %in% c("all", "htmt", "htmt2"))) {
    # HTMT 
    if(inherits(.object, "cSEMResults_default")) {
      if(any(.quality_criterion %in% c("all", "htmt"))){
       out[["HTMT"]]  <- calculateHTMT(
        .object,
        .only_common_factors  = .only_common_factors,
        .type_htmt = "htmt",
        ...
       )}
      
      if(any(.quality_criterion %in% c("all", "htmt2"))){
        out[["HTMT2"]]  <- calculateHTMT(
          .object,
          .only_common_factors  = .only_common_factors,
          .type_htmt = "htmt2",
          ...
        )}
    
      # Get argument values
      args_htmt <- list(...)
      if(any(names(args_htmt) == ".inference")) {
        out$Information[[".inference"]] <- args_htmt[[".inference"]]
      } else {
        out$Information[[".inference"]] <- formals(calculateHTMT)[[".inference"]]
      }
      
      if(any(names(args_htmt) == ".alpha")) {
        out$Information[[".alpha"]] <- args_htmt[[".alpha"]]
      } else {
        out$Information[[".alpha"]] <- formals(calculateHTMT)[[".alpha"]]
      } 
      
      if(any(names(args_htmt) == ".ci")) {
        out$Information[[".ci"]] <- args_htmt[[".ci"]]
      } else {
        out$Information[[".ci"]] <- "CI_percentile"
      } 
      
      if(any(names(args_htmt) == ".type_htmt")) {
        out$Information[[".type_htmt"]] <- args_htmt[[".type_htmt"]]
      } else {
        # If .type_htmt is not set in the function
        out$Information[[".type_htmt"]] <- "htmt"
      } 
    } else { # 2nd_order
      warning("Computation of the HTMT",
              " not supported for models containing second-order constructs:\n",
              "Argument 'htmt' is ignored.", call. = FALSE) 
    }
  }
  if(any(.quality_criterion %in% c("all", "r2")) && !all(x22$Model$structural == 0)) {
    # R2
    out[["R2"]]  <- x21$R2
    if(inherits(.object, "cSEMResults_2ndorder")) {
      # Remove temp from the R2s of the second stage
      names(out$R2) <- gsub("_temp", "", names(x21$R2))
      # out$R2 <- out$R2[intersect(second_order, cfs)]
    }
  }
  if(any(.quality_criterion %in% c("all", "r2_adj")) && !all(x22$Model$structural == 0)) {
    # Adjusted R2
    out[["R2_adj"]]  <- x21$R2adj
    if(inherits(.object, "cSEMResults_2ndorder")) {
      # Remove temp from the R2s of the second stage
      names(out$R2_adj) <- gsub("_temp", "", names(x21$R2adj))

    }
  }
  if(any(.quality_criterion %in% c("all", "reliability"))) {
    # RhoT
    out[["Reliability"]]  <- list()
    
    # Cronbachs alpha (rho_T)
    out$Reliability[["Cronbachs_alpha"]] <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      ...
    )
    
    # Joereskogs rho (rho_C_mm)
    out$Reliability[["Joereskogs_rho"]] <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors,
      .weighted = FALSE,
      .model_implied = TRUE
    )
    
    # Dijkstra-Henselers rho A (rho_C_weighted_mm)
    out$Reliability[["Dijkstra-Henselers_rho_A"]] <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors, 
      .weighted = TRUE,
      .model_implied = FALSE
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_T"))) {
    # RhoT
    out[["RhoT"]]  <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      ...
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_T_weighted"))) {
    # RhoC weighted
    out[["RhoT_weighted"]]  <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      .weighted            = TRUE, 
      ...
    )
  }
  if(any(.quality_criterion %in% c("all", "vif")) && !all(x22$Model$structural == 0)) {
    # VIF
    out[["VIF"]]  <- x21$VIF
    
    # Make output a matrix:
    # Note: this is necessary to be able to bootstrap the VIFs
    #       via the .user_funs argument. Currently, .user_funs functions 
    #       need to return a vector or a matrix. I may change that in the future.
    m <- matrix(0, nrow = length(names(out$VIF)), ncol = length(unique(unlist(sapply(out$VIF, names)))),
                dimnames = list(names(out$VIF), unique(unlist(sapply(out$VIF, names)))))
    
    for(i in names(out$VIF)) {
      m[i, match(names(out$VIF[[i]]), colnames(m))] <- out$VIF[[i]]
    }
    
    out$VIF <- m
    
  }
  if(any(.quality_criterion %in% c("all", "vifmodeB"))) {
      out[["VIF_modeB"]] <- calculateVIFModeB(.object)
   }
  
  out$Information[["All"]] <- FALSE
  
  if(any(.quality_criterion == "all")) {
    out$Information$All <- TRUE
  }
  
  class(out) <- "cSEMAssess"
  return(out)
}