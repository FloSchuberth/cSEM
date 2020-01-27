#' Assess model
#'
#' Assess a model using common quality criteria.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article on the
#' \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} for details.
#' 
#' The function is essentially a wrapper around a number of internal functions
#' that perform an "assessment task" (called a **quality criterion** in \pkg{cSEM}
#' parlance) like computing reliability estimates,
#' the effect size, the heterotrait-monotrait ratio of correlations (HTMT) etc.
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
#' \item{Congeneric reliability; "rho_C", "rho_C_mm,", "rho_C_weighted", "rho_C_weighted_mm"}{
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
#' \item{Cronbach's alpha; "cronbachs_alpha"}{An estimate of the
#'   reliability assuming a tau-equivalent measurement model (i.e., a measurement
#'   model with equal loadings) and a test score (proxy) based on unit weights. 
#'   To compute Cronbach's alpha based on a score that uses the weights of the
#'   weight approach used to obtain `.object`, use `"cronbachs_alpha_weighted"` instead.
#'   Cronbach's alpha is an alias for `"rho_T"` the tau-equivalent
#'   reliability which is
#'   the prefered name for this kind of reliability in \pkg{cSEM}, as it clearly states what
#'   it actually estimates (the tau-equivalent reliability as opposed to
#'   the congeneric reliability). "rho_T" and "cronbachs_alpha" are therefore
#'   always identical. 
#'   The tau-equivalent
#'   reliability (Cronbach's alpha) is inherently
#'   tied to the common factor model. It is therefore unclear how to meaningfully 
#'   interpret tau-equivalent
#'   reliability estimates for constructs modeled as composites. 
#'   It is possible to report tau-equivalent
#'   reliability estimates for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateRhoT()]}
#' \item{Distance measures; "dg", "dl", "dml"}{Measures of the distance
#'   between the model-implied and the empirical indicator correlation matrix.
#'   Currently, the geodesic distance (`"dg"`), the squared Euclidian distance
#'   (`"dl"`) and the the maximum likelihood-based distance function are implemented 
#'   (`"dml"`). Calculation is done by [calculateDL()], [calculateDG()], 
#'   and [calculateDML()].}
#' \item{Degrees of freedom, "df"}{
#'   Returns the degrees of freedom. Calculation is done by [calculateDf()].
#'   }
#' \item{Effect size; "esize"}{An index of the effect size of an independent
#'   variable in a structural regression equation. The effect size of the k'th
#'   independent variable in this case
#'   is definied as the ratio (R2_included - R2_excluded)/(1 - R2_included), where 
#'   R2_included and R2_excluded are the R squares of the 
#'   original structural model regression equation (R2_included) and the
#'   alternative specification with the k'th variable dropped (R2_excluded).
#'   This measure is commonly known as Cohen's f^2.
#'   Calculation is done by [calculateEffectSize()].}
#' \item{Fit indices; "chi_square", "chi_square_df", "cfi", "gfi", "ifi", "nfi", 
#'       "nnfi",  "rmsea", "rms_theta", "srmr"}{
#'   Several absolute and incremental fit indices. Note that their suitability
#'   for models containing constructs modeled as composites is still an
#'   open research question. Also note that fit indices are not tests in a 
#'   hypothesis testing sense and
#'   decisions based on common cut-offs proposed in the literature should be
#'   considered with caution!. Calculation is done by [calculateChiSquare()],
#'   [calculateChiSquareDf()], [calculateCFI()], 
#'   [calculateGFI()], [calculateIFI()], [calculateNFI()], [calculateNNFI()], 
#'   [calculateRMSEA()], [calculateRMSTheta()] and [calculateSRMR()].}
#' \item{Fornell-Larcker criterion; "fl_criterion"}{A rule suggested by \insertCite{Fornell1981;textual}{cSEM}
#'   to assess discriminant validity. The Fornell-Larcker
#'   criterion is a decision rule based on a comparison between the squared
#'   construct correlations and the average variance extracted. FL returns
#'   a matrix with the squared construct correlations on the off-diagonal and 
#'   the AVE's on the main diagonal. Calculation is done by `assess()`.}
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
#'   An estimate of the correlation between latent variables. The HTMT is used 
#'   to assess convergent and/or discriminant validity of a construct. 
#'   The HTMT is inherently tied to the common factor model. If the model contains
#'   less than two constructs modeled as common factors and 
#'   `.only_common_factors = TRUE`, `NA` is returned.
#'   It is possible to report the HTMT for constructs modeled as 
#'   composites by setting `.only_common_factors = FALSE`, however, result should be 
#'   interpreted with caution as they may not have a conceptual meaning.
#'   Calculation is done by [calculateHTMT()].}
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
#'   that assume a tau-equivalent measurment model such as Cronbach's alpha.
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
#' we explicitly warn to interpret quality criteria in this case with caution, 
#' as they may not even have a conceptual meaning. 
#'
#' \subsection{Resampling}{
#' To resample a given quality criterion supply the name of the function
#' that calculates the desired quality criterion to [csem()]'s `.user_funs` argument.
#' See [resamplecSEMResults()] for details.
#' }
#' 
#' @usage assess(
#'   .object              = NULL, 
#'   .only_common_factors = TRUE, 
#'   .quality_criterion   = c("all", "ave", "rho_C", "rho_C_mm", "rho_C_weighted", 
#'                            "rho_C_weighted_mm", "cronbachs_alpha", 
#'                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
#'                           "esize", "chi_square", "chi_square_df",
#'                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
#'                           "reliability", 
#'                           "rmsea", "rms_theta", "srmr",
#'                           "gof", "htmt", "r2", "r2_adj",
#'                           "rho_T", "rho_T_weighted", "vif", 
#'                           "vifmodeB",  "fl_criterion"),
#'   ...
#' )
#' 
#' @inheritParams csem_arguments
#' @param ... Further arguments passed to functions called by [assess()].
#'   See [args_assess_dotdotdot] for a complete list of available arguments.
#'
#' @seealso [csem()], [resamplecSEMResults()]
#'
#' @return A named list of quality criteria. Note that if only a single quality
#'   criteria is computed the return value is still a list!
#'   
#' @example inst/examples/example_assess.R
#' 
#' @export

assess <- function(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "chi_square", "chi_square_df",
                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB",  "fl_criterion"),
  ...
){
  UseMethod("assess")
}

#' @export

assess.cSEMResults_default <- function(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "chi_square", "chi_square_df",
                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB",  "fl_criterion"),
  ...
){
  
  ## Check arguments
  match.arg(.quality_criterion, 
            args_default(.choices = TRUE)$.quality_criterion, several.ok = TRUE)
  
  ## Set up empty list
  out <- list()
  
  ## Select quality criteria
  if(any(.quality_criterion %in% c("all", "ave"))) {
    # AVE
    out[["AVE"]] <- calculateAVE(
      .object, 
      .only_common_factors = .only_common_factors
    )
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
  if(any(.quality_criterion %in% c("all", "esize"))) {
    # Effect size
    out[["Effect_size"]] <- calculateEffectSize(.object)
  }
  if(any(.quality_criterion %in% c("all", "chi_square"))) {
    # Effect size
    out[["Chi_square"]] <- calculateChiSquare(.object)
  }
  if(any(.quality_criterion %in% c("all", "chi_square_df"))) {
    # Effect size
    out[["Chi_square_df"]] <- calculateChiSquareDf(.object)
  }
  if(any(.quality_criterion %in% c("all", "cfi"))) {
    # Effect size
    out[["CFI"]] <- calculateCFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "gfi"))) {
    # Effect size
    out[["GFI"]] <- calculateGFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "ifi"))) {
    # Effect size
    out[["IFI"]] <- calculateIFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "nfi"))) {
    # Effect size
    out[["NFI"]] <- calculateNFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "nnfi"))) {
    # Effect size
    out[["NNFI"]] <- calculateNNFI(.object)
  }
  if(any(.quality_criterion %in% c("all", "rmsea"))) {
    # Effect size
    out[["RMSEA"]] <- calculateRMSEA(.object)
  }
  if(any(.quality_criterion %in% c("all", "rms_theta"))) {
    # Effect size
    out[["RMS_theta"]] <- calculateRMSTheta(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "srmr"))) {
    # Effect size
    out[["SRMR"]] <- calculateSRMR(.object, ...)
  }
  if(any(.quality_criterion %in% c("all", "fl_criterion"))) {
    # Fornell-Larcker
    ## Get relevant objects
    con_types <-.object$Information$Model$construct_type
    names_cf  <- names(con_types[con_types == "Common factor"])
    P         <- .object$Estimates$Construct_VCV
    
    if(.only_common_factors) {
      P <- P[names_cf, names_cf]
    }
    
    if(sum(dim(P)) > 0) {
      FL_matrix <- cov2cor(P)^2
      diag(FL_matrix) <- calculateAVE(.object, 
                                      .only_common_factors = .only_common_factors)
      out[["Fornell-Larcker"]] <- FL_matrix
    }
    
    
  }
  if(any(.quality_criterion %in% c("all", "gof"))) {
    # GoF
    out[["GoF"]]   <- calculateGoF(
      .object, 
      .only_common_factors = .only_common_factors
    )
  }
  if(any(.quality_criterion %in% c("all", "htmt"))) {
    # HTMT 
    out[["HTMT"]]  <- calculateHTMT(
      .object, 
      .only_common_factors = .only_common_factors
    )
  }
  if(any(.quality_criterion %in% c("all", "r2"))) {
    # R2
    out[["R2"]]  <- .object$Estimates$R2
  }
  if(any(.quality_criterion %in% c("all", "r2_adj"))) {
    # Adjusted R2
    out[["R2_adj"]]  <- .object$Estimates$R2adj
  }
  if(any(.quality_criterion %in% c("all", "reliability"))) {
    # RhoT
    out[["reliability"]]  <- list()
    
    # Cronbachs alpha (rho_T)
    out$reliability[["Cronbachs_alpha"]] <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      ...
    )
    
    # Joereskogs rho (rho_C_mm)
    out$reliability[["Joereskogs_rho"]] <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors,
      .weighted = FALSE,
      .model_implied = TRUE
    )
    
    # Dijkstra-Henselers rho A (rho_C_weighted_mm)
    out$reliability[["Dijkstra-Henselers_rho_A"]] <- calculateRhoC(
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
  if(any(.quality_criterion %in% c("all", "vif"))) {
    # VIF
    out[["VIF"]]  <- .object$Estimates$VIF
    
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
    # VIFModeB
    out[["VIF_modeB"]] <- calculateVIFModeB(.object)
  }
  
  out[["Information"]] <- list()
  out$Information[["All"]] <- FALSE
  
  if(any(.quality_criterion == "all")) {
    out$Information$All <- TRUE
  }
  
  class(out) <- "cSEMAssess"
  return(out)
}

#' @export

assess.cSEMResults_multi <- function(
  .object              = NULL,
  .only_common_factors = TRUE,
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "cronbachs_alpha",  
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "chi_square", "chi_square_df",
                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB",  "fl_criterion"),
  ...
){
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    lapply(.object, assess.cSEMResults_2ndorder, 
           .only_common_factors = .only_common_factors, 
           .quality_criterion = .quality_criterion, 
           ...
    )
  } else {
    lapply(.object, assess.cSEMResults_default, 
           .only_common_factors = .only_common_factors,
           .quality_criterion = .quality_criterion,
           ...
    )
  }
}

#' @export

assess.cSEMResults_2ndorder <- function(
  .object              = NULL,
  .only_common_factors = TRUE,
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_mm", "rho_C_weighted", 
                           "rho_C_weighted_mm", "cronbachs_alpha",  
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "chi_square", "chi_square_df",
                           "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "reliability",
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB",  "fl_criterion"),
  ...
  ){
  
  stop2("Currently, models containing second-order constructs are not supported by `assess()`.")
}