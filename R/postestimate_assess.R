#' Assess model
#'
#' Assess a model using common quality criteria.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html}{Postestimation: Assessing a model} 
#' article on the
#' \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website} for details.
#' 
#' The function is essentially a wrapper around a number of internal functions
#' that perform an "assessment task" (called a **quality criterion** in \pkg{cSEM}
#' parlance) like computing the (congeneric) reliability,
#' the effect size, the heterotrait-monotrait ratio of correlations (HTMT) etc.
#' 
#' By default every possible quality criterion is calculated (`.quality_criterion = "all"`). 
#' If only a subset of quality criteria are needed a single character string
#' or a vector of character strings naming the criteria to be computed may be 
#' supplied to [assess()] via the `.quality_criterion` argument. Currently, the
#' following quality criteria are implemented (in alphabetical order):
#' \describe{
#' \item{Average variance extracted (AVE); "ave"}{An estimate of the 
#'   amount of variation in the indicators that is due to the underlying latent variable. 
#'   Practically, it is calculated as the ratio of the (indicator) true score variances 
#'   (i.e., the sum of the squared loadings)
#'   relative to the sum of the total indicator variances. Calculation is done
#'   by [calculateAVE()].}
#' \item{Congeneric reliability; "rho_C"}{An estimate of the 
#'   reliability assuming a congeneric measurement model (i.e., loadings are
#'   allowed to differ) and a test score (proxy) based on unit weights.
#'   To compute the congeneric reliability based on a score that uses the weights of the
#'   weight approach used to obtain `.object`, use `"rho_C_weighted"` instead.
#'   Congeneric reliability is the unified name for 
#'   reliability estimates that assume a congeneric measurement model. 
#'   Alternative but synonemmous names for `"rho_C"` are: 
#'   composite reliability, construct reliablity, reliability coefficient, 
#'   Jöreskog's rho, coefficient omega, or Dillon-Goldstein's rho. 
#'   For `"rho_C_weighted"`: rho_A, or rho_B. Calculation is done by [calculateRhoC()].}
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
#'   always identical. Calculation is done by [calculateRhoT()]}
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
#'   Calculation is done by [calculateEffectSize()].}
#' \item{Fit indices; "cfi", "gfi", "ifi", "nfi", "nnfi",  "rmsea", "rms_theta"
#'   "srmr"}{
#'   Several absolute and incremental fit indices. Note that their suitability
#'   for models containing constructs modeled as common factors is still an
#'   open research question. Also note that fit indices are not tests in a 
#'   hypothesis testing sense and
#'   decisions based on common cut-offs proposed in the literature should be
#'   considered with caution!. Calculation is done by [calculateCFI()], 
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
#'   An estimate of the latent variable correlation used to assess
#'   convergent and/or discriminant validity of a construct. Calculation is done
#'   by [calculateHTMT()].}
#' \item{R square and R square adjusted; "r2", "r2_adj"}{The R square and the adjusted
#'   R square for each structural regression equation.
#'   Calculated when running [csem()].}
#' \item{Redundancy analysis (RA); "ra"}{The process of regressing the scores 
#'   of a reflectively measured construct on the scores of a formatively measured 
#'   construct in order to gain empirical evidence for convergent validity of a 
#'   formatively measured construct. 
#'   RA is therefore confined to PLS, specifically PLS with at least one construct
#'   whose mode is Mode B. This is the case if the construct is modeled as a 
#'   composite or if the construct was explicitly given Mode B.
#'   Hence RA is only done if `.object` was obtained using 
#'   `.approach_weights = "PLS-PM"` and if at least one constructs mode is Mode B.
#'   Performed by [doRedundancyAnalysis()].}
#' \item{Tau-equivalent reliability; "rho_T"}{An estimate of the
#'   reliability assuming a tau-equivalent measurement model (i.e. a measurement
#'   model with equal loadings) and a test score (proxy) based on unit weights.
#'   Tau-equivalent reliability is the preferred name for reliability estimates
#'   that assume a tau-equivalent measurment model such as Cronbach's alpha.
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
#' For details on all quality criteria see the \href{https://m-e-rademaker.github.io/cSEM/articles/Using-assess.html#methods}{Methods and Formulae} section
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
#'   .quality_criterion   = c("all", "ave", "rho_C", "rho_C_weighted", "cronbachs_alpha", 
#'                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
#'                           "esize", "cfi", "gfi", "ifi", "nfi", "nnfi", 
#'                           "rmsea", "rms_theta", "srmr",
#'                           "gof", "htmt", "r2", "r2_adj", "ra",
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
#' @examples 
#' # ===========================================================================
#' # Using the threecommonfactors dataset
#' # ===========================================================================
#' model <- "
#' # Structural model
#' eta2 ~ eta1
#' eta3 ~ eta1 + eta2
#' 
#' # Each concept os measured by 3 indicators, i.e., modeled as latent variable
#' eta1 =~ y11 + y12 + y13
#' eta2 =~ y21 + y22 + y23
#' eta3 =~ y31 + y32 + y33
#' "
#' 
#' res <- csem(threecommonfactors, model)
#' a   <- assess(res) # computes all quality criteria (.quality_criterion = "all")
#' a
#' 
#' ## The return value is a named list
#' str(a)
#' a$HTMT
#' 
#' # You may also just compute a subset of quality criteria
#' assess(res, .quality_criterion = c("ave", "rho_C", "htmt"))
#' 
#' ## Resampling ---------------------------------------------------------------
#' # To resample a given quality criterion use csem()'s .user_funs argument
#' 
#' res <- csem(threecommonfactors, model, 
#'             .resample_method = "bootstrap", 
#'             .user_funs       = cSEM:::calculateHTMT,
#'             .R               = 80 
#' )
#' 
#' ## Look at the resamples
#' res$Estimates$Estimates_resample$Estimates1$User_fun$Resampled[1:4, ]
#' 
#' ## Use infer() to compute e.g. the 95% percentile confidence interval
#' res_infer <- infer(res, .quantity = "CI_percentile")
#' res_infer$User_fun
#' 
#'   
#' # Note that .user_funs expects a function that returns a vector or a matrix!
#' 
#' @export

assess <- function(
  .object              = NULL, 
  .only_common_factors = TRUE, 
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_weighted", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj", "ra",
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
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_weighted", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj", "ra",
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
      .only_common_factors = .only_common_factors
    )
  }
  if(any(.quality_criterion %in% c("all", "rho_C_weighted"))) {
    # RhoC weighted
    out[["RhoC_weighted"]]  <- calculateRhoC(
      .object, 
      .only_common_factors = .only_common_factors, 
      .weighted = TRUE
    )
  }
  if(any(.quality_criterion %in% c("all", "cronbachs_alpha"))) {
    # Cronbach's alpha aka RhoT
    out[["Cronbachs_alpha"]]  <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      ...
    )
  }
  if(any(.quality_criterion %in% c("all", "cronbachs_alpha_weighted"))) {
    # Cronbach's alpha weighted aka RhoT weighted
    out[["Cronbachs_alpha_weighted"]]  <- calculateRhoT(
      .object, 
      .only_common_factors = .only_common_factors, 
      .output_type         = "vector",
      .weighted            = TRUE,
      ...
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
  if(any(.quality_criterion %in% c("all", "ra"))) {
    # Redundancy analysis (RA)
    out[["RA"]] <- doRedundancyAnalysis(.object)
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
  if(any(.quality_criterion %in% c("all", "rho_C_weighted"))) {
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
  
  class(out) <- "cSEMAssess"
  return(out)
}

#' @export

assess.cSEMResults_multi <- function(
  .object              = NULL,
  .only_common_factors = TRUE,
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_weighted", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj", "ra",
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
  .quality_criterion   = c("all", "ave", "rho_C", "rho_C_weighted", "cronbachs_alpha", 
                           "cronbachs_alpha_weighted", "dg", "dl", "dml", "df",
                           "esize", "cfi", "gfi", "ifi", "nfi", "nnfi", 
                           "rmsea", "rms_theta", "srmr",
                           "gof", "htmt", "r2", "r2_adj", "ra",
                           "rho_T", "rho_T_weighted", "vif", 
                           "vifmodeB",  "fl_criterion"),
  ...
  ){
  
  stop2("Currently, models containing second-order constructs are not supported by `assess()`.")
}