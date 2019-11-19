#' Inference
#'
#' Calculate common inferential quantities. For ssers interested in the
#' estimated standard errors and/or confidences intervals [summarize()] will usually
#' be more helpful as it has a much more user-friendly print method.
#'
#' Calculate common inferential quantities (e.g., estimated standard errors, estimated bias,
#' confidence intervals) based on a `cSEMResults_resampled` object as obtained
#' by calling [resamplecSEMResults()] or by setting `.resample_method = "bootstrap"`
#' or `"jackknife"` when calling [csem()]. Currently, the following quantities are
#' returned by default (`.quantity = "all"`):
#' \describe{
#' \item{`"mean"`, `"sd"` and `"bias"`}{The mean, the standard 
#'   deviation and the estimated bias (defined as the difference between the resample mean
#'   and the original estimate).}
#' \item{`"CI_standard_z"` and `"CI_standard_t"`}{The standard confidence interval 
#'   with standard errors estimated by the resample standard deviation. 
#'   While `"CI_standard_z"` assumes a standard normally distributed statistic,
#'   `"CI_standard_t"` assumes a t-statistic with N - 1 degrees of freedom.}
#' \item{`"CI_percentile"`}{The percentile confidence interval}
#' \item{`"CI_basic"`}{The basic confidence interval}
#' \item{`"CI_bc"`}{The bias corrected confidence interval}
#' }
#' 
#' In addition, the bias-corrected and accelerated (`.quantity = "CI_bca"`) and/or the "studentized"
#' confidence interval (`".quantity = CI_t_interval"`) may be computed. The former requires 
#' jackknife estimates to compute influence values and the latter requires 
#' double bootstrap. Both can potentially be time consuming. Hence, computation is triggerd 
#' only if explicitly chosen.
#' 
#' @usage infer(
#'  .object            = NULL,
#'  .alpha             = 0.05,
#'  .bias_corrected    = TRUE,
#'  .quantity          = c("all", "mean", "sd", "bias", "CI_standard_z", 
#'                         "CI_standard_t", "CI_percentile", "CI_basic", 
#'                         "CI_bc", "CI_bca", "CI_t_interval")
#' )
#'
#' @inheritParams csem_arguments
#' 
#' @return A list of class `cSEMInfer`.
#' 
#' @seealso [csem()], [resamplecSEMResults()], [summarize()] [cSEMResults]
#' 
#' @example inst/examples/example_infer.R
#' @export
#'

infer <- function(
  .object = NULL,
  .alpha           = 0.05,
  .bias_corrected  = TRUE,
  .quantity        = c("all", "mean", "sd", "bias", "CI_standard_z", "CI_standard_t",
                       "CI_percentile", "CI_basic", "CI_bc", "CI_bca", "CI_t_interval")
) {
  
  ## Check arguments
  match.arg(.quantity, args_default(.choices = TRUE)$.quantity, several.ok = TRUE)
  
  ## Check if "all" is part of .quantity. If yes, set .quantity = "all"
  if(any(.quantity == "all")) {
    .quantity <- "all"
  }
  
  if(!inherits(.object, "cSEMResults")) {
    stop2("The following error occured in the `infer()` function:\n",
          "Object must be of class `cSEMResults`")
  }
  
  ## If multi object, do recursive call
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, function(x) {
      infer(
        .object = x,
        .alpha = .alpha,
        .bias_corrected = .bias_corrected,
        .quantity = .quantity
      )
    })
    
    ## Add/ set class
    class(out) <- c("cSEMInfer", "cSEMInfer_multi")
    return(out)
  }
  
  if(!inherits(.object, "cSEMResults_resampled")) {
    stop2("The following error occured in the `infer()` function:\n",
          "Object must contain resamples.", 
          " Use `resamplecSEMResults(.object = .object, ...)` first."
    )
  }
  
  if(any(class(.object) == "cSEMResults_2ndorder")) {
    first_resample  <- .object$Second_stage$Information$Resamples$Estimates$Estimates1
    second_resample <- .object$Second_stage$Information$Resamples$Estimates$Estimates2
    info            <- .object$Second_stage$Information$Resamples$Information_resample
  } else {
    first_resample  <- .object$Estimates$Estimates_resample$Estimates1
    second_resample <- .object$Estimates$Estimates_resample$Estimates2
    info            <- .object$Information$Information_resample
  }
  
  ## Compute quantiles/critical values -----------------------------------------
  probs  <- c()
  .alpha <- .alpha[order(.alpha)]
  for(i in seq_along(.alpha)) { 
    if(.alpha[i] < 0 | .alpha[i] > 1) {
      stop2("The following error occured in the `infer()` function:\n",
            "`.alpha` must be between 0 and 1.")
    }
    # Both two sided and one sided confidence intervalls may be needed.
    # Therefore for every alpha four values will be put in a vector 
    # 1. alpha
    # 2. 1 - alpha
    # 3. alpha/2
    # 4. 1 - alpha/2
    # to make sure the corresponding quantile is available later on. Values are
    # round to four digits.
    probs <- c(probs, 
               # round(.alpha[i], 4), round(1 - .alpha[i], 4),
               .alpha[i]/2, 1 - .alpha[i]/2) 
  }
  
  ## Compute statistics and quantities
  out <- list()
  
  if(any(.quantity %in% c("all", "mean"))) {
    out[["mean"]] <- MeanResample(first_resample)
  }
  
  if(any(.quantity %in% c("all", "sd"))) {
    out[["sd"]] <- SdResample(
      .first_resample  = first_resample, 
      .resample_method = info$Method, 
      .n               = info$Number_of_observations
    )
  }
  
  if(any(.quantity %in% c("all", "bias"))) {
    out[["bias"]] <- BiasResample(
      .first_resample  = first_resample, 
      .resample_method = info$Method, 
      .n               = info$Number_of_observations
    )
  }
  
  if(any(.quantity %in% c("all", "CI_standard_z"))) {
    out[["CI_standard_z"]] <- StandardCIResample(
      .first_resample = first_resample, 
      .bias_corrected = .bias_corrected,
      .df             = NULL,
      .dist           = "z",
      .resample_method= info$Method, 
      .n              = info$Number_of_observations,
      .probs          = probs
    )
  }
  
  if(any(.quantity %in% c("all", "CI_standard_t"))) {
    out[["CI_standard_t"]] <- StandardCIResample(
      .first_resample = first_resample, 
      .bias_corrected = .bias_corrected,
      .df             = "type1",
      .dist           = "t",
      .resample_method= info$Method, 
      .n              = info$Number_of_observations,
      .probs          = probs
    )
  }
  
  if(any(.quantity %in% c("all", "CI_percentile"))) {
    out[["CI_percentile"]] <- PercentilCIResample(
      .first_resample = first_resample, 
      .probs = probs
    )
  }
  
  if(any(.quantity %in% c("all", "CI_basic"))) {
    out[["CI_basic"]] <- BasicCIResample(
      .first_resample = first_resample, 
      .bias_corrected = .bias_corrected,
      .probs          = probs
    )
  }
  
  if(any(.quantity %in% c("all", "CI_bc"))) {
    out[["CI_bc"]] <- BcCIResample(first_resample, probs)
  }
  
  if(any(.quantity == "CI_bca")) {
    out[["CI_bca"]] <- BcaCIResample(.object = .object, 
                                     first_resample, probs)
  }
  
  if(any(.quantity %in% c("all", "CI_t_interval"))) {
    if(!anyNA(second_resample)) {
      out[["CI_t_interval"]] <- TStatCIResample(
        .first_resample     = first_resample, 
        .second_resample    = second_resample, 
        .bias_corrected     = .bias_corrected,
        .resample_method    = info$Method, 
        .resample_method2   = info$Method2, 
        .n                  = info$Number_of_observations, 
        .probs              = probs
      ) 
    } else if(any(.quantity == "CI_t_interval")) {
      stop2("The following error occured in the `infer()` function:\n",
            "`CI_t_interval` requires (jackknife) resamples for each resample.",
            " Rerun your original estimation using .resample_method2 = 'jackknife'.",
            " and try again.")
    }
  }
  
  out <- purrr::transpose(out)
  ## Add/ set class
  class(out) <- c("cSEMInfer")
  out
}