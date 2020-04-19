#' Inference
#'
#' Calculate common inferential quantities. For users interested in the
#' estimated standard errors, t-values, p-values and/or confidences
#' intervals of the path, weight or loading estimates, calling [summarize()]
#' directly will usually be more convenient as it has a much more 
#' user-friendly print method. [infer()] is useful for comparing 
#' different confidence interval estimates.
#'
#' [infer()] is a convenience wrapper around a 
#' number of internal functions that compute a particular inferential
#' quantity, i.e., a value or set of values to be used in statistical inference.
#' 
#' \pkg{cSEM} relies on resampling (bootstrap and jackknife) as the basis for 
#' the computation of e.g., standard errors or confidence intervals.
#' Consequently, [infer()] requires resamples to work. Technically, 
#' the [cSEMResults] object used in the call to [infer()] must 
#' therefore also have class attribute `cSEMResults_resampled`. If 
#' the object provided by the user does not contain resamples yet,
#' [infer()] will obtain bootstrap resamples first. 
#' Naturally, computation will take longer in this case.
#' 
#' [infer()] does as much as possible in the  background. Hence, every time 
#' [infer()] is called on a [cSEMResults] object the quantities chosen by 
#' the user are automatically computed for every estimated parameter 
#' contained in the object. By default all possible quantities are 
#' computed (`.quantity = all`). The following table list the available 
#' inferential quantities alongside a brief description. Implementation and 
#' terminology of the confidence intervals is based on 
#' \insertCite{Hesterberg2015;textual}{cSEM} and 
#' \insertCite{Davison1997;textual}{cSEM}.
#' \describe{
#' \item{`"mean"`, `"sd"`}{The mean or the standard deviation 
#'   over all `M` resample estimates of a generic statistic or parameter.}
#' \item{`"bias"`}{The difference between the resample mean and the original 
#'   estimate of a generic statistic or parameter.}
#' \item{`"CI_standard_z"` and `"CI_standard_t"`}{The standard confidence interval 
#'   for a generic statistic or parameter with standard errors estimated by 
#'   the resample standard deviation. While `"CI_standard_z"` assumes a 
#'   standard normally distributed statistic,
#'   `"CI_standard_t"` assumes a t-statistic with N - 1 degrees of freedom.}
#' \item{`"CI_percentile"`}{The percentile confidence interval. The lower and 
#'   upper bounds of the confidence interval are estimated as the alpha and 
#'   1-alpha quantiles of the distribution of the resample estimates.}
#' \item{`"CI_basic"`}{The basic confidence interval also called the reverse 
#'   bootstrap percentile confidence interval. See \insertCite{Hesterberg2015;textual}{cSEM}
#'   for details.}
#' \item{`"CI_bc"`}{The bias corrected (Bc) confidence interval. See 
#'   \insertCite{Davison1997;textual}{cSEM} for details.}
#' \item{`"CI_bca"`}{The bias-corrected and accelerated (Bca) confidence interval.
#'   Requires additional jackknife resampling to compute the influence values. 
#'   See \insertCite{Davison1997;textual}{cSEM} for details.}
#' \item{`"CI_t_interval"`}{The "studentized" t-confidence interval. If based on bootstrap 
#'   resamples the interval is also called the bootstrap t-interval 
#'   confidence interval. See \insertCite{Hesterberg2015;textual}{cSEM} on page 381. 
#'   Requires resamples of resamples. See [resamplecSEMResults()].}
#' }
#' 
#' By default, all but the studendized t-interval confidence interval and the
#' bias-corrected and accelerated confidence interval are calculated. The 
#' reason for excluding these quantities by default are that both require 
#' an additional resampling step. The former requires 
#' jackknife estimates to compute influence values and the latter requires 
#' double bootstrap. Both can potentially be time consuming. 
#' Hence, computation is triggered only if explicitly chosen.
#' 
#' @usage infer(
#'  .object            = NULL,
#'  .quantity          = c("all", "mean", "sd", "bias", "CI_standard_z", 
#'                         "CI_standard_t", "CI_percentile", "CI_basic", 
#'                         "CI_bc", "CI_bca", "CI_t_interval"),
#'  .alpha             = 0.05,
#'  .bias_corrected    = TRUE
#' )
#'
#' @inheritParams csem_arguments
#' 
#' @return A list of class `cSEMInfer`.
#' 
#' @references
#'   \insertAllCited{} 
#' 
#' @seealso [csem()], [resamplecSEMResults()], [summarize()] [cSEMResults]
#' 
#' @example inst/examples/example_infer.R
#' @export

infer <- function(
  .object          = NULL,
  .quantity        = c("all", "mean", "sd", "bias", "CI_standard_z", 
                       "CI_standard_t", "CI_percentile", "CI_basic", 
                       "CI_bc", "CI_bca", "CI_t_interval"),
  .alpha           = 0.05,
  .bias_corrected  = TRUE
) {
  
  ## Check arguments
  match.arg(.quantity, args_default(.choices = TRUE)$.quantity, several.ok = TRUE)
  
  ## Check if "all" is part of .quantity. If yes, set .quantity = "all"
  if(any(.quantity == "all")) {
    .quantity <- "all"
  }
  
  
  if(!inherits(.object, "cSEMResults_resampled")) {
    # Bootstrap if necessary
    .object <- resamplecSEMResults(.object)
  }
  
  if(inherits(.object, "cSEMResults_multi")) {
    ## If multi object, do recursive call
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
  } else if(inherits(.object, "cSEMResults_2ndorder")) {
    first_resample  <- .object$Second_stage$Information$Resamples$Estimates$Estimates1
    second_resample <- .object$Second_stage$Information$Resamples$Estimates$Estimates2
    info            <- .object$Second_stage$Information$Resamples$Information_resample
  } else if(inherits(.object, "cSEMResults_default")) {
    first_resample  <- .object$Estimates$Estimates_resample$Estimates1
    second_resample <- .object$Estimates$Estimates_resample$Estimates2
    info            <- .object$Information$Information_resample
  } else {
    stop2("The following error occured in the `infer()` function:\n",
          "Object must be of class `cSEMResults`")
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