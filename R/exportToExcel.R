#' Export to Excel (.xlsx)
#'
#' Export results from a postestimation function such as [assess()] or 
#' [summarize()] to an .xlsx (Excel) file. The function uses the openxlsx
#' package which does not depend on Java!
#' 
#' The function is deliberately kept simple: all it does is to take all the 
#' relevant elements in `.object` and write them (worksheet by worksheet) into 
#' an .xlsx file named `.filename` in the directory given by `.path` (the current 
#' working directory by default).
#' 
#' Note: rerunning [exportToExcel()] without changing `.filename` and `.path` 
#' overwrites the file!
#' 
#' @usage exportToExcel(
#'   .object    = NULL, 
#'   .filename  = "results.xlsx",
#'   .path      = NULL
#'   )
#'
#' @inheritParams csem_arguments
#'   
#' @seealso [assess()], [summarize()], [predict()], [verify()]
#'
#' @export
exportToExcel <- function(
  .object    = NULL, 
  .filename = "results.xlsx",
  .path      = NULL) {
  
  ## Install openxlsx if not already installed
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop2(
      "Package `openxlsx` required. Use `install.packages(\"openxlsx\")` and rerun.")
  }
  
  wb <- openxlsx::createWorkbook()
  
  ## Check class
  if(inherits(.object, "cSEMSummarize_default")) {
    est      <- .object$Estimates
    elements <- c("Path coefficients", "Loadings", "Weights", "Inner weights",
                  "Residual correlation", "Construct scores", 
                  "Indicator correlation matrix", "Composite correlation matrix",
                  "Construct correlation matrix", "Total Effects", "Indirect Effects")
    
    for(element in elements) {
      ## Add worksheets
      openxlsx::addWorksheet(wb, element)
    }

    ## Write data to worksheets
    openxlsx::writeData(wb, sheet = "Path coefficients", est$Path_estimates)
    openxlsx::writeData(wb, sheet = "Loadings", est$Loading_estimates)
    openxlsx::writeData(wb, sheet = "Weights", est$Weight_estimates)
    openxlsx::writeData(wb, sheet = "Inner weights", est$Inner_weight_estimates)
    openxlsx::writeData(wb, sheet = "Residual correlation", est$Residual_correlation)
    openxlsx::writeData(wb, sheet = "Construct scores", est$Construct_scores)
    openxlsx::writeData(wb, sheet = "Indicator correlation matrix", est$Indicator_VCV)
    openxlsx::writeData(wb, sheet = "Composite correlation matrix", est$Proxy_VCV)
    openxlsx::writeData(wb, sheet = "Construct correlation matrix", est$Construct_VCV)
    openxlsx::writeData(wb, sheet = "Total Effects", est$Effect_estimates$Total_effect)
    openxlsx::writeData(wb, sheet = "Indirect Effects", est$Effect_estimates$Indirect_effect)
    
  } else if(inherits(.object, "cSEMAssess")) {
    elements <- c("AVE", "R2", "R2_adj", "Reliability",
                  "Distance and Fit measures", "Model selection criteria", 
                  "VIFs", "Effect sizes", "HTMT", "Fornell-Larcker matrix")
    
    for(element in elements) {
      ## Add worksheets
      openxlsx::addWorksheet(wb, element)
    }
    
    ## Write data to worksheets
    openxlsx::writeData(wb, sheet = "AVE", data.frame("Name" = names(.object$AVE), "AVE" = .object$AVE))
    openxlsx::writeData(wb, sheet = "R2", data.frame("Name" = names(.object$R2), "R2" = .object$R2))
    openxlsx::writeData(wb, sheet = "R2_adj", data.frame("Name" = names(.object$R2_adj), "R2_adj" = .object$R2_adj))
    
    # Reliability
    d <- data.frame(
      "Name" = names(.object$Reliability$Cronbachs_alpha),
      "Cronbachs_alpha" = .object$Reliability$Cronbachs_alpha,
      "Joereskogs_rho"  = .object$Reliability$Joereskogs_rho,
      "Dijkstra-Henselers_rho_A" = .object$Reliability$`Dijkstra-Henselers_rho_A`,
      "RhoC" = .object$RhoC,
      "RhoC_mm" = .object$RhoC_mm,
      "RhoC_weighted" = .object$RhoC_weighted,
      "RhoC_weighted_mm" = .object$RhoC_weighted_mm,
      "RhoT" = .object$RhoT,
      "RhoT_weighted" = .object$RhoT_weighted
    )
    openxlsx::writeData(wb, sheet = "Reliability", d)
    
    # Distance and fit
    d <- data.frame(
      "Geodesic distance"          = .object$DG,
      "Squared Euclidian distance" = .object$DL,
      "ML distance"                = .object$DML,
      "Chi_square"                 = .object$Chi_square,
      "Chi_square_df"              = .object$Chi_square_df,
      "CFI"                        = .object$CFI,
      "CN"                         = .object$CN,
      "GFI"                        = .object$GFI,
      "IFI"                        = .object$IFI,
      "NFI"                        = .object$NFI,
      "NNFI"                       = .object$NNFI,
      "RMSEA"                      = .object$RMSEA,
      "RMS_theta"                  = .object$RMS_theta,
      "SRMR"                       = .object$SRMR,
      "Degrees of freedom"         = .object$Df
    )
    openxlsx::writeData(wb, sheet = "Distance and Fit measures", d)
    
    # Model selection criteria
    d <- data.frame(
      "Name" = names(.object$AIC),
      "AIC"  = .object$AIC,
      "AICc" = .object$AICc,
      "AICu" = .object$AICu,
      "BIC"  = .object$BIC,
      "FPE"  = .object$FPE,
      "GM"   = .object$GM,
      "HQ"   = .object$HQ,
      "HQc"  = .object$HQc,
      "Mallows_Cp" = .object$Mallows_Cp
    )
    openxlsx::writeData(wb, sheet = "Model selection criteria", d)
    openxlsx::writeData(wb, sheet = "VIFs", data.frame("Name" = rownames(.object$VIF), .object$VIF))
    openxlsx::writeData(wb, sheet = "Effect sizes", data.frame("Name" = rownames(.object$F2), .object$F2))
    openxlsx::writeData(wb, sheet = "HTMT", .object$HTMT)
    openxlsx::writeData(wb, sheet = "Fornell-Larcker matrix", .object$`Fornell-Larcker`)
    
    
  } else {
    stop2(
      "The following error occured in the exportToExcel() function:\n",
      paste0("`.object` of class '", class(.object), "' not supported.")
    )
  }
  
  ## Save workbook
  if(is.null(.path)) {
    .path = getwd()
  } else {
    .path = normalizePath(.path)
  }
  
  file = file.path(.path, .filename)
  
  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  print(paste0("File saved as: '", .filename, "' in: ", .path))
}