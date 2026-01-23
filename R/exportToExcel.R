#' Export to Excel (.xlsx)
#'
#' \lifecycle{experimental}
#' 
#' Export results from postestimation functions [assess()], [predict()],
#' [summarize()] and [testOMF()] to an .xlsx (Excel) file. The function uses the \link[openxlsx]{openxlsx}
#' package which does not depend on Java!
#' 
#' The function is deliberately kept simple: it takes all the 
#' relevant elements in `.postestimation_object` and writes them (worksheet by worksheet) into 
#' an .xlsx file named `.filename` in the directory given by `.path` (the current 
#' working directory by default).
#' 
#' If `.postestimation_object` has class attribute `_2ndorder` two .xlsx files
#' named `".filename_first_stage.xlsx"` and `".filename_second_stage.xlsx"`
#' are created. If `.postestimation_object` is a list of appropriate objects,
#' one file for each list elements is created.
#' 
#' Note: rerunning [exportToExcel()] without changing `.filename` and `.path` 
#' overwrites the file!
#' 
#' @usage exportToExcel(
#'   .postestimation_object = NULL, 
#'   .filename              = "results.xlsx",
#'   .path                  = NULL
#'   )
#'
#' @inheritParams csem_arguments
#'   
#' @seealso [assess()], [summarize()], [predict()], [testOMF()]
#'
#' @export
exportToExcel <- function(
  .postestimation_object  = NULL, 
  .filename               = "results.xlsx",
  .path                   = NULL) {
  
  ## Install openxlsx if not already installed
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop2(
      "Package `openxlsx` required. Use `install.packages(\"openxlsx\")` and rerun.")
  }
  
  if(inherits(.postestimation_object, c("list", "cSEMPredict_multi"))) {
    for(i in seq_along(.postestimation_object)) {
      filename = paste0(gsub(".xlsx", "", .filename), "_", i, ".xlsx")
      exportToExcel(.postestimation_object[[i]], .filename = filename, .path = .path)
    }
    return(invisible())
  }
  
  wb <- openxlsx::createWorkbook()
  
  ### summarize-----------------------------------------------------------------
  ## Check class
  if(inherits(.postestimation_object, "cSEMSummarize")) {
    
    if(inherits(.postestimation_object, "cSEMSummarize_2ndorder")) {
      for(name in names(.postestimation_object)) {
        filename = paste0(gsub(".xlsx", "", .filename), "_", name, ".xlsx")
        exportToExcel(.postestimation_object[[name]], .filename = filename, .path = .path)
      }
      return(invisible())
    }

    elements <- c("Path coefficients", "Loadings", "Weights", "Inner weights",
                  "Residual correlation", "Construct scores", 
                  "Indicator correlation matrix", "Composite correlation matrix",
                  "Construct correlation matrix", "Total Effects", "Indirect Effects")
    
    est      <- .postestimation_object$Estimates
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
    
  } else if(inherits(.postestimation_object, "cSEMAssess")) {
    elements <- c("AVE", "R2", "R2_adj", "Reliability",
                  "Distance and Fit measures", "Model selection criteria", 
                  "VIFs", "Effect sizes", "HTMT", "HTMT2", "Fornell-Larcker matrix")
    
    for(element in elements) {
      ## Add worksheets
      openxlsx::addWorksheet(wb, element)
    }
    # For ease of implementation: 
    #  If the user selected only a subset of quality criteria via .quality_criterion
    #  some quality criteria are not in .postestimate_object. This would require
    #  to write an "if present then use, else dont" check, which is just a lot
    #  of code for probably a rather rare instance. Therefore, if .quality_criterion != "all",
    #  the user is asked to run assess() again with .quality_criterion = "all" to make sure all
    #  elements are present.
    if(!.postestimation_object$Information$All) {
      stop2("Rerun assess() with .quality_criterion = 'all' and try again.")
    }
    
    ## Write data to worksheets
    openxlsx::writeData(wb, sheet = "AVE", data.frame("Name" = names(.postestimation_object$AVE), "AVE" = .postestimation_object$AVE))
    openxlsx::writeData(wb, sheet = "R2", data.frame("Name" = names(.postestimation_object$R2), "R2" = .postestimation_object$R2))
    openxlsx::writeData(wb, sheet = "R2_adj", data.frame("Name" = names(.postestimation_object$R2_adj), "R2_adj" = .postestimation_object$R2_adj))
    
    # Reliability
    d <- data.frame(
      "Name" = names(.postestimation_object$Reliability$Cronbachs_alpha),
      "Cronbachs_alpha" = .postestimation_object$Reliability$Cronbachs_alpha,
      "Joereskogs_rho"  = .postestimation_object$Reliability$Joereskogs_rho,
      "Dijkstra-Henselers_rho_A" = .postestimation_object$Reliability$`Dijkstra-Henselers_rho_A`,
      "RhoC" = .postestimation_object$RhoC,
      "RhoC_mm" = .postestimation_object$RhoC_mm,
      "RhoC_weighted" = .postestimation_object$RhoC_weighted,
      "RhoC_weighted_mm" = .postestimation_object$RhoC_weighted_mm,
      "RhoT" = .postestimation_object$RhoT,
      "RhoT_weighted" = .postestimation_object$RhoT_weighted
    )
    openxlsx::writeData(wb, sheet = "Reliability", d)
    
    # Distance and fit
    d <- data.frame(
      "Geodesic distance"          = .postestimation_object$DG,
      "Squared Euclidean distance" = .postestimation_object$DL,
      "ML distance"                = .postestimation_object$DML,
      "Chi_square"                 = .postestimation_object$Chi_square,
      "Chi_square_df"              = .postestimation_object$Chi_square_df,
      "CFI"                        = .postestimation_object$CFI,
      "CN"                         = .postestimation_object$CN,
      "GFI"                        = .postestimation_object$GFI,
      "IFI"                        = .postestimation_object$IFI,
      "NFI"                        = .postestimation_object$NFI,
      "NNFI"                       = .postestimation_object$NNFI,
      "RMSEA"                      = .postestimation_object$RMSEA,
      "RMS_theta"                  = if(is.null(.postestimation_object$RMS_theta)) {
        NA
      } else {
        .postestimation_object$RMS_theta
      },
      "SRMR"                       = .postestimation_object$SRMR,
      "Degrees of freedom"         = .postestimation_object$Df
    )
    openxlsx::writeData(wb, sheet = "Distance and Fit measures", d)
    
    # Model selection criteria
    d <- data.frame(
      "Name" = names(.postestimation_object$AIC),
      "AIC"  = .postestimation_object$AIC,
      "AICc" = .postestimation_object$AICc,
      "AICu" = .postestimation_object$AICu,
      "BIC"  = .postestimation_object$BIC,
      "FPE"  = .postestimation_object$FPE,
      "GM"   = .postestimation_object$GM,
      "HQ"   = .postestimation_object$HQ,
      "HQc"  = .postestimation_object$HQc,
      "Mallows_Cp" = .postestimation_object$Mallows_Cp
    )
    openxlsx::writeData(wb, sheet = "Model selection criteria", d)
    openxlsx::writeData(wb, sheet = "VIFs", data.frame("Name" = rownames(.postestimation_object$VIF), .postestimation_object$VIF))
    openxlsx::writeData(wb, sheet = "Effect sizes", data.frame("Name" = rownames(.postestimation_object$F2), .postestimation_object$F2))
    openxlsx::writeData(wb, sheet = "HTMT", if(is.null(.postestimation_object$HTMT)) {
      NA
    } else {
      .postestimation_object$HTMT$htmts
    })
    openxlsx::writeData(wb, sheet = "HTMT2", if(is.null(.postestimation_object$HTMT2)) {
      NA
    } else {
      .postestimation_object$HTMT2$htmts
    })
    openxlsx::writeData(wb, sheet = "Fornell-Larcker matrix", if(is.null(.postestimation_object$`Fornell-Larcker`)) {
      NA
    } else {
      .postestimation_object$`Fornell-Larcker`
    })
  
  } else if(inherits(.postestimation_object, "cSEMPredict")) {
    for(element in names(.postestimation_object)) {
      ## Add worksheets
      openxlsx::addWorksheet(wb, element)
    }
    openxlsx::writeData(wb, sheet = "Information", data.frame(
      "Info"  = names(.postestimation_object$Information),
      "Value" = unname(unlist(.postestimation_object$Information)) 
    ))
    
    for(x in setdiff(names(.postestimation_object), "Information")) {
      openxlsx::writeData(wb, sheet = x, .postestimation_object[[x]])
    }
  } else if(inherits(.postestimation_object, "cSEMTestOMF")) {
    for(element in c(names(.postestimation_object), "Bootstrap_values")) {
      ## Add worksheets
      openxlsx::addWorksheet(wb, element)
    }
    # Remove NA's in bootstrap values to be able bind rows
    boot_values <- Filter(Negate(anyNA), .postestimation_object$Information$Bootstrap_values)
    openxlsx::writeData(wb, sheet = "Bootstrap_values", dplyr::bind_rows(boot_values))
    openxlsx::writeData(wb, sheet = "Information", data.frame(
      "Info"  = setdiff(names(.postestimation_object$Information), "Bootstrap_values"),
      "Value" = unname(unlist(.postestimation_object$Information[setdiff(names(.postestimation_object$Information), "Bootstrap_values")])) 
    ))
    openxlsx::writeData(wb, sheet = "Test_statistic", data.frame(
      "Distance measure"  = names(.postestimation_object$Test_statistic),
      "Value" = .postestimation_object$Test_statistic,
      check.names = FALSE
    ))
    openxlsx::writeData(wb, sheet = "Critical_value", 
                        data.frame("Distance measure" = rownames(.postestimation_object$Critical_value), 
                                   .postestimation_object$Critical_value, check.names = FALSE))
    openxlsx::writeData(wb, sheet = "Decision", 
                        data.frame("Distance measure" = rownames(.postestimation_object$Decision), 
                                   .postestimation_object$Decision, check.names = FALSE))
  } else {
    stop2(
      "The following error occured in the exportToExcel() function:\n",
      paste0("`.postestimation_object` of class '", class(.postestimation_object), "' not supported.")
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