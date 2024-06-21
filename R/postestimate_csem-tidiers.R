#' tidy Method for cSEMResults Objects
#' 
#' This method and documentation is based on the `tidy` method for `lavaan`, as found in `broom:::tidy.lavaan()`.
#' 
#' **Cautionary Note**: As with all `tidy()` methods, mis-specified arguments may be silently ignored.
#' 
#' @param x Fitted `csem` model with class `cSEMResults` as returned from [cSEM::csem()]
#' @param parameters Character string. Selects the desired parameter estimates as outputted from [cSEM::summarize()]. Defaults to `all`.
#' @inheritDotParams summarize
#' @return A `data.frame` with one row for each estimated parameter. The
#'   description of the column names are as follows: 
#'   
#'   \item{term}{The result of `paste(lhs, op, rhs)`}
#'   \item{op}{The operator in the model syntax (e.g., `~~` for correlation, `<~` for weights, `=~` for loadings, and `~` for path estimates)}
#'   \item{group}{The group as specified by `.id` in the [cSEM::csem()].}
#'   \item{estimate}{The parameter estimate}
#'   \item{construct.type}{Whether the modelled construct is a common factor or composite variable}
#'   \item{std.error}{}
#'   \item{statistic}{The t-statistic as returned from [cSEM::summarize()] and `cSEM:::addInfer()`}
#'   \item{p_value}{}
#' 
#' @author Michael s. truong
#' 
#' @importFrom generics tidy
#' @export
#'
tidy.cSEMResults <- function(x,
                             parameters = c(
                               "all",
                               "Path_estimates",
                               "Loadings_estimates",
                               "Weight_estimates",
                               "D2",
                               "Effect_estimates",
                               "Exo_construct_correlation"
                             ),
                             ...) {

  parameters <- match.arg(
    parameters,
    c(
      "all",
      "Path_estimates",
      "Loadings_estimates",
      "Weight_estimates",
      "D2",
      "Effect_estimates",
      "Exo_construct_correlation"
    )
  )
  
  if (identical(names(x), c("Estimates", "Information"))) {
    
    if (identical(names(x), c("Estimates", "Information")) &
        identical(names(x[[1]]), c("Estimates", "Information"))) {
      stop(
        'This tidy method does not support group levels called "Estimates" and "Information", as this conflicts with internal cSEM functionality.'
      )
      
    }

    # Single Group Summary
    summarized_cSEMResults_list <- list(summarize(x, ...))
  } else {
    # Multi Group Summary
    summarized_cSEMResults_list <- lapply(x, summarize, ...)
    
  }
  
 out <- lapply(summarized_cSEMResults_list, function(summarized_cSEMResults, parameters) {
   
   ## Selection ---------------------------------------------------------------
   if (parameters == "all") {
     
     
     if (is.null(summarized_cSEMResults$Estimates$D2)) {
       general_estimates <- summarized_cSEMResults[["Estimates"]][c("Path_estimates", "Loading_estimates", "Weight_estimates", "Exo_construct_correlation")]
     } else {
       general_estimates <- summarized_cSEMResults[["Estimates"]][c("Path_estimates", "Loading_estimates", "Weight_estimates", "D2", "Exo_construct_correlation")]
     }
     
     
     effect_estimates <- summarized_cSEMResults[["Estimates"]][["Effect_estimates"]][c("Direct_effect", "Indirect_effect", "Total_effect")]
     
     # Add type of effect as column
     # See user1317221_G on April 26/2015 https://stackoverflow.com/a/29878732
     effect_estimates <- mapply(function(df, df_name) {
       df$op <- rep(df_name, nrow(df))
       return(df)
     },
     effect_estimates,
     names(effect_estimates),
     SIMPLIFY = FALSE)
     
     
     out <- c(general_estimates, effect_estimates)
     
   } else if (parameters != "Effect_estimates") {
     out <- summarized_cSEMResults[["Estimates"]][parameters] 
     
   } else if (parameters == "Effect_estimates") {
     effect_estimates <- summarized_cSEMResults[["Estimates"]][["Effect_estimates"]][c("Direct_effect", "Indirect_effect", "Total_effect")]
     
     # Add type of effect as column
     # See user1317221_G on April 26/2015 https://stackoverflow.com/a/29878732
     effect_estimates <- mapply(function(df, df_name) {
       df$op <- rep(df_name, nrow(df))
       return(df)
     },
     effect_estimates,
     names(effect_estimates),
     SIMPLIFY = FALSE)
     
     out <- effect_estimates
     
   } 
   
   # Take out the empty data.frames to facilitate the later concatenation
   out <- out[sapply(out, function(x)
     nrow(x) > 0)]
   # Add operand as op column
   out <- mapply(function(df, df_name) {
     if ((df_name %in% c("Direct_effect", "Indirect_effect", "Total_effect"))) {
       return(df)
     } else {
       
       # Test for what operand is consistent for the iterated dataframe
       ## This needs spaces so that we don't conflate ~~ with ~
       operands <- c(" =~ ", " <~ ", " ~~ ", " ~ ") 
       names(operands) <- c("=~", "<~", "~~", "~")
       df_type <- sapply(operands, function(x) grepl(x, df$Name[1]))
       
       if (length(names(df_type)[df_type]) == 1) {
         df$op <- rep(names(df_type)[df_type], nrow(df))
         
       } else {
         stop("There is more than one type of allowed operand in the iterated data.frame. Cannot determine what to assign to op column.")
         
       }
       
       return(df)
     }
   },
   out,
   names(out),
   SIMPLIFY = FALSE) 
   
   ## Conversion to broom Style -------------------------------------------------------------
   # Merge list of data.frames into one
   # See response by Charles on November 11/2011 https://stackoverflow.com/a/8097519
   out <- Reduce(function(...) merge(..., all = TRUE, sort = FALSE), out)
   
   # Rename columns to be consistent with tidy glossary standard
   # See Zon Shi Wu Jie from May 5/2014 https://stackoverflow.com/a/23475492
   colnames(out)[colnames(out) == "Name"] <- 'term'
   colnames(out)[colnames(out) == "Construct_type"] <- 'type'
   colnames(out)[colnames(out) == "Estimate"] <- 'estimate'
   colnames(out)[colnames(out) == "Std_err"] <- 'std.error'
   colnames(out)[colnames(out) == "t_stat"] <- 'statistic'
   colnames(out)[colnames(out) == "p_value"] <- 'p.value'
   
   # Column reordering
   
   out <- out[, c("term", "op", "estimate", "std.error", "statistic", "p.value")]
   
   return(out)
 },
 parameters = parameters
 )
  
 if(length(out) == 1) {
   out <- out[[1]]
 } else if (length(out) > 1) {
   out <- mapply(FUN = function(summary, name) {
     summary$group <- name
     summary<-summary[, c("term", "op", "group", "estimate", "std.error", "statistic", "p.value")]
     return(summary)
     
   }, 
   summary = out,
   name = names(out), SIMPLIFY = FALSE
   )
   
   out <- Reduce(function(...) merge(..., all = TRUE, sort = FALSE), out)
   
 }
  
   
  return(out)
  
}

#' glance Method for cSEMResults Objects
#' 
#' Placeholder function for future development. 
#' 
#' TODO: This method should be used to obtain a one row `data.frame` with a different column for every fit-statistic/estimation information.
#' 
#' This method and documentation is based on the `glance` method for `lavaan`, as found in `broom:::glance.lavaan()`.
#'
#' @inheritParams tidy.cSEMResults
#' @inheritDotParams assess
#'
#' @return A one-row `data.frame` with columns: 
#' @importFrom generics glance
#' @author Michael s. truong
glance.cSEMResults <- function(x, ...) {
  # TODO: May want to look at assess
  # cSEM:::assess()
  print("The glance method is currently not implemented in cSEM.")
}

#' augment Method for cSEMResults Objects
#'
#' Placeholder function for future development
#' 
#' TODO: This method should be used to obtain the construct scores and other information provided by [cSEM::predict()]
#' 
#' 
#' @inheritParams tidy.cSEMResults
#' @inheritDotParams predict
#' 
#' @return A `data.frame` of the original dataset alongside its construct scores and...
#' @importFrom generics augment
#' @author Michael s. truong
augment.cSEMResults <- function(x, ...) {
  # TODO: Look at predict and maybe also get the construct scores
  # cSEM::predict()
  print("The augment method is currently not implemented in cSEM.")
}