#' Takes GSCAPro input and Creates a Lavaan-style Table
#'
#' TODO: Rename the functions here
#' 
#' Assumes that every indicator loads onto only one latent variable (composite/factor)
#' 
#' Expects output from parse_GSCAPro_FullResults.
#' 
#' @param gscapro_in 
#' @param model 
#'
#' @return Table of GSCA Pro results of Weights, Loadings, Path Coefficients and Uniqueness Terms in lavaan style. 
#'
get_lavaan_table_igsca_gscapro <- function(gscapro_in, model) {
  
  table <- lavaan::lavaanify(model = model)[, c("lhs", "op", "rhs")]
  # Remove unnecessary rows
  table <- table[table$op %in% c("=~", "<~", "~"),]
  # Pre-allocate Columns
  table <-
    cbind(table, list(
      "weights" = 0,
      "loadings" = 0,
      "uniqueD" = 0,
      "paths" = 0
    ))
  
  
  # Correct names of latent variables
  gscapro_in$Weights$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  gscapro_in$Weights$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Weights$V1,
         replacement = "")
  
  
  
  lv_idx <- list()
  span <- list() 
  
  lv_idx$weights <- which(is.na(gscapro_in$Weights$Estimate))
  # Number of rows after the row of the LV that correspond to that LV's indicators
  span$weights <- diff(lv_idx$weights)
  gscapro_in$Weights$lv <- ""
  gscapro_in$Weights$lv <-
    c(with(gscapro_in$Weights, rep(V1[is.na(Estimate) &
                                        (nchar(V1) > 0)], times = span$weights)), "")
  
  
  # Retrieving indicators from gscapro
  gscapro_indicators <-
    with(gscapro_in$Weights, V1[!(V1 %in% lv) & (nchar(V1) > 0)])
  gscapro_lv <-
    with(gscapro_in$Weights, unique(lv)[nchar(unique(lv)) > 0])
  
  rownames(gscapro_in$Weights)<-gscapro_in$Weights$V1
  
  for (indicator in gscapro_indicators) {
    for (lv in gscapro_lv) {
      table[((table$lhs == lv &
                table$rhs == indicator) &
               table$op %in% c("<~", "=~")), "weights"] <- gscapro_in$Weights[indicator, "Estimate"]
    }
  }
  
  
  # Loadings Parser Shortcut
  gscapro_in$Loadings$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  gscapro_in$Loadings$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$Loadings$V1,
         replacement = "")
  
  if (identical(gscapro_in$Weights$V1, gscapro_in$Loadings$V1)) {
    rownames(gscapro_in$Loadings) <- gscapro_in$Loadings$V1
    
    for (indicator in gscapro_indicators) {
      for (lv in gscapro_lv) {
        table[((table$lhs == lv &
                  table$rhs == indicator) &
                 table$op %in% c("<~", "=~")), "loadings"] <-
          gscapro_in$Loadings[indicator, "Estimate"]
      }
    }
  } else {
    stop("Cannot take the shortcut of weight indicators for loadings")
  }
  
  # Paths need their own parser 
  
  gscapro_in$`Path coefficients` <-
    gscapro_in$`Path coefficients`[!is.na(gscapro_in$`Path coefficients`$Estimate), ]
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = " ",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "->",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = " ")
  
  gscapro_in$`Path coefficients`$V1 <-
    gsub(pattern = "-",
         x = gscapro_in$`Path coefficients`$V1,
         replacement = "")
  
  # ephpostfacto on November 6, 2009; editted by Jilber Urbina on Jan 1/2014 https://stackoverflow.com/a/1690753
  paths <-strsplit(gscapro_in$`Path coefficients`$V1, " ")
  gscapro_in$`Path coefficients`$lvfrom <-sapply(paths, FUN = \(.) .[1])
  gscapro_in$`Path coefficients`$lvto <-sapply(paths, FUN = \(.) .[2])
  
  for (lv_to_iter in gscapro_in$`Path coefficients`$lvto) {
    for (lv_from_iter in gscapro_in$`Path coefficients`$lvfrom) {
      table[((table$rhs == lv_from_iter &
                table$lhs == lv_to_iter) &
               table$op == "~"), "paths"] <- with(gscapro_in$`Path coefficients`, Estimate[(lvfrom == lv_from_iter) & (lvto == lv_to_iter)])
      
    }
  }
  
  # Parsing for Uniqueness Terms
  for (indicator in names(gscapro_in$`Unique D^2`)) {
    table[((table$rhs == indicator) &
             (table$op == "=~")), "uniqueD"] <- gscapro_in$`Unique D^2`[indicator]
  }
  
  # Remove zeros for cells that shouldn't have values
  table[!(table$op %in% c("<~", "=~")), "weights"] <- NA
  table[!(table$op %in% c("<~", "=~")), "loadings"] <- NA
  table[!(table$op %in% c("=~")), "uniqueD"] <- NA
  table[!(table$op %in% c("~")), "paths"] <- NA
  
  return(table)
}