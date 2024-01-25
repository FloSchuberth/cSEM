#' Writes the extracted matrices to run with igsca_sim_test.m
#'
#' @param extracted_matrices Object returned by extract_parseModel
#'
#' @return Writes the matrices from extract_parseModel() to the appropriate test directory for the Matlab implementation of igsca_sim to run.
#' @importFrom here here
#' @examples
#' 
#' 
#' tutorial_igsca_model <- "
#' # Composite Model
#' NetworkingBehavior <~ Behavior1 + Behavior2 + Behavior3 + Behavior5 + Behavior7 + Behavior8 + Behavior9
#' Numberofjobinterviews <~ Interview1 + Interview2
#' Numberofjoboffers <~ Offer1 + Offer2 
#' 
#' # Reflective Measurement Model
#' HonestyHumility =~ Honesty1 + Honesty2 + Honesty3 + Honesty4 + Honesty5 + Honesty6 + Honesty7 + Honesty8 + Honesty9 + Honesty10
#' Emotionality =~ Emotion1 + Emotion2 + Emotion3 + Emotion4 + Emotion5 + Emotion6 + Emotion8 + Emotion10
#' Extraversion =~ Extraver2 + Extraver3 + Extraver4 + Extraver5 + Extraver6 + Extraver7 + Extraver8 + Extraver9 + Extraver10
#' Agreeableness =~ Agreeable1 + Agreeable3 + Agreeable4 + Agreeable5 + Agreeable7 + Agreeable8 + Agreeable9 + Agreeable10
#' Conscientiousness =~ Conscientious1 + Conscientious3 + Conscientious4 + Conscientious6 + Conscientious7 + Conscientious8 + Conscientious9 + Conscientious10
#' OpennesstoExperience =~ Openness1 + Openness2 + Openness3 + Openness5 + Openness7 + Openness8 + Openness9 + Openness10
#' 
#' # Structural Model
#' NetworkingBehavior ~ HonestyHumility + Emotionality + Extraversion + Agreeableness + Conscientiousness + OpennesstoExperience
#' Numberofjobinterviews ~ NetworkingBehavior
#' Numberofjoboffers ~ NetworkingBehavior
#' "
#' 
#' data("LeDang2022")
#' 
#' write_for_matlab(extract_parseModel(model = tutorial_igsca_model, data = LeDang2022, ind_domi_as_first = TRUE))
write_for_matlab <- function(extracted_matrices) {
  indir <- list("tests", "comparisons", "igsca_translation", "matlab_in")
  extracted_matrices$lv_type <-
    as.numeric(extracted_matrices$lv_type)
  extracted_matrices$ov_type <-
    as.numeric(extracted_matrices$ov_type)
  
  mapply(
    write.csv,
    x = extracted_matrices,
    file = paste0(here::here(indir), "/", names(extracted_matrices), ".csv"),
    row.names = FALSE
  )
  
  invisible()
}

#' Helper Function for parsing GSCA Pro V1.2.1 Full Results Files 
#' 
#' The commented headers require more work to get working
#' 
#' @return List of extracted matrices from GSCA Pro Full results
#' @export
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom here here
#' 
#' @examples
#' parse_GSCAPro_FullResults()
parse_GSCAPro_FullResults <-
  function(file_path = list(
    "tests",
    "comparisons",
    "igsca_translation",
    "GSCAPro_1_2_1Output",
    "Tutorial_IGSCA_model1_full_result_ver0.csv"
  ),
  headers = c(# "Model fit measures",
    "Weights",
    "Loadings",
    "Path coefficients",
    # "Component/Factor correlations",
    # "HTMT",
    # "Construct quality measures",
    # "Fornell-Larcker criterion values",
    # "R squared values of indicators in measurement model",
    # "R squared values of component/correlations in structural model",
    # "VIFs",
    # "F squared values",
    # "Sample correlations (lower diagonal) & Residual correlations (upper diagonal)",
    # "Correlations between indicators and component/correlations",
    "Unique D^2")) {
    file <- here::here(file_path)
    
    to_parse <-
      data.table::fread(file = file, fill = TRUE, sep = ",")
    
    col_1_strings <- unlist(to_parse[, 1])
    idx_of_headers <- which((col_1_strings %in% headers) == TRUE)
    
    idx_of_blank_rows <- which((col_1_strings %in% "") == TRUE)
    
    # Gets the idx of idx_of_blank_rows of the first blank that occurs after the header
    list_of_blanks <-
      lapply(idx_of_headers, FUN = \(x) which.max(x < idx_of_blank_rows))
    idx_first_blank_after_header <-
      lapply(list_of_blanks, FUN = \(x) idx_of_blank_rows[x]) |>
      unlist()
    
    # Add one to headers because the first row with the column names is one row after the header name
    nrows_of_tables <-
      idx_first_blank_after_header - (idx_of_headers + 1)
    
    if (length(idx_of_headers) != length(nrows_of_tables)) {
      stop("The number of rows and the the number of indices of headers does not match")
    }
    
    tables <- mapply(
      FUN = data.table::fread,
      file = file,
      fill = TRUE,
      skip = idx_of_headers + 1,
      # Add one because in Version 1.2.1 of GSCA Pro the Header was one row ahead of the column labels
      nrows = nrows_of_tables,
      # Minus one for adjusting for the idx_of_headers and another for the empty-gaps between tables
      sep = ","
    )
    names(tables) <- headers
    
    # Write Tables to Directory
    mapply(
      FUN = data.table::fwrite,
      x = tables,
      file = here::here(file_path[1:(length(file_path) - 1)], paste0(names(tables), ".csv"))
    )
    
    # FIXME: get_lavaan_table_igsca_gscapro was giving me some trouble for reasons I don't quite understand
    # The strangeness is that this issue didn't seem to occur before I changed data.table to be an explicit dependency?
    tables <- lapply(tables, as.data.frame)
    
    return(tables)
    
  }