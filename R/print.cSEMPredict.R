#' `cSEMPredict` method for `print()`
#'
#' The `cSEMPredict` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [predict()]
#'
#' @export
#' @keywords internal
print.cSEMPredict <- function(x, ...) {
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview"), 
    "\n"
  )
  
  if(inherits(x, "cSEMPredict_multi")) {
    l <- max(nchar(names(x)), nchar("Dataset")) + 2
    
    cat2(
      "\n\t", "Prediction results have been saved for each data set.\n",
      "\tUse `<object_name>$...` to print the results."
    )
  } else {
    
    x1 <- x$Prediction_metrics
    x2 <- x$Information
    
    # cat2("\n\tGeneral information:\n\t","------------------------")
    cat2(
      col_align("\n\tNumber of obs. training", 35), "= ", x2$Number_of_observations_training,
      col_align("\n\tNumber of obs. test", 35), "= ", x2$Number_of_observations_test,
      col_align("\n\tNumber of cv folds", 35), "= ", x2$Number_of_folds,
      col_align("\n\tNumber of repetitions", 35), "= ", x2$Number_of_repetitions,
      col_align("\n\tHandle inadmissibles", 35), "= ", x2$Handle_inadmissibles,
      col_align("\n\tBenchmark", 35), "= ", paste0("'", x2$Benchmark, "'")
    )
    
    ### Prediction metricts-------------------------------------------------------
    cat2("\n\n", rule2("Prediction metrics"), "\n\n")
    
    l <- max(nchar(x1[, "Name"]))
    
    cat2(
      "\n  ", 
      col_align("Name", l + 2), 
      col_align("MAE target", 13, align = "right"), 
      col_align("MAE benchmark", 15, align = "right"), 
      col_align("RMSE target", 13, align = "right"),
      col_align("RMSE benchmark", 15, align = "right"),
      col_align("Q2_predict", 13, align = "right")
    )
    
    for(i in 1:nrow(x1)) {
      cat2(
        "\n  ", 
        col_align(x1[i, "Name"], l + 2), 
        col_align(sprintf("%.4f", x1[i, "MAE_target"]), 13, align = "right"),
        col_align(sprintf("%.4f", x1[i, "MAE_benchmark"]), 15, align = "right"),
        col_align(sprintf("%.4f", x1[i, "RMSE_target"]), 13, align = "right"),
        col_align(sprintf("%.4f", x1[i, "RMSE_benchmark"]), 15, align = "right"),
        col_align(sprintf("%.4f", x1[i, "Q2_predict"]), 13, align = "right")
      )
    } 
  }
  
  cat2("\n", rule2(type = 2))
}