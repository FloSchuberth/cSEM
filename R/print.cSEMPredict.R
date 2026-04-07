#' `cSEMPredict` method for `print()`
#'
#' The `cSEMPredict` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#' 
#' @param .metrics Character string or a vector of character strings. 
#'   Which prediction metrics should be displayed? One of: "*MAE*", "*RMSE*", "*Q2*", 
#'   "*MER*", "*MAPE*, "*MSE2*", "*U1*", "*U2*", "*UM*", "*UR*", or "*UD*". 
#'   Default to c("*MAE*", "*RMSE*", "*Q2*").
#'
#' @seealso [csem()], [cSEMResults], [predict()]
#'
#' @export
#' @keywords internal
print.cSEMPredict <- function(x, 
                              .metrics = c("MAE", "RMSE", "Q2"),...) {
  
  diff <- setdiff(.metrics, args_default(TRUE)$.metrics)
  
  if(length(diff) != 0) {
    stop2(
      "The following error occured in the predict() function:\n",
      "Unknown approach: ", paste0(diff, collapse = ", "), ".",
      " Possible choices are: ",
      paste0(args_default(TRUE)$.metrics, collapse = ", "))
  }
  
  
  metrics_to_display <- .metrics
  
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
    x3 <- x$Information$Estimator_benchmark
    
    # cat2("\n\tGeneral information:\n\t","------------------------")
    cat2(
      ansi_align("\n\tNumber of obs. training", 35), "= ", x2$Number_of_observations_training,
      ansi_align("\n\tNumber of obs. test", 35), "= ", x2$Number_of_observations_test,
      ansi_align("\n\tNumber of cv folds", 35), "= ", x2$Number_of_folds,
      ansi_align("\n\tNumber of repetitions", 35), "= ", x2$Number_of_repetitions,
      ansi_align("\n\tHandle inadmissibles", 35), "= ", x2$Handle_inadmissibles,
      ansi_align("\n\tEstimator target", 35), "= ", paste0("'", x2$Estimator_target, "'"),
      ansi_align("\n\tEstimator benchmark", 35), "= ", paste0("'", x2$Estimator_benchmark, "'"),
      ansi_align("\n\tDisattenuation target", 35), "= ", paste0("'", x2$Disattenuation_target, "'"),
      ansi_align("\n\tDisattenuation benchmark", 35), "= ", paste0("'", x2$Disattenuation_benchmark, "'"),
      ansi_align("\n\tApproach to predict", 35), "= ", paste0("'", x2$Approach_to_predict, "'")
    )
    
    ### Prediction metricts-------------------------------------------------------
    cat2("\n\n", rule2("Prediction metrics"), "\n\n")
    
    l <- max(nchar(x1[, "Name"]))
    
    cat2(
      "\n  ", 
      ansi_align("Name", l + 2)
    )
    
    if(x3 != "NA"){
    if(any(metrics_to_display == "MAE")){
      cat2(
        ansi_align("MAE target", 13, align = "right"), 
        ansi_align("MAE benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "RMSE")){
      cat2(
        ansi_align("RMSE target", 13, align = "right"), 
        ansi_align("RMSE benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "Q2")){
      cat2(
        ansi_align("Q2_predict", 13, align = "right")
      )
    }  
    if(any(metrics_to_display == "MER")){
      cat2(
        ansi_align("MER target", 13, align = "right"), 
        ansi_align("MER benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "MAPE")){
      cat2(
        ansi_align("MAPE target", 13, align = "right"), 
        ansi_align("MAPE benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "MSE2")){
      cat2(
        ansi_align("MSE2 target", 13, align = "right"), 
        ansi_align("MSE2 benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "U1")){
      cat2(
        ansi_align("U1 target", 13, align = "right"), 
        ansi_align("U1 benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "U2")){
      cat2(
        ansi_align("U2 target", 13, align = "right"), 
        ansi_align("U2 benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "UM")){
      cat2(
        ansi_align("UM target", 13, align = "right"), 
        ansi_align("UM benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "UR")){
      cat2(
        ansi_align("UR target", 13, align = "right"), 
        ansi_align("UR benchmark", 15, align = "right")
      )
    }
    if(any(metrics_to_display == "UD")){
      cat2(
        ansi_align("UD target", 13, align = "right"), 
        ansi_align("UD benchmark", 15, align = "right")
      )
    }
    }else if(x3 == "NA"){
      if(any(metrics_to_display == "MAE")){
        cat2(
          ansi_align("MAE target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "RMSE")){
        cat2(
          ansi_align("RMSE target", 13, align = "right")
        )
      }

      if(any(metrics_to_display == "MER")){
        cat2(
          ansi_align("MER target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "MAPE")){
        cat2(
          ansi_align("MAPE target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "MSE2")){
        cat2(
          ansi_align("MSE2 target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "U1")){
        cat2(
          ansi_align("U1 target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "U2")){
        cat2(
          ansi_align("U2 target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "UM")){
        cat2(
          ansi_align("UM target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "UR")){
        cat2(
          ansi_align("UR target", 13, align = "right")
        )
      }
      if(any(metrics_to_display == "UD")){
        cat2(
          ansi_align("UD target", 13, align = "right")
        )
      }
    }


    if(x3 != "NA"){
      for(i in 1:nrow(x1)){
        cat2(
          "\n  ", 
          ansi_align(x1[i, "Name"], l + 2))
          if(any(metrics_to_display == "MAE")){ 
            cat2(ansi_align(sprintf("%.4f", x1[i, "MAE_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "MAE_benchmark"]), 15, align = "right"))
            }
          if(any(metrics_to_display == "RMSE")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "RMSE_target"]), 13, align = "right"),
            ansi_align(sprintf("%.4f",x1[i, "RMSE_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "Q2")){
            cat2(ansi_align(sprintf("%.4f",x1[i, "Q2_predict"]), 13, align = "right"))
          }
          if(any(metrics_to_display == "MER")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "MER_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "MER_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "MAPE")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "MAPE_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "MAPE_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "MSE2")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "MSE2_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "MSE2_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "U1")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "U1_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "U1_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "U2")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "U2_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "U2_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "UM")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "UM_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "UM_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "UR")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "UR_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "UR_benchmark"]), 15, align = "right"))
          }
          if(any(metrics_to_display == "UD")){
            cat2(ansi_align(sprintf("%.4f", x1[i, "UD_target"]), 13, align = "right"),
          ansi_align(sprintf("%.4f",x1[i, "UD_benchmark"]), 15, align = "right"))
          }
      }
    }else if(x3 == "NA"){
      for(i in 1:nrow(x1)){
        cat2(
          "\n  ", 
          ansi_align(x1[i, "Name"], l + 2))
        if(any(metrics_to_display == "MAE")){ 
          cat2(ansi_align(sprintf("%.4f", x1[i, "MAE_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "RMSE")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "RMSE_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "MER")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "MER_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "MAPE")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "MAPE_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "MSE2")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "MSE2_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "U1")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "U1_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "U2")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "U2_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "UM")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "UM_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "UR")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "UR_target"]), 13, align = "right"))
        }
        if(any(metrics_to_display == "UD")){
          cat2(ansi_align(sprintf("%.4f", x1[i, "UD_target"]), 13, align = "right"))
        }
      }
    }    
  }
  
  cat2("\n", rule2(type = 2), "\n")
}