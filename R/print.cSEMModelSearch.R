#' `cSEMModelSearch` method for `print()`
#'
#' The `cSEMModelSearch` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [doModelSearch()]
#'
#' @export
#' @keywords internal
print.cSEMModelSearch <- function(x, ...) {
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview")
  )
  
  if(inherits(x, "cSEMModelSearch_multi")) {
    for(j in seq_along(x)){
      cat2("\nGroup:", names(x)[j], "\n")
      printModelSearchDetails(.x = x[[j]])
    }
    
  } else {
    printModelSearchDetails(.x = x)
  }
  

  cat(rule2(type = 2), "\n")
}




printModelSearchDetails <- function(.x) {

  cat2("\n\tAGAS-PLS information:\n\t","------------------------")
  cat2(
    col_align("\n\tPopulation size", width = 36), "= ", 
    .x$Information$pop_size,
    col_align("\n\tNumber of generations", width = 36), "= ", 
    .x$Information$n_generation,
    col_align("\n\tMutation probability", width = 36), "= ", 
    .x$Information$prob_mutation,
    col_align("\n\tCrossover probability", width = 36), "= ", 
    .x$Information$prob_crossover,
    col_align("\n\tPenalize value (fbar)", width = 36), "= ", 
    format(.x$Information$fbar, scientific = FALSE),
    col_align("\n\tModel selection criterion used", width = 36), "= ", 
    .x$Information$model_selection_criterion,
    col_align("\n\tChosen exogenous construct(s)", width = 36), "= ", 
    paste(.x$Information$cons_exo, collapse = " "),
    col_align("\n\tSeed", width = 36), "= ", 
    if(is.null(.x$Information$seed)){
      'none'
    }else{
      .x$Information$seed
    }
  )
  
  cat2('\n')
  
  cat2("\n\tAGAS-PLS results:\n\t","------------------------")
  cat2(
    col_align("\n\tFitness value original model", width = 36), "= ", 
    .x$Results$original_fitness,
    col_align("\n\tBest fitness value", width = 36), "= ",
    .x$Results$best_fitness,
    col_align("\n\tNumber of best fitting models", width = 36), "= ",
    length(.x$Results$best_matrix),
    col_align("\n\n\tStructural model structure:")
  )
  
  mat_output <- capture.output(print(.x$Results$best_matrix[[1]]))
  indented_mat <- paste0("\t", mat_output, collapse = "\n")
  cat2("\n", indented_mat, "\n")
  
  if(length(.x$Results$best_matrix)>1){
    cat2("\n\tNote: ", length(.x$Results$best_matrix), " models with the same fitness value found.\n")
  }
  
  if(.x$Results$original_fitness < .x$Results$best_fitness){
    cat2(col_align(red("\n The best fitting model shows a worse fit than the original one.")),
         red("\n Try to increase `.pop_size` and/or `.n_generations`.\n"))
  }

}

