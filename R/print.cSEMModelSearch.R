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
  cat2("\n\tAGAS-PLS information:\n\t","------------------------")
  cat2(
    col_align("\n\tPopulation size", width = 36), "= ", 
    x$Information$pop_size,
    col_align("\n\tNumber of generations", width = 36), "= ", 
    x$Information$n_generation,
    col_align("\n\tMutation probability", width = 36), "= ", 
    x$Information$prob_mutation,
    col_align("\n\tCrossover probability", width = 36), "= ", 
    x$Information$prob_crossover,
    col_align("\n\tPenalize value (fbar)", width = 36), "= ", 
    format(x$Information$fbar, scientific = FALSE),
    col_align("\n\tModel selection criterion used", width = 36), "= ", 
    x$Information$model_selection_criterion
  )
  
  cat2('\n')
  
  cat2("\n\tOriginal model:\n\t","------------------------")
  cat2(col_align("\n\tFitness value", width = 36), "= ", 
       x$original_fitness)
  
  cat2('\n')
  
  cat2("\n\tAGAS-PLS results:\n\t","------------------------")
  cat2(
    col_align("\n\tFitness value", width = 36), "= ", 
    x$best_fitness,
    col_align("\n\tNumber of models found", width = 36), "= ", 
    length(x$best_matrix),
    col_align("\n\tStructure of the model found", width = 36), "=\n ", 
    x$bext_matrix[[1]]
  )
  
  
}
